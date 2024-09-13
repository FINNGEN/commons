from collections import defaultdict
import imp, logging, subprocess, shlex, re, sys, time, gzip, pysam, threading
from subprocess import CalledProcessError
import requests

import os
import tempfile
from typing import Dict
from collections.abc import Callable
from abc import abstractmethod
from functools import partial
import argparse

LD_DISPL_DECIMALS = 3

class VariantNotFoundException(Exception):
    pass


def get_region_mapping(chrom:str,pos:str,ref:str,alt:str, mapfile:str, window:int):
    """
    Gets tomahawk position for the query variant and a tomahawk_position-to-variant mapping for variants in the given panel within the given window
    Returns a tuple (tomahawk chr:pos, dict from tomahawk chr:pos to chr:pos:ref:alt)
    Raises VariantNotFoundException if variant not found
    """
    twk2cpra = {} # mapping from tomahawk positions chr:pos to actual variants chr:pos:ref:alt
    twk = None # tomahawk position of query variant
    
    tb = pysam.TabixFile(mapfile, parser=None)

    cpra = "{}:{}:{}:{}".format(chrom, pos, ref, alt)
    tabix_iter = tb.fetch(chrom, max(1,int(pos)-window-5), int(pos)+window+2, parser=None)
    for row in tabix_iter:
        s = row.split('\t')
        if s[2] == cpra:
            twk = s[3]
        twk2cpra[s[3]] = s[2]
    if twk is None:
        raise VariantNotFoundException({'variant not found: ' + cpra})
    
    return (twk, twk2cpra)

def parse_ld(data:str, cpra:str, r2_thresh:float, twk2cpra:list[Dict[str,str]], limit_to_vars:dict[str,str]=None):
    """
    Parses the given tomahawk LD output using the given query variant and tomahawk_position-to-variant mapping
    Returns a list of dicts with keys variation1,variation2,r2,d_prime where variation1 is the query variant
    """
    data_iter = iter(data.split('\n'))
    for line in data_iter:
        if not line.startswith('#') and line != '':
            break
    hdr = {h:i for i,h in enumerate(line.strip().split('\t'))}
    res = []
    used = {}
    for line in data_iter:
        if line=="":
            continue

        s = line.strip().split('\t')
        if float(s[hdr['R2']]) < r2_thresh:
            continue
        
        var1 = s[hdr['ridA']] + ':' + s[hdr['posA']]
        var2 = s[hdr['ridB']] + ':' + s[hdr['posB']]


        if var1 not in twk2cpra:
            print(var1 + ' tomahawk position not in given mapping (query variant ' + cpra + '), this is an issue only if the position is not at the boundary of the window. Ignoring position', file=sys.stderr)
        elif var2 not in twk2cpra:
            print(var2 + ' tomahawk position not in given mapping (query variant ' + cpra + '), this is an issue only if the position is not at the boundary of the window. Ignoring position', file=sys.stderr)
        else:
            var1 = twk2cpra[var1]
            var2 = twk2cpra[var2]
            if var2 == cpra:
                temp = var1
                var1 = var2
                var2 = temp
                
            if var1 == cpra and var2 not in used and (not limit_to_vars or var2 in limit_to_vars):
                res.append({'variation1': var1, 'variation2': var2, 'r2': round(float(s[hdr['R2']]), LD_DISPL_DECIMALS), 'd_prime': round(float(s[hdr['Dprime']]), LD_DISPL_DECIMALS)})
                used[var2] = True
    return res

@abstractmethod
def get_ld_vars(chrom:str, pos:int, ref:str, alt:str, r2:float, ld_w:int) -> Dict[str,str]:
    pass

def get_ld_vars_tomahawk(chrom, pos, ref, alt, r2:float, ld_w:int,tomafile:str, mapfile:str,
                         tomahawk_threads=None ,limit_to_vars:dict[str,str]=None) -> Dict[str,str]:
    """
        Gets LD data for the given variant using tomahawk. 
        Throws VariantNotFoundException if variant not found
    """

    if(chrom.startswith("chr")):
        chrom = chrom.replace("chr","")


    map = get_region_mapping(chrom,pos, ref, alt,mapfile,ld_w)      

    tomahawk_chrpos = map[0].split(':')
    tf = tempfile.NamedTemporaryFile(suffix=".two")
    
    tomafile = tomafile.replace("{CHR}",chrom)

    searchwidth=str(ld_w)
            
    var = tomahawk_chrpos[0] + ':' + str(int(tomahawk_chrpos[1])-1) 
    print("searching from: " + var + " width "+ searchwidth )
    
    opts = ""
    if tomahawk_threads:
        opts = " -t " + str(tomahawk_threads)

    out =tf.file.name
    print(f'tomahawk scalc {opts} -i {tomafile} -o {out} -I {tomahawk_chrpos[0]}:{int(tomahawk_chrpos[1])-1}-{tomahawk_chrpos[1]} -w {searchwidth}', file=sys.stderr )
    cmd_scalc = 'tomahawk scalc ' + opts + ' -i ' + tomafile + ' -o ' + out + ' -I ' + tomahawk_chrpos[0] + ':' + str(int(tomahawk_chrpos[1])-1) + '-' + tomahawk_chrpos[1] + ' -w ' + searchwidth

    pr = subprocess.run(
                shlex.split(cmd_scalc),
                stdout=subprocess.PIPE,
                encoding="ascii",
            )
      
    if pr.returncode!=0:
        raise Exception(f'Error running tomahawk scalc: {pr.stderr} {pr.stdout}')
    cmd_view = 'tomahawk view -i ' + out
    pr = subprocess.run(shlex.split(cmd_view), stderr=subprocess.PIPE, 
                             stdout=subprocess.PIPE, encoding="ascii")
    if pr.returncode!=0:
        raise Exception(f'Error running tomahawk view: {pr.stderr} {pr.stdout}')

    cpra = "{}:{}:{}:{}".format(chrom, pos, ref, alt)
    return parse_ld(pr.stdout, cpra, r2, map[1])


def get_ld_vars_api( chrom:str, pos:int, ref:str, alt:str, r2:float, ld_w:int, ld_source,
                    retries:int=5) -> Dict[str,str]:
    snooze=2
    snooze_multiplier = 2
    max_snooze=256

    ld_w=max(min(5000000,ld_w), 500000)

    url = f'http://api.finngen.fi/api/ld?variant={chrom}:{pos}:{ref}:{alt}&panel={ld_source}&window={ld_w}&r2_thresh={r2}'
    print(f'requesting LD {url}')
    r = requests.get(url)
    retries_left = retries
    while r.status_code!=200 and retries_left>0:
        retries_left-=1
        if r.status_code!=200:
            print("Error requesting ld for url {}. Error code: {}".format(url, r.status_code) ,file=sys.stderr)
            time.sleep(snooze)
            snooze = min(snooze * snooze_multiplier, max_snooze)

        r = requests.get(url)

    if r.status_code!=200:
        raise Exception("LD server response failed for {} attempts".format(retries) )

    return(r.json()["ld"])



def get_common_LD_API_arg_parser(parser:argparse.ArgumentParser):
    '''
        Adds common arguments for LD API to the given parser
    '''

    parser.add_argument('-ld', type=float, default=0.2)
    parser.add_argument('-ld_source', default="sisu42")
    parser.add_argument('-local_tomahawk_LD', action='store_true', help=' Tomahawk must be locally installed and FinnGen produced tomhawk LD  and variant mapping must be available. See: https://github.com/FINNGEN/ld_server for generating')
    parser.add_argument('-local_tomahawk_threads', type=int)
    parser.add_argument('-max_ld_width',type=int, help="restrict max lad width search to this.",)
    parser.add_argument('-fixed_ld_search_width', type=int, help="Don't try to optimize the search width based on cluster min/max bp position but use fixed width. LD server can be very innacurate in width")
    parser.add_argument('-clump_expected_chisq', type=float,
        help="Use only for single phenotype file pruning hits caused by residual LD from stronger signal !!! Clumps variants if based on LD the observed variant chisq is this much attributable to excepted LD from stronger hit.")
    parser.add_argument('-clump_expected_chisq_af', type=float, default=0.01,help="AF limit where rarer than this variants are subject to clump_expected_chisq")
    parser.add_argument('-clump_expected_chisq_filter_af_col', type=str, help="AF column if clump_expected_chisq applies only to low freq variants ")
    parser.add_argument('-ld_w', default=500000, type=int, help='How close hits are first clustered together for LD based pruning. Cluster region can become larger as this is between adjacent hits.')
    parser.add_argument('-min_region', help="column with chr:start-stop of minimum region to include in ld search for each hit ")
    parser.add_argument('-n_retries_ld', type=int, default=5)
    
    main_args, _ = parser.parse_known_args()
    if main_args.local_tomahawk_LD:
        parser.add_argument("-tomahawk_template", type=str, required=main_args.local_tomahawk_LD, help="Template to tomahawk files where chromosome number is replaced with {CHR}")
        parser.add_argument("-tomahawk_mapfile", type=str, required=main_args.local_tomahawk_LD, help="Path to location map file")
        args = parser.parse_args()
       
    else:
        args = parser.parse_args()
        

    return parser


def get_ld_api_by_cmdargs(args) -> Callable[[str, int, str, str, float, int], Dict[str,str]]:
    '''
        Returns a function that can be used to get LD data based on the given command line arguments.
        Function implements the signature of get_ld_vars
    '''

    if args.local_tomahawk_LD:
        if not os.path.exists(args.tomahawk_mapfile):
            raise Exception("tomahawk map file not found")
        ld_interface = partial(get_ld_vars_tomahawk, tomafile=args.tomahawk_template, mapfile=args.tomahawk_mapfile, 
                               tomahawk_threads=args.local_tomahawk_threads)
    else:
        ld_interface = partial(get_ld_vars_api, ld_source=args.ld_source, retries=args.n_retries_ld)
    
    return ld_interface






