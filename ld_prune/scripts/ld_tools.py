from collections import defaultdict
import imp, logging, subprocess, shlex, re, sys, time, gzip, pysam, threading
from subprocess import CalledProcessError
import requests

import tempfile
from typing import Dict
from abc import abstractmethod

LD_DISPL_DECIMALS = 3

def get_region_mapping(chrom:str,pos:str,ref:str,alt:str, mapfile:str, window:int):
    """
    Gets tomahawk position for the query variant and a tomahawk_position-to-variant mapping for variants in the given panel within the given window
    Returns a tuple (tomahawk chr:pos, dict from tomahawk chr:pos to chr:pos:ref:alt)
    Raises if variant not found
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
        raise Exception({'variant not found: ' + cpra})
    
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
                         tomahawk_threads=None, limit_to_vars:dict[str,str]=None) -> Dict[str,str]:

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


    cmd_scalc = 'tomahawk scalc ' + opts + ' -i ' + tomafile + ' -o ' +  tf.file.name + ' -I ' + tomahawk_chrpos[0] + ':' + str(int(tomahawk_chrpos[1])-1) + '-' + tomahawk_chrpos[1] + ' -w ' + searchwidth
    out = subprocess.check_output(shlex.split(cmd_scalc)).decode(sys.stdout.encoding).strip()
    cmd_view = 'tomahawk view -i ' + tf.file.name
    try:
        out = subprocess.check_output(shlex.split(cmd_view)).decode(sys.stdout.encoding).strip()
    except:
        print("ERROR viewing results." + out, file=sys.stderr )

    cpra = "{}:{}:{}:{}".format(chrom, pos, ref, alt)
    return parse_ld(out, cpra, r2, map[1])


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









