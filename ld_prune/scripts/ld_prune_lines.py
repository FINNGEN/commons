#!/usr/bin/env python3
import argparse
import gzip
import sys
import time
import requests


if sys.version_info[0] != 3 or sys.version_info[1]<7 or (sys.version_info[1]==7 and sys.version_info[2]<2):
    raise Exception(" Python Version >= 3.7.2 required")

from functools import partial
from collections import namedtuple
from collections import defaultdict
from collections import OrderedDict
from queue import PriorityQueue
from io import TextIOBase
from typing import Tuple, List, OrderedDict, Dict, Iterator
from scipy.stats import chi2
import math

def get_ld_vars( chrom:str, pos:int, ref:str, alt:str, r2:float, ld_w:int, ld_source, retries:int=5) -> Dict[str,str]:
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

    return(r.json())


class Hit:

    def __init__(self, chrom:str, pos:int, ref:str, alt:str, priority:float, data: OrderedDict):
        '''
            Args: data row data with keys are the column names and values. Dictionary will be stored and data returned in original order
        '''
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.variant = f'{chrom}_{pos}_{ref}_{alt}'
        self.priority = priority

        gettable = [ t for t in data.items() ]
        gettable.extend([ ("var",self.variant), ("chrom", chrom), ("pos",pos), ("ref", ref), ("alt", alt),("priority",priority) ] )
        self.gettable = { k:v for (k, v) in  gettable }

        self.linedata = [v for v in data.values()]

    @property
    def varid(self):
        return self.variant

    @property
    def data(self):
        return self.linedata

    def __str__(self):
        return self.variant

    def __getitem__(self,key):
        return self.gettable[key]

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return other.data==self.data

    def __lt__(self, other):
        self.priority < other.priority

    def __gt__(self, other):
        self.priority > other.priority

def read_hit(f, h:OrderedDict, args) -> Iterator[Hit]:
    '''
        Args:
        h: header name to index dict in the order they are in the file
    '''
    line=2
    while True:
        l = f.readline()
        if l =="":
            break

        dat = l.strip("\n").split("\t")

        chrom = dat[h[args.chromcol]]
        chrom = chrom if chrom.startswith("chr") else "chr"+chrom
        pos = int(dat[h[args.poscol]])
        priority = float(dat[h[args.pcol]])

        ref = dat[h[args.refcol]]
        alt = dat[h[args.altcol]]
        if(len(dat)!=len(h.keys())):
            raise Exception(f'Not enought elements in line {line}, excepted {len(h.keys())} observed {len(dat)}')

        data = OrderedDict([ (k,dat[v]) for (k,v) in h.items()])
        line+=1
        yield Hit(chrom, pos,ref, alt, priority, data)


class Cluster(object):

    def __init__(self):
        self.hits = PriorityQueue()
        self.hit_cache=defaultdict(list)

        self.start = float('inf')
        self.end = -float('inf')
        self.n_hits = 0
        self.removed = {}
        self.peak_priority = float('inf')


    def add_hit(self, hit:Hit):

        if hit.priority<self.peak_priority:
            self.peak_priority=hit.priority

        self.hit_cache[ hit.varid ].append(hit)
        self.hits.put( (hit.priority, hit) )

        if hit.pos > self.end:
            self.end = hit.pos

        if hit.pos < self.start:
            self.start = hit.pos

        self.n_hits += 1




    @property
    def boundaries(self):
        return (self.start,self.end)

    def width(self):
        return self.end - self.start

    def remove_hits(self, var_ids:List[Tuple[str,float]], r2:float, chisqtop:float, obs_exp_chi2_thr=None,
                        clump_expected_chisq_filter_af_col:str=None, af_threshold_chi_prune:float=None ) -> List[Hit]:
        """
            Removes variants in cluster if r2 is higher than threshold or if observed chisq is attributable to
            stronger variants chisq
            var_ids: tuple of variant_id, r2
            r2: r2 threshold
            chisqtop: strongest variant chisq in the cluster
        """
        rm_hits = []
        for vi in var_ids:

            if vi[0] in self.hit_cache:

                eligible_for_chi_prune=True
                delvar_af=0
                if obs_exp_chi2_thr is not None:
                    ## this option should only be used for single pheno files so only one of each variant can exist.
                    if len(self.hit_cache[vi[0]])>1:
                        raise Exception("You should only use the --clump_expected_chisq to single pheno files (no replicated variants)")
                    delvar = self.hit_cache[vi[0]][0]

                    if clump_expected_chisq_filter_af_col is not None:
                        delvar_af=float(delvar[clump_expected_chisq_filter_af_col])
                        if delvar_af>af_threshold_chi_prune:
                            eligible_for_chi_prune=False

                    chisq = chi2.isf(delvar.priority, df=2)

                if vi[1]>r2 or (eligible_for_chi_prune and obs_exp_chi2_thr is not None and (chisq/chisqtop)>=obs_exp_chi2_thr):
                    self.removed[vi[0]]=""
                    rm_hits.extend(self.hit_cache[vi[0]])
                    del self.hit_cache[vi[0]]

        self.n_hits-=len(rm_hits)

        return rm_hits

    def size(self)->int:
        return self.n_hits

    def empty(self) -> bool:
        return self.n_hits==0

    def get_best(self)-> (Hit, List[Hit]):
        """
            returns best hit in first element and all other hits in the same variant in second element and removes
            all of them from cluster
        """
        if self.n_hits==0:
            return None

        b = self.hits.get(block=False)[1]

        while self.hits.qsize()>0 and b.varid in self.removed:
            b = self.hits.get(block=False)[1]

        self.n_hits-=1

        merged = [ h for h in self.hit_cache[b.varid] if h != b]

        self.n_hits-=len(merged)
        self.hit_cache.pop(b.varid,None)

        # mark this variant as removed from the cluster so we know to skip if they arrive
        # from priorityqueu
        self.removed[b.varid]=""

        return (b,merged)

    def get_all(self) -> List[Hit]:
        hits = []
        while not self.empty():
            h = self.get_best()
            hits.append(h[0])
            hits.extend(h[1])
            self.hit_cache.pop(h[0].varid,None)

        return hits

    def clear(self):
        self.hits = PriorityQueue()
        self.start = float('inf')
        self.end = -float('inf')
        self.n_hits = 0
        self.hit_cache=defaultdict(list)
        self.removed = {}
        self.peak_priority = float('inf')

        #get_ld_vars(chrom, var[1], var[2], var[3], args.ld, args.ld_w)

    def __str__(self):
        cls = "\n".join([ v.varid + ":" + format(v.priority) for varhits in self.hit_cache.values() for v in varhits ])
        return (f'Cluster size {self.size()}, width {self.width()}. \n'
                f'Clusterhits:\n { cls }'
                )


def prune_cluster(cl:Cluster, r2:float, n_retry:int, ld_source,clump_expected_chisq:float=None,
                    clump_expected_chisq_filter_af_col:str=None, af_threshold_chi_prune:float=None, fixed_ld_search_width=None) ->List[Tuple[Hit,List[Hit]]]:
    '''
        Takes cluster and prunes by LD of Hits in cluster, retaining the top in in the clustered
        Args:
            cl:
            r2: threshold above which to cluster Hits
            n_retry: number of times to retry error in ld server
        Returns: list of tuples where first element is the lead Hit of the cluster and second is list of clustered Hits
    '''

    if cl.size()==1:
        print("size 0")
        return [(cl.get_best()[0],[])]
    ### self.queue[0]
    outdat=[]

    redhit=0

    while not cl.empty():
        print(cluster)
        best = cl.get_best()
        worse_same = best[1]
        top = best[0]

        print(f'Current {top}')
        if top is None:
            return outdat

        if cl.width()==0 or cl.size()==0:
            ## all same lead. cluster to best
            worse_same.extend([ h for h in cl.get_all() ])
            outdat.append( (top, worse_same ))
            return outdat

        if args.min_region:
            range = hit[arg.min_region].split(":").split("-")
            min_width = 2 * max( abs(pos-range[0]), abs(pos-range[1]) ) + 100

        ## ld server is remarkably inexact on the range so add a lot of buffer.
        safety_buffer=50000

        if fixed_ld_search_width:
            width=fixed_ld_search_width
        else:
            width = max(abs(top.pos-cl.boundaries[0]-safety_buffer)*2, abs(top.pos-cl.boundaries[1]+safety_buffer)*2 )

        search_r2=r2
        chisqtop=chi2.isf(cl.peak_priority, df=1)

        if clump_expected_chisq:
            ## r2 needed for nominal significance
            #search_r2 = 5/chisqtop
            # we need all ld variants for chisq clumping.
            search_r2=0

        ld = get_ld_vars(top.chrom, top.pos, top.ref, top.alt, search_r2,
            width, ld_source=ld_source,retries=n_retry )

        varid = ":".join( [top.chrom,str(top.pos),top.ref,top.alt] )
        ldvars = [ ("chr" + ldv["variation2"].replace(":","_"), float(ldv["r2"])) for ldv in ld["ld"] if ldv["variation2"]!=varid]

        removed = cl.remove_hits( ldvars, r2, chisqtop, clump_expected_chisq, clump_expected_chisq_filter_af_col, af_threshold_chi_prune=af_threshold_chi_prune )

        print(f"removed {removed}")
        removed.extend(worse_same)
        print(f'merged: {len(removed)} hits to top. Top: {top.varid}, removed {[ h.varid for h in removed]}')
        redhit+=1
        outdat.append( (top, removed) )

    return outdat

def write_cluster(pruned:list[tuple[Hit,List[Hit]]], outcols:list[str] ,out:TextIOBase):
    for h in pruned:
        out.write("\t".join(h[0].data) + "\t" + ",".join([ ";".join([pr[c] for c in outcols]) for pr in h[1] ]) + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', action='store', type=str, help='')
    parser.add_argument('outfile', action='store', type=str, help='')
    parser.add_argument('-pcol', default="lead_pval", help="Priority columns, smaller is more important. Data needs to be convertable to float")
    parser.add_argument('-ld_source', default="sisu4")
    parser.add_argument('-chromcol', default="chrom")
    parser.add_argument('-poscol', default="pos")
    parser.add_argument('-refcol', default="ref")
    parser.add_argument('-altcol', default="alt")
    parser.add_argument('-n_retries_ld', type=int, default=5)
    parser.add_argument('-ld', type=float, default=0.5)
    parser.add_argument('-fixed_ld_search_width', type=int, help="Don't try to optimize the search width based on cluster min/max bp position but use fixed width. LD server can be very innacurate in width")
    parser.add_argument('-clump_expected_chisq', type=float,
        help="Use only for single phenotype file pruning hits caused by residual LD from stronger signal !!! Clumps variants if based on LD the observed variant chisq is this much attributable to excepted LD from stronger hit.")
    parser.add_argument('-clump_expected_chisq_af', type=float, default=0.01,help="AF limit where rarer than this variants are subject to clump_expected_chisq")
    parser.add_argument('-clump_expected_chisq_filter_af_col', type=str, help="AF column if clump_expected_chisq applies only to low freq variants ")
    parser.add_argument('-ld_w', default=500000, type=int, help='How close hits are first clustered together for LD based pruning. Cluster region can become larger as this is between adjacent hits.')
    parser.add_argument('-prune_column_list', default="phenocode", help="Comma separated list of column names that will be printed out in the last column where pruned endpoints are listed")
    parser.add_argument('-min_region', help="column with chr:start-stop of minimum region to include in ld search for each hit ")

    args = parser.parse_args()
    of = gzip.open if args.file.endswith(".gz") else open

    with of(args.file, 'rt') as infile:

        h = OrderedDict([ (hd,i) for i,hd in enumerate(infile.readline().strip().split("\t")) ])

        reqs = [args.pcol,
            args.poscol, args.chromcol, args.refcol, args.altcol]

        pruned_cols = args.prune_column_list.split(",")
        reqs.extend(pruned_cols)

        if args.min_region:
            reqs.extend(args.min_region)

        missing = [r for r in reqs if not r in h]
        if len(missing)>0:
            raise Exception(f"Given columns {missing} not in file ")

        seen_chroms = {}
        progress_rep=100

        with open(args.outfile,"wt") as out:
            out.write("\t".join( list(h.keys())+ ["ld_pruned_phenos"] ) + "\n")
            n_vars=0
            progress_rep_step=10
            cluster = Cluster()
            nexthit = read_hit(infile, h, args)
            hit = next(nexthit)
            cluster.add_hit( hit )
            last_hit = hit
            hits_proc=0
            for hit in nexthit:
                ## read cluster of hits and take top one and remove all in LD.....
                ## spit to top and removed and continue until cluster empty

                if n_vars % progress_rep_step == 0 and n_vars>0:
                   print("Retrieved LD for {} variants".format( n_vars) ,file=sys.stderr)
                if (hit.chrom != last_hit.chrom and hit.chrom in seen_chroms) or (hit.chrom == last_hit.chrom and hit.pos < last_hit.pos):
                    raise Exception(f"Data must be sorted by chromosome and column before running. Offending row{hit.data}")

                if hit.chrom != last_hit.chrom or hit.pos-last_hit.pos > args.ld_w:
                    print("#### starting prune cluster")
                    pruned = prune_cluster(cluster, args.ld, n_retry=args.n_retries_ld,ld_source=args.ld_source,
                        clump_expected_chisq=args.clump_expected_chisq, clump_expected_chisq_filter_af_col=args.clump_expected_chisq_filter_af_col,
                        af_threshold_chi_prune=args.clump_expected_chisq_af, fixed_ld_search_width=args.fixed_ld_search_width)
                    print(f'#### Cluster pruned to {len(pruned)} hits')
                    write_cluster(pruned,pruned_cols, out)
                    cluster.clear()

                hits_proc+=1
                if hits_proc % progress_rep == 0:
                    print(f'{hits_proc} hits processed!')
                cluster.add_hit(hit)
                last_hit = hit


            if cluster.size()>0:
                # write last cluster not pruned yet
                print(f"#### starting prune cluster {cluster}")
                pruned = prune_cluster(cluster, args.ld, n_retry=args.n_retries_ld, clump_expected_chisq=args.clump_expected_chisq)
                print(f'#### Cluster pruned to {len(pruned)} hits')
                write_cluster(pruned,pruned_cols, out)
