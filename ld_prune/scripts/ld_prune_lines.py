#!/usr/bin/env python3
import argparse
import gzip
import requests
import sys
import time

from functools import partial
from collections import namedtuple
from collections import defaultdict
from collections import OrderedDict

from queue import PriorityQueue

from io import TextIOBase

def get_ld_vars( chrom, pos, ref, alt, r2, ld_w, retries=5):
    snooze=2
    print("requesting LD")
    api="http://api.finngen.fi/api/ld?variant={}:{}:{}:{}&panel=sisu3&window={}&r2_thresh={}"
    url = api.format(chrom,pos,ref,alt,ld_w,r2)
    r = requests.get(url)
    while r.status_code!=200 or retries>0:
        if r.status_code!=200:
            print("Error requesting ld for url {}. Error code: {}".format(url, r.status_code) ,file=sys.stderr)
            time.sleep(snooze)
        retries-=1
        r = requests.get(url)

    if r.status_code!=200:
        raise Exception("LD server response failed for {} attempts".format(retries) )

    print("got LD")
    return(r.json())

Hit = namedtuple('Hit', 'var, chrom, pos, ref, alt, priority, data')

class Hit:

    def __init__(self, chrom, pos, ref, alt, priority, data: OrderedDict):
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

    def __getitem__(self,key):
        return self.gettable[key]


    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        return other.data==self.data


def read_hit(f, h:OrderedDict, args):
    '''
        Args:
        h: header name to index dict in the order they are in the file
    '''

    while True:
        l = f.readline()
        if l =="":
            break

        dat = l.strip().split("\t")

        chrom = dat[h[args.chromcol]]
        chrom = chrom if chrom.startswith("chr") else "chr"+chrom
        pos = int(dat[h[args.poscol]])
        priority = float(dat[h[args.pcol]])

        ref = dat[h[args.refcol]]
        alt = dat[h[args.altcol]]

        data = OrderedDict([ (k,dat[v]) for (k,v) in h.items()])

        yield Hit(chrom, pos,ref, alt, priority, data)



class Cluster(object):

    def __init__(self):
        self.hits = PriorityQueue()
        self.hit_cache=defaultdict(list)

        self.start =0
        self.end =0
        self.n_hits = 0
        self.removed = {}


    def add_hit(self, hit:Hit):
        if self.n_hits == 0:
            self.start = hit.pos

        print(hit.varid)
        self.hit_cache[ hit.varid ].append(hit)
        self.end = hit.pos
        self.hits.put( (hit.priority, hit) )
        self.n_hits += 1

    @property
    def boundaries(self):
        return (self.start,self.end)

    def width(self):
        return self.end - self.start

    def remove_hits(self, var_ids:list[str]) -> list[Hit]:
        rm_hits = []
        for vi in var_ids:
            if vi in self.hit_cache:
                self.removed[vi]=""
                rm_hits.extend(self.hit_cache[vi])
                del self.hit_cache.pop[var_id]
        self.n_hits-=len(rm_hits)

        return rm_hits

    def size(self)->int:
        return self.n_hits

    def empty(self) -> bool:
        return self.n_hits==0

    def get_best(self)-> (Hit, list[Hit]):
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

    def get_all(self) -> list[Hit]:
        hits = []
        while not self.empty():
            h = self.get_best()
            hits.append(h[0])
            hits.extend(h[1])
            self.hit_cache.pop(h[0].varid,None)

        return hits

    def clear(self):
        self.hits = PriorityQueue()
        self.start = 0
        self.end = 0
        self.n_hits = 0
        self.hit_cache=defaultdict(list)
        self.removed = {}

        #get_ld_vars(chrom, var[1], var[2], var[3], args.ld, args.ld_w)

    def __str__(self):
        cls = "\n".join([ v.varid + ":" + format(v.priority) for varhits in self.hit_cache.values() for v in varhits ])
        return (f'Cluster size {self.size()}, width {self.width()}. \n'
                f'Clusterhits:\n { cls }'
                )


def prune_cluster(cl:Cluster, r2:float) ->list[tuple[Hit,list[Hit]]]:
    '''
        Takes cluster and prunes by LD of Hits in cluster, retaining the top in in the clustered
        Args:
            cl:
            r2: threshold above which to cluster Hits
        Returns: list of tuples where first element is the lead Hit of the cluster and second is list of clustered Hits
    '''

    if cl.size()==1:
        return [(cl.get_best()[0],[])]
    ### self.queue[0]
    outdat=[]

    redhit=0
    while not cl.empty():
        print(cluster)
        best = cl.get_best()
        worse_same = best[1]
        top = best[0]

        if top is None:
            return outdat

        if cl.width()==0 or cl.size()==0:
            ## all same lead. cluster to best
            worse_same.extend([ h for h in cl.get_all() ])
            outdat.append( (top, worse_same ))
            return outdat

        ld = get_ld_vars(top.chrom, top.pos, top.ref, top.alt, r2, max(100000,  max(100000, abs(top.pos-cl.boundaries[0]) + 100, abs(top.pos-cl.boundaries[1])+ 100 ) )  )
        varid = ":".join( [top.chrom,str(top.pos),top.ref,top.alt] )
        ldvars = [ ldv["variation2"].replace(":","_") for ldv in ld["ld"] if ldv["variation2"]!=varid]

        removed = cl.remove_hits( ldvars )
        removed.extend(worse_same)
        print(f'merged: {len(removed)} hits to top')
        redhit+=1
        outdat.append( (top, removed) )

    return outdat

def write_cluster(pruned:list[tuple[Hit,list[Hit]]], outcols:list[str] ,out:TextIOBase):
    for h in pruned:
        out.write("\t".join(h[0].data) + "\t" + ",".join([ "_".join([pr[c] for c in outcols]) for pr in h[1] ]) + "\n")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', action='store', type=str, help='')
    parser.add_argument('outfile', action='store', type=str, help='')
    parser.add_argument('-pcol', default="lead_pval", help="Priority columns, smaller is more important. Data needs to be convertable to float")
    parser.add_argument('-chromcol', default="chrom")
    parser.add_argument('-poscol', default="pos")
    parser.add_argument('-refcol', default="ref")
    parser.add_argument('-altcol', default="alt")
    parser.add_argument('-ld', default=0.5)
    parser.add_argument('-ld_w', default=500000, help='How close hits are first clustered together for LD based pruning. Cluster region can become larger as this is between adjacent hits.')
    parser.add_argument('-prune_column_list', default="phenocode", help="Comma separated list of column names that will be printed out in the last column where pruned endpoints are listed")
    args = parser.parse_args()
    of = gzip.open if args.file.endswith(".gz") else open

    with of(args.file, 'rt') as infile:

        h = OrderedDict([ (hd,i) for i,hd in enumerate(infile.readline().strip().split("\t")) ])

        reqs = [args.pcol,
            args.poscol, args.chromcol, args.refcol, args.altcol]

        pruned_cols = args.prune_column_list.split(",")

        reqs.extend(pruned_cols)
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
                    pruned = prune_cluster(cluster, args.ld)
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
                print("#### starting prune cluster")
                pruned = prune_cluster(cluster, args.ld)
                print(f'#### Cluster pruned to {len(pruned)} hits')
                write_cluster(pruned,pruned_cols, out)
