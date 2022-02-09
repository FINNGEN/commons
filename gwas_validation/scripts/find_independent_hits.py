from typing import NamedTuple
from collections import defaultdict

RANGE=1500000

class CPRA(NamedTuple):
    c:str
    p:int
    r:str
    a:str
    pval:float

class Region(NamedTuple):
    c:str
    start:int
    end:int

def region_region_overlap(r1: Region, r2: Region) -> bool:
    if r1.c != r2.c or (r1.end < r2.start or r2.end < r1.start):
        return False
    return True

def region_variant_overlap(region: Region, variant: CPRA) -> bool:
    if region.c != variant.c or (variant.p < region.start or variant.p > region.end):
        return False
    return True

with open("~{phenotype}_5e8_previous.tsv" ) as f:
    header = f.readline().strip("\n").split("\t")
    h_idx = {a:i for i,a in enumerate(header)}
    c_idx = h_idx["chrom"]
    p_idx = h_idx["pos"]
    regions =[]
    for line in f:
        a = line.strip("\n").split("\t")
        regions.append(Region(str(a[c_idx]),int(a[p_idx])-RANGE,int(a[p_idx])+RANGE) )
    #merge
    region_chrdict = defaultdict(list)
    for r in regions:
        region_chrdict[r.c].append(r)
    for chrom in region_chrdict.keys():
        rs = region_chrdict[chrom]
        rs = sorted(rs,key=lambda x: x.start)
        stack = []
        stack.append(rs[0])
        for reg in rs[1:]:
            top = stack.pop()
            if region_region_overlap(top,reg):
                stack.append(Region(top.c,min(top.start,reg.start),max(top.end,reg.end)) )
            else:
                stack.append(top)
                stack.append(reg)
        region_chrdict[chrom]=stack

with open("~{phenotype}_1e5_current.tsv") as f:
    with open("~{phenotype}_new_hits.tsv","w") as of:
        hline=f.readline()
        header = hline.strip("\n").split("\t")
        h_idx = {a:i for i,a in enumerate(header)}
        (c,p,r,a,pval) = (
            h_idx["chrom"],
            h_idx["pos"],
            h_idx["ref"],
            h_idx["alt"],
            h_idx["p2"],
        )

        #write of header
        of.write(hline)

        for line in f:
            data = line.strip("\n").split("\t")
            var = CPRA(
                data[c],
                int(data[p]),
                data[r],
                data[a],
                float(data[pval])
            )
            if var.pval <= 5e-8:
                if not any([region_variant_overlap(a,var) for a in region_chrdict[var.c]]):
                    of.write(line)