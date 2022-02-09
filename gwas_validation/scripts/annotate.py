import os, sys
tsv = sys.argv[1]
annotation = sys.argv[2]
out = sys.argv[3]
#read in the data
with open(tsv,"r") as tsvfile:
    header = tsvfile.readline().strip("\n").split("\t")
    h_idx = {a:i for i,a in enumerate(header)}
    data = [a.strip("\n").split("\t") for a in tsvfile]
#read annotations to dict
anndict = {}
with open(annotation) as annfile:
    
    
    a_header = annfile.readline().strip("\n").split("\t")
    a_idx = {a:i for i,a in enumerate(a_header)}
    for line in annfile:
        splitline = line.strip("\n").split("\t")
        cpra = (
            splitline[a_idx["chr"]],
            splitline[a_idx["pos"]],
            splitline[a_idx["ref"]],
            splitline[a_idx["alt"]]
        )
        anndict[cpra] = [
            splitline[a_idx["INFO"]],
            splitline[a_idx["AF"]]
        ]

with open(out,"w") as outf:
    outheader = header + ["INFO","AF"]
    outf.write("\t".join(outheader)+"\n")
    for d in data:
        d_idx = (
            d[h_idx["chrom"]],
            d[h_idx["pos"]],
            d[h_idx["ref"]],
            d[h_idx["alt"]]
        )
        ann = ["NA","NA"]
        if d_idx in anndict:
            ann = anndict[d_idx]
        outline = "\t".join(d+ann) + "\n"
        outf.write(outline)