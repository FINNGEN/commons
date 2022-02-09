#args
# tsv with variants
# summstat to copy vars from
# output name
# direction: the copied variants end up "left" or "right" 
import gzip, os, sys
var_tsv = sys.argv[1]
summstat = sys.argv[2]
out_name = sys.argv[3]
direction = sys.argv[4]
chrom = "~{colnames[0]}"
pos = "~{colnames[1]}"
ref = "~{colnames[2]}"
alt = "~{colnames[3]}"
pval = "~{colnames[4]}"
beta = "~{colnames[5]}"
mlogp = "~{colnames[6]}"
data = {}
with open(var_tsv) as f:
    header = f.readline().strip("\n").split("\t")
    h_idx = {a:i for (i,a) in enumerate(header)}
    for line in f:
        splitline = line.strip("\n").split("\t")
        idx = (
            splitline[h_idx[chrom]],
            splitline[h_idx[pos]],
            splitline[h_idx[ref]],
            splitline[h_idx[alt]]
        )
        values = (
            splitline[h_idx[pval]],
            splitline[h_idx[beta]],
            splitline[h_idx[mlogp]]
        )
        data[idx] = values

with gzip.open(summstat,"rt") as prevfile:
    with open(out_name,"wt") as outf:
        outheader=(
            "chrom",
            "pos",
            "ref",
            "alt",
            "p1",
            "p2",
            "b1",
            "b2",
            "mlogp1",
            "mlogp2",
        )
        outf.write("\t".join(outheader)+"\n")

        header = prevfile.readline().strip("\n").split("\t")
        h_idx = {a:i for (i,a) in enumerate(header)}
        for line in prevfile:
            splitline = line.strip("\n").split("\t")
            idx = (
                splitline[h_idx[chrom]],
                splitline[h_idx[pos]],
                splitline[h_idx[ref]],
                splitline[h_idx[alt]]
            )

            if idx in data:
                c_stats = data[idx]
                if direction == "left":
                    outline = [
                        idx[0],
                        idx[1],
                        idx[2],
                        idx[3],
                        c_stats[0],
                        splitline[h_idx[pval]],
                        c_stats[1],
                        splitline[h_idx[beta]],
                        c_stats[2],
                        splitline[h_idx[mlogp]]
                    ]
                else:
                    outline = [
                        idx[0],
                        idx[1],
                        idx[2],
                        idx[3],
                        splitline[h_idx[pval]],
                        c_stats[0],
                        splitline[h_idx[beta]],
                        c_stats[1],
                        splitline[h_idx[mlogp]],
                        c_stats[2]
                    ]
                outf.write("\t".join(outline)+"\n")