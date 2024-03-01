#!/usr/bin/env python3

import argparse
import ld_tools
import gzip
from collections import OrderedDict
from multiprocessing import Pool

CODING_CONSEQUENCES = ["transcript_ablation",
    "splice_donor_variant",
    "stop_gained",
    "splice_acceptor_variant",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant"] 

def read_variant_annot(variant_annot:str, consequence_col:str, other_annot:list[str]):
    """
    Reads the given variant annotation file and returns a dict from chr:pos:ref:alt to a dict with the given columns
    """

    of = gzip.open if variant_annot.endswith(".gz") else open
    with of(variant_annot, 'rt') as f:
        h = f.readline().strip().split('\t')
        reqs = ['chr', 'pos', 'ref', 'alt',consequence_col]
        reqs.extend(other_annot)

        missing = [r for r in reqs if not r in h]
        if len(missing)>0:
            raise Exception(f"Given columns {missing} not in file observed headers {h}")
        h = OrderedDict([ (hd,i) for i,hd in enumerate(h) ])

        res = {}
        for line in f:
            s = line.strip().split('\t')
            chrom = s[h['chr']]
            pos = s[h['pos']]
            ref = s[h['ref']]
            alt = s[h['alt']]

            cpra = ":".join([chrom,pos,ref,alt])
            
            res[cpra] = {hd:s[i] for i,hd in enumerate(h)}
    return res 

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', action='store', type=str, help='')
    parser.add_argument('outfile', action='store', type=str, help='')

    parser.add_argument('variant_annot', action='store', type=str, help='')
    parser.add_argument('-consequence_col', default="most_severe", action='store', type=str, help='column containing consequence information in the variant annotation file')
    parser.add_argument('-other_annot', default="gene_most_severe", action='store', type=str, help='comma separated list of other annotations to add to the output file.')

    parser.add_argument('-chromcol', default="chrom")
    parser.add_argument('-poscol', default="pos")
    parser.add_argument('-refcol', default="ref")
    parser.add_argument('-altcol', default="alt")
    ## 
    
    parser = ld_tools.get_common_LD_API_arg_parser(parser)
    args = parser.parse_args()
    ld_interface = ld_tools.get_ld_api_by_cmdargs(args)
    other_annot = [a.strip() for a in args.other_annot.split(",")]
    annots = read_variant_annot(args.variant_annot, args.consequence_col, other_annot)

    of = gzip.open if args.file.endswith(".gz") else open
    with of(args.file, 'rt', encoding="ISO-8859-1") as infile:
        h = OrderedDict([ (hd,i) for i,hd in enumerate(infile.readline().strip("\n").split("\t")) ])
        reqs = [args.poscol, args.chromcol, args.refcol, args.altcol]
    
        missing = [r for r in reqs if not r in h]
        if len(missing)>0:
            raise Exception(f"Given columns {missing} not in file observed headers {h}")
        
        with open(args.outfile,"wt") as out:
            
            out.write("\t".join( list(h.keys())+ ["functional_variants_in_ld"] ) + "\n") 

            for line in infile:
                s = line.strip("\n").split("\t")
                chrom = s[h[args.chromcol]]
                pos = s[h[args.poscol]]
                ref = s[h[args.refcol]]
                alt = s[h[args.altcol]]
                cpra = ":".join([chrom,pos,ref,alt])

                ld = ld_interface(chrom, pos, ref, alt, r2=args.ld, ld_w=args.ld_w)

                variants_in_ld = []
                
                print("checking", cpra, "for coding consequence.")
                if cpra in annots and annots[cpra][args.consequence_col] in CODING_CONSEQUENCES:
                    print("adding", cpra, "to list of variants in LD.")
                    cons = annots[cpra][args.consequence_col] 
                    other_annot_vals = [annots[cpra][oa] for oa in other_annot]
                    dat = [cpra,cons, "1.0"]
                    dat.extend(other_annot_vals)
                    variants_in_ld.append(dat)
            
                for ldvar in ld:
                    var = ldvar["variation2"]
                    
                    r2= ldvar["r2"]
                    if r2 >= args.ld and var in annots and annots[var][args.consequence_col] in CODING_CONSEQUENCES:
                        cons = annots[var][args.consequence_col] 
                        other_annot_vals = [annots[var][oa] for oa in other_annot]

                        dat = [var,cons, str(r2)]
                        dat.extend(other_annot_vals)
                        variants_in_ld.append(dat)
                funcstring =  ";".join([ ",".join(v) for v in variants_in_ld]) 
                out.write("\t".join(s) + "\t" + funcstring + "\n")        

            

                
