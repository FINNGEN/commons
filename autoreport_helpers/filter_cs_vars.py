#!/usr/bin/env python3

import gzip

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tophits', type=str)
    parser.add_argument('--ld_clumped', action='store_true', help="In case autoreporting was ran with ld clumping. get highest r2 functional variant") 
    parser.add_argument('--allow_no_cs', action='store_true')
    args = parser.parse_args()

    ofunc = open
    if args.tophits.endswith(".gz"):
        ofunc = gzip.open
        

    with ofunc(args.tophits,'rt', encoding="ISO-8859-1") as inf:
        hi = { h:i for i,h in enumerate(inf.readline().rstrip("\n").split("\t")) }
        clusters = []
        clust_idx_map = {}

        
        func_var_field = "functional_variants_strict"

        if args.ld_clumped:
            func_var_field = "functional_variants_relaxed"


        with open(args.tophits + ".cs_filter.tsv","wt") as out:
            out.write( "\t".join( list(hi.keys()) + ["best_func_var","best_func_var_pip","best_func_var_cons","best_func_var_gene"]  ) + "\n")
            for l in inf:
                ld = l.rstrip("\n").split("\t")
                cs_vars = { fd.split("|")[0]:fd.split("|") for fd in ld[hi["credible_set_variants"] ].split(";") }
                cs_f_vars = [ fd.split("|") for fd in ld[hi[func_var_field]].split(";") if fd!="" and fd!="NA" ]

                best_f_var=["",0.0,"",""]
                for fvar in cs_f_vars:
                    if not args.ld_clumped and fvar[0] not in cs_vars:
                            msg=f'Missing credible set PIP annotation for functional variant {fvar[0]} in locus {ld[hi["locus_id"]]} for pheno {{ld[hi["locus_id"]]}}'
                            if ld[hi["credible_set_variants"] ]=="NA" and args.allow_no_cs:
                                print(msg)
                                continue
                            else:
                                raise Exception(msg)
                    score=0
                    if args.ld_clumped:
                        score = float(fvar[3])
                    else:
                        csdat = cs_vars[fvar[0]]
                        score = float(csdat[1])

                    if score>best_f_var[1]:
                        best_f_var[0]=fvar[0]
                        best_f_var[1]=score
                        best_f_var[2]=fvar[1]
                        best_f_var[3]=fvar[2]

                best_f_var[1]=str(best_f_var[1])
                out.write("\t".join( ld + list(best_f_var) ) + "\n")

