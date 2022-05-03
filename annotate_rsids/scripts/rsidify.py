#!/usr/bin/env python3
import gzip
from typing import NamedTuple
from argparse import ArgumentParser
import sys

class CPRA(NamedTuple):
    c: str
    p: str
    r: str
    a: str

CHR_MAP={
    "1": 1,
    "2": 2,
    "3": 3,
    "4": 4,
    "5": 5,
    "6": 6,
    "7": 7,
    "8": 8,
    "9": 9,
    "10": 10,
    "11": 11,
    "12": 12,
    "13": 13,
    "14": 14,
    "15": 15,
    "16": 16,
    "17": 17,
    "18": 18,
    "19": 19,
    "20": 20,
    "21": 21,
    "22": 22,
    "23": 23,
    "24": 24,
    "25": 25,
    "X": 23,
    "Y": 24,
    "M": 25,
    "MT": 25,
    "x":23,
    "y":24,
    "m":25,
    "mt":25
}

if __name__=="__main__":
    ap = ArgumentParser(prog="RSIDify",description="Annotate a tsv/csv file with rsids from a dbSNP reference. Input file must be sorted by chromosome and position.")
    ap.add_argument("file")
    ap.add_argument("--reference",help="dbSNP reference VCF",required=True)
    ap.add_argument("--out",help="Output file path (not compressed). If output is not specified, output will be written to stdout.")
    ap.add_argument("-c","--chrom",required=True,help="Input file chromosome column")
    ap.add_argument("-p","--pos",required=True,help="Input file position column")
    ap.add_argument("-r","--ref",required=True,help="Input file reference allele column")
    ap.add_argument("-a","--alt",required=True,help="Input file alternate allele column")
    ap.add_argument("-s","--sep",help="File separator in input file",default="\t")
    args = ap.parse_args()

    fopen = gzip.open if args.file.endswith("gz") else open

    if args.out:
        OUTPUT = open(args.out,"wt")
        out_write = OUTPUT.write
    else:
        out_write = sys.stdout.write

    cpra = CPRA(args.chrom,args.pos,args.ref,args.alt)
    SEP=args.sep
    OUT = args.out
    REF = args.reference
    IN = args.file
    last_chrom = -1
    last_pos = -1
    ####
    # Skip comments in reference VCF and get header
    ####
    fp_ref = gzip.open(REF, 'rt')
    ref_has_lines = True
    ref_line = fp_ref.readline()
    while ref_line.startswith("##"):
        ref_line = fp_ref.readline()
    if ref_line.startswith('#'):
        ####
        # Check that all required columns are present
        ####
        req_cols = ["#CHROM","POS","REF","ALT","ID"]
        headercols = ref_line.strip("\r\n").split("\t")
        if not all([a in headercols for a in req_cols]):
            raise ValueError(f"Error on reference file header: Not all required columns present. Found columns {[a for a in req_cols if a in headercols]}, did not find columns {[a for a in req_cols if a not in headercols]}")
    ref_h_idx = {h:i for i,h in enumerate(ref_line.rstrip('\r\n').split('\t'))}
    ref_pos = -1
    ref_chr = -1
    #open input file
    with fopen(IN, 'rt') as f:
        header = f.readline().strip("\r\n")
        h_idx = {h:i for i,h in enumerate(header.split(SEP))}
        line_len = len(h_idx.keys())
        ####
        # Check if header contains all columns 
        ####
        keys_not_in_header =  [a not in h_idx.keys() for a in [cpra.c,cpra.p,cpra.r,cpra.a]]
        if any(keys_not_in_header):
            msg = f"ERROR: columns {keys_not_in_header} not in file header"
            print(msg,file=sys.stderr)
            raise KeyError(msg)

        out_write(f"{header}{SEP}rsid\n")

        #variants that are on same c:p
        ref_vars = []

        for line_idx, line in enumerate(f):
            ####
            # Get variant cpra
            ####
            line = line.strip("\r\n")
            s = line.split(SEP)
            if line_len != len(s):
                raise ValueError(f"Error on input line {line_idx+1}:Input file line had different amount of columns ({len(s)}) than the header ({line_len}). This is very likely an error on the input gile (e.g. fields containing line separator). Line: '{line}'. Header:'{header}'")
            try:
                chrom = int(CHR_MAP[s[h_idx[cpra.c]]])
                pos = int(s[h_idx[cpra.p]])
                ref = s[h_idx[cpra.r]]
                alt = s[h_idx[cpra.a]]
            except ValueError:
                print(f"Error on input line {line_idx+1}:Chromosome or position was not a number or coercable to number. Chromosome:{s[h_idx[cpra.c]]} Position:{s[h_idx[cpra.p]]} ",file= sys.stderr)
                raise
            except KeyError:
                print(f"Error on input line {line_idx+1}:Chromosome was not recognized. Chromosome: {s[h_idx[cpra.c]]}  Allowed values are numerical chromosomes (1-25) with or without 'chr' prefix, and X,Y,M,MT lower or uppercase. ",file=sys.stderr)
                raise
            ####
            # Check that chrom,pos are in order
            ####
            if chrom < last_chrom or (chrom == last_chrom and pos < last_pos):
                if chrom < last_chrom:
                    raise ValueError(f"Error on input line {line_idx+1}: Values not in sort order. Last chromosome value was parsed as {last_chrom}, but current line's chromosome {s[h_idx[cpra.c]]} was parsed as {chrom}, which is smaller than last line's chromosome value.")
                else:
                    raise ValueError(f"Error on input line {line_idx+1}: Values not in sort order. Last chrom:pos was parsed as {last_chrom}:{last_pos}, but current chrom:pos ({s[h_idx[cpra.c]]}:{s[h_idx[cpra.p]]}) was parsed as {chrom}:{pos}, which is lower than last chrom:pos.")
            last_pos = pos
            last_chrom = chrom
            ####
            # If the variant is in same chrompos as the stored variants from dbsnp reference, then keep them. Otherwise, we need to reset them.
            ####
            if any([(int(a[ref_h_idx['#CHROM']]) != chrom) or (int(a[ref_h_idx['POS']]) != pos) for a in ref_vars ]):
                ref_vars = []

            ####
            # Spool reference to current position
            ####
            while ref_has_lines and (ref_chr < chrom or (ref_chr == chrom and ref_pos < pos) ):
                ref_line = fp_ref.readline().rstrip('\r\n').split('\t')
                if ref_line == [""]:
                    ref_has_lines = False
                    break
                try:
                    ref_chr = int(ref_line[ref_h_idx['#CHROM']])
                    ref_pos = int(ref_line[ref_h_idx['POS']])
                except ValueError:
                    print(f"Error: Reference file chrom or pos could not be evaluated. Offending line: {ref_line}",file=sys.stderr)
                    raise
            ####
            # Store all reference lines on same chrom:pos as our variant
            ####
            while ref_has_lines and ref_chr == chrom and ref_pos == pos:
                ref_vars.append(ref_line)
                ref_line = fp_ref.readline().strip().split('\t')
                if ref_line == [""]:
                    ref_has_lines = False
                    break
                try:
                    ref_chr = int(ref_line[ref_h_idx['#CHROM']])
                    ref_pos = int(ref_line[ref_h_idx['POS']])
                except ValueError:
                    print(f"Error: Reference file chrom or pos could not be evaluated. Offending line: {ref_line}",file=sys.stderr)
                    raise

            
            ####
            # Set RSID for this variant
            ####
            #if chrom and pos match and there are items in the ref vars, we try to match.
            # If chrom and pos don't match, we reset ref vars
            rsid = 'NA'
            
            for r in ref_vars:
                if r[ref_h_idx['REF']] == ref and alt in r[ref_h_idx['ALT']].split(','):
                    rsid = r[ref_h_idx['ID']]
                    break
            
            out_write(f"{line}{SEP}{rsid}\n")
    if args.out:
        OUTPUT.close()