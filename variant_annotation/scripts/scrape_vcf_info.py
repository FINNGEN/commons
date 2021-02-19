#! /usr/bin/env python3
import argparse
import gzip
import re

def read_info_fields(file):
    fields = []
    for line in file:
        if line.startswith("#CHROM"):
            break
        if line.startswith("##INFO=<"):
            line = line.rstrip("\n")
            dat = { elem[0]:elem[1] for elem in map(lambda x: x.split("="), re.sub("^##INFO=<|>$","",line).split(",")) }
            fields.append(dat)

    return fields

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Export info fields as TSV file')
    parser.add_argument('vcf', type=str)
    parser.add_argument('outfile', type=str)

    args = parser.parse_args()

    with open(args.outfile, 'w') as out:
        with gzip.open(args.vcf,mode='rt') as infile:
            fields = read_info_fields(infile)
            fields.sort( key=lambda x: x["ID"]  )

            items = map(lambda x: x["ID"], fields )

            out.write( "variant\tchr\tpos\t" + "\t".join(items) + "\n" )

            for line in infile:
                dat = line.split("\t")[0:8]
                var = dat[0].replace('chr', '').replace('X', '23').replace('Y', '24').replace('MT', '25').replace('M', '25') + ":" + dat[1] + ":" + dat[3] + ":" + dat[4]
                info = { d[0]:(d[1] if len(d)>1 else "") for d in list(map(lambda x: x.split("="), dat[7].split(";")) ) }

                def getdat( elem, info):
                    if elem["Type"]=="Flag":
                        return "1" if elem["ID"] in info else "0"
                    else:
                        return info[elem["ID"]] if elem["ID"] in info else "NA"

                out.write( var + "\t" + dat[0] + "\t" + dat[1] + "\t" + "\t".join([  getdat(elem, info) for elem in fields  ] ) + "\n")
