
import hail as hl
from hail_functions import annot_most_severe

import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create annotation file from sites vcf.')
    parser.add_argument('vcf', type=str, help="A single vcf or glob with *.")
    parser.add_argument('outfile', type=str)
    parser.add_argument('--ref_genome', type=str, default="GRCh38")
    parser.add_argument('--vep_conf', type=str, default="gs://hail-eu-vep/vep95-GRCh38-loftee-gcloud.json")
    parser.add_argument('--overwrite', action='store_true')
    parser.add_argument('--no_maf', action='store_true')

    args = parser.parse_args()
    print("Starting annotations")
    print("input: {} and ref {} ".format(args.vcf,args.ref_genome ))

    ## uncomment this if getting this: IllegalArgumentException: requirement failed
    hl._set_flags(no_whole_stage_codegen='1')
    
    ## assume already annot
    if args.vcf.endswith(".mt"):
        mt = hl.read_matrix_table(args.outfile + ".mt")
    else:
        mt = hl.import_vcf(args.vcf,force_bgz=True, reference_genome=args.ref_genome, drop_samples=True)
        mt = mt.filter_rows(mt.alleles[1] != '*')
        mt = hl.vep(mt,args.vep_conf)
        mt.write(args.outfile + ".mt", overwrite=args.overwrite )

    mt= annot_most_severe(mt)
    rows= mt.rows()
    if args.no_maf:
        rows.select('rsid',
                    variant=rows.locus.contig.replace("^chr","").replace("X","23") + ":" + hl.format('%d',rows.locus.position)
                    + ":" + rows.alleles[0] + ":" + rows.alleles[1],
                    gene_most_severe=rows.gene_most_severe,most_severe=rows.vep.most_severe_consequence,
                    genes_most_severe=rows.genes_most_severe).export(args.outfile + "_annot.tsv.bgz")
    else:
        rows.select('rsid',maf=rows.info.MAF[0],
                    variant=rows.locus.contig.replace("^chr","").replace("X","23") + ":" + hl.format('%d',rows.locus.position)
                    + ":" + rows.alleles[0] + ":" + rows.alleles[1],
                    gene_most_severe=rows.gene_most_severe,most_severe=rows.vep.most_severe_consequence,
                    genes_most_severe=rows.genes_most_severe).export(args.outfile + "_annot.tsv.bgz")
