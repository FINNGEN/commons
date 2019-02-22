from hail import *
from hail.utils import hadoop_read
import argparse


def annotate_variants(vds):
    if not vds.was_split():
        raise Exception("Split multi-allelics before running!")

    return vds.annotate_variants_expr(
        'va.LOF = if ( !isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding" ).canonical) && ( let cons = va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding"  && ( isMissing(c.variant_allele) || c.variant_allele==v.alt() ) ) in !isMissing(cons.lof) && cons.lof=="HC" && (isMissing(cons.lof_flags) || cons.lof_flags=="" ) )   || '
        '( isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding" && ( isMissing(c.variant_allele) || c.variant_allele==v.alt() ) ).canonical) '
        '&& ( let cons=va.vep.transcript_consequences.find(c =>c.consequence_terms.toSet.contains(va.vep.most_severe_consequence) ) in !isMissing(cons.lof) && cons.lof == "HC" && (isMissing(cons.lof_flags) || cons.lof_flags=="" ) ) ) )    '
        'true else false, '
        'va.missense = if( (va.vep.transcript_consequences.find(c => c.canonical == 1 && ( isMissing(c.variant_allele) || c.variant_allele==v.alt() ) && c.biotype=="protein_coding").consequence_terms.toSet.contains("missense_variant") ) || (isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding").canonical) && va.vep.most_severe_consequence == "missense_variant") || va.vep.transcript_consequences.find(c => c.canonical == 1 && ( isMissing(c.variant_allele) || c.variant_allele==v.alt() ) && c.biotype=="protein_coding").consequence_terms.toSet.contains("inframe_deletion") || (isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding").canonical) && va.vep.most_severe_consequence == "inframe_deletion") ) true else false, '
        'va.synonymous=va.vep.most_severe_consequence=="synonymous_variant", va.gene = if (isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding").canonical)) va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).gene_symbol else va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding").gene_symbol,  va.genenames = if (isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding").canonical)) va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).gene_symbol else str(va.vep.transcript_consequences.filter(c => c.canonical == 1 && c.biotype=="protein_coding").map(c => c.gene_symbol)), va.LOF_flag = if( ! isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding" && c.lof=="HC") ) ) va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding" && c.lof=="HC").lof_flags  else if (  isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding").canonical) && !isMissing( va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence) && c.lof == "HC" ) ) ) va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence) && c.lof == "HC" ).lof_flags else "", '
        'va.gene = if (isMissing(va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding").canonical))'
        'va.vep.transcript_consequences.find(c => c.consequence_terms.toSet.contains(va.vep.most_severe_consequence)).gene_symbol else'
        'va.vep.transcript_consequences.find(c => c.canonical == 1 && c.biotype=="protein_coding").gene_symbol')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Annotate variants to separate file")

    parser.add_argument('inputfiles', action='store', type=str,
                        help='vcf or glob to vcfs in Google Bucket')

    parser.add_argument('output_prefix', action='store', type=str,
                        help='Prefix to google bucket')

    parser.add_argument('-vep_conf', action='store', type=str,
                        help='vep config file', default="/vep/vep-gcloud.properties")

    hc = HailContext()

    args = parser.parse_args()
    output_prefix = args.output_prefix
    files = [  f.rstrip("\n") for f in hadoop_read(args.inputfiles) ]

    sites = hc.import_vcf(files, drop_samples=True, force_bgz=True)
    sites.write(output_prefix + "_sites.vds")
    sites = sites.split_multi()
    ## assumes you have ran the cluster up with cloudtools

    sites = sites.vep(config="/vep/vep-gcloud-grch38.properties")

    sites.write(output_prefix +  "_annot.vds" )
    sites.persist()
    sites.export_variants( output_prefix, "variant=v, chr=v.contig, pos=v.start, ref=v.ref, alt=v.alt, hc_lof=va.LOF, gene=va.gene, most_severe=va.vep.most_severe_consequence" )
