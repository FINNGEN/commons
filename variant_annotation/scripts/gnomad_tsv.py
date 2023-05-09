# adapted from
# https://github.com/FINNGEN/commons/blob/master/variant_annotation/scripts/hail_functions.py

import hail as hl


def filter_table(table):
    # MUC5B:
    # return table.filter(
    #    (table.locus == hl.Locus.parse("chr11:1219991", "GRCh38"))
    #    & (table.alleles[1] == "T")
    # )
    # chr21:
    # return table.filter(table.locus.contig == "chr21")
    # passing:
    # return table.filter(hl.len(table.filters) == 0)
    # non-passing:
    # return table.filter(hl.len(table.filters) > 0)
    return table


def annotate_table(table):

    canon_pc = table.vep.transcript_consequences.filter(
        lambda x: (x.canonical == 1)
        & (x.biotype == "protein_coding")
        & (x.gene_symbol_source != "Clone_based_ensembl_gene")
        & (x.consequence_terms.contains(table.vep.most_severe_consequence))
    )

    most_severe = table.vep.transcript_consequences.filter(
        lambda x: (x.biotype == "protein_coding")
        & (x.gene_symbol_source != "Clone_based_ensembl_gene")
        & (x.consequence_terms.contains(table.vep.most_severe_consequence))
    )

    most_severe_others = table.vep.transcript_consequences.filter(
        lambda x: (x.consequence_terms.contains(table.vep.most_severe_consequence))
    )

    return table.annotate(
        chr=table.locus.contig.replace("^chr", ""),
        pos=hl.int32(table.locus.position),
        ref=table.alleles[0],
        alt=table.alleles[1],
        rsids=hl.str(",").join(table.rsid),
        filters=hl.str(",").join(table.filters),
        AN=table.freq[table.freq_index_dict["adj"]].AN,
        AF=table.freq[table.freq_index_dict["adj"]].AF,
        AF_afr=table.freq[table.freq_index_dict["afr-adj"]].AF,
        AF_amr=table.freq[table.freq_index_dict["amr-adj"]].AF,
        AF_asj=table.freq[table.freq_index_dict["asj-adj"]].AF,
        AF_eas=table.freq[table.freq_index_dict["eas-adj"]].AF,
        AF_fin=table.freq[table.freq_index_dict["fin-adj"]].AF,
        AF_mid=table.freq[table.freq_index_dict["mid-adj"]].AF,
        AF_nfe=table.freq[table.freq_index_dict["nfe-adj"]].AF,
        AF_oth=table.freq[table.freq_index_dict["oth-adj"]].AF,
        AF_sas=table.freq[table.freq_index_dict["sas-adj"]].AF,
        consequences=hl.array(
            hl.set(
                table.vep.transcript_consequences.map(
                    lambda x: hl.struct(
                        gene_symbol=x.gene_symbol,
                        gene_id=x.gene_id,
                        consequences=x.consequence_terms,
                    )
                )
            )
        ),
        most_severe=table.vep.most_severe_consequence,
        gene_most_severe=hl.if_else(
            hl.any(
                lambda x: (x.canonical == 1)
                & (x.biotype == "protein_coding")
                & (x.gene_symbol_source != "Clone_based_ensembl_gene")
                & (x.consequence_terms.contains(table.vep.most_severe_consequence)),
                table.vep.transcript_consequences,
            ),
            canon_pc.first().gene_symbol,
            hl.if_else(
                hl.any(
                    lambda x: (x.biotype == "protein_coding")
                    & (x.gene_symbol_source != "Clone_based_ensembl_gene")
                    & (x.consequence_terms.contains(table.vep.most_severe_consequence)),
                    table.vep.transcript_consequences,
                ),
                most_severe.first().gene_symbol,
                most_severe_others.first().gene_symbol,
                missing_false=True,
            ),
            missing_false=True,
        ),
        genes_most_severe=hl.if_else(
            hl.any(
                lambda x: (x.canonical == 1)
                & (x.biotype == "protein_coding")
                & (x.gene_symbol_source != "Clone_based_ensembl_gene")
                & (x.consequence_terms.contains(table.vep.most_severe_consequence)),
                table.vep.transcript_consequences,
            ),
            canon_pc.map(
                lambda x: hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id)
            ),
            hl.if_else(
                hl.any(
                    lambda x: (x.biotype == "protein_coding")
                    & (x.gene_symbol_source != "Clone_based_ensembl_gene")
                    & (x.consequence_terms.contains(table.vep.most_severe_consequence)),
                    table.vep.transcript_consequences,
                ),
                most_severe.map(
                    lambda x: hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id)
                ),
                most_severe_others.map(
                    lambda x: hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id)
                ),
                missing_false=True,
            ),
            missing_false=True,
        ),
    )


def export(table, outfile):
    table.select(
        "rsids",
        "filters",
        "AN",
        "AF",
        "AF_afr",
        "AF_amr",
        "AF_asj",
        "AF_eas",
        "AF_fin",
        "AF_mid",
        "AF_nfe",
        "AF_oth",
        "AF_sas",
        "most_severe",
        "gene_most_severe",
        "consequences",
    ).export(outfile)


table = hl.read_table(
    "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht"
)

table = filter_table(table)
# rerunning VEP, gnomAD uses a different transcript set with a lot more non-coding transcripts and we want more coding genes
table = hl.vep(table, "gs://hail-eu-vep/vep95-GRCh38-loftee-gcloud.json")
table = (
    annotate_table(table).rename({"chr": "#chr"}).key_by("#chr", "pos", "ref", "alt")
)
export(
    table,
    "gs://gnomad3/genomes_3.0/gnomad.genomes.v3.1.2.sites.all.vep95.coding.tsv.bgz",
)

# to start a cluster:
# (only use a high number of workers if sure that this works)
# hailctl dataproc start gnomad --region europe-west1 --zone europe-west1-b --num-workers 10 --max-idle 30m --subnet projects/finngen-refinery-dev/regions/europe-west1/subnetworks/default
#
# to run the job:
# hailctl dataproc submit gnomad --region europe-west1 gnomad_tsv.py --overwrite
#
# to stop the cluster:
# (don't rely on the automatic idle shutdown)
# hailctl dataproc stop gnomad --region europe-west1
