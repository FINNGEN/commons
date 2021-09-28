import hail as hl


def annot_most_severe(mt):

    canon_pc = mt.row.vep.transcript_consequences.filter(
        lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding') &
                  (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                  (x.consequence_terms.contains(mt.row.vep.most_severe_consequence))
    )

    most_severe = mt.row.vep.transcript_consequences.filter(
        lambda x: (x.biotype == "protein_coding") & (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                    (x.consequence_terms.contains(mt.row.vep.most_severe_consequence))
    )

    most_severe_others = mt.row.vep.transcript_consequences.filter(
        lambda x: (x.consequence_terms.contains(mt.row.vep.most_severe_consequence))
    )

    mt = mt.annotate_rows(
        gene_most_severe=
            hl.if_else(hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding') &
                      (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                      (x.consequence_terms.contains(mt.row.vep.most_severe_consequence)),
            mt.row.vep.transcript_consequences),
                canon_pc.first().gene_symbol,
            hl.if_else(hl.any( lambda x: (x.biotype == "protein_coding") & (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                      (x.consequence_terms.contains(mt.row.vep.most_severe_consequence)),
                       mt.row.vep.transcript_consequences),
            most_severe.first().gene_symbol,
            most_severe_others.first().gene_symbol, missing_false=True)
            , missing_false=True),
        genes_most_severe=
            hl.if_else(hl.any(lambda x: (x.canonical == 1) & (x.biotype == 'protein_coding') &
                      (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                      (x.consequence_terms.contains(mt.row.vep.most_severe_consequence)),
            mt.row.vep.transcript_consequences),
                canon_pc.map(lambda x:
                            hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id)),
            hl.if_else(hl.any( lambda x: (x.biotype == "protein_coding") & (x.gene_symbol_source != "Clone_based_ensembl_gene") &
                      (x.consequence_terms.contains(mt.row.vep.most_severe_consequence)),
                       mt.row.vep.transcript_consequences),
                most_severe.map(lambda x:
                            hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id)),
                most_severe_others.map(lambda x:
                            hl.struct(gene_symbol=x.gene_symbol, gene_id=x.gene_id)),
                       missing_false=True)
            , missing_false=True)
    )
    return mt
