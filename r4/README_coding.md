# Coding variant munging

Create annotation and association result tables of coding variants

## Required files

```
gs://fg-datateam-analysisteam-share/gnomad/fin_enriched_exomes_select_columns.txt.bgz
gs://r4_data_west1/annotations/R4_annotated_variants_v1.gz
gs://r4_data_west1/annotations/R4_annotated_variants_v1.gz.tbi
gs://finngen-production-library-red/finngen_R4/finngen_R4_phenotype_2.0/finngen_R4_phenotype_documentation/finngen_R4_endpoint_definition.tsv
matrix.tsv.gz and matrix.tsv.gz from pheweb
```

Use e.g. `gcr.io/finngen-refinery-dev/bioinformatics:0.4`

## Create table of coding variant annotations in gnomAD and FinnGen

This takes about 2 hours like this

### Filter coding variants from gnomAD

```
gunzip -c fin_enriched_exomes_select_columns.txt.bgz | awk '
BEGIN{FS=OFS="\t"}
NR==1{for (i=1;i<=NF;i++) {sub("^locus$", "grch37_locus", $i); a[$i]=i} print "#variant",$0}
NR>1{if ($a["consequence"] == "pLoF" ||
         $a["consequence"] == "LC" ||
         $a["consequence"] == "start_lost" ||
         $a["consequence"] == "stop_lost" ||
         $a["consequence"] == "inframe_indel" ||
         $a["consequence"] == "missense_variant") {
        sub("chr", "", $a["chrom"]);
        sub("X", "23", $a["chrom"]);
        if ($a["chrom"] ~ /^[0-9]+$/) {
            print $a["chrom"]":"$a["pos"]":"$a["ref"]":"$a["alt"],$0
        }
    }
}' | bgzip -@ 4 > gnomad_func.txt.gz
```

### Get FinnGen annotations for coding variants

```
gunzip -c gnomad_func.txt.gz | awk '
BEGIN{FS=OFS="\t"}
NR==1{for (i=1;i<NF;i++) {a[$i]=i}}
NR>1{print $a["chrom"],$a["pos"],$a["pos"]}
' | uniq > gnomad_func.tabix

tabix -h -R gnomad_func.tabix R4_annotated_variants_v1.gz | bgzip -@ 4 > R4_annotated_variants_v1_gnomad_func.txt.gz
```

### Join FinnGen annotations to gnomAD and filter on INFO

```
# bioinformatics crane kick
n_join=$((`gunzip -c gnomad_func.txt.gz | head -1 | awk 'BEGIN{FS="\t"} {print NF}'` + `gunzip -c R4_annotated_variants_v1_gnomad_func.txt.gz | head -1 | awk 'BEGIN{FS="\t"} {print NF-1}'`))
join -a 1 -t $'\t' -1 2 -2 1 \
<(gunzip -c gnomad_func.txt.gz | nl -nln | sort -b -k 2,2) \
<(gunzip -c R4_annotated_variants_v1_gnomad_func.txt.gz | sort -b -k 1,1) | sort -g -k 2,2 | cut -f1,3- | \
awk -v nf=$n_join 'BEGIN{FS=OFS="\t"} {printf $0; if(NF<nf) {for(i=NF+1;i<=nf;i++) printf "\tNA"} printf "\n"}' | \
awk '
BEGIN{FS=OFS="\t"}
NR==1 {for (i=1;i<=NF;i++) {a[$i]=i} print $0}
NR>1 {if ($a["INFO"] != "NA" && $a["INFO"] > 0.6) print $0}' | \
uniq | bgzip -@ 4 > gnomad_func_info_0.6_r4_annotation.txt.gz
```

`gnomad_func_info_0.6_r4_annotation.txt.gz` now has coding variant INFO > 0.6 annotations from gnomAD exomes and FinnGen

## Add variant associations

### Get association results for the coding variants

This takes an hour

```
gunzip -c gnomad_func_info_0.6_r4_annotation.txt.gz | tail -n+2 | \
awk 'BEGIN{FS=OFS="\t"}{split($2,s,":"); print(s[1],s[2],s[2])}' | \
uniq > gnomad_func_info_0.6_r4_annotation.tabix

tabix -h -R <(awk '{sub("chr", "", $1); print $0}' gnomad_func_info_0.6_r4_annotation.tabix) \
matrix.tsv.gz | bgzip -@ 4 > matrix.func.txt.gz

# get a table format filtering for p-value
python3 scripts/significant_from_matrix.py matrix.func.txt.gz 1e-4 | bgzip -@4 > func_1e-4.txt.gz
```

`func_1e-4.txt.gz` now has p < 1e-4 associations for the coding variants

### Join association results with annotations

The files are small at this point.. too much awk

```
# join phenotype names
join -a 1 -t$'\t' -1 3 -2 1 \
<(gunzip -c func_1e-4.txt.gz | nl -nln | sort -k 3,3) \
<(cut -f4,5 finngen_R4_endpoint_definition.tsv | sed 's/NAME/pheno/' | sed 's/LONGNAME/phenoname/' | sort -k1,1) \
| sort -g -k2,2 | awk '
BEGIN{FS=OFS="\t"}
NR==1{n=NF; print $0}
NR>1{if(NF<n) {printf $0; for(i=NF;i<n;i++) printf "\t"; printf "\n"} else print $0}
' | bgzip > func_1e-4_phenoname.txt.gz
```

```
# join annotations to associations and select columns
join -t$'\t' -1 3 -2 1 \
<(gunzip -c func_1e-4_phenoname.txt.gz | sort -b -k 3,3) \
<(gunzip -c gnomad_func_info_0.6_r4_annotation.txt.gz | sort -b -k 1,1) \
| sort -g -k3,3 | awk '
BEGIN{FS=OFS="\t"}
NR==1 {
    j=1;
    for (i=1;i<=NF;i++) {
        sub("^#variant$", "variant", $i);
        sub("^consequence$", "variant_category", $i);
        a[$i] = i;
        if ($i ~ "^AF_" || $i ~ "^INFO_" || $i ~ "^CHIP" || $i ~ "^HWE") {
            extra_cols[j] = i;
            j++;
        }
    }
}
{
    sub("NA", -1, $a["enrichment_nfsee"]); sub("NaN", -2, $a["enrichment_nfsee"]); sub("Infinity", 1e6, $a["enrichment_nfsee"]);
    printf($a["pheno"]"\t"$a["phenoname"]"\t"$a["variant"]"\t"$a["rsid"]"\t"$a["chrom"]"\t"$a["pos"]"\t"$a["ref"]"\t"$a["alt"]"\t"\
           $a["pval"]"\t"$a["beta"]"\t"$a["sebeta"]"\t"$a["maf"]"\t"$a["maf_cases"]"\t"$a["maf_controls"]"\t"$a["INFO"]"\t"\
           $a["AF"]"\t"$a["AC"]"\t"$a["AC_Hom"]"\t"$a["AC_Het"]"\t"$a["AN"]"\t"\
           $a["variant_category"]"\t"$a["most_severe"]"\t"$a["gene_most_severe"]"\t"$a["grch37_locus"]"\t"\
           $a["fin.AF"]"\t"$a["fin.AC"]"\t"$a["fin.AN"]"\t"$a["fin.homozygote_count"]"\t"$a["nfsee.AF"]"\t"$a["nfsee.AC"]"\t"$a["nfsee.AN"]"\t"$a["nfsee.homozygote_count"]"\t"\
           $a["enrichment_nfsee"]"\t"$a["fet_nfsee.p_value"]"\t"$a["fet_nfsee.odds_ratio"])
    for (i=1;i<=length(extra_cols);i++) {
        printf("\t"$extra_cols[i])
    }
    printf("\n")
}' | iconv -f ISO-8859-1 -t UTF-8 | bgzip > r4_coding_variant_associations_1e-4.txt.gz
```

`r4_coding_variant_associations_1e-4.txt.gz` now has annotations and associations

### Filter Finnish-enriched variants

```
gunzip -c r4_coding_variant_associations_1e-4.txt.gz | \
awk '
BEGIN{FS=OFS="\t"}
NR==1 {for(i=1;i<=NF;i++) a[$i]=i; print $0}
NR>1 {if ($a["fet_nfsee.p_value"] < 0.01 &&
(($a["fin.AF"] > 0.001 && ($a["enrichment_nfsee"] > 7.5 || $a["nfsee.AC"] == 0)) ||
($a["fin.AF"] > 0.0005 && $a["fin.AF"] <= 0.001 && ($a["enrichment_nfsee"] > 10 || $a["nfsee.AC"] == 0))))
print $0
}' | bgzip > r4_coding_variant_associations_1e-4_fin_enriched.txt.gz
```

## List functional variants regardless of association

```
gunzip -c gnomad_func_info_0.6_r4_annotation.txt.gz | \
awk 'BEGIN{FS=OFS="\t"}
NR==1 {
    j=1;
    for (i=1;i<=NF;i++) {
        sub("^#variant$", "variant", $i);
        sub("^consequence$", "variant_category", $i);
	a[$i] = i;
        if ($i ~ "^AF_" || $i ~ "^INFO_" || $i ~ "^CHIP" || $i ~ "^HWE") {
            extra_cols[j] = i;
            j++;
        }
    }
}
NR>=1 {
    sub("NA", -1, $a["enrichment_nfsee"]); sub("NaN", -2, $a["enrichment_nfsee"]); sub("Infinity", 1e6, $a["enrichment_nfsee"]);
    printf($a["variant"]"\t"$a["rsid"]"\t"$a["chrom"]"\t"$a["pos"]"\t"$a["ref"]"\t"$a["alt"]"\t"\
           $a["INFO"]"\t"$a["AF"]"\t"$a["AC"]"\t"$a["AC_Hom"]"\t"$a["AC_Het"]"\t"$a["AN"]"\t"\
           $a["variant_category"]"\t"$a["most_severe"]"\t"$a["gene_most_severe"]"\t"$a["grch37_locus"]"\t" \
           $a["fin.AF"]"\t"$a["fin.AC"]"\t"$a["fin.AN"]"\t"$a["fin.homozygote_count"]"\t"$a["nfsee.AF"]"\t"$a["nfsee.AC"]"\t"$a["nfsee.AN"]"\t"$a["nfsee.homozygote_count"]"\t" \
           $a["enrichment_nfsee"]"\t"$a["fet_nfsee.p_value"]"\t"$a["fet_nfsee.odds_ratio"])
    for (i=1;i<=length(extra_cols);i++) {
        printf("\t"$extra_cols[i])
    }
    printf("\n")
}' | bgzip -@ 4 > r4_coding_variants.txt.gz
```

## Add recessive and dominant results and PIPs

Get `full.formatted.txt.gz` files from recessive and dominant GWAS pipeline

```
gsutil cp gs://fg-cromwell/saige_tests/92ceb874-8d7e-4697-b545-afdd1738f1ac/call-combine/full.formatted.txt.gz r4_coding_recessive_full.formatted.txt.gz
gsutil cp gs://fg-cromwell/saige_tests/e7a42906-3dd5-4fba-bd4f-16597f170cae/call-combine/full.formatted.txt.gz r4_coding_dominant_full.formatted.txt.gz
```

Get SuSiE `.snp` files from finemapping pipeline

```
mkdir -p finemapping/susie/snp
gsutil -mq cp gs://r4_data_west1/finemap/*.susie.snp finemapping/susie/snp/
```

Add recessive/dominant p-values and PIPs:

```
python3 scripts/add_rec_dom_pip.py r4_coding_variant_associations_1e-4.txt.gz r4_coding_recessive_full.formatted.txt.gz r4_coding_dominant_full.formatted.txt.gz "finemapping/susie/snp/*.susie.snp"
```

Writes `coding_web.txt` for the browser
