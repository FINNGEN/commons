task join_annot {

	Array[File] files
	File? external_annot
	Int local_disk=500
  	String docker

	runtime {
		docker: "${docker}"
		memory:"8G"
		disks: "local-disk ${local_disk} HDD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 0
		noAddress: true
	}

	command <<<
		set -eux
		#save headercheck to file
		cat > headercheck.awk <<'CODE'
		{
		for(i = 1; i <= NF; i++) {
			h[$i] = i
			if($i==value){print i}
		}
		exists=value in h
		if(!exists)
		{
			system("echo ERROR header did not contain column "value"  >&2")
			exit 1
		}
		exit 0
		}
		CODE

		#if there are no external annots, then just simply concat files
		if [[ -z "${external_annot}" ]]; then
			cat <(head -n 1 ${files[0]}) <(awk 'FNR>1' ${sep=" " files}) | bgzip > annotated_variants.gz
		else
			#else sort VEP and annotation files in similar order, then join
			
			#headercheck: check that the column is in header, else abort with error 1
			vep_varcol=$(zcat ${external_annot} |head -n1|awk -v value="variant" -F "\t" -f headercheck.awk)
			vep_conscol=$(zcat ${external_annot}|head -n1|awk -v value="most_severe" -F "\t" -f headercheck.awk )
			vep_genecol=$(zcat ${external_annot}|head -n1|awk -v value="gene_most_severe" -F "\t" -f headercheck.awk)
			cat <(zcat ${external_annot}|head -n1|cut -f $vep_varcol,$vep_conscol,$vep_genecol) <(zcat ${external_annot}|cut -f $vep_varcol,$vep_conscol,$vep_genecol|sort -V -k 1,1)|bgzip > sorted_vep.gz
			#sort annotation
			cat <(head -n 1 ${files[0]}) <(awk 'FNR>1' ${sep=" " files}) | bgzip > annotated_variants_1.gz
			cat <(zcat annotated_variants_1.gz|head -n1)  <(zcat annotated_variants_1.gz|tail -n+2|sort -V -k 1,1 -T ./) |bgzip  > sorted_annotated.gz
			#join
			join -t $'\t' -1 1 -2 1 -e "NA" --header -a 1 <(zcat sorted_annotated.gz) <(zcat sorted_vep.gz) |bgzip  > annotated_variants.gz
		fi

		tabix -b 3 -e 3 -s 2 annotated_variants.gz

	>>>



	output {
		File out="annotated_variants.gz"
		File index="annotated_variants.gz.tbi"
	}


}


task small {
	File annotation
	String wanted_columns
	Int local_disk = ceil(size(annotation, "G")) + 10
	String docker
	runtime {
		docker: "${docker}"
		memory:"8G"
		disks: "local-disk ${local_disk} HDD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 0
		noAddress: true
	}

	command <<<
		set -eux
		cat > cols.awk <<'AWK'
		BEGIN {
		FS=OFS="\t"
		len_head=split("${wanted_columns}", wanted, " ")
		}
		NR==1 {
		for(i=1;i<=NF;i++) h[$i]=i
		for(i in wanted) if (!(wanted[i] in h)) {print wanted[i]" expected in header">>"/dev/stderr"; exit 1}
		for(i=1;i<=len_head;i++) { printf $h[wanted[i]]"\t" } printf "index\n"
		}
		NR >1 {
		for(i=1;i<=len_head;i++) { printf $h[wanted[i]]"\t" } printf NR-2"\n"
		}
		AWK

		zcat ${annotation}|awk -f cols.awk |bgzip  > annotated_variants.small.gz
		tabix -b 3 -e 3 -s 2 annotated_variants.small.gz
	>>>

	output {
		File small_out="annotated_variants.small.gz"
		File small_index="annotated_variants.small.gz.tbi"
	}
}

task extract {

	File vcf
	Int local_disk = ceil(size(vcf, "G")) + 100
	String docker
	String outfile=basename(vcf) + ".annot.gz"

	runtime {
		docker: "${docker}"
		memory:"2G"
		disks: "local-disk ${local_disk} HDD"
		cpu:"1"
		zones: "europe-west1-b europe-west1-c europe-west1-d"
		preemptible: 0
		noAddress: true
	}

	command <<<

        python3 <<EOF

        import gzip
        import re

        def read_info_fields(file):
            fields = []
            for line in file:
                if line.startswith("#CHROM"):
                    break
                if line.startswith("##INFO=<"):
                    line = line.rstrip("\n").replace(', ', '-') # next line splits on commas also within Description so avoid breaking
                    dat = { elem[0]:elem[1] for elem in map(lambda x: x.split("="), re.sub("^##INFO=<|>$","",line).split(",")) }
                    fields.append(dat)

            return fields

        with open('${outfile}', 'w') as out:
            with gzip.open('${vcf}',mode='rt') as infile:
                fields = read_info_fields(infile)
                fields.sort( key=lambda x: x["ID"]  )

                items = map(lambda x: x["ID"], fields )

                out.write( "#variant\tchr\tpos\tref\talt\t" + "\t".join(items) + "\n" )
                chromsanitize = lambda c: c.replace("chr","").replace("X","23").replace("Y","24").replace("MT","25").replace("M","25")
                for line in infile:
                    dat = line.split("\t")[0:8]
                    var = chromsanitize(dat[0]) + ":" + dat[1] + ":" + dat[3] + ":" + dat[4]
                    info = { d[0]:(d[1] if len(d)>1 else "") for d in list(map(lambda x: x.split("="), dat[7].split(";")) ) }

                    def getdat( elem, info):
                        if elem["Type"]=="Flag":
                            return "1" if elem["ID"] in info else "0"
                        else:
                            return info[elem["ID"]] if elem["ID"] in info else "NA"

                    out.write( var + "\t" + chromsanitize(dat[0]) + "\t" + dat[1] + "\t" + dat[3] + "\t" + dat[4] + "\t" + "\t".join([  getdat(elem, info) for elem in fields  ] ) + "\n")

        EOF

	>>>

	output {
		File out="${outfile}"
	}

}

# Add rsids
task add_rsids {

	File annotation_file

	File ref_file
	String docker

    Int local_disk = (ceil(size(annotation_file, "G")) + ceil(size(ref_file, "G"))) * 2

    command <<<

        set -euxo pipefail

        echo "`date` Adding rsids"

        python3 <<EOF | bgzip > annotation.gz

        import gzip
        fp_ref = gzip.open('${ref_file}', 'rt')
        ref_has_lines = True
        ref_chr = 1
        ref_pos = 0
        ref_line = fp_ref.readline()
        while ref_line.startswith("##"):
            ref_line = fp_ref.readline()
        if ref_line.startswith('#'):
            assert ref_line.rstrip('\r\n').split('\t') == '#CHROM POS ID REF ALT QUAL FILTER INFO'.split(), repr(ref_line)
        ref_h_idx = {h:i for i,h in enumerate(ref_line.rstrip('\r\n').split('\t'))}

        with gzip.open('${annotation_file}', 'rt') as f:
            header = f.readline().strip()
            h_idx = {h:i for i,h in enumerate(header.split('\t'))}
            print(header + '\trsid')
            #variants that are on same c:p
            ref_vars = []
            for line in f:
                #get variant values
                line = line.strip()
                s = line.split('\t')
                chrom = int(s[h_idx['chr']])
                pos = int(s[h_idx['pos']])
                ref = s[h_idx['ref']]
                alt = s[h_idx['alt']]
                #check if ref vars c:p has moved forward
                if ref_vars:
                    if int(ref_vars[0][ref_h_idx['#CHROM']]) != chrom or int(ref_vars[0][ref_h_idx['POS']]) != pos:
                        ref_vars = []
                while ref_has_lines and ref_chr < chrom or (ref_chr == chrom and ref_pos < pos):
                    ref_line = fp_ref.readline().rstrip('\r\n').split('\t')
                    try:
                        ref_chr = int(ref_line[ref_h_idx['#CHROM']])
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False
                while ref_has_lines and int(ref_chr) == chrom and ref_pos == pos:
                    ref_vars.append(ref_line)
                    ref_line = fp_ref.readline().strip().split('\t')
                    try:
                        ref_chr = int(ref_line[ref_h_idx['#CHROM']])
                        ref_pos = int(ref_line[ref_h_idx['POS']])
                    except ValueError:
                        ref_has_lines = False

                rsid = 'NA'
                #if chrom and pos match and there are items in the ref vars, we try to match.
                # If chrom and pos don't match, we reset ref vars
                if ref_vars:
                    if (chrom == int(ref_vars[0][ref_h_idx["#CHROM"]]) ) and (pos == int(ref_vars[0][ref_h_idx["POS"]]) ):
                        for r in ref_vars:
                            if r[ref_h_idx['REF']] == ref and alt in r[ref_h_idx['ALT']].split(','):
                                rsid = r[ref_h_idx['ID']]
                                break
                    else:
                        #chrom and pos have advanced further than ref vars, clear ref vars
                        ref_vars = []
                
                print(line + '\t' + rsid)

        EOF

        echo "`date` done"

    >>>

    output {
        File rsid_out = "annotation.gz"
    }

    runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "4 GB"
        disks: "local-disk ${local_disk} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }

}

task add_gnomad {
	File annotation_file

	File gnomad_genomes
	File gnomad_exomes
    Array[String] exome_cols
    Array[String] genome_cols

    Int local_disk = (ceil(size(annotation_file, "G")) + ceil(size(gnomad_genomes, "G")) + ceil(size(gnomad_exomes, "G"))) * 2

	String docker

	command <<<
        set -euxo pipefail

        echo "`date` Adding gnomad annotations"
        
        python3 <<EOF | bgzip > annotation_gnomad.gz
        
        import gzip
        class Reference:
            def __init__(self, fp):
                self.fp = fp
                self.chrom = 1
                self.pos = 0
                self.has_lines=True
                self.line = self.fp.readline()
                self.header = self.line.rstrip('\r\n').split('\t')
                self.h_idx = {h:i for i,h in enumerate(self.line.rstrip('\r\n').split('\t'))}
                self.has_lines = True

        fp_genome = gzip.open('${gnomad_genomes}', 'rt')
        fp_exome = gzip.open('${gnomad_exomes}', 'rt')

        gen_ref = Reference(fp_genome)
        exo_ref = Reference(fp_exome)

        exome_cols = ["${sep='\",\"' exome_cols}"]
        genome_cols = ["${sep='\",\"' genome_cols}"]

        exome_cols_headernames = [f"EXOME_{val}" for val in exome_cols]
        genome_cols_headernames = [f"GENOME_{val}" for val in genome_cols]

        exome_rename_dict = {a:b for a,b in zip(exome_cols,exome_cols_headernames)}
        genome_rename_dict = {a:b for a,b in zip(genome_cols,genome_cols_headernames)}

        with gzip.open('${annotation_file}', 'rt') as f:
            header = f.readline().strip()
            h_idx = {h:i for i,h in enumerate(header.split('\t'))}
            ex_head = '\t'.join(exome_cols_headernames)
            gen_head = '\t'.join(genome_cols_headernames)
            out_header = f"{header}\tb37_coord\t{ex_head}\t{gen_head}"
            print(out_header)

            #create vars for each of the gen/exomes
            gen_vars = []
            exo_vars=[]
            for line in f:
                oline = line.rstrip("\r\n")
                s = oline.split('\t')
                chrom = int(s[h_idx['chr']])
                pos = int(s[h_idx['pos']])
                ref = s[h_idx['ref']]
                alt = s[h_idx['alt']]
                #check if gen and exo vars have different c:p than variant. If yes, empty them.
                if gen_vars:
                    if int(gen_vars[0][gen_ref.h_idx['#CHROM']]) != chrom or int(gen_vars[0][gen_ref.h_idx['POS']]) != pos:
                        gen_vars = []

                if exo_vars:
                    if int(exo_vars[0][exo_ref.h_idx['#CHROM']]) != chrom or int(exo_vars[0][exo_ref.h_idx['POS']]) != pos:
                        exo_vars = []

                #genome
                while gen_ref.has_lines and gen_ref.chrom < chrom or ( gen_ref.chrom == chrom and gen_ref.pos < pos):
                    gen_ref.line = gen_ref.fp.readline().rstrip('\r\n').split('\t')
                    try:
                        gen_ref.pos = int(gen_ref.line[gen_ref.h_idx["POS"]])
                        gen_ref.chrom = int(gen_ref.line[gen_ref.h_idx["#CHROM"]])
                    except ValueError:
                        gen_ref.has_lines = False
                    
                #exome
                while exo_ref.has_lines and exo_ref.chrom < chrom or (exo_ref.chrom == chrom and exo_ref.pos < pos):
                    exo_ref.line = exo_ref.fp.readline().rstrip('\r\n').split('\t')
                    try:
                        exo_ref.pos = int(exo_ref.line[exo_ref.h_idx["POS"]])
                        exo_ref.chrom = int(exo_ref.line[exo_ref.h_idx["#CHROM"]])
                    except ValueError:
                        exo_ref.has_lines = False
                
                #drain all lines with same chrom & pos to gen_vars
                while gen_ref.has_lines and ( gen_ref.chrom == chrom and gen_ref.pos == pos ):
                    gen_vars.append(gen_ref.line)
                    gen_ref.line = gen_ref.fp.readline().rstrip("\r\n").split('\t')
                    try:
                        gen_ref.chrom = int(gen_ref.line[gen_ref.h_idx['#CHROM']])
                        gen_ref.pos = int(gen_ref.line[gen_ref.h_idx['POS']])
                    except ValueError:
                        gen_ref.has_lines = False
                
                while exo_ref.has_lines and ( exo_ref.chrom == chrom and exo_ref.pos == pos ):
                    exo_vars.append(exo_ref.line)
                    exo_ref.line = exo_ref.fp.readline().rstrip("\r\n").split('\t')
                    try:
                        exo_ref.chrom = int(exo_ref.line[exo_ref.h_idx['#CHROM']])
                        exo_ref.pos = int(exo_ref.line[exo_ref.h_idx['POS']])
                    except ValueError:
                        exo_ref.has_lines = False

                #if the chrom and pos are the same, we try to match the variants. Else, we reset the respective {gen,exo}_vars
                gen_values = ["NA"]*len(genome_cols)
                if gen_vars:
                    if chrom == int(gen_vars[0][gen_ref.h_idx['#CHROM']]) and pos == int(gen_vars[0][gen_ref.h_idx['POS']]):
                        for genvar in gen_vars:
                            match = (ref == genvar[gen_ref.h_idx['REF']]) and (alt == genvar[gen_ref.h_idx['ALT']])
                            if match:
                                gen_values = [genvar[gen_ref.h_idx[col]] for col in genome_cols]
                                break
                    else: #reset
                        gen_vars = []

                exo_values = ["NA"]*len(exome_cols)
                if exo_vars:
                    if chrom == int(exo_vars[0][exo_ref.h_idx['#CHROM']]) and pos == int(exo_vars[0][exo_ref.h_idx['POS']]):
                        for exvar in exo_vars:
                            match = (ref == exvar[exo_ref.h_idx['REF']]) and (alt == exvar[exo_ref.h_idx['ALT']])
                            if match:
                                exo_values = [exvar[exo_ref.h_idx[col]] for col in exome_cols]
                                break
                    else: #reset
                        exo_vars=[]
                
                exvals = '\t'.join(exo_values)
                gevals = '\t'.join(gen_values)
                #b37 coord
                b37_coord = "NA"
                for exo in exo_vars:
                    match = (chrom == int(exo[exo_ref.h_idx['#CHROM']]) ) and (pos == int(exo[exo_ref.h_idx['POS']]) ) and (ref == exo[exo_ref.h_idx['REF']]) and (alt == exo[exo_ref.h_idx['ALT']])
                    if match:
                        #CHROM_37	POS_37	REF_37	ALT_37
                        b37_coord = f"{exo[exo_ref.h_idx['CHROM_37']]}:{exo[exo_ref.h_idx['POS_37']]}:{exo[exo_ref.h_idx['REF_37']]}:{exo[exo_ref.h_idx['ALT_37']]}"
                        break
                for geno in gen_vars:
                    if b37_coord != "NA":
                        break
                    match = (chrom == int(geno[gen_ref.h_idx['#CHROM']]) ) and (pos == int(geno[gen_ref.h_idx['POS']]) ) and (ref == geno[gen_ref.h_idx['REF']]) and (alt == geno[gen_ref.h_idx['ALT']])
                    if match:
                        b37_coord = f"{geno[gen_ref.h_idx['CHROM_37']]}:{geno[gen_ref.h_idx['POS_37']]}:{geno[gen_ref.h_idx['REF_37']]}:{geno[gen_ref.h_idx['ALT_37']]}"
                        break
                outline = f"{oline}\t{b37_coord}\t{exvals}\t{gevals}"
                print(outline)
        gen_ref.fp.close()
        exo_ref.fp.close()
        EOF

        echo "`date` tabixing"
        tabix -s 2 -b 3 -e 3 annotation_gnomad.gz
	>>>
	output {
		File gnomad_joined_out = "annotation_gnomad.gz"
		File gnomad_joined_out_tbi = "annotation_gnomad.gz.tbi"
	}
	runtime {
        docker: "${docker}"
        cpu: "1"
        memory: "4 GB"
        disks: "local-disk ${local_disk} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: true
    }
}

task download_genes {

    String docker
    String genome_build
    String genes_out_fname = "genes-b" + genome_build + ".bed"

    command <<<
    
    python3 << CODE

    import io
    import os
    import re
    import csv
    import wget
    import gzip
    import boltons.iterutils
    from contextlib import contextmanager

    good_genetypes = set('''
    protein_coding
    IG_C_gene
    IG_D_gene
    IG_J_gene
    IG_V_gene
    TR_C_gene
    TR_D_gene
    TR_J_gene
    TR_V_gene
    '''.split())
    nonpseudo_genetypes = set('''
    3prime_overlapping_ncRNA
    antisense
    bidirectional_promoter_lncRNA
    lincRNA
    macro_lncRNA
    miRNA
    misc_RNA
    Mt_rRNA
    Mt_tRNA
    non_coding
    processed_transcript
    ribozyme
    rRNA
    sRNA
    scRNA
    scaRNA
    sense_intronic
    sense_overlapping
    snRNA
    snoRNA
    TEC
    vaultRNA
    '''.split()).union(good_genetypes)

    chrom_order_list = [str(i) for i in range(1, 22 + 1)] + ["X", "Y", "MT"]
    chrom_order = {chromosome: index for index, chromosome in enumerate(chrom_order_list)}
    chrom_order["23"] = 22
    chrom_order["24"] = 23
    chrom_order["25"] = 24
    chrom_aliases = {"23": "X", "24": "Y", "25": "MT", "M": "MT"}

    @contextmanager
    def read_gzip(filepath):
        # handles buffering # TODO: profile whether this is fast.
        with gzip.open(
            filepath, "rb"
        ) as f:  # leave in binary mode (default), let TextIOWrapper decode
            with io.BufferedReader(f, buffer_size=2 ** 18) as g:  # 256KB buffer
                with io.TextIOWrapper(g) as h:  # bytes -> unicode
                    yield h

    def get_all_genes(gencode_filepath):
        with read_gzip(gencode_filepath) as f:
            for l in f:
                if l.startswith('#'): continue
                r = l.split('\t')
                if r[2] != 'gene': continue

                try:
                    # Remove pseudogenes and other unwanted types of genes.
                    genetype = re.search(r'gene_type "(.+?)"', r[8]).group(1)
                    assert 'pseudogene' in genetype or genetype in nonpseudo_genetypes
                    if genetype not in good_genetypes: continue

                    assert r[0].startswith('chr')
                    chrom = r[0][3:]
                    if chrom in chrom_aliases: chrom = chrom_aliases[chrom]
                    elif chrom not in chrom_order: continue
                    pos1, pos2 = int(r[3]), int(r[4])
                    assert pos1 < pos2
                    symbol = re.search(r'gene_name "(.+?)"', r[8]).group(1)
                    full_ensg = re.search(r'gene_id "(ENSGR?[0-9\._A-Z]+?)"', r[8]).group(1)
                    ensg = full_ensg.split('.')[0]
                except Exception:
                    print('ERROR on line:', r)
                    raise

                yield {
                    'chrom': chrom,
                    'start': pos1,
                    'end': pos2,
                    'symbol': symbol,
                    'ensg': ensg,
                    'full_ensg': full_ensg,
                }

    def dedup_ensg(genes):
        # If two genes share the same "ENSGXXXX" (before period), then use their "ENSGXXXX.XXX" instead.
        for ensg_group in boltons.iterutils.bucketize(genes, key=lambda g:g['ensg']).values():
            if len(ensg_group) == 1:
                del ensg_group[0]['full_ensg']
                yield ensg_group[0]
            else:
                # These are all psuedo-autosomals across X/Y
                assert sorted(g['chrom'] for g in ensg_group) == ['X', 'Y']
                for g in ensg_group:
                    g['ensg'] = g['symbol'] = g.pop('full_ensg')
                    yield g

    def dedup_symbol(genes):
        # If genes share the same SYMBOL, check that they are adjacent and then merge them
        for symbol_group in boltons.iterutils.bucketize(genes, key=lambda g:g['symbol']).values():
            if len(symbol_group) == 1:
                yield symbol_group[0]
            elif (boltons.iterutils.same(g['chrom'] for g in symbol_group) and
                all(g1['end'] + 600e3 > g2['start'] for g1,g2 in boltons.iterutils.pairwise(sorted(symbol_group, key=lambda g:g['start'])))):
                # 600kb is long enough to resolve all problems.
                yield {
                    'chrom': symbol_group[0]['chrom'],
                    'start': min(g['start'] for g in symbol_group),
                    'end': min(g['end'] for g in symbol_group),
                    'symbol': symbol_group[0]['symbol'],
                    'ensg': ','.join(g['ensg'] for g in symbol_group),
                }
            else:
                print('broken symbol_group:')
                for g in symbol_group:
                    print('- {:12,}\t{:12,}\t{}'.format(g['start'], g['end'], g))
                raise
    
    if ${genome_build} == '37':
        url="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz"
    else:
        url="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz"
    
    wget.download(url=url, out="gencode.annotation.gtf.gz")
    genes = get_all_genes("gencode.annotation.gtf.gz")
    genes = dedup_ensg(genes)
    genes = dedup_symbol(genes)

    with open("${genes_out_fname}", 'w') as f:
        writer = csv.DictWriter(f, delimiter='\t', fieldnames='chrom start end symbol ensg'.split(), lineterminator='\n')
        writer.writerows(genes)

    CODE

    >>>    

    output {
        File genes_out = "${genes_out_fname}"
    }

    runtime {
        docker: "${docker}"
        memory: "2 GB"
        cpu: "1"
        disks: "local-disk 10 SSD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
        noAddress: false
    }

}

task add_nearest_genes {

    File genes
    String docker
    File annotation
    String fout = basename(annotation)
    Int disk_size = ceil(size(annotation, "GB")*3) + 20

    command <<<

    touch ${fout}

    python3 << CODE

    import csv
    import sys
    import gzip
    import bgzip
    import bisect
    import argparse
    import collections
    if sys.version_info.major == 3 and sys.version_info.minor >= 10:
        from collections.abc import MutableSet
        collections.MutableSet = collections.abc.MutableSet
    else: 
        from collections import MutableSet
    import intervaltree
    import boltons.iterutils

    chrom_order_list = [str(i) for i in range(1, 22 + 1)] + ["X", "Y", "MT"]
    chrom_order = {chromosome: index for index, chromosome in enumerate(chrom_order_list)}
    chrom_order["23"] = 22
    chrom_order["24"] = 23
    chrom_order["25"] = 24

    def get_gene_tuples(genes_filepath=None, include_ensg=False):
        """
        Get gene tuples.

        Very unsure what this is about.

        @param include_ensg:
        @return: 4-tuple
        """
        with open(genes_filepath, encoding="utf-8") as file:
            for row in csv.reader(file, delimiter="\t"):
                assert row[0] in chrom_order, row[0]
                if include_ensg:
                    yield row[0], int(row[1]), int(row[2]), row[3], row[4]
                else:
                    yield row[0], int(row[1]), int(row[2]), row[3]

    class BisectFinder(object):
        '''Given a list like [(123, 'foo'), (125, 'bar')...], BisectFinder helps you find the things before and after 124.'''
        def __init__(self, tuples):
            '''tuples is like [(123, 'foo'),...]'''
            tuples = sorted(tuples, key=lambda t:t[0])
            self._nums, self._values = list(zip(*tuples))
        def get_item_before(self, pos):
            '''If we get an exact match, let's return it'''
            idx = bisect.bisect_right(self._nums, pos) - 1 # note: bisect_{left,right} just deals with ties.
            if idx < 0: return None # It's fallen off the beginning!
            return (self._nums[idx], self._values[idx])
        def get_item_after(self, pos):
            if pos > self._nums[-1]: return None # it's fallen off the end!
            idx = bisect.bisect_left(self._nums, pos)
            return (self._nums[idx], self._values[idx])

    class GeneAnnotator(object):
        def __init__(self, interval_tuples):
            '''intervals is like [('22', 12321, 12345, 'APOL1'), ...]'''
            self._its = {}
            self._gene_starts = {}
            self._gene_ends = {}
            for interval_tuple in interval_tuples:
                chrom, pos_start, pos_end, gene_name = interval_tuple
                assert isinstance(pos_start, int)
                assert isinstance(pos_end, int)
                if chrom not in self._its:
                    self._its[chrom] = intervaltree.IntervalTree()
                    self._gene_starts[chrom] = []
                    self._gene_ends[chrom] = []
                self._its[chrom].add(intervaltree.Interval(pos_start, pos_end, gene_name))
                self._gene_starts[chrom].append((pos_start, gene_name))
                self._gene_ends[chrom].append((pos_end, gene_name))
            for chrom in self._its:
                self._gene_starts[chrom] = BisectFinder(self._gene_starts[chrom])
                self._gene_ends[chrom] = BisectFinder(self._gene_ends[chrom])
        
        def annotate_position(self, chrom, pos):
            if chrom == 'MT': chrom = 'M'
            if chrom not in self._its:
                return ''
            overlapping_genes = self._its[chrom].search(pos)
            if overlapping_genes:
                return ','.join(sorted(boltons.iterutils.unique_iter(og.data for og in overlapping_genes)))
            nearest_gene_end = self._gene_ends[chrom].get_item_before(pos)
            nearest_gene_start = self._gene_starts[chrom].get_item_after(pos)        
            if nearest_gene_end is None or nearest_gene_start is None:
                if nearest_gene_end is not None: return nearest_gene_end[1]
                if nearest_gene_start is not None: return nearest_gene_start[1]
                print('This is very surprising - {!r} {!r}'.format(chrom, pos))
                return ''
            dist_to_nearest_gene_end = abs(nearest_gene_end[0] - pos)
            dist_to_nearest_gene_start = abs(nearest_gene_start[0] - pos)
            if dist_to_nearest_gene_end < dist_to_nearest_gene_start:
                return nearest_gene_end[1]
            return nearest_gene_start[1]

    def annotate_genes(in_filepath, genes_filepath, out_filepath):
        '''Both args are filepaths'''
        ga = GeneAnnotator(get_gene_tuples(genes_filepath=genes_filepath))
        with open(out_filepath, 'ab') as raw, gzip.open(in_filepath, 'rt') as in_f:
            reader = csv.reader(in_f, delimiter='\t', lineterminator="\n")
            for i, v in enumerate(reader):
                if i == 0:
                    h = {f:j for j,f in enumerate(v)}
                    v.append('nearest_genes')
                else:
                    ng = ga.annotate_position(v[h['chr']], int(v[h['pos']]))
                    v.append(ng)
                line = '\t'.join(v) + '\n'
                with bgzip.BGZipWriter(raw) as out_f:
                    out_f.write(line.encode())

    annotate_genes("${annotation}", "${genes}", "${fout}")

    CODE

    tabix -b 3 -e 3 -s 2 "${fout}"

    >>>    

    output {
        File anno_ng_out = "${fout}"
		File anno_ng_out_tbi="${fout}.tbi"
    }

    runtime {
        docker: "${docker}"
        memory: "16 GB"
        cpu: "4"
        disks: "local-disk ${disk_size} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 2
    }

}


workflow scrape_annots {

	File vcfs

	Array[String] files = read_lines(vcfs)

	String docker

    File? genes_bed

    String genome_build

	## write readme and specify annot: output will be in the same order as input files.... external annot will be matched by first col. id chr:pos:ref:alt
	## join external annotation already in the extract for speed!!

	File? external_annot

	scatter(file in files) {
		call extract {
			input: vcf=file,
      		docker=docker
		}
	}
	call join_annot {
		input: files=extract.out,
		external_annot=external_annot,
		docker=docker
	}

	call add_rsids {
		input: docker=docker,
		annotation_file = join_annot.out
	}

    call add_gnomad {
        input: docker=docker,
        annotation_file = add_rsids.rsid_out,
    }

    if (!defined(genes_bed)) {
        call download_genes {
            input: docker = docker, 
            genome_build = genome_build
        }
    }

    # Optional types can be coalesced by using the select_first function
    File genes = select_first([genes_bed, download_genes.genes_out])

    call add_nearest_genes {
        input: annotation  = add_gnomad.gnomad_joined_out,
        genes = genes,
        docker = docker
    }

	call small {
		input: annotation =  add_nearest_genes.anno_ng_out,
		docker=docker
	}

	output {
		File annotation = add_nearest_genes.anno_ng_out
		File annotatation_index = add_nearest_genes.anno_ng_out_tbi
		File annotation_small = small.small_out
		File annotation_small_index = small.small_index
	}

}
