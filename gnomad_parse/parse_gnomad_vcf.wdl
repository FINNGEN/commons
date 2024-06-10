version 1.0

workflow parse_gnomad_vcf {

    input {
        Array[File] vcfs
        String docker
    }

    scatter(vcf in vcfs) {
        call parse {
            input:
                vcf = vcf,
                docker = docker
        }
    }

    call concat {
        input:
            parsed = parse.parsed,
            docker = docker
    }

    output {
        Array[File] parsed = parse.parsed
        File merged = concat.merged
    }
}

task parse {

    input {
        File vcf
        String docker
    }

    String base = sub(sub(basename(vcf, ".bgz"), ".gz$", ""), ".vcf$", "")
    Int local_disk = ceil(size(vcf, "G") * 2 + 2)

    command <<<
        
        set -euxo pipefail

        bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\t%INFO/AF_fin\t%INFO/AF_nfe\t%INFO/AF_afr\t%INFO/AF_amr\t%INFO/AF_asj\t%INFO/AF_eas\t%INFO/AF_mid\t%INFO/AF_sas\t%INFO/AC_fin\t%INFO/AC_nfe\t%INFO/AC_afr\t%INFO/AC_amr\t%INFO/AC_asj\t%INFO/AC_eas\t%INFO/AC_mid\t%INFO/AC_sas\t%INFO/AN_fin\t%INFO/AN_nfe\t%INFO/AN_afr\t%INFO/AN_amr\t%INFO/AN_asj\t%INFO/AN_eas\t%INFO/AN_mid\t%INFO/AN_sas\n' ~{vcf} | \
        awk '
        BEGIN {
            FS=OFS="\t"
        } NR == 1 {
            for(i=1; i<=NF; i++) {
                sub(/\[.*\]/, "", $i);
                sub(/ /, "", $i)
                h[$i] = i
            }
            print $0, "enrichment_nfe"
        } NR > 1 {
            for(i=1; i<=NF; i++) {
                if ($i == ".") {
                    $i = "NA"
                }
            }
            enr = "NA"
            if ($h["AF_nfe"] > 0 && $h["AF_fin"] != "NA" && $h["AF_nfe"] != "NA") {
                enr = $h["AF_fin"] / $h["AF_nfe"]
            }
            print $0, enr
        }' | bgzip > ~{base}.tsv.gz

        tabix -s 1 -b 2 -e 2 ~{base}.tsv.gz

    >>>

    output {
        File parsed = "~{base}.tsv.gz"
        File parsed_tbi = "~{base}.tsv.gz.tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "4"
        memory: "32 GB"
        disks: "local-disk ~{local_disk} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 1
        noAddress: true
    }
}

task concat {

    input {
        Array[File] parsed
        String docker
    }

    String base = sub(basename(parsed[0], ".tsv.gz"), ".chr[0-9]*$", "")
    Int local_disk = ceil(size(parsed, "G") * 2 + 2)

    command <<<
        
        set -euxo pipefail

        cat <(zcat -f ~{parsed[0]} | head -1) \
            <(for file in ~{sep=" " parsed}; do
                zcat -f $file | tail -n +2;
            done) \
        | bgzip > ~{base}.merged.tsv.gz

        echo "`date` tabixing"
        tabix -s 1 -b 2 -e 2 ~{base}.merged.tsv.gz
        echo "`date` done"

    >>>

    output {
        File merged = "~{base}.merged.tsv.gz"
        File merged_tbi = "~{base}.merged.tsv.gz.tbi"
    }

    runtime {
        docker: "~{docker}"
        cpu: "2"
        memory: "8 GB"
        disks: "local-disk ~{local_disk} HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: 1
        noAddress: true
    }
}