version 1.0

workflow long_range_prune {

    input {
        String fmapfile
        Array[String] fmaps = read_lines(fmapfile)

        String tomahawk_pattern
        String path_to_tw_map 
    }

    scatter (fmap in fmaps) {

        call split_fmap {
            input:
                fmap = fmap
        }
        
        scatter (fmap_chunk in split_fmap.fmap_chunks) {
            call ld_prune {
                input:
                    fmap = fmap_chunk,
                    tomahawk_pattern = tomahawk_pattern, 
                    tw_map = path_to_tw_map
            }
        }
    }

  }

  task split_fmap {

    input {
        File fmap
        String docker
        String chromcol = "chrom"
    }

    String base = basename(fmap)

    command <<<

        zcat -f ~{fmap} | awk 'BEGIN{FS=OFS="\t"} FNR==1 {header=$0; for(i=1;i<=NF;i++) {h[$i]=i} next} {out="~{base}.chr"$h["~{chromcol}"]".tsv"; if (!seen[$h["~{chromcol}"]]++) print header > out; print > out}'

    >>>

    runtime {
        docker: "~{docker}"
        memory: "2 GB"
        cpu: "1"
        disks: "local-disk 20 HDD"
        preemptible: 1
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        noAddress: true
    }

    output {
        Array[File] fmap_chunks = glob("~{base}.chr*")
    }
  }

  task ld_prune {
    
    input {
        File fmap 
        File tw_map
        String tomahawk_pattern
        Int ld_width=20000000
        Float clump_ld
        Boolean enable_fuse = true
        String docker
        String chromcol = "chrom"
        String poscol = "pos"
        String refcol = "ref"
        String altcol = "alt"
        String pcol = "pval"
    }
    
    File tw_map_index=tw_map+".tbi"

    String toma_bucket = sub(sub(tomahawk_pattern, "^gs://", ""), "/.+$", "")
    String mountpoint = "/cromwell_root/gcsfuse/" + toma_bucket
    String new_toma_pattern = sub(tomahawk_pattern,"^gs://" + toma_bucket, mountpoint)

    String base = basename(fmap)

    command<<<

        #!/usr/bin/env bash
        # mount bgen bucket
        mkdir -p ~{mountpoint}
        echo "created mount point bucket"
        if [[ ~{enable_fuse} == "true" ]]
        then
            echo "mounting bucket ~{toma_bucket} at ~{mountpoint}"
            gcsfuse --debug_fuse --debug_fuse_errors --implicit-dirs ~{toma_bucket} ~{mountpoint}
        else
            dir=$(dirname ~{new_toma_pattern})
            mkdir -p $dir
            echo gsutil -q cp $(echo ~{tomahawk_pattern}| sed 's/{CHR}/*/') $dir/
            gsutil -q cp $(echo ~{tomahawk_pattern}| sed 's/{CHR}/*/') $dir/
        fi

        echo "mounted bucket ~{toma_bucket} at ~{mountpoint}"
        ls -latr ~{mountpoint}

        ls -latr /cromwell_root/gcsfuse/finngen-imputation-panel/sisu4.2
        ls -latr /cromwell_root/gcsfuse/finngen-imputation-panel/sisu4.2/tomahawk/

        echo ld_prune_lines.py ~{fmap} ~{base}.dynclump -ld_w ~{ld_width} \
            -ld ~{clump_ld} -pcol ~{pcol} \
            -ld_source sisu42 -chromcol ~{chromcol} \
            -poscol ~{poscol} -refcol ~{refcol}  -altcol ~{altcol} -prune_column_list locus_id,phenotype,lead_beta,pval \
            -local_tomahawk_LD -tomahawk_template ~{new_toma_pattern}  \
            -tomahawk_mapfile ~{tw_map}
        
        ld_prune_lines.py ~{fmap} ~{base}.dynclump -ld_w ~{ld_width} \
        -ld ~{clump_ld} -pcol ~{pcol} \
        -ld_source sisu42 -chromcol ~{chromcol} \
        -poscol ~{poscol} -refcol ~{refcol}  -altcol ~{altcol} -prune_column_list locus_id,phenotype,lead_beta,pval \
        -local_tomahawk_LD -tomahawk_template ~{new_toma_pattern}  \
        -tomahawk_mapfile ~{tw_map} -sort_first

    >>>    

    runtime {
        docker: "~{docker}"
        memory: "8 GB"
        cpu: "16"
        disks: "local-disk 30 HDD"
        preemptible: 1
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        noAddress: true
    }

   output {
        File dynclump = "~{base}.dynclump"
    }

  }




