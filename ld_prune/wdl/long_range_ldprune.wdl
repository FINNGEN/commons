workflow long_range_prune {
    String fmapfile
    Array[String] fmaps = read_lines(fmapfile)

    String tomahawk_pattern
    String path_to_tw_map 


    scatter (fmap in fmaps) {
        call ld_prune{input: fmap = fmap, tomahawk_pattern = tomahawk_pattern, 
        tw_map = path_to_tw_map}
    }


  }

  task ld_prune {
    
    File fmap 
    File tw_map
    File tw_map_index=tw_map+".tbi"
    
    String tomahawk_pattern

    String toma_bucket = sub(sub(tomahawk_pattern, "^gs://", ""), "/.+$", "")
    String mountpoint = "/cromwell_root/gcsfuse/" + toma_bucket
    String new_toma_pattern_fuse = sub(tomahawk_pattern,"^gs://" + toma_bucket, mountpoint)
    String new_toma_pattern_nofuse = sub(tomahawk_pattern, "^gs://.*/", mountpoint + "/")

    Int ld_width=20000000
    Int clump_expected_chisq=5
    Float clump_expected_chisq_af=1.0

    Boolean enable_fuse = true

    String base = basename(fmap)

    command<<<

        #!/usr/bin/env bash
        # mount bgen bucket
        mkdir -p ${mountpoint}
        echo "created mount point bucket"
        if [[ ${enable_fuse} == "true" ]]
        then
            echo "mounting bucket ${toma_bucket} at ${mountpoint}"
            gcsfuse --debug_fuse --debug_fuse_errors --implicit-dirs ${toma_bucket} ${mountpoint}
            new_toma_pattern_local="${new_toma_pattern_fuse}"
        else
            dir=$(dirname ${mountpoint})
            mkdir -p $dir
            echo gsutil -q cp $(echo ${tomahawk_pattern}| sed 's/{CHR}/*/') ${mountpoint}/
            gsutil -q cp $(echo ${tomahawk_pattern}| sed 's/{CHR}/*/') ${mountpoint}/
            new_toma_pattern_local="${new_toma_pattern_nofuse}"
        fi

        echo "mounted bucket ${toma_bucket} at ${mountpoint}"
        ls -latr ${mountpoint}

        ls -latr /cromwell_root/gcsfuse/finngen-imputation-panel/sisu4.2
        ls -latr /cromwell_root/gcsfuse/finngen-imputation-panel/sisu4.2/tomahawk/

        echo ld_prune_lines.py ${fmap} ${base}.dynclump -ld_w ${ld_width} \
            -clump_expected_chisq ${clump_expected_chisq}  -clump_expected_chisq_af ${clump_expected_chisq_af} \
            -ld_source sisu42 -chromcol 'chrom' -poscol pos -refcol ref -altcol alt -pcol pval \
            -prune_column_list locus_id,phenotype,lead_beta,pval \
            -local_tomahawk_LD -tomahawk_template $new_toma_pattern_local  \
            -tomahawk_mapfile ${tw_map} -sort_first
        
        ld_prune_lines.py ${fmap} ${base}.dynclump -ld_w ${ld_width} \
            -clump_expected_chisq ${clump_expected_chisq}  -clump_expected_chisq_af ${clump_expected_chisq_af} \
            -ld_source sisu42 -chromcol 'chrom' -poscol pos -refcol ref -altcol alt -pcol pval \
            -prune_column_list locus_id,phenotype,lead_beta,pval \
            -local_tomahawk_LD -tomahawk_template $new_toma_pattern_local  \
            -tomahawk_mapfile ${tw_map} -sort_first

    >>>    

    runtime {
        docker: "europe-west1-docker.pkg.dev/finngen-refinery-dev/fg-refinery-registry/ld_prune:latest"
        memory: "8 GB"
        cpu: "16"
        disks: "local-disk 20 HDD"
        zones: "europe-west1-b"
    }
   output {
        File dynclump = "${base}.dynclump"
    }

  }




