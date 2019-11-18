import liftover.wdl as lift


workflow copyfiles {
    File filestocopy
    String docker

    Array[Array[String]] files = read_tsv(filestocopy)

    scatter (file in files) {
        call docopy {
            input: filename=file[0],
                   url=file[1],
                   docker=docker
        }

        call lift.dolift {
            input: resfile=docopy.resfile,
            docker=docker
        }
    }

}


task docopy {
    String url
    String filename
    String docker

    runtime {
          docker: "${docker}"
          memory:"1G"
          disks: "local-disk 100 SSD"
          cpu:"1"
    }

    command <<<
        wget ${url} -O ${filename}
    >>>

    output {
        File resfile = "${filename}"
    }
}


task dolift {
  File resfile
  Int memory
  String docker
  Int local_disk

  command <<<
    zcat ${resfile} | awk 'NR>1{ split($1,a,":"); len=length(a[4]); if(length(a[3])>len) len=length(a[3]); print "chr"a[1],a[2]-1,a[2]+len-1,$1  }'  > variants.bed

    liftOver variants.bed /hg19ToHg38.over.chain.gz variants_lifted errors

    export LCTYPE=C
    export LANG=C

    cols=$(zcat ${resfile} | head -n 1 | awk '{ print NF}')

    join -1 1 -2 3 -t$'\t' <(  zcat ${resfile} | awk 'NR==1{ print "#"$0 } NR>1{ print $0}' | sort -b -k 1,1 -t $'\t'  ) \
    <( awk 'BEGIN{ OFS="\t"; print "achr38","apos38","#variant" }{ print $1,$2+1,$4}'  variants_lifted | sort -t$'\t' -b -k 3,3 ) \
    | sort -V -k $((cols+1)),$((cols+1))  -k $((cols+2)),$((cols+2)) | awk  'BEGIN{ FS="\t"; OFS="\t"} NR==1{ print $0,"REF","ALT"} NR>1{ split($1,a,":"); print $0,a[3],a[4] } '|  bgzip > ${basename(resfile)}".lifted.gz"
    tabix -s $((cols+1)) -b $((cols+2)) -e $((cols+2)) ${basename(resfile)}".lifted.gz"
  >>>

  output {
        File lifted="${basename(resfile)}.lifted.gz"
        File index="${basename(resfile)}.lifted.gz.tbi"
    }

  runtime {
        docker: "${docker}"
        memory:"${memory}G"
        disks: "local-disk ${local_disk} SSD"
        cpu:"1"
    }

}
