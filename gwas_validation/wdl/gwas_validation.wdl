version 1.0

workflow gwas_validation {
    input{
        File phenolist
    }
    Array[String] phenos = read_lines(phenolist)
    
    
    scatter( p in phenos ) {
        call validate {
            input: 
                phenotype=p
        }
    }

    output {
        Array[File] current_1e5_filtered = validate.table
        Array[File] previous_5e8_filtered = validate.table2
        Array[File] curent_new_hits = validate.table3
        Array[File] plots = validate.plots

    }
}

task validate {
    input{
        String docker
        String phenotype
        String current_path
        String previous_path
        String current_label
        String prev_label
        Array[String] colnames
        File annotation
        String zones
    }
    File current_file = sub(current_path,"\\{PHENO\\}",phenotype)
    File previous_file = sub(previous_path,"\\{PHENO\\}",phenotype)
    File annotation_tbi = annotation + ".tbi"

    Float pval_threshold =  0.00001
    Float gwsig_threshold = 0.00000005
    Float new_hit_range = 1500000
    
    command <<<
        set -eux
        # filter significant values in first file, by level 1e-5
        CURRENT="~{current_file}"
        PREVIOUS="~{previous_file}"
        PVAL=$(zcat $CURRENT|head -n1|tr "\t" "\n"|grep -n "~{colnames[4]}"|cut -d ":" -f1)
        cat <(zcat ~{current_file}|head -n1) <(zcat ~{current_file}|tail -n+2|awk -v pval=$PVAL '{if($pval<~{pval_threshold}){print $0}}') > current_sig.tsv
        cat <(zcat ~{previous_file}|head -n1) <(zcat ~{previous_file}|tail -n+2|awk -v pval=$PVAL '{if($pval<~{gwsig_threshold}){print $0}}') > previous_gwsig.tsv
        #python code to get the corresponding values in other file
        #args
        # tsv with variants
        # summstat to copy vars from
        # output name
        # direction: the copied variants end up "left" or "right" 
        cat > extract_others.py <<PYCODE
        import gzip, os, sys
        var_tsv = sys.argv[1]
        summstat = sys.argv[2]
        out_name = sys.argv[3]
        direction = sys.argv[4]
        chrom = "~{colnames[0]}"
        pos = "~{colnames[1]}"
        ref = "~{colnames[2]}"
        alt = "~{colnames[3]}"
        pval = "~{colnames[4]}"
        beta = "~{colnames[5]}"
        mlogp = "~{colnames[6]}"
        data = {}
        with open(var_tsv) as f:
            header = f.readline().strip("\n").split("\t")
            h_idx = {a:i for (i,a) in enumerate(header)}
            for line in f:
                splitline = line.strip("\n").split("\t")
                idx = (
                    splitline[h_idx[chrom]],
                    splitline[h_idx[pos]],
                    splitline[h_idx[ref]],
                    splitline[h_idx[alt]]
                )
                values = (
                    splitline[h_idx[pval]],
                    splitline[h_idx[beta]],
                    splitline[h_idx[mlogp]]
                )
                data[idx] = values

        with gzip.open(summstat,"rt") as prevfile:
            with open(out_name,"wt") as outf:
                outheader=(
                    "chrom",
                    "pos",
                    "ref",
                    "alt",
                    "p1",
                    "p2",
                    "b1",
                    "b2",
                    "mlogp1",
                    "mlogp2",
                )
                outf.write("\t".join(outheader)+"\n")

                header = prevfile.readline().strip("\n").split("\t")
                h_idx = {a:i for (i,a) in enumerate(header)}
                for line in prevfile:
                    splitline = line.strip("\n").split("\t")
                    idx = (
                        splitline[h_idx[chrom]],
                        splitline[h_idx[pos]],
                        splitline[h_idx[ref]],
                        splitline[h_idx[alt]]
                    )

                    if idx in data:
                        c_stats = data[idx]
                        if direction == "left":
                            outline = [
                                idx[0],
                                idx[1],
                                idx[2],
                                idx[3],
                                c_stats[0],
                                splitline[h_idx[pval]],
                                c_stats[1],
                                splitline[h_idx[beta]],
                                c_stats[2],
                                splitline[h_idx[mlogp]]
                            ]
                        else:
                            outline = [
                                idx[0],
                                idx[1],
                                idx[2],
                                idx[3],
                                splitline[h_idx[pval]],
                                c_stats[0],
                                splitline[h_idx[beta]],
                                c_stats[1],
                                splitline[h_idx[mlogp]],
                                c_stats[2]
                            ]
                        outf.write("\t".join(outline)+"\n")
        PYCODE
        #get current up to 1e-5, and previous that match
        python3 extract_others.py current_sig.tsv $PREVIOUS current_1e5.tsv right
        #get previous up to 5e-8, and current that match
        python3 extract_others.py previous_gwsig.tsv $CURRENT previous_gwsig.tsv left

        #annotate files with INFO and AFs
        #TODO
        tail -n+2 current_1e5.tsv | cut -f 1-2|awk   'OFS="\t"{print $1,$2,$2}'|sort -V|uniq > current.tabix
        tail -n+2 previous_gwsig.tsv | cut -f 1-2|awk   'OFS="\t"{print $1,$2,$2}'|sort -V|uniq > previous.tabix
        #get annotations
        tabix -hR current.tabix ~{annotation}|uniq > current.annotation
        tabix -hR previous.tabix ~{annotation}|uniq > previous.annotation

        #combine annotations with tsvs
        #args: tsv to annotate, annotation,output filename
        cat > annotate.py <<PYCODE
        import os, sys
        tsv = sys.argv[1]
        annotation = sys.argv[2]
        out = sys.argv[3]
        #read in the data
        with open(tsv,"r") as tsvfile:
            header = tsvfile.readline().strip("\n").split("\t")
            h_idx = {a:i for i,a in enumerate(header)}
            data = [a.strip("\n").split("\t") for a in tsvfile]
        #read annotations to dict
        anndict = {}
        with open(annotation) as annfile:
            
            
            a_header = annfile.readline().strip("\n").split("\t")
            a_idx = {a:i for i,a in enumerate(a_header)}
            for line in annfile:
                splitline = line.strip("\n").split("\t")
                cpra = (
                    splitline[a_idx["chr"]],
                    splitline[a_idx["pos"]],
                    splitline[a_idx["ref"]],
                    splitline[a_idx["alt"]]
                )
                anndict[cpra] = [
                    splitline[a_idx["INFO"]],
                    splitline[a_idx["AF"]]
                ]

        with open(out,"w") as outf:
            outheader = header + ["INFO","AF"]
            outf.write("\t".join(outheader)+"\n")
            for d in data:
                d_idx = (
                    d[h_idx["chrom"]],
                    d[h_idx["pos"]],
                    d[h_idx["ref"]],
                    d[h_idx["alt"]]
                )
                ann = ["NA","NA"]
                if d_idx in anndict:
                    ann = anndict[d_idx]
                outline = "\t".join(d+ann) + "\n"
                outf.write(outline)
        PYCODE

        python3 annotate.py current_1e5.tsv current.annotation ~{phenotype}_1e5_current.tsv
        python3 annotate.py previous_gwsig.tsv previous.annotation ~{phenotype}_5e8_previous.tsv

        cat > find_independent_hits.py <<PYCODE
        from typing import NamedTuple
        from collections import defaultdict

        RANGE=~{new_hit_range}

        class CPRA(NamedTuple):
            c:str
            p:int
            r:str
            a:str
            pval:float

        class Region(NamedTuple):
            c:str
            start:int
            end:int

        def region_region_overlap(r1: Region, r2: Region) -> bool:
            if r1.c != r2.c or (r1.end < r2.start or r2.end < r1.start):
                return False
            return True

        def region_variant_overlap(region: Region, variant: CPRA) -> bool:
            if region.c != variant.c or (variant.p < region.start or variant.p > region.end):
                return False
            return True

        with open("~{phenotype}_5e8_previous.tsv" ) as f:
            header = f.readline().strip("\n").split("\t")
            h_idx = {a:i for i,a in enumerate(header)}
            c_idx = h_idx["chrom"]
            p_idx = h_idx["pos"]
            regions =[]
            for line in f:
                a = line.strip("\n").split("\t")
                regions.append(Region(str(a[c_idx]),int(a[p_idx])-RANGE,int(a[p_idx])+RANGE) )
            #merge
            region_chrdict = defaultdict(list)
            for r in regions:
                region_chrdict[r.c].append(r)
            for chrom in region_chrdict.keys():
                rs = region_chrdict[chrom]
                rs = sorted(rs,key=lambda x: x.start)
                stack = []
                stack.append(rs[0])
                for reg in rs[1:]:
                    top = stack.pop()
                    if region_region_overlap(top,reg):
                        stack.append(Region(top.c,min(top.start,reg.start),max(top.end,reg.end)) )
                    else:
                        stack.append(top)
                        stack.append(reg)
                region_chrdict[chrom]=stack

        with open("~{phenotype}_1e5_current.tsv") as f:
            with open("~{phenotype}_new_hits.tsv","w") as of:
                hline=f.readline()
                header = hline.strip("\n").split("\t")
                h_idx = {a:i for i,a in enumerate(header)}
                (c,p,r,a,pval) = (
                    h_idx["chrom"],
                    h_idx["pos"],
                    h_idx["ref"],
                    h_idx["alt"],
                    h_idx["p2"],
                )

                #write of header
                of.write(hline)

                for line in f:
                    data = line.strip("\n").split("\t")
                    var = CPRA(
                        data[c],
                        int(data[p]),
                        data[r],
                        data[a],
                        float(data[pval])
                    )
                    if var.pval <= 5e-8:
                        if not any([region_variant_overlap(a,var) for a in region_chrdict[var.c]]):
                            of.write(line)
        PYCODE
        python3 find_independent_hits.py
        # start up R, calculate r2, slope of beta + logpval
        cat > script.R << 'CODE'
        library(data.table)
        library(ggplot2)
        library(tidyverse)
        #define function for plotting
        model_and_plot <- function(data,plot_title,pheno,x_label,y_label) {

        #if b1 < 0, b1*-1, b2*-1
        sig = sign(data$b1)
        data$b1 = sig*data$b1
        data$b2 = sig*data$b2
        #models
        if (nrow(data)>=3){ 
        betamodel = lm(b2~0+b1,data=data)
        summary(betamodel)
        mlogpmodel = lm(mlogp2~0+mlogp1,data=data)
        summary(mlogpmodel)
        #extract params
        br2 = summary(betamodel)$r.squared
        bslope = summary(betamodel)$coefficients[1]

        mr2 = summary(mlogpmodel)$r.squared
        mslope = summary(mlogpmodel)$coefficients[1]
        }
        else {
        bslope=0
        br2=0
        mslope=0
        mr2=0
        }
        #plot limits
        b_min = min(c( min(data$b1),min(data$b2)))
        b_max = max(c( max(data$b1),max(data$b2)))
        mlogp_min = min(c( min(data$mlogp1),min(data$mlogp2)))
        mlogp_max = max(c( max(data$mlogp1),max(data$mlogp2)))

        bplot <- ggplot(data,aes(x=b1,y=b2))+geom_point()+
        geom_abline(intercept=0,slope=bslope,linetype=5)+
        geom_abline(intercept=0,slope=1,linetype=2,color="gray")+
        annotate("text",x=b_min+0.2*(b_max-b_min),y=b_min+0.8*(b_max-b_min),label=paste("r2:",format(round(br2,2),nsmall=2)),size=8)+
        xlim(b_min,b_max)+ylim(b_min,b_max)+ggtitle(label=paste(plot_title,"Betas",pheno))+
        xlab(x_label)+ylab(y_label)

        mplot <- ggplot(data,aes(x=mlogp1,y=mlogp2))+geom_point()+
        geom_abline(intercept=0,slope=mslope,linetype=5)+
        geom_abline(intercept=0,slope=1,linetype=2,color="gray")+
        annotate("text",x=mlogp_min +0.2*(mlogp_max-mlogp_min),y=mlogp_min + 0.8*(mlogp_max-mlogp_min),label=paste("r2:",format(round(mr2,2),nsmall=2)),size=8)+
        xlim(mlogp_min,mlogp_max)+ylim(mlogp_min,mlogp_max)+ggtitle(label=paste(plot_title,"mlogp",pheno))+
        xlab(x_label)+ylab(y_label)

        print(bplot)
        print(mplot)
        }
        pheno <- "~{phenotype}"
        prev_label <- "~{prev_label}"
        current_label <- "~{current_label}"

        data_1e5 <- fread("~{phenotype}_1e5_current.tsv")
        data_old_gws <- fread("~{phenotype}_5e8_previous.tsv")
        data_new_hits <- fread("~{phenotype}_new_hits.tsv")
        data_new_gws_not_in_old <- data_1e5 %>%filter(p2<5e-8 & p1>5e-8)

        pdf("~{phenotype}_plot.pdf")
        model_and_plot(data_1e5,"new 1e-5 vars",pheno,prev_label,current_label)
        model_and_plot(data_old_gws,"Old GWS vars",pheno,prev_label,current_label)
        model_and_plot(data_new_gws_not_in_old,"Vars GWS in new but not old",pheno,prev_label,current_label)
        model_and_plot(data_new_hits,"New hits",pheno,prev_label,current_label)

        dev.off()
        CODE
        Rscript script.R

    >>>

    output {
        File table = phenotype + "_1e5_current.tsv"
        File table2 = phenotype +"_5e8_previous.tsv"
        File table3 = phenotype +"_new_hits.tsv"
        File plots = phenotype + "_plot.pdf"
    }

    runtime {
        docker: "${docker}"
        cpu: 1
        memory: "6 GB"
        disks: "local-disk 20 HDD"
        zones: "${zones}"
        preemptible: 2
        noAddress: true
    }
}