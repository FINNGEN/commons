#!/usr/bin/env bash
##
##     @script.name [option] ARGUMENTS...
##
## Options:
##     -h, --help       All client scripts have this, it can be omitted.
##     --var=VALUE      Index for variant column chr:pos:ref:alt. This or next four must be specified
##     --chr=VALUE      Columnn indexes
##     --pos=VALUE
##     --ref=VALUE
##     --alt=VALUE
##     --tmp=VALUE      Specify tmp directory location for sorting and joining...
##     --chr_as_is      Don't convert chromosome names to numeric
##     --no_duplicates  Don't report duplicates



export LCTYPE=C
export LANG=C
my_dir="$(dirname "$0")"
source "${my_dir}/easyoptions.sh" || exit

inputfile=$1
liftedfile=$2

cat_cmd="cat"
file $inputfile | grep gzip > /dev/null 2>&1
if [[ $? -eq 0 ]];
then
    cat_cmd='zcat'
fi

cols=$( $cat_cmd < "$inputfile" | head -n 1 | awk 'BEGIN {FS="\t"}{ print NF}')

temploc="-T /tmp/"

if [[ -n "$tmp" ]]
then
	temploc="-T $tmp"
fi

if [[ ! -n "$var" ]]
then
    if [[ ! -n "$chr" ]] || [[ ! -n "$pos" ]] || [[ ! -n "$ref" ]] || [[ ! -n "$alt" ]]
    then
        echo "chr pos ref alt columns must be specified if var not given"
        exit 1
    fi
    
    cols=$((cols+1))
    anew_chr=$((cols+1))

    join -1 1 -2 3 -t$'\t' \
    <( $cat_cmd < "$inputfile" | awk -v chr=$chr -v pos=$pos -v ref=$ref -v alt=$alt 'BEGIN{OFS="\t"} NR==1{ print "#variant",$0 } NR>1{ print $chr":"$pos":"$ref":"$alt,$0 }' | sort -b -k 1,1 -t $'\t' $temploc ) \
    <( awk 'BEGIN{ OFS="\t"; print "anew_chr","anew_pos","#variant" }{ print $1,$2+1,$4}' $liftedfile | sort $temploc -t$'\t' -b -k 3,3 ) \
    | sort $temploc -t$'\t' -V -k $((cols+1)),$((cols+1)) -k $((cols+2)),$((cols+2)) \
    | awk -v anew_chr=$anew_chr -v chr=$chr_as_is 'BEGIN{ FS="\t"; OFS="\t"} NR==1{ print $0,"REF","ALT"} NR>1{ split($1,a,":"); if(chr!="yes"){ gsub("chr", "", $anew_chr); gsub("X", "23", $anew_chr); gsub("Y", "24", $anew_chr); gsub("M", "25", $anew_chr) } print $0,a[3],a[4] } ' \
    | bgzip > $(basename $inputfile)".lifted.gz"
else
    anew_chr=$((cols+1))

    join -1 1 -2 3 -t$'\t' \
    <( $cat_cmd < "$inputfile" | awk -v var=$var 'NR==1{ printf "#variant"; for(i=1;i<=NF; i++) if(i!=var) printf "\t"$i; printf "\n"; } NR>1{ printf $var; for(i=1;i<=NF; i++) if(i!=var) printf "\t"$i; printf "\n"; }' | sort -b -k 1,1 -t $'\t' $temploc ) \
    <( awk 'BEGIN{ OFS="\t"; print "anew_chr","anew_pos","#variant" }{ print $1,$2+1,$4}' $liftedfile | sort $temploc -t$'\t' -b -k 3,3 ) \
    | sort $temploc -t$'\t' -V -k $((cols+1)),$((cols+1)) -k $((cols+2)),$((cols+2)) \
    | awk -v anew_chr=$anew_chr -v chr=$chr_as_is 'BEGIN{ FS="\t"; OFS="\t"} NR==1{ print $0,"REF","ALT"} NR>1{ split($1,a,":"); if(chr!="yes"){ gsub("chr", "", $anew_chr); gsub("X", "23", $anew_chr); gsub("Y", "24", $anew_chr); gsub("M", "25", $anew_chr) } print $0,a[3],a[4] } ' \
    | bgzip > $(basename $inputfile)".lifted.gz"
fi

if [[ ! -n "$no_duplicates" ]]
then
    echo "Finding duplicated variants"
    zcat $(basename $inputfile)".lifted.gz" \
    | awk -v anew_chr=$anew_chr '
    {
        variant=$anew_chr"_"$(anew_chr+1)"_"$(anew_chr+2)"_"$(anew_chr+3)
        variant_line=$0
        if (variant==prev_variant) {
            printf "%s\n%s\n", prev_variant_line, variant_line
        }
        prev_variant=variant
        prev_variant_line=variant_line
    }' > duplicated_lifted.tsv
fi

tabix -s $((cols+1)) -b $((cols+2)) -e $((cols+2)) $(basename $inputfile)".lifted.gz"