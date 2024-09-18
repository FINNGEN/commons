
## on .top.out file.... headers taken from this and rest merged 
one=$1
out=$2
base=$(dirname $one)

cat <(head -n 1 $one | awk 'BEGIN{ FS=OFS="\t"} { print "pheno",$0}') <(find $base/*.top.out | while read f; do p=$(basename $f); p=${p%%.top.out};  awk -v p=$p 'BEGIN{ FS=OFS="\t" } NR>1{ if(NF<19) {print p,$0,"NA","NA"} else { print p,$0}  }' $f; done) > $out


