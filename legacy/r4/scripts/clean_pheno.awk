BEGIN {OFS=FS="\t"}
{
    if (NR==1) {
	printf($1)
	for (i=2; i<=NF; i++) {
	    if (!($i ~ /^PREVAL_/ || $i ~ /^INCIDENT/ || $i ~ /_AGEDIFF$/ || $i ~ /_AGE$/ || $i ~ /_YEAR$/ || $i ~ /_NEVT$/)) {
		ix[i] = 1
		printf("\t"$i)
	    }
	    if ($i ~ /^DEATH_AGE$/ || $i ~ /^BL_AGE$/ || $i ~ /^BL_YEAR$/) {
		ix[i] = 1
		printf("\t"$i)
	    }
	}
	printf("\n")
    }
    if(NR>1) {
	printf($1)
	for (i=2; i<=NF; i++) {
	    if (ix[i]) printf("\t"$i)
	}
	printf("\n")
    }
}
