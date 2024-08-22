#!/bin/sh

awk '
BEGIN { FS=OFS="\t"; headers[1]="marker"; headers[2]="allele1"; headers[3]="allele2"}
FNR==NR { headers[NR+3]=$1; next }
        { sep=""
        for (j=1; j<=length(headers); j++ ){
	    for (i=1; i<=NF; i++) {
		if (FNR==1 && ($i == headers[j])) {
		    k++
		    cols[k] = i
		}

	    }
	    }
	    for (n=1; n<=length(cols); n++) {
		    printf "%s%s",sep,$cols[n]
		    sep=OFS
	    }
	    print ""
	}
 ' $1 $2

