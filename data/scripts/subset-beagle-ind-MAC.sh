#!/bin/sh

############################################
# Cut beagle file individuals into k partitions
# Matt DeSaix
# August 23, 2024
# Description: User provides beagle file and it will be split into 
# k new files by outputting the header and every kth line (marker) to a new file
############################################

############################################
# Define help function to help in use and debugging
############################################

helpfunction()
{
  echo ""
  echo "Usage: $0 -b beagle -i individuals -o outname"
  echo "\t-b [string] Name of Beagle file"
  echo "\t-i [string] Name of individuals file"
  echo "\t-o [string] Full out name"
}

while getopts "b:i:o:" opt
do
  case "$opt" in
    b ) beagle="$OPTARG";;
    i ) individuals="$OPTARG";;
    o ) outname="$OPTARG";;
    ? ) helpfunction ;; # Print help function in case parameter is non-existent
  esac
done

if [ -z "$beagle" ] || [ -z "$individuals" ] || [ -z "$outname" ]
then
  echo "Some or all of the parameters are empty"
  helpfunction
fi

set +o posix

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
 ' ${individuals} <(zcat -f < ${beagle}) | gzip > ${outname}

