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
  echo "Usage: $0 -b beagle -k partitions -o outname"
  echo "\t-b [string] Name of Beagle file"
  echo "\t-k [integer] K partitions to split"
  echo "\t-o [string] outname prefix to name "
}

while getopts "b:k:o:" opt
do
  case "$opt" in
    b ) beagle="$OPTARG";;
    k ) k="$OPTARG";;
    o ) out="$OPTARG";;
    ? ) helpfunction ;; # Print help function in case parameter is non-existent
  esac
done

if [ -z "$beagle" ] || [ -z "$k" ] || [ -z "$out" ]
then
  echo "Some or all of the parameters are empty"
  helpfunction
fi

for i in $(seq 1 "$k")
do
  outname=${out}.${i}.beagle.gz
  index=$((i-1))
  zcat -f ${beagle} | awk -v N=$index -v K=$k 'NR==1 || NR%K==N' | gzip > ${outname}
done
