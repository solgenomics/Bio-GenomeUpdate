#!/bin/sh

# Surya Saha
# BTI
# Purpose: Uses mummer tools. Makes dot plots for full chr length alignment with 1-1 mapping between ref and query

usage(){
	echo "usage:
	$0 <ref.fa> <query.fa> <prefix for outfiles> <min tile length> <min identity perc>"
	exit 1
}

if [ "$#" -ne 5 ]
then
	usage
fi

printf "Ref file: %s \n" "$1"
printf "Query file: %s \n" "$2"
printf "Prefix: %s \n" "$3"
printf "Min tile length : %d \n" "$4"
printf "Min id perc: %d \n\n" "$5"

time nucmer -p "$3" -o --noextend "$1" "$2" 
#only 1-1 mapping between ref and query
delta-filter -1 -l "$4" -i "$5" "${3}".delta > filtered."${3}".delta
mummerplot --png --prefix="${3}" filtered."${3}".delta

