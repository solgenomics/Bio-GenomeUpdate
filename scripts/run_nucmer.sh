#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Run nucmer pipeline to align BACs (param 1) to chr ref (param 2). NOT USING delta-filter -u 99 as diff BACs can contain same repeat and that does not make the alignment to ref chr invalid. Also BAC can map to chr00 and chr1-12 which is important to identify cases where contigs from chr00 can be integrated into chr1-12.

usage(){
	echo "usage:
	$0 <BAC fasta> <REF fasta>"
	exit 1
}

if [ "$#" -ne 2 ]
then
	usage
fi

printf "Query : %s \n" "$1"
printf "Reference : %s \n" "$2"

nucmer -l 100 -c 500 --noextend -p 500bp_qry_"${1}"__ref_"${2}" "$2" "$1"
#delta-filter -l 500 -u 99 500bp_qry_"${1}"__ref_"${2}".delta > 500bp_qry_"${1}"__ref_"${2}".delta.filtered
delta-filter -l 500 500bp_qry_"${1}"__ref_"${2}".delta > 500bp_qry_"${1}"__ref_"${2}".delta.filtered
show-coords -c -d -l -q -T -o 500bp_qry_"${1}"__ref_"${2}".delta.filtered > 500bp_qry_"${1}"__ref_"${2}".delta.filtered.coords


