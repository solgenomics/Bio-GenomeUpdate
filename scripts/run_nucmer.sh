#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Run nucmer pipeline to align BACs (param 1) to chr ref (param 2). NOT USING delta-filter -u 99 as diff BACs can contain same repeat and that does not make the alignment to ref chr invalid. Also BAC can map to chr00 and chr1-12 which is important to identify cases where contigs from chr00 can be integrated into chr1-12.

usage(){
	echo "usage:
	$0 <BAC fasta> <REF fasta> <word size> <min cluster size>
	Recommend 100 bp word size and 500 bp min cluster size. Default word size is 20 amd default min cluster size is 65"
	exit 1
}

if [ "$#" -ne 4 ]
then
	usage
fi

printf "Query : %s \n" "$1"
printf "Reference : %s \n" "$2"
printf "Word size : %s \n" "$3"
printf "Min cluster size : %s \n" "$4"
printf "Searching with word size of %i and min cluster size of %i\n\n" "$3" "$4"

nucmer -l "$3" -c "$4" --noextend -p "${3}"bp_qry_"${1}"__ref_"${2}" "$2" "$1"
#delta-filter -l 500 -u 99 500bp_qry_"${1}"__ref_"${2}".delta > 500bp_qry_"${1}"__ref_"${2}".delta.filtered
delta-filter -l "$3" "${3}"bp_qry_"${1}"__ref_"${2}".delta > "${3}"bp_qry_"${1}"__ref_"${2}".delta.filtered
show-coords -c -d -l -q -T -o "${3}"bp_qry_"${1}"__ref_"${2}".delta.filtered > "${3}"bp_qry_"${1}"__ref_"${2}".delta.filtered.coords

printf "\n\nOutput in %s \n" "${3}"bp_qry_"${1}"__ref_"${2}".delta.filtered.coords
