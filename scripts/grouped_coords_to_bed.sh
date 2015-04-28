#!/bin/sh

# Surya Saha
# BTI/PPath@Cornell
# Purpose: Convert STDOUT from group_coords.pl to BED file
# http://bedtools.readthedocs.org/en/latest/content/general-usage.html

usage(){
	echo "usage:
	$0 <group_coords stdout>"
	exit 1
}

if [ "$#" -ne 1 ]
then
	usage
fi

# print file to STDERR
printf "group_coords file : %s \n" "$1" >&2

#skipping header line 
tail -n +2 "$1"| awk '{print $2"\t"$3-1"\t"$4"\t"$1}'
