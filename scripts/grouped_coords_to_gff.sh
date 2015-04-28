#!/bin/sh

# Surya Saha
# BTI/PPath@Cornell
# Purpose: Convert STDOUT from group_coords.pl to GFF file

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

