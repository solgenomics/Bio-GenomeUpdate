#!/bin/sh

# Surya Saha
# BTI/PPath@Cornell
# Purpose: Make dot plots for full chr length alignment with all and 1-1 mapping between ref and query
# Use with v4 mummer https://github.com/mummer4/mummer/releases
# All matches are 100% identical and only ATGC regions are reported

usage(){
	echo "usage:
	$0 <ref.fa> <query.fa> <prefix for outfiles> <min tile length>"
	exit 1
}

if [ "$#" -ne 4 ]
then
	usage
fi

printf "Ref file: %s \n" "$1"
printf "Query file: %s \n" "$2"
printf "Prefix: %s \n" "$3"
printf "Min tile length : %d \n" "$4"

# unique matches
mummer -mum -l "$4" -b -c -n -qthreads 6 "$1" "$2" > "$3".mummer.uniq.both.out
mummerplot --png --prefix "$3".uniq  --title "$3".uniq "$3".mummer.uniq.both.out

# all matches
mummer -maxmatch -l "$4" -b -c -n -threads 6 -qthreads 6 "$1" "$2" > "$3".mummer.all.both.out
mummerplot --png --prefix "$3".all --title "$3".all "$3".mummer.all.both.out

#cleanup
/bin/rm -f "${3}".uniq.rplot
/bin/rm -f "${3}".uniq.gp
/bin/rm -f "${3}".uniq.fplot
/bin/rm -f "${3}".all.rplot
/bin/rm -f "${3}".all.gp
/bin/rm -f "${3}".all.fplot

