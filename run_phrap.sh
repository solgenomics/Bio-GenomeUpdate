#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Run phrap assemblies (def, -shatter_greedy -penalty -4 -minmatch 20, bandwidth 0) on a chr, phrap should be in $PATH

set -u #exit if uninit var
set -e #exit if non-zero return value (error), use command || {echo 'command failed'; exit 1;}
set -o nounset
set -o errexit

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(readlink -m $(dirname $0))

usage() {
    echo "usage:
    $PROGNAME <chr name> <fasta file>
    
    Example:
    $PROGNAME 0 0.fas"
    exit 1
}

if [ "$#" -ne 2 ]
then
	usage
fi

#default
mkdir "${1}_def"
cp "$2" "${1}_def"
cd "${1}_def"
phrap -new_ace "$2" 1> phrap.out 2> phrap.stderr

for F in ${1}.fas.*
	do
		pigz $F
	done

cd $PROGDIR

#mix (-shatter_greedy -penalty -4 -minmatch 20)
mkdir "${1}_mix"
cp "$2" "${1}_mix"
cd "${1}_mix"
phrap -new_ace -shatter_greedy -penalty -4 -minmatch 20 "$2" 1> phrap.out 2> phrap.stderr

for F in ${1}.fas.*
	do
		pigz $F
	done

cd $PROGDIR

#gapless (-bandwidth 0)
mkdir "${1}_bw0"
cp "$2" "${1}_bw0"
cd "${1}_bw0"
phrap -new_ace -bandwidth 0 "$2" 1> phrap.out 2> phrap.stderr

for F in ${1}.fas.*
	do
		pigz $F
	done

cd $PROGDIR

