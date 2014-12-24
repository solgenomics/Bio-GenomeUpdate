#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Run phrap assemblies (def, -shatter_greedy -penalty -4 -minmatch 20, bandwidth 0) on a chr

set -u #exit if uninit var
set -e #exit if non-zero return value (error), use command || {echo 'command failed'; exit 1;}
set -o nounset
set -o errexit

readonly PROGNAME=$(basename $0)
readonly PROGDIR=$(readlink -m $(dirname $0))
readonly ARGS="$@"
readonly WDIR=`pwd`

usage() {
    echo "usage:
    $PROGNAME <chr num>
    
    Example:
    $PROGNAME 0"
    exit 1
}

if [ "$#" -ne 1 ]
then
	usage
fi

#default
mkdir "${1}_def"
ln -s $WDIR/data/${1}.fas "${1}_def"
cd "${1}_def"
phrap -new_ace ${1}.fas 1> phrap.out 2> phrap.stderr

for F in ${1}.fas.*
	do
	pigz $F
	done

cd $WDIR

#mix (-shatter_greedy -penalty -4 -minmatch 20)
mkdir "${1}_mix"
ln -s $WDIR/data/${1}.fas "${1}_mix"
cd "${1}_mix"
phrap -new_ace -shatter_greedy -penalty -4 -minmatch 20 ${1}.fas 1> phrap.out 2> phrap.stderr

for F in ${1}.fas.*
	do
	pigz $F
	done

cd $WDIR

#gapless (-bandwidth 0)
mkdir "${1}_bw0"
ln -s $WDIR/data/${1}.fas "${1}_bw0"
cd "${1}_bw0"
phrap -new_ace -bandwidth 0 ${1}.fas 1> phrap.out 2> phrap.stderr
for F in ${1}.fas.*
	do
	pigz $F
	done
cd $WDIR

