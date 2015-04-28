#!/bin/bash

# count number of unique chromosomes in file
CHROMS=($(cut -f1 $1 | uniq | grep -v "^#"))

# split file into separate files for each individual chromosome
for chrom in ${CHROMS[*]}; do
   # echo chrom = $chrom
    grep $chrom $1 > temp_update_files/$1.$chrom
    LINES=($(wc -l temp_update_files/$1.$chrom))
   # echo LINES = $LINES
    SPLITNUMBER=$(($LINES/2000))
    ((SPLITNUMBER+=1)) 
   # echo SPLITNUMBER = $SPLITNUMBER
# split large individual chromsome files into smaller files in order to prevent bottleneck at ./update_coordinates_gff.pl step 
    if [[ "$chrom" == *00 ]]; then
            :
            else
            for ((f=1; f<=$SPLITNUMBER; f++)); do
		head -n 2000 temp_update_files/$1.$chrom > temp_update_files/$1.$chrom.$f
		echo temp_update_files/$1.$chrom.$f
		if [ "$f" -lt  "$SPLITNUMBER" ]; then
	            tail -n +2001 temp_update_files/$1.$chrom > temp_update_files/$1.$chrom.temp
		    rm temp_update_files/$1.$chrom
		    mv temp_update_files/$1.$chrom.temp temp_update_files/$1.$chrom
		else 
	            rm temp_update_files/$1.$chrom
		fi
	    done
    fi
done

