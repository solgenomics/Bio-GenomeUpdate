#!/bin/bash
#-------------------------------------------------------------------------------------------------------------------------------
# NAME
#
# auto_vcf_update.sh
#
# SYNOPSIS
# Shell script for updating VCF file coordinates. AGP files and new ref seq must be given as command line arguments.
# Script will automatically update all vcf files in the current directory unless glob at step 3 is modified.
# Step 5 and 6 can also be modified to change speed / CPU core usage. This script should be run from the
# scripts directory
#
# ./auto_vcf_update.sh -o [old AGP file] -n [new AGP file] -f [new ref seq FASTA file]
#
# To run, this script requires:
#                      samtools
#                      create_pseudo_gff_from_vcf.pl
#                      split_by_lines.sh
#                      update_coordinates_gff.pl
#                      copy_updated_coordinates_to_vcf.pl
#-------------------------------------------------------------------------------------------------------------------------------


function split_by_chromosome {
    # count number of unique chromosomes in file

    CHROMS=($(cut -f1 < $1 | uniq | grep -v '^#'))

    LABEL=($(echo ${CHROMS[0]} | sed -r 's/([0-9]+)$//'))

    NUMBERS=($(for chrom in ${CHROMS[*]}; do echo $chrom; done | sed -r 's/^(.*)([0-9][0-9])$/\2/' | grep -v '00'))

    for num in ${NUMBERS[*]}; do
	echo one chrom number = "$num"
    done

    for chrom in ${CHROMS[*]}; do
	grep $chrom $1 > $1.$chrom
    done

    echo "$1 split into ${#CHROMS[@]} individual chromosome files"
}

function combine_by_chromosome {
    #take split files with updated coordinates and recombine them in the right order

    EXTENSION="_mapped"
    FIRSTCHROM=($( head -n 2 $1 | awk '{ print $1}' ))
    echo FIRSTCHROM = ${FIRSTCHROM[1]}

    if [ "${FIRSTCHROM[1]}" == "${OLD_CHROM_LABEL}00" ]; then

	cat temp_update_files/$1.${OLD_CHROM_LABEL}00 >> temp_update_files/$1.${NEW_CHROM_LABEL}

	for num in ${NUMBERS[*]}; do
            NUMFILES=($(ls temp_update_files/$1.${OLD_CHROM_LABEL}$num* | wc -l))
            echo NUMFILES = $NUMFILES
            NUMFILES=$((NUMFILES / 2))
            echo now NUMFILES = $NUMFILES
            for ((i=1; i<=$NUMFILES; i+=1)); do
		cat temp_update_files/$1.${OLD_CHROM_LABEL}$num.$i$EXTENSION >> temp_update_files/$1.${NEW_CHROM_LABEL}
            done
	done

    else

	for num in ${NUMBERS[*]}; do
            NUMFILES=($(ls temp_update_files/$1.${OLD_CHROM_LABEL}$num* | wc -l))
            echo NUMFILES = $NUMFILES
            NUMFILES=$((NUMFILES / 2))
            echo now NUMFILES = $NUMFILES
            for ((i=1; i<=$NUMFILES; i+=1)); do
		cat temp_update_files/$1.${OLD_CHROM_LABEL}$num.$i$EXTENSION >> temp_update_files/$1.${NEW_CHROM_LABEL}
            done
	done

	cat temp_update_files/$1.${OLD_CHROM_LABEL}00 >> temp_update_files/$1.${NEW_CHROM_LABEL}
    fi
}

#--------------------------------------------------------------------------------
# 1 Parse command line arguments:
#-------------------------------------------------------------------------------

while [[ $# > 1 ]]
do
key="$1"

case $key in
  -f|--fasta)
  FASTA="$2"
  shift
  ;;
  -o|--old.agp)
  OLD_AGP="$2"
  shift
  ;;
  -n|--new.agp)
  NEW_AGP="$2"
  ;;
esac
shift
done
echo FASTA  = "$FASTA"
echo OLD_AGP  = "$OLD_AGP"
echo NEW_AGP  = "$NEW_AGP"
if [ -z "$FASTA" ] || [ -z "$OLD_AGP" ] || [ -z "$NEW_AGP" ]
then
    echo "Trouble reading command line arguments, make sure
    -f [new ref seq FASTA] and
    -o [old AGP] and
    -n [new AGP] are all specified";
    exit
fi

#----------------------------------------------------------------------------------
# 2 split AGP files by chromosome. Store chrom names in shell variables
#----------------------------------------------------------------------------------

echo Splitting agp files...

split_by_chromosome ${OLD_AGP}
OLD_CHROM_LABEL=$LABEL
echo OLD_CHROM_LABEL = "$OLD_CHROM_LABEL"

split_by_chromosome ${NEW_AGP}
NEW_CHROM_LABEL=$LABEL
echo NEW_CHROM_LABEL = "$NEW_CHROM_LABEL"

#--------------------------------------------------------------------------------
# 3 create pseudo gff files from each vcf file in parallel
#-------------------------------------------------------------------------------

ls *.vcf | parallel ./create_pseudo_gff_from_vcf.pl


#--------------------------------------------------------------------------------------------------
# 4 make directories for temporary and final output
# 5 split pseudo gff files by chromosome and line,   *split_by_lines script can be modified
#                                                     to determine line cutoff. Default is 2000
#
# 6 then update their coordinates in parallel        * parallel command is restricted to use 12
#                                                      CPUs as default. Feel free to change.
#--------------------------------------------------------------------------------------------------

mkdir updated_vcfs

echo Output directory updated_vcfs created
echo Starting parallelized coordinate update tasks

for file in *pseudo.gff* ; do

mkdir temp_update_files

./split_by_lines.sh $file | parallel -j 60 --verbose --rpl '{..} s/^.*([0-9][0-9])\.[0-9]+$/\1/;' "../update_coordinates_gff.pl -o ${OLD_AGP}.${OLD_CHROM_LABEL}{..} -n ${NEW_AGP}.${NEW_CHROM_LABEL}{..} -g {} -m {}_mapped -c 0" ;

#------------------------------------------------------------------------------------------------
# 7 compile  all split chromosome files for current accession back into a single updated gff file
#------------------------------------------------------------------------------------------------

combine_by_chromosome $file ;

#-----------------------------------------------------------------------------------------------
# 8 take updated coordinates and info from original vcf file to create updated vcf file
#-----------------------------------------------------------------------------------------------

samtools faidx $FASTA;
./copy_updated_coordinates_to_vcf.pl temp_update_files/$file.$NEW_CHROM_LABEL $file $FASTA $NEW_CHROM_LABEL;

#-----------------------------------------------------------------------------------------------
# 9 remove files temporarily created during update process
#------------------------------------------------------------------------------------------------

rm -r temp_update_files
rm $file

done

#--------------------------------------------------------------------------------------------------
# 10 report successful update
#--------------------------------------------------------------------------------------------------

echo Removing temp files...

function finish {
    rm *.agp.$OLD_CHROM_LABEL*
    rm *.agp.$NEW_CHROM_LABEL*
}
trap finish EXIT

echo Finished. Updated vcf files stored in updated_vcfs directory
