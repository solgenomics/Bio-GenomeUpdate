#!/bin/bash

#-------------------------------------------------
# Shell commands for updating VCF file coordinates
#-------------------------------------------------

# split old.agp file by chromosome

./split_by_chromosome.sh  S_lycopersicum_chromosomes_from_scaffolds.2.40.agp 

# split new.agp file by chromosome

./split_by_chromosome.sh  SL2.50ch_from_sc.agp 

# create pseudo gff files from each vcf file in parallel

parallel ./create_pseudo_gff_from_vcf.pl TS-{}* ::: {1..360}

# split pseudo gff files by chromosome in parallel, then use these files along with old and new agp files to create new pseudo gffs with updated coordinates (in parallel)

ls *.pseudo.gff | parallel -I @@ -j 10 './split_by_chromosome.sh @@; parallel -X ./update_coordinates_gff.pl -o S_lycopersicum_chromosomes_from_scaffolds.2.40.agp.SL2.40ch{} -n SL2.50ch_from_sc.agp.SL2.50ch{} -g @@.SL2.40ch{} -m @@.SL2.50ch{} -c 0 ::: {01..12}'

# find all individual updated pseudo gff chromosome files for each specific accession and combine them 

./combine_updated_gffs.sh

# take full gff file with updated coordinates and old vcf file; use them to create new vcf file with updated coordinates 

mkdir updated_vcfs

parallel ./copy_updated_coordinates_to_vcf.pl RF_{}*.pseudo.gff.SL2.50 RF_{}*.vcf ::: {001..105} 

# remove all files temporarily created during process

parallel rm *.SL2.40ch{} ::: {00..12}
parallel rm *.SL2.50ch{} ::: {00..12}
rm *.pseudo.gff.SL2.50
rm *.pseudo.gff
