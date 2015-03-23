#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Find overlaps of filled gaps on ref chr (that were filled by BACs/assembled BACs) with ITAG2.4 genes
# Requires bedtools
# Please change Bio-GenomeUpdate directory to reflect your installation
# Please modify group_coords.pl output filenames according to your runs

BG_DIR="/home/surya/work/Eclipse/Bio-GenomeUpdate/"

printf "Stats for proper alignments\n"

ch=0
while [ $ch -le 12 ]
	do printf "Processing $ch\n"
	cd $ch

	${BG_DIR}scripts/grouped_coords_to_bed.sh 500bp.mixedoutoforder.agp.group_coords.stdout > 500bp.mixedoutoforder.agp.group_coords.stdout.bed
	
	if [ $ch -le 9 ]
	then
		bedtools intersect -a ch0${ch}.genes.ITAG2.4_gene_models.gff3 -b 500bp.mixedoutoforder.agp.group_coords.stdout.bed > 500bp.mixedoutoforder.agp.group_coords.stdout.bed.genes.gff3
	else
		bedtools intersect -a ch${ch}.genes.ITAG2.4_gene_models.gff3 -b 500bp.mixedoutoforder.agp.group_coords.stdout.bed > 500bp.mixedoutoforder.agp.group_coords.stdout.bed.genes.gff3

	fi

	#total region covered
	covered=`awk '{covered+=$3-$2} END {print covered}' 500bp.mixedoutoforder.agp.group_coords.stdout.bed`
	printf "Chr coverage: $covered bp\n"

	fixed=`wc -l 500bp.mixedoutoforder.agp.group_coords.stdout.bed.genes.gff3 | awk '{print $1}'`
	#if any genes were covered
	if [ "$fixed" -gt "0" ]
		then
		printf "Genes overlapping filled gaps: $fixed\n"
	else
		unlink 500bp.mixedoutoforder.agp.group_coords.stdout.bed.genes.gff3
	fi
	
	cd ..
	ch=$(($ch+1)); 
done

