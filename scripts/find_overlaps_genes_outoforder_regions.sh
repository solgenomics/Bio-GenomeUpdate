#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Find overlaps of ref chr regions ()that align out of order with BACs) with ITAG2.4 genes
# Requires bedtools
# Please change Bio-GenomeUpdate directory to reflect your installation
# Please modify group_coords.pl output filenames according to your runs

BG_DIR="/home/surya/work/Eclipse/Bio-GenomeUpdate/"

printf "Stats for out of order alignments\n"

ch=0
while [ "$ch" -le 12 ]
	do printf "Processing %s\n" "$ch"
	cd "$ch" || exit
	
	if [ "$ch" -le 9 ]
	then
		if [ -f noncolinear_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out ]
		then
			"${BG_DIR}"scripts/grouped_coords_to_bed.sh noncolinear_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out > noncolinear_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed
			bedtools intersect -a data/ch0"${ch}".genes.ITAG2.4_gene_models.gff3 -b noncolinear_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed > noncolinear_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3
			#total region covered
			covered=$(awk '{covered+=$3-$2} END {print covered}' noncolinear_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed)
			printf "Chr coverage: %s bp\n" "$covered"

			errors=$(wc -l noncolinear_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3 | awk '{print $1}')

			#if any genes were covered
			if [ "$errors" -gt "0" ]
				then
				printf "Genes overlapping noncolinear regions: %s\n" "$errors"
				else
				unlink noncolinear_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3
			fi
		fi
	else
		if [ -f noncolinear_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out ]
		then
			"${BG_DIR}"scripts/grouped_coords_to_bed.sh noncolinear_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out > noncolinear_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed
			bedtools intersect -a ch"${ch}".genes.ITAG2.4_gene_models.gff3 -b noncolinear_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed > noncolinear_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3

			#total region covered
			covered=$(awk '{covered+=$3-$2} END {print covered}' noncolinear_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed)
			printf "Chr coverage: %s bp\n" "$covered"

			errors=$(wc -l noncolinear_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3 | awk '{print $1}')

			#if any genes were covered
			if [ "$errors" -gt "0" ]
				then
				printf "Genes overlapping noncolinear regions: %s\n" "$errors"
				else
				unlink noncolinear_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3
			fi
		fi
	fi
	
	cd ..
	ch=$(("$ch"+1)); 
done

