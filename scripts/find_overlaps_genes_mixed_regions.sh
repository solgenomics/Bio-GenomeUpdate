#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Find overlaps of ref chr regions (that align in mixed orientation with BACs) with ITAG2.4 genes
# Requires bedtools
# Please change Bio-GenomeUpdate directory to reflect your installation
# Please modify group_coords.pl output filenames according to your runs

BG_DIR="/home/surya/source_code_repos/Bio-GenomeUpdate/"

printf "Stats for mixed orientation alignments\n"

ch=0
while [ "$ch" -le 12 ]
	do printf "Processing %s\n" "$ch"
	cd "$ch" || exit
	
	if [ "$ch" -le 9 ]
	then
		if [ -f mixed_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out ]
		then
			"${BG_DIR}"scripts/grouped_coords_to_bed.sh mixed_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out > mixed_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed
			bedtools intersect -a ch0"${ch}".genes.ITAG2.4_gene_models.gff3 -b mixed_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed > mixed_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3
			#total region covered
			covered=$(awk '{covered+=$3-$2} END {print covered}' mixed_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed)
			printf "Chr coverage: %s bp\n" "$covered"

			errors=$(wc -l mixed_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3 | awk '{print $1}')

			#if any genes were covered
			if  [ "$errors" -gt "0" ]
				then
				printf "Genes overlapping mixed regions: %s\n" "$errors"
				else
				unlink mixed_500bp.ch0"${ch}"_asm_BACs__SL2.50ch0"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3
			fi
		fi
	else
		if [ -f mixed_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out ]
		then
			"${BG_DIR}"scripts/grouped_coords_to_bed.sh mixed_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out > mixed_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed
			bedtools intersect -a ch"${ch}".genes.ITAG2.4_gene_models.gff3 -b mixed_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed > mixed_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3

			#total region covered
			covered=$(awk '{covered+=$3-$2} END {print covered}' mixed_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed)
			printf "Chr coverage: %s bp\n" "$covered"

			errors=$(wc -l mixed_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3 | awk '{print $1}')

			#if any genes were covered
			if  [ "$errors" -gt "0" ]
				then
				printf "Genes overlapping mixed regions: %s\n" "$errors"
				else
				unlink mixed_500bp.ch"${ch}"_asm_BACs__SL2.50ch"${ch}".delta.filtered.coords_group_coords.out.bed.genes.gff3
			fi

		fi
	fi
	
	cd ..
	ch=$(("$ch"+1)); 
done

