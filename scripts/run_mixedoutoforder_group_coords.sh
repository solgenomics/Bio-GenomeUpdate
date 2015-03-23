#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Run group_coords.pl that creates separate files for valid coords, out of order coords and mixed coords for 0-12 chromosomes
# Please change Bio-GenomeUpdate directory to reflect your installation
# Please modify group_coords.pl output filenames according to your runs

BG_DIR="/home/surya/work/Eclipse/Bio-GenomeUpdate/"

ch=0
while [ $ch -le 12 ]
	do printf "Processing $ch\n"
	cd $ch
	
	if [ $ch -le 9 ]
		then
		perl -I ${BG_DIR}lib ${BG_DIR}group_coords.pl -i 500bp.ch0${ch}_asm_BACs__SL2.50ch0${ch}.delta.filtered.coords -g 100000 -u "dummy" -r SL2.50ch0${ch}.fa -q ch0${ch}_asm_BACs.fas -c noscafgaps.chr0${ch}_fish2_gaps.comp.agp -s chr0${ch}_fish2_gaps.chr.agp  > 500bp.mixedoutoforder.agp.group_coords.stdout 2> 500bp.mixedoutoforder.agp.group_coords.stderr
	
	else
		perl -I ${BG_DIR}lib ${BG_DIR}group_coords.pl -i 500bp.ch${ch}_asm_BACs__SL2.50ch${ch}.delta.filtered.coords -g 100000 -u "dummy" -r SL2.50ch${ch}.fa -q ch${ch}_asm_BACs.fas -c noscafgaps.chr${ch}_fish2_gaps.comp.agp -s chr${ch}_fish2_gaps.chr.agp  > 500bp.mixedoutoforder.agp.group_coords.stdout 2> 500bp.mixedoutoforder.agp.group_coords.stderr
	
	fi
	cd ..
	
	ch=$(($ch+1)); 
done
