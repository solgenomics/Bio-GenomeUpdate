#!/bin/sh

# Surya Saha
# PPath@Cornell/BTI
# Purpose: Run nucmer pipeline to align BACs (param 1) to chr ref (param 2). NOT USING delta-filter -u 99 as diff BACs can contain same repeat and that does not make the alignment to ref chr invalid. 

nucmer -l 100 -c 500 --noextend -p 500bp_qry_${1}__ref_${2} "$2" "$1"
#delta-filter -l 500 -u 99 500bp_qry_${1}__ref_${2}.delta > 500bp_qry_${1}__ref_${2}.delta.filtered
#delta-filter -l 500 -u 99 500bp_qry_${1}__ref_${2}.delta > 500bp_qry_${1}__ref_${2}.delta.filtered
delta-filter -l 500 500bp_qry_${1}__ref_${2}.delta > 500bp_qry_${1}__ref_${2}.delta.filtered
show-coords -c -d -l -q -T -o 500bp_qry_${1}__ref_${2}.delta.filtered > 500bp_qry_${1}__ref_${2}.delta.filtered.coords

#perl -I /lib group_coords.pl -i 500bp_qry_${1}__ref_${2}.delta.filtered.coords -g 100000 -u "dummy" -f "$2" -q "$1"  > 500bp.group_coords.out

