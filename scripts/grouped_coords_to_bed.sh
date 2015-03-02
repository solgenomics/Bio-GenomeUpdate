#!/bin/sh

# Surya Saha
# BTI/PPath@Cornell
# Purpose: Convert STDOUT from group_coords.pl to BED file
# http://bedtools.readthedocs.org/en/latest/content/general-usage.html

#skipping header line 
tail -n +2 "$1"| awk '{print $2"\t"$3-1"\t"$4"\t"$1}'
