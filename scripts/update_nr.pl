#!/usr/bin/perl -w
# Solgenomics@BTI
# Surya Saha 06/10/2016 
# Purpose: For updating coords for Neelima's BIL's (from NR files)


unless (@ARGV == 2){
	print "USAGE: $0 <updated gff> <nr file> \n";
	exit;
}

use strict;
use warnings;

my ($gfffname,$nrfname);

$gfffname=$ARGV[0];
$nrfname=$ARGV[1];

unless(open(INGFF,$gfffname)){print "not able to open ".$gfffname."\n\n";exit;}
unless(open(INNR,$nrfname)){print "not able to open ".$nrfname."\n\n";exit;}

unless(open(OUT,">updated.${nrfname}")){print "not able to open updated".$nrfname."results\n\n";exit;}

my $nr_rec=<INNR>;
print OUT $nr_rec;

while (my $gff_rec=<INGFF>){
	my @gff_rec_arr = split("\t", $gff_rec);
	$nr_rec=<INNR>;
	my @nr_rec_arr = split("\t", $nr_rec);
	
	if ($gff_rec_arr[6] eq '+'){
		print OUT $gff_rec_arr[0],"\t",$gff_rec_arr[4],"\t",$nr_rec_arr[2],"\t",$nr_rec_arr[3],"\t",$nr_rec_arr[4],"\t",$nr_rec_arr[5],"\t",$nr_rec_arr[6];
	}
	elsif ($gff_rec_arr[6] eq '-'){
		#flipping both REF and ALT allele
		$nr_rec_arr[2]=~ tr/ATGC/TACG/;
		$nr_rec_arr[3]=~ tr/ATGC/TACG/;
		print OUT $gff_rec_arr[0],"\t",$gff_rec_arr[4],"\t",$nr_rec_arr[2],"\t",$nr_rec_arr[3],"\t",$nr_rec_arr[4],"\t",$nr_rec_arr[5],"\t",$nr_rec_arr[6];
	}

}

close (INGFF);
close (INNR);
close (OUT);

