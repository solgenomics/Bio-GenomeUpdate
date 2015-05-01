#!/usr/bin/perl

=head1 NAME
copy_updated_coordinates_to_vcf.pl
=head1 SYNOPSIS
copy_updated_coordinates_to_vcf.pl [updated pseudo GFF file] [old VCF file]
=cut

use warnings;
use strict;

my $GFF = $ARGV[0];
my $VCF = $ARGV[1];
my $chrom_label = $ARGV[2];

my $build_name = $chrom_label;
$build_name =~ s/(.*)ch$/$1/;

$VCF =~ s/.pseudo.gff(.*)/.vcf$1/;

my $new_vcf = $VCF;
$new_vcf =~ s/(.*)/updated_vcfs\/$1/;

open(NEWCOORD, "<", $GFF) || die "Can't open gff file $GFF"; 

open(OLDVCF, "<", $VCF) || die "Can't open vcf file $VCF"; 

open(my $NEWVCF, ">", $new_vcf) || die "Can't create new VCF file $new_vcf";

my $data_line_counter = 0;

for (<OLDVCF>) {
    chomp;
    if (m/^##/) {
	print $NEWVCF $_ . "\n";
    } elsif (m/^#/) {
	print $NEWVCF "##mapped to $build_name " . localtime . "\n";
	print $NEWVCF $_ . "\n";
    } else {
	my ($chrom, $pos, @extra) = split /\t/;
	chomp (my $new_coord = <NEWCOORD>);
       	my ($ignore1, $ignore2, $ignore3, $new_pos, @ignore) = split /\t/, $new_coord;
	$pos = $new_pos;
	$chrom =~ s/^(.*)([0-9][0-9])$/$chrom_label$2/;
	print $NEWVCF join "\t", ($chrom, $new_pos, @extra); 
	print $NEWVCF "\n";
    }
}

# sort new vcf file to ensure newly mapped coordinates are in order

system ("grep '^#' $new_vcf > $new_vcf.sorted");
system ("grep -v '^#' $new_vcf | sort -V >> $new_vcf.sorted");
system ("mv $new_vcf.sorted $new_vcf");

print "New vcf file $new_vcf with updated coordinates successfully created.\n";
   















