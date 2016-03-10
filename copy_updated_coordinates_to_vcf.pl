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

open(MAPPEDGFF, "<", $GFF) || die "Can't open gff file $GFF"; 

open(OLDVCF, "<", $VCF) || die "Can't open vcf file $VCF"; 

open(my $NEWVCF, ">", $new_vcf) || die "Can't create new VCF file $new_vcf";

# for troubleshooting
#open(my $FLIPLOG, '>', 'flip_log.txt') || die "Can't create flip_log.txt: $!";

my $data_line_counter = 0;

for (<OLDVCF>) {
    chomp;
    if (m/^##/) {
	print $NEWVCF $_ . "\n";
    }
    elsif (m/^#/) {
	print $NEWVCF "##mapped to $build_name " . localtime . "\n";
	print $NEWVCF $_ . "\n";
    }
    else {
	my ($chrom, $position, $id, $ref, $alt, @extra) = split /\t/;
	chomp (my $gff_string = <MAPPEDGFF>);
       	my @values = split /\t/, $gff_string;
	my $new_position = $values[3];
	my $new_orientation = $values[6];
	$chrom =~ s/^(.*)([0-9][0-9])$/$chrom_label$2/;
	
	# before printing, correct ref and alt alleles if scaffold has been flipped in new genome build
	if ($new_orientation eq '-') {
	    #print $FLIPLOG join "\t", ($chrom, $new_orientation, $new_position);
	    #print $FLIPLOG "\n";
	    my $new_ref = &replace_with_complementary_base($ref);
	    my $new_alt = &replace_with_complementary_base($alt);
	    print $NEWVCF join "\t", ($chrom, $new_position, $id, $new_ref, $new_alt, @extra);
	    print $NEWVCF "\n";
	}
	else {
	    print $NEWVCF join "\t", ($chrom, $new_position, $id, $ref, $alt, @extra); 
	    print $NEWVCF "\n";
	}
    }
}

# sort new vcf file to ensure newly mapped coordinates are in order

system ("grep '^#' $new_vcf > $new_vcf.sorted");
system ("grep -v '^#' $new_vcf | sort -V >> $new_vcf.sorted");
system ("mv $new_vcf.sorted $new_vcf");

print "New vcf file $new_vcf with updated coordinates successfully created.\n";
   
sub replace_with_complementary_base {
    if ($_[0] eq 'A') {
	return 'T';
    }
    elsif ($_[0] eq 'T') {
	return 'A';
    }
    elsif ($_[0] eq 'C') {
	return 'G';
    }
    elsif ($_[0] eq 'G') {
	return 'C';
    }
}















