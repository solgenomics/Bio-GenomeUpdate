#!/usr/bin/perl

=head1 NAME
create_pseudo_gff_from_vcf.pl
=head1 SYNOPSIS
create_pseudo_gff_from_vcf.pl [VCF file to be processed]
=cut

use warnings;                                                                                      
use strict;

my $file  = shift;

open(VCF, "<", $file) || die "Can't open file $file";

my $pseudo_gff = $file;
$pseudo_gff =~ s/.vcf/.pseudo.gff/;

open(my $G, ">", $pseudo_gff) || die "Can't open gff file $pseudo_gff";

#print $G "##gff-version 3 pseudo\n";

my $INC = 00000;

for (<VCF>) {
    shift;
    if ($_ !~ m/^#/) {
        chomp;
	$INC++;
        my ($CHROM, $POS, @extra) = split /\t/;
        print $G "$CHROM\t.\texon\t$POS\t$POS\t.\t+\t.\tID=ex$INC\n";
    } else {
    }
}
print "$pseudo_gff successfully created.\n";
