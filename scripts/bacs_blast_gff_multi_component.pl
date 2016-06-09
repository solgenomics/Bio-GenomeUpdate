#!/usr/bin/perl

=head1 NAME

bacs_blast_gff_multi_component.pl

=head1 SYNOPSIS

bacs_blast_gff_multi_component.pl -a [ AGP file] -f [Fasta of BAC accessions]

=head1 COMMAND-LINE OPTIONS

 -g  GFF file from blast output with one record per BAC (required, use https://github.com/suryasaha/Utils/blob/master/Blast2HitGFF3.pl, Mirella's script)
 -a  AGP file from GRC output (required)
 -l  minimum length of feature or HSP (default 5000bp)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;

use Getopt::Std;

our ( $opt_g, $opt_a, $opt_l, $opt_d, $opt_h );
getopts('g:a:l:d:h');
if ($opt_h) {
	help();
	exit;
}

if ( !$opt_a ) {
	print "\n AGP file from GRC output file is required. See help below\n\n\n";
	help();
}
if ( !$opt_g ) {
	print "\n GFF file from blast output of BACs to reference is required. See help below\n\n\n";
	help();
}
my $minimum_length = defined $opt_l ? $opt_l : 5000;
unless (open (GFF, "<${opt_g}")) {die "Not able to open ${opt_g}\n\n";}
unless (open (AGP, "<${opt_a}")) {die "Not able to open ${opt_a}\n\n";}

CHR01_FISH2_GAPS        1       3610    1       W       AEKE02013047.1  1       3610    -
CHR01_FISH2_GAPS        3611    3863    2       N       253     clone   yes
CHR01_FISH2_GAPS        3864    12848   3       W       AEKE02013046.1  1       8985    -

while (my $rec = <AGP>){

}

while (my $rec = <GFF>){
	if( $rec =~ /^#/ ){ next;}
	
	my @rec_arr = split ("\t", $rec );

	
}


#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

    This script will print a list of BACs that cover multiple components in the AGP file. These BACs can be used to cover gaps in the reference genome.
     
    Usage:
      agp_to_bacs_fasta.pl -a [ AGP file] -f [Fasta of BAC accessions]
      
    Flags:

           -g  GFF file from blast output (required, use https://github.com/suryasaha/Utils/blob/master/Blast2HitGFF3.pl)
           -a  AGP file from GRC output (required)
           -d  debugging messages (1 or 0)
           -h  Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha at cornell.edu , @SahaSurya>

=cut

