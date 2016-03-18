#!/usr/bin/perl

=head1 NAME

bacs_blast_gff_multi_component.pl

=head1 SYNOPSIS

bacs_blast_gff_multi_component.pl -a [ AGP file] -f [Fasta of BAC accessions]

=head1 COMMAND-LINE OPTIONS

 -g  GFF file from blast output (required, use https://github.com/suryasaha/Utils/blob/master/Blast2HitGFF3.pl)
 -a  AGP file from GRC output (required)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;

use Getopt::Std;

our ( $opt_g, $opt_a, $opt_d, $opt_h );
getopts('f:a:d:h');
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

