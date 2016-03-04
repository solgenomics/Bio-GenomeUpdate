#!/usr/bin/perl

=head1 NAME

agp_to_bacs_fasta.pl

=head1 SYNOPSIS

agp_to_bacs_fasta.pl -a [ AGP file] -f [Fasta of BAC accessions]

=head1 COMMAND-LINE OPTIONS

 -a  AGP file from GRC output (required)
 -f  Fasta file of BAC accessions. The sequence names should be in Genbank format gi|118344470|gb|AC171726.2| (required)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;

use Getopt::Std;

our ( $opt_f, $opt_a, $opt_d, $opt_h );
getopts('f:a:d:h');
if ($opt_h) {
	help();
	exit;
}

if ( !$opt_a ) {
	print "\n AGP file from GRC output file is required. See help below\n\n\n";
	help();
}
if ( !$opt_f ) {
	print "\n Fasta file of BAC accessions is required. See help below\n\n\n";
	help();
}




#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

    This script will create a pseudomolecule of only the sequences in the Fasta file. This is helpful for creating a dummy fasta file with selected accessions for doing dot-plots and other analyses.
     
    Usage:
      agp_to_bacs_fasta.pl -a [ AGP file] -f [Fasta of BAC accessions]
      
    Flags:

          -a  AGP file from GRC output (required)
          -f  Fasta file of BAC accessions. The sequence names should be in Genbank format gi|118344470|gb|AC171726.2| (required)
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

