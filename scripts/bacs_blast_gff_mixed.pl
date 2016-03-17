#!/usr/bin/perl

=head1 NAME

bacs_blast_gff_mixed.pl

=head1 SYNOPSIS

bacs_blast_gff_mixed.pl -g [ GFF file] -l [ 1000 ]

=head1 COMMAND-LINE OPTIONS

 -g  GFF file from blast output (required, use Blast2HitGFF3.pl)
 -l  minimum length of feature or HSP (default 1000bp)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;

use Getopt::Std;

our ( $opt_g, $opt_l, $opt_d, $opt_h );
getopts('g:l:d:h');
if ($opt_h) {
	help();
	exit;
}

if ( !$opt_g ) {
	print "\n GFF3 of BAC alignments to reference is required. Use blastn and Blast2HitGFF3.pl. See help below\n\n\n";
	help();
}




#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

    This script will print a list of BACs that align in mixed orientation to the reference. These regions could be errors in the reference assembly. This does NOT check if the alignments are OUT OF ORDER.
     
    Usage:
      bacs_blast_gff_mixed.pl -g [ GFF file] -l [ 1000 ]
      
    Flags:

          -g  GFF file from blast output (required, use Blast2HitGFF3.pl)
          -l  minimum length of feature or HSP (default 1000bp)
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

