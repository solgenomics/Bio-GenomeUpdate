#!/usr/bin/perl

=head1 NAME

add_OGSids_maker-apollo_GFF.pl

=head1 SYNOPSIS

add_OGSids_maker-apollo_GFF.pl -g [old Merged Maker and Apollo GFF file] -a [AHRD file for mRNA with Maker and Apollo ids] -p [species acronym or prefix] -c [Prefix for chromosome in GFF]  -s [starting value for naming] -o [output formatted gff file with OGS ids]

=head1 COMMAND-LINE OPTIONS

 -g  Merged Maker and Apollo GFF file (required)
 -a  AHRD file for mRNA with Maker and Apollo ids (required)
 -p  Prefix for name, e.g DcitrP (required)
 -c  Prefix for chromosome in GFF e.g. Dc3.0sc (required)
 -s  Starting seed, e.g. 1 (required)
 -o  output GFF file with OGS ids
 -h  Help

=cut

use strict;
use warnings;

use File::Slurp;
use Getopt::Std;

our ( $opt_g, $opt_a, $opt_p, $opt_s, $opt_c, $opt_o, $opt_h );
getopts('g:a:p:s:c:g:h');
if ($opt_h) {
  help();
  exit;
}
if ( !$opt_g || !$opt_a || !$opt_p || !$opt_c || !$opt_s || !$opt_o) {
  print
"\nOld GFF file, AHRD file, name prefix, chr prefix, starting seed, output GFF file is required.
See help below\n\n\n";
  help();
}





#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     Renames maker and Apollo assigned gene names to gene name with version numbers, e.g.  maker-ScVcwli_1-pred_gff_maker-gene-0.0-mRNA-1 becomes DcitrP00001.1.1. Creates an index file with old and new mRNA ids. This is hard coded for <1,000,000 mRNAs. Counter skips over 10 gene models so manually curated genes can be added.
     
    NOTE:


    Usage:
      cmdline_perldoc.pl -g [?? file] <other params>
      
    Flags:

		  -g  Merged Maker and Apollo GFF file (required)
      -a  AHRD file for mRNA with Maker and Apollo ids (required)
      -p  Prefix for name, e.g DcitrP (required)
      -c  Prefix for chromosome in GFF e.g. Dc3.0sc (required)
      -s  Starting seed, e.g. 1 (required)
      -o  output GFF file with OGS ids
      -h  Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

