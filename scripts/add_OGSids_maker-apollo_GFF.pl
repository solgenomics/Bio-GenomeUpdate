#!/usr/bin/perl

=head1 NAME

add_OGSids_maker-apollo_GFF.pl

=head1 SYNOPSIS

add_OGSids_maker-apollo_GFF.pl -i [old Merged Maker and Apollo GFF file] -p [species acronym or prefix] -c [Prefix for chromosome in GFF]  -s [starting value for naming] -g [gff file with chromosomal location]

=head1 COMMAND-LINE OPTIONS

 -i  Merged Maker and Apollo GFF file (required)
 -p  Prefix for name, e.g DcitrP (required)
 -c  Prefix for chromosome in GFF e.g. Dc3.0sc (required)
 -s  Starting seed, e.g. 1 (required)
 -g  GFF file with chromosomal location
 -h  Help

=cut

use strict;
use warnings;







#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     ??????????????
     
    NOTE:


    Usage:
      cmdline_perldoc.pl -g [?? file] <other params>
      
    Flags:

		 -g  old GFF file (required)
		 -d  debugging messages (1 or 0)
		 -h  Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl. Change??????????

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

