#!/usr/bin/perl
=head1

filter_delta.pl

=head1 SYNOPSIS

    filter_delta.pl -c [coords file] -d [delta file] -o [delta file]

=head1 COMMAND-LINE OPTIONS

 -c  COORDS file created by show-coords (required)
 -d  DELTA file created by nucmer (required) 
 -o  Output filtered DELTA file 
 -h  Help

=head1 TODO

  Print delta output file

=cut

use strict;
use warnings;

use Getopt::Std;
use File::Slurp;
use Bio::GenomeUpdate::AlignmentCoords;
use Bio::GenomeUpdate::AlignmentCoordsGroup;











sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script filters the DELTA file produced by nucmer by removing BAC alignments that align in mixed orientation or in non co-linear order.

    Usage:
      filter_delta.pl -c [coords file] -d [delta file] -o [delta file]

    Flags:

      -c  COORDS file created by show-coords (required)
      -d  DELTA file created by nucmer (required) 
      -o  Output filtered DELTA file 
      -h  Help
 
EOF
exit (1);
}

__END__