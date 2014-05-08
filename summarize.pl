#!/usr/bin/perl

=head1 NAME

summarize.pl

=head1 SYNOPSIS

summarize.pl -g [old GFF file] -f [old Fasta file] -u [updated GFF file] -n [new Fasta file] <other params>

=head1 COMMAND-LINE OPTIONS

 -g  old GFF file (required)
 -f  old Fasta file based on old scaffold AGP (required)
 -u  Updated GFF file output by update_coordinates_gff.pl (required)
 -n  New Fasta file based on new scaffold AGP (required)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;
