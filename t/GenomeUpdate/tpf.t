#!/usr/bin/perl

=head1 NAME

genome_update.t
A test for Bio::GenomeUpdate::TPF class

=cut

=head1 SYNOPSIS

perl tpf.t


=head1 DESCRIPTION



=head2 Author

Jeremy Edwards <jde22@cornell.edu>
=cut

use strict;
use warnings;
use autodie;


use Test::More tests => 3;
BEGIN {use_ok( 'Bio::GenomeUpdate::TPF' ); }
require_ok( 'Bio::GenomeUpdate::TPF::TPFSequenceLine' );
require_ok( 'Bio::GenomeUpdate::TPF::TPFGapLine' );
