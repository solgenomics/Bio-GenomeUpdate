#!/usr/bin/perl

=head1 NAME

genome_update.t
A test for Bio::GenomeUpdate::AGP class

=cut

=head1 SYNOPSIS

perl agp.t


=head1 DESCRIPTION



=head2 Author

Jeremy Edwards <jde22@cornell.edu>
=cut

use strict;
use warnings;
use autodie;


use Test::More tests => 4;
BEGIN {use_ok( 'Bio::GenomeUpdate::AGP' ); }
require_ok( 'Bio::GenomeUpdate::AGP::AGPSequenceLine' );
require_ok( 'Bio::GenomeUpdate::AGP::AGPGapLine' );
require_ok( 'Bio::GenomeUpdate::AGP::AGPConvert' );
