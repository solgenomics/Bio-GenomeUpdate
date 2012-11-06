#!/usr/bin/perl

=head1 NAME

genome_update.t
A test for Bio::GenomeUpdate class

=cut

=head1 SYNOPSIS

perl genome_update.t


=head1 DESCRIPTION



=head2 Author

Jeremy Edwards <jde22@cornell.edu>
=cut

use strict;
use warnings;
use autodie;


use Test::More tests => 5;
BEGIN {use_ok( 'Bio::GenomeUpdate' ); }
require_ok( 'Bio::GenomeUpdate::TPF' );
require_ok( 'Bio::GenomeUpdate::AGP' );
require_ok( 'Bio::GenomeUpdate::AlignmentCoords' );
require_ok( 'Bio::GenomeUpdate::AlignmentCoordsGroup' );


#ok (my $pedigree = Bio::GeneticRelationships::Pedigree->new());
#is ($pedigree->get_name(),'test_name');
