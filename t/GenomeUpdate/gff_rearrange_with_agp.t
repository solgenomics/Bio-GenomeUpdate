#!/usr/bin/perl

=head1 NAME

gff_rearrange_with_agp.t


=cut

=head1 SYNOPSIS

perl gff_rearrange_with_agp.t


=head1 DESCRIPTION



=head2 Author

Jeremy Edwards <jde22@cornell.edu>
=cut

use strict;
use warnings;
use autodie;

use Test::More tests => 10;
BEGIN {use_ok( 'Bio::GenomeUpdate::AGP' ); }
require_ok( 'Bio::GenomeUpdate::AGP::AGPSequenceLine' );
require_ok( 'Bio::GenomeUpdate::AGP::AGPGapLine' );

#create AGP and add lines and gaps
ok(my $agp = Bio::GenomeUpdate::AGP->new());
ok(my $agp_sequence_line = Bio::GenomeUpdate::AGP::AGPSequenceLine->new());
ok(my $agp_gap_line = Bio::GenomeUpdate::AGP::AGPGapLine->new());

#my $data_str = <DATA>;


my $old_agp_str = q(# ORGANISM: Organism
# TAX_ID: 1234
# ASSEMBLY NAME: name of assembly
# ASSEMBLY DATE: 01-January-2012
# GENOME CENTER: genome center
# DESCRIPTION: testing AGP
# COMMENTS:
# first comment
Chr01	1	100	1	F	SC01	1	100	+
Chr01	101	150	2	N	50	contig	no
Chr01	151	260	3	F	SC02	1	110	-
Chr01	261	300	4	N	40	contig	no
Chr01	301	330	5	F	SC03	1	30	+
);

my $new_agp_str = q(# ORGANISM: Organism
# TAX_ID: 1234
# ASSEMBLY NAME: name of assembly
# ASSEMBLY DATE: 01-January-2012
# GENOME CENTER: genome center
# DESCRIPTION: testing AGP
# COMMENTS:
# first comment
Chr01	1	100	1	F	SC01	1	100	+
Chr01	101	150	2	N	50	contig	no
Chr01	151	260	3	F	SC02	1	110	+
Chr01	261	300	4	N	40	contig	no
Chr01	301	330	5	F	SC03	1	30	+
);

ok($agp->parse_agp($old_agp_str));


ok(my $component_result = $agp->get_component_from_coordinates('Chr01','1'));
ok(my $component_result_2 = $agp->get_component_from_coordinates('Chr01','151'));
ok(my $coord_result_3 = $agp->get_obj_coordinates_from_component('SC03','2'));

if ($component_result) {
  my %result = %{$component_result};
  print STDERR "Component:  ".$result{'component_id'}."\n";
  print STDERR "Coordinate:  ".$result{'component_coordinate'}."\n";
}
if ($component_result_2) {
  my %result2 = %{$component_result_2};
  print STDERR "Component:  ".$result2{'component_id'}."\n";
  print STDERR "Coordinate:  ".$result2{'component_coordinate'}."\n";
}
if ($coord_result_3) {
  my %result = %{$coord_result_3};
  print STDERR "Component:  ".$result{'obj_id'}."\n";
  print STDERR "Coordinate:  ".$result{'obj_coordinate'}."\n";
}
