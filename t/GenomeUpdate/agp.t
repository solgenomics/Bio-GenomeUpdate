#!/usr/bin/perl

=head1 NAME

agp.t
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

use Test::More tests => 51;
BEGIN { use_ok('Bio::GenomeUpdate::AGP'); }
require_ok('Bio::GenomeUpdate::AGP::AGPSequenceLine');
require_ok('Bio::GenomeUpdate::AGP::AGPGapLine');
require_ok('Bio::GenomeUpdate::AGP::AGPConvert');

#create AGP and add lines and gaps
ok( my $agp               = Bio::GenomeUpdate::AGP->new() );
ok( my $agp_sequence_line = Bio::GenomeUpdate::AGP::AGPSequenceLine->new() );
ok( my $agp_gap_line      = Bio::GenomeUpdate::AGP::AGPGapLine->new() );
ok( $agp->set_organism("Organism") );
ok( $agp->set_tax_id("1234") );
ok( $agp->set_assembly_name("name of assembly") );

#my $assembly_date = DateTime->new(year => 2012, month => 1, day => 1);
my $assembly_date = "01-January-2012";
ok( $agp->set_assembly_date($assembly_date) );
ok( $agp->set_genome_center("genome center") );
ok( $agp->set_description("testing AGP") );
ok( $agp->add_comment_line("first comment") );
ok( $agp->add_comment_line("next comment") );
ok( $agp->add_comment_line("last comment") );

#create AGP sequence line
ok( $agp_sequence_line->set_object_being_assembled("Chromosome1") );
ok( $agp_sequence_line->set_object_begin("1001") );
ok( $agp_sequence_line->set_object_end("1400") );
ok( $agp_sequence_line->set_component_type("A") );
ok( $agp_sequence_line->set_component_id("component id") );
ok( $agp_sequence_line->set_component_begin("1") );
ok( $agp_sequence_line->set_component_end("400") );
ok( $agp_sequence_line->set_orientation("+") );

#create AGP gap line
ok( $agp_gap_line->set_object_being_assembled("Chromosome1") );
ok( $agp_gap_line->set_object_begin("1001") );
ok( $agp_gap_line->set_object_end("1400") );
ok( $agp_gap_line->set_component_type("N") );
ok( $agp_gap_line->set_gap_length("400") );
ok( $agp_gap_line->set_gap_type("scaffold") );
ok( $agp_gap_line->set_linkage("yes") );
ok( $agp_gap_line->add_linkage_evidence("paired-ends") );
ok( $agp_gap_line->add_linkage_evidence("map") );

#add lines to AGP file
ok( $agp->add_line_to_end($agp_gap_line) );
ok( $agp->add_line_to_end($agp_sequence_line) );
ok( $agp->add_line_to_beginning($agp_sequence_line) ); #commented out previously
ok( $agp->add_line_to_beginning($agp_gap_line) );      #commented out previoulsy
ok( $agp->insert_line_before( 2, $agp_gap_line ) );
ok( $agp->insert_line_after( 3, $agp_gap_line ) );     #commented out previously
ok( $agp->delete_line( 2, $agp_gap_line ) );           #commented out previously

#adding chr2 lines
ok( $agp_sequence_line->set_object_being_assembled("Chromosome2") );
ok( $agp_gap_line->set_object_being_assembled("Chromosome2") );
ok( $agp->add_line_to_end($agp_gap_line) );
ok( $agp->add_line_to_end($agp_sequence_line) );
ok( $agp->add_line_to_beginning($agp_sequence_line) );
ok( $agp->add_line_to_beginning($agp_gap_line) );
ok( $agp->insert_line_before( 2, $agp_gap_line ) );
ok( $agp->insert_line_after( 3, $agp_gap_line ) );
ok( $agp->delete_line( 2, $agp_gap_line ) );


#get formatted AGP string and compare to expected output
ok( my $out_str = $agp->get_formatted_agp() );
my $compare_str = q(##agp-version	2.0
# ORGANISM: Organism
# TAX_ID: 1234
# ASSEMBLY NAME: name of assembly
# ASSEMBLY DATE: 01-January-2012
# GENOME CENTER: genome center
# DESCRIPTION: testing AGP
# COMMENTS:
# first comment
# next comment
# last comment
Chromosome1	1001	1400	1	N	400	scaffold	yes	paired-ends;map
Chromosome1	1001	1400	2	A	component id	1	400	+
Chromosome1	1001	1400	3	N	400	scaffold	yes	paired-ends;map
Chromosome1	1001	1400	4	N	400	scaffold	yes	paired-ends;map
Chromosome1	1001	1400	5	A	component id	1	400	+
Chromosome2	1001	1400	1	N	400	scaffold	yes	paired-ends;map
Chromosome2	1001	1400	2	A	component id	1	400	+
Chromosome2	1001	1400	3	N	400	scaffold	yes	paired-ends;map
Chromosome2	1001	1400	4	N	400	scaffold	yes	paired-ends;map
Chromosome2	1001	1400	5	A	component id	1	400	+
);
is( $out_str, $compare_str, 'AGP output string is as expected' );

#parse formatted AGP string into an AGP object and compare to expected output
ok( $agp->parse_agp($compare_str) );
ok( my $out_str_from_agp_parsed = $agp->get_formatted_agp() );
is( $out_str_from_agp_parsed, $compare_str,
	'AGP output from parsed is as expected' );

#get lines for AGP
ok( $agp->get_description());
is ($agp->get_description(),'testing AGP','Desc is testing as expected');
#is($agp->get_next_line(), 'Bio::GenomeUpdate::AGP::AGPSequenceLine | Bio::GenomeUpdate::AGP::AGPGapLine', 'Got object');
is($agp->get_current_agp_line_number(), 1, 'Got 1 for line number');
ok(my $formatted_agp_line=$agp->get_next_formatted_agp_line());
print $formatted_agp_line;
is($agp->get_current_agp_line_number(), 2, 'Got 2 for line number');
ok(my $formatted_agp_line=$agp->get_next_formatted_agp_line());
print $formatted_agp_line;
