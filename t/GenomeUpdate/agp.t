#!/usr/bin/perl

=head1 NAME

agp.t
A test for Bio::GenomeUpdate::AGP class

=cut

=head1 SYNOPSIS

perl agp.t


=head1 DESCRIPTION

Test script for Bio::GenomeUpdate::AGP class. 
TODO: Should edit coordinates in lines so that it looks like a real AGP file

=head2 Author

Surya Saha <suryasaha at cornell.edu>
Jeremy Edwards <jde22@cornell.edu>

=cut

use strict;
use warnings;
use autodie;

use Test::More tests => 100;
BEGIN { use_ok('Bio::GenomeUpdate::AGP'); }
require_ok('Bio::GenomeUpdate::AGP::AGPSequenceLine');
require_ok('Bio::GenomeUpdate::AGP::AGPGapLine');
require_ok('Bio::GenomeUpdate::AGP::AGPConvert');

#create AGP and add lines and gaps
ok( my $agp               = Bio::GenomeUpdate::AGP->new() );
ok( my $agp_sequence_line_1 = Bio::GenomeUpdate::AGP::AGPSequenceLine->new() );
ok( my $agp_sequence_line_2 = Bio::GenomeUpdate::AGP::AGPSequenceLine->new() );
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

#create AGP sequence line 1
ok( $agp_sequence_line_1->set_object_being_assembled("Chromosome1") );
ok( $agp_sequence_line_1->set_object_begin("1") );
ok( $agp_sequence_line_1->set_object_end("1400") );
#ok( $agp_sequence_line_1->set_part_number("1") );
ok( $agp_sequence_line_1->set_component_type("A") );
ok( $agp_sequence_line_1->set_component_id("component_1") );
ok( $agp_sequence_line_1->set_component_begin("1") );
ok( $agp_sequence_line_1->set_component_end("1400") );
ok( $agp_sequence_line_1->set_orientation("+") );

#create AGP sequence line 2
ok( $agp_sequence_line_2->set_object_being_assembled("Chromosome1") );
ok( $agp_sequence_line_2->set_object_begin("1601") );
ok( $agp_sequence_line_2->set_object_end("2000") );
ok( $agp_sequence_line_2->set_component_type("A") );
ok( $agp_sequence_line_2->set_component_id("component_2") );
ok( $agp_sequence_line_2->set_component_begin("1") );
ok( $agp_sequence_line_2->set_component_end("400") );
ok( $agp_sequence_line_2->set_orientation("+") );

#create AGP gap line
ok( $agp_gap_line->set_object_being_assembled("Chromosome1") );
ok( $agp_gap_line->set_object_begin("1401") );
ok( $agp_gap_line->set_object_end("1600") );
ok( $agp_gap_line->set_component_type("N") );
ok( $agp_gap_line->set_gap_length("200") );
ok( $agp_gap_line->set_gap_type("scaffold") );
ok( $agp_gap_line->set_linkage("yes") );
ok( $agp_gap_line->add_linkage_evidence("paired-ends") );
ok( $agp_gap_line->add_linkage_evidence("map") );

#add lines to AGP file
ok( $agp->add_line_to_end($agp_sequence_line_1) );
ok( $agp->add_line_to_end($agp_gap_line) );
ok( $agp->add_line_to_beginning($agp_sequence_line_1) );
ok( $agp->delete_line( 1, $agp_sequence_line_1 ) );
ok( $agp->insert_line_before( 2, $agp_gap_line ) );
ok( $agp->delete_line( 2, $agp_gap_line ) );
ok( $agp->insert_line_after( 2, $agp_gap_line ) );
ok( $agp->delete_line( 3, $agp_gap_line ) );
ok( $agp->add_line_to_end( $agp_sequence_line_2 ) );

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
Chromosome1	1	1400	1	A	component_1	1	1400	+
Chromosome1	1401	1600	2	N	200	scaffold	yes	paired-ends;map
Chromosome1	1601	2000	3	A	component_2	1	400	+
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
ok($formatted_agp_line=$agp->get_next_formatted_agp_line());
print $formatted_agp_line;

#test summary methods
is($agp->get_number_of_gap_lines(),1,'Gap line count is as expected');
is($agp->get_number_of_sequence_lines(),2,'Sequence line count is as expected');
my @compare_gap_lengths=(200);
ok(my @gap_lengths_from_agp=$agp->get_gap_lengths());
is_deeply(\@gap_lengths_from_agp,\@compare_gap_lengths,'Gap lengths as expected');
my @compare_sequence_lengths=(1400,400);
ok(my @sequence_lengths_from_agp=$agp->get_sequence_lengths());
is_deeply(\@sequence_lengths_from_agp,\@compare_sequence_lengths,'Sequence lengths as expected');

#test overlap methods
print STDERR "\ntesting overlap methods..\n";
ok(my ($cov_seq_count, $cov_seq_length, $par_cov_seq_count, $par_cov_seq_length) = $agp->get_sequence_overlap(1201,1400), 'testing seq region 1201-1400bp'); 
is($cov_seq_count,0,'Covered sequence component count as expected for 1201-1400bp');
is($cov_seq_length,0,'Covered sequence component length as expected for 1201-1400bp');
is($par_cov_seq_count,1,'Partially covered sequence component count as expected for 1201-1400bp');
is($par_cov_seq_length,200,'Partially covered sequence component length as expected for 1201-1400bp');

ok(($cov_seq_count, $cov_seq_length, $par_cov_seq_count, $par_cov_seq_length) = $agp->get_sequence_overlap(1,1400), 'testing seq region 1-1400bp');
is($cov_seq_count,1,'Covered sequence component count as expected for 1-1400bp');
is($cov_seq_length,1400,'Covered sequence component length as expected for 1-1400bp');
is($par_cov_seq_count,0,'Partially covered sequence component count as expected for 1-1400bp');
is($par_cov_seq_length,0,'Partially covered sequence component length as expected for 1-1400bp');

#region covers both gap and seq region  
ok(($cov_seq_count, $cov_seq_length, $par_cov_seq_count, $par_cov_seq_length) = $agp->get_sequence_overlap(1201,1500), 'testing seq region 1201-1500bp');
is($par_cov_seq_count,1,'Partially covered sequence component count as expected for 1201-1500bp');
is($par_cov_seq_length,200,'Partially covered sequence component length as expected for 1201-1500bp');

ok(my ($cov_gap_count, $cov_gap_length, $par_cov_gap_count, $par_cov_gap_length) = $agp->get_gap_overlap(1401,1500), 'testing gap region 1401-1500bp'); 
is($cov_gap_count,0,'Covered gap component count as expected for 1401-1500bp');
is($cov_gap_length,0,'Covered gap component length as expected for 1401-1500bp');
is($par_cov_gap_count,1,'Partially covered gap component count as expected for 1401-1500bp');
is($par_cov_gap_length,100,'Partially covered gap component length as expected for 1401-1500bp');

ok(($cov_gap_count, $cov_gap_length, $par_cov_gap_count, $par_cov_gap_length) = $agp->get_gap_overlap(1401,1600), 'testing gap region 1401-1600bp');
is($cov_gap_count,1,'Covered gap component count as expected for 1401-1600bp');
is($cov_gap_length,200,'Covered gap component length as expected for 1401-1600bp');
is($par_cov_gap_count,0,'Partially covered gap component count as expected for 1401-1600bp');
is($par_cov_gap_length,0,'Partially covered gap component length as expected for 1401-1600bp');

#region covers both gap and seq region  
ok(($cov_gap_count, $cov_gap_length, $par_cov_gap_count, $par_cov_gap_length) = $agp->get_gap_overlap(1301,1500), 'testing gap region 1301-1500bp');
is($par_cov_gap_count,1,'Partially covered gap component count as expected for 1301-1500bp');
is($par_cov_gap_length,100,'Partially covered gap component length as expected for 1301-1500bp');

#region covers both gap and seq region for seq_id in AGP  
ok(($cov_gap_count, $cov_gap_length, $par_cov_gap_count, $par_cov_gap_length) = $agp->get_gap_overlap(1301,1500,'Chromosome1'), 'testing gap region 1301-1500bp for Chromosome1');
is($par_cov_gap_count,1,'Partially covered gap component count as expected for 1301-1500bp');
is($par_cov_gap_length,100,'Partially covered gap component length as expected for 1301-1500bp');

#region covers both gap and seq region for seq_id NOT in AGP  
ok(($cov_gap_count, $cov_gap_length, $par_cov_gap_count, $par_cov_gap_length) = $agp->get_gap_overlap(1301,1500,'Chromosome2'), 'testing gap region 1301-1500bp for dummy Chromosome2');
is($par_cov_gap_count,0,'Partially covered gap component count as expected for 1301-1500bp');
is($par_cov_gap_length,0,'Partially covered gap component length as expected for 1301-1500bp');
