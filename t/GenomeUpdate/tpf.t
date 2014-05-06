#!/usr/bin/perl

=head1 NAME

tpf.t
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

use Test::More tests => 34;
BEGIN {use_ok( 'Bio::GenomeUpdate::TPF' ); }
require_ok( 'Bio::GenomeUpdate::TPF::TPFSequenceLine' );
require_ok( 'Bio::GenomeUpdate::TPF::TPFGapLine' );

#create a sequence line
ok(my $sequence_line = Bio::GenomeUpdate::TPF::TPFSequenceLine->new());
#ok(my$sequence_line->set_contains("CONTAINED_TURNOUT"));
ok($sequence_line->set_local_contig_identifier("local_contig"));
ok($sequence_line->set_containing_accession("Accession"));
ok($sequence_line->set_containing_clone_name("Clone"));
ok ($sequence_line->set_orientation("PLUS"));

#create a gap line
ok(my $gap_line = Bio::GenomeUpdate::TPF::TPFGapLine->new());
ok($gap_line->set_gap_size("100"));
ok($gap_line->set_gap_type("TYPE-3"));
ok($gap_line->add_gap_method("PAIRED ENDS"));
ok($gap_line->add_gap_method("PCR"));

#create a TPF and add lines and gaps
ok(my $tpf = Bio::GenomeUpdate::TPF->new());
ok($tpf->add_line_to_end($gap_line));
ok($tpf->add_line_to_end($sequence_line));
ok($tpf->add_line_to_beginning($sequence_line));
ok($tpf->add_line_to_beginning($gap_line));
ok($tpf->insert_line_before(2,$gap_line));
ok($tpf->insert_line_after(3,$gap_line));
ok($tpf->delete_line(2,$gap_line));
ok($tpf->set_organism("An organism"));
ok($tpf->set_assembly_name("Assembly name"));
ok($tpf->set_chromosome("1"));
ok($tpf->set_strain_haplotype_cultivar("cultivar"));

#get formatted TPF string and compare to expected output
ok(my $out_str = $tpf->get_formatted_tpf());
#print STDERR $out_str;
my $compare_str = q(##Organism: An organism
##Chromosome: 1
##Assembly Name: Assembly name
##Strain/Haplotype/Cultivar: cultivar
##Type: Complete Chromosome

##=== Beginning of TPF Data ===

GAP	TYPE-3	100	PAIRED ENDS;PCR
??	?	local_contig	PLUS
GAP	TYPE-3	100	PAIRED ENDS;PCR
GAP	TYPE-3	100	PAIRED ENDS;PCR
??	?	local_contig	PLUS
##=== End of TPF Data ===
);
is ($out_str,$compare_str,'TPF output string is as expected');

#parse formatted TPF string into an TPF object and compare to expected output
ok($tpf->parse_tpf($compare_str));
ok(my $out_str_from_tpf_parsed = $tpf->get_formatted_tpf());
is ($out_str_from_tpf_parsed,$compare_str,'TPF output from parsed is as expected');

#summary functions
is($tpf->get_number_of_gap_lines(),3,'Gap line count is as expected');
is($tpf->get_number_of_sequence_lines(),2,'Sequence line count is as expected');
my @compare_gap_lengths=(100,100,100);
ok(my @gap_lengths_from_tpf=$tpf->get_gap_lengths());
is(@gap_lengths_from_tpf,@compare_gap_lengths,'Gap lengths as expected');