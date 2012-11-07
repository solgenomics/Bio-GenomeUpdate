#!/usr/bin/perl

=head1 NAME

genome_update.t
A test for Bio::GenomcdeUpdate::AGP class

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
use File::Slurp;

use Test::More tests => 71;
BEGIN {use_ok( 'Bio::GenomeUpdate::AGP' ); }
require_ok( 'Bio::GenomeUpdate::TPF::TPFSequenceLine' );
require_ok( 'Bio::GenomeUpdate::TPF::TPFGapLine' );
require_ok( 'Bio::GenomeUpdate::TPF' );
require_ok( 'Bio::GenomeUpdate::AGP::AGPSequenceLine' );
require_ok( 'Bio::GenomeUpdate::AGP::AGPGapLine' );
require_ok( 'Bio::GenomeUpdate::AGP' );
require_ok( 'Bio::GenomeUpdate::AGP::AGPConvert' );

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

#create AGP and add lines and gaps
ok(my $agp = Bio::GenomeUpdate::AGP->new());
ok(my $agp_sequence_line = Bio::GenomeUpdate::AGP::AGPSequenceLine->new());
ok(my $agp_gap_line = Bio::GenomeUpdate::AGP::AGPGapLine->new());
ok($agp->set_organism("Organism"));
ok($agp->set_tax_id("1234"));
ok($agp->set_assembly_name("name of assembly"));
#my $assembly_date = DateTime->new(year => 2012, month => 1, day => 1);
my $assembly_date = "01-January-2012";
ok($agp->set_assembly_date($assembly_date));
ok($agp->set_genome_center("genome center"));
ok($agp->set_description("testing AGP"));
ok($agp->add_comment_line("first comment"));
ok($agp->add_comment_line("next comment"));
ok($agp->add_comment_line("last comment"));

#create AGP sequence line
ok($agp_sequence_line->set_object_being_assembled("Chromosome 1"));
ok($agp_sequence_line->set_object_begin("1001"));
ok($agp_sequence_line->set_object_end("1400"));
ok($agp_sequence_line->set_component_type("A"));
ok($agp_sequence_line->set_component_id("component id"));
ok($agp_sequence_line->set_component_begin("1"));
ok($agp_sequence_line->set_component_end("400"));
ok($agp_sequence_line->set_orientation("+"));

#create AGP gap line
ok($agp_gap_line->set_object_being_assembled("Chromosome 1"));
ok($agp_gap_line->set_object_begin("1001"));
ok($agp_gap_line->set_object_end("1400"));
ok($agp_gap_line->set_component_type("N"));
ok($agp_gap_line->set_gap_length("400"));
ok($agp_gap_line->set_gap_type("scaffold"));
ok($agp_gap_line->set_linkage("yes"));
ok($agp_gap_line->add_linkage_evidence("paired-ends"));
ok($agp_gap_line->add_linkage_evidence("map"));

#add lines to AGP file
ok($agp->add_line_to_end($agp_gap_line));
ok($agp->add_line_to_end($agp_sequence_line));
ok($agp->add_line_to_beginning($agp_sequence_line)); #commented out previously
ok($agp->add_line_to_beginning($agp_gap_line)); #commented out previoulsy
ok($agp->insert_line_before(2,$agp_gap_line));
ok($agp->insert_line_after(3,$agp_gap_line)); #commented out previously
ok($agp->delete_line(2,$agp_gap_line)); #commented out previously

#get formatted AGP string and compare to expected output
ok(my $out_str = $agp->get_formatted_agp());
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
Chromosome 1	1001	1400	1	N	400	scaffold	yes	paired-ends;map
Chromosome 1	1001	1400	2	A	component id	1	400	+
Chromosome 1	1001	1400	3	N	400	scaffold	yes	paired-ends;map
Chromosome 1	1001	1400	4	N	400	scaffold	yes	paired-ends;map
Chromosome 1	1001	1400	5	A	component id	1	400	+
);
is ($out_str,$compare_str,'AGP output string is as expected');

#parse formatted AGP string into an AGP object
ok($agp->parse_agp($compare_str));
ok(my $out_str_from_agp_parsed = $agp->get_formatted_agp());
is ($out_str_from_agp_parsed,$compare_str,'AGP output from parsed is as expected');



#my %lines = %{$tpf->get_tpf_lines()};
#my $out = $lines{1}->get_line_type();
#$tpf->print_formatted_tpf();
#print "\n\n";
#$agp->print_formatted_agp();
#print "\n\n";
#$AGP->set_version('1.1');
#$agp->parse_agp("/home/jeremy/Code/short.agp");
#print "\n\n";



my %ncbi_to_local;
my %local_to_ncbi;
my $input_file = "/home/jeremy/Code/component_localID2acc";
my @lines = read_file($input_file);
foreach my $line (@lines){
    chomp($line);
    my @tab_parsed_line = split(/\t/, $line);
    $ncbi_to_local{$tab_parsed_line[1]} = $tab_parsed_line[0];
    $local_to_ncbi{$tab_parsed_line[0]} = $tab_parsed_line[1];
}
	


my $agp_convert = Bio::GenomeUpdate::AGP::AGPConvert->new();
$agp_convert->set_agp_to_convert($agp);
$agp_convert->set_component_id_is_local(1);
$agp_convert->set_ncbi_to_local_conversion(\%ncbi_to_local);
$agp_convert->set_local_to_ncbi_conversion(\%local_to_ncbi);
my @converted_tpfs;
@converted_tpfs = $agp_convert->to_tpf();
#$converted_tpfs[0]->print_formatted_tpf("formatted_tpf_file.tpf");


#gap
#sequence
#gap
#gap
#sequence

