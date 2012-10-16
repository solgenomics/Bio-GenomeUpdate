#!/usr/bin/perl
use strict;

use Bio::GenomeUpdate::TPF::TPFSequenceLine;
use Bio::GenomeUpdate::TPF::TPFGapLine;
use Bio::GenomeUpdate::TPF;
use Bio::GenomeUpdate::AGP::AGPSequenceLine;
use Bio::GenomeUpdate::AGP::AGPGapLine;
use Bio::GenomeUpdate::AGP;
use Bio::GenomeUpdate::AGP::AGPConvert;
use File::Slurp;

my $sequence_line = TPFSequenceLine->new();
#$sequence_line->set_contains("CONTAINED_TURNOUT");
$sequence_line->set_local_contig_identifier("local_contig");
$sequence_line->set_containing_accession("Accession");
$sequence_line->set_containing_clone_name("Clone");
$sequence_line->set_orientation("PLUS");

my $gap_line = TPFGapLine->new();
$gap_line->set_gap_size("100");
$gap_line->set_gap_type("TYPE-3");
$gap_line->add_gap_method("PAIRED ENDS");
$gap_line->add_gap_method("PCR");
my $tpf = TPF->new();
$tpf->add_line_to_end($gap_line);
$tpf->add_line_to_end($sequence_line);
$tpf->add_line_to_beginning($sequence_line);
$tpf->add_line_to_beginning($gap_line);
$tpf->insert_line_before(2,$gap_line);
$tpf->insert_line_after(3,$gap_line);
$tpf->delete_line(2,$gap_line);
$tpf->set_organism("An organism");
$tpf->set_assembly_name("Assembly name");
$tpf->set_chromosome("1");
$tpf->set_strain_haplotype_cultivar("cultivar");
my $agp = AGP->new();

my $agp_sequence_line = AGPSequenceLine->new();
my $agp_gap_line = AGPGapLine->new();

$agp->set_organism("Organism");
$agp->set_tax_id("1234");
$agp->set_assembly_name("name of assembly");
#my $assembly_date = DateTime->new(year => 2012, month => 1, day => 1);
my $assembly_date = "01-January-2012";

$agp->set_assembly_date($assembly_date);
$agp->set_genome_center("genome center");
$agp->set_description("testing AGP");
$agp->add_comment_line("first comment");
$agp->add_comment_line("next comment");
$agp->add_comment_line("last comment");

$agp_sequence_line->set_object_being_assembled("Chromosome 1");
$agp_sequence_line->set_object_begin("1001");
$agp_sequence_line->set_object_end("1400");
$agp_sequence_line->set_component_type("A");
$agp_sequence_line->set_component_id("component id");
$agp_sequence_line->set_component_begin("1");
$agp_sequence_line->set_component_end("400");
$agp_sequence_line->set_orientation("+");

$agp_gap_line->set_object_being_assembled("Chromosome 1");
$agp_gap_line->set_object_begin("1001");
$agp_gap_line->set_object_end("1400");
$agp_gap_line->set_component_type("N");
$agp_gap_line->set_gap_length("500");
$agp_gap_line->set_gap_type("scaffold");
$agp_gap_line->set_linkage("yes");
$agp_gap_line->add_linkage_evidence("paired-ends");
$agp_gap_line->add_linkage_evidence("map");

$agp->add_line_to_end($agp_gap_line);
$agp->add_line_to_end($agp_sequence_line);
#$agp->add_line_to_beginning($agp_sequence_line);
#$agp->add_line_to_beginning($agp_gap_line);
$agp->insert_line_before(2,$agp_gap_line);
#$agp->insert_line_after(3,$agp_gap_line);
#$agp->delete_line(2,$agp_gap_line);



#my %lines = %{$tpf->get_tpf_lines()};
#my $out = $lines{1}->get_line_type();

#$tpf->print_formatted_tpf();
#print "\n\n";
#$agp->print_formatted_agp();
#print "\n\n";
#$agp->set_version('1.1');
$agp->parse_agp("short.agp");
#print "\n\n";
$agp->print_formatted_agp('agp_outfile_test.agp');



my %ncbi_to_local;
my %local_to_ncbi;
my $input_file = "component_localID2acc";
my @lines = read_file($input_file);
foreach my $line (@lines){
    chomp($line);
    my @tab_parsed_line = split(/\t/, $line);
    $ncbi_to_local{$tab_parsed_line[1]} = $tab_parsed_line[0];
    $local_to_ncbi{$tab_parsed_line[0]} = $tab_parsed_line[1];
}
	


my $agp_convert = AGPConvert->new();
$agp_convert->set_agp_to_convert($agp);
$agp_convert->set_component_id_is_local(1);
$agp_convert->set_ncbi_to_local_conversion(\%ncbi_to_local);
$agp_convert->set_local_to_ncbi_conversion(\%local_to_ncbi);
my @converted_tpfs;
@converted_tpfs = $agp_convert->to_tpf();
$converted_tpfs[0]->print_formatted_tpf("formatted_tpf_file.tpf");


#gap
#sequence
#gap
#gap
#sequence

