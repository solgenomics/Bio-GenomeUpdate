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
use Test::More tests => 1;

my %lines;
my $previous_accession;

my $input_file = '/home/jeremy/test_files/chr1_tpf.txt';
my $input_tpf = read_file($input_file);
my $tpf = Bio::GenomeUpdate::TPF->new();
$tpf->parse_tpf($input_tpf);
$tpf->move_scaffold_before("SL2.40sc03594","SL2.40sc03666");
$tpf->flip_scaffold("SL2.40sc03666");



my $out_str_from_tpf_parsed = $tpf->get_formatted_tpf();
write_file('/home/jeremy/test_files/processed_tpf.txt',$out_str_from_tpf_parsed);





# my $scaffold_to_modify="SL2.40sc03666";

my $tpf_str = q(##Organism: An organism
##Chromosome: 1
##Assembly Name: Assembly name
##Strain/Haplotype/Cultivar: cultivar
##Type: Complete Chromosome

##=== Beginning of TPF Data ===

ACCESSION1	?	Scaffold1	MINUS
GAP	TYPE-2	253	PAIRED ENDS
ACCESSION2	?	Scaffold1	MINUS
GAP	TYPE-2	20	PAIRED ENDS
ACCESSION3	?	Scaffold1	MINUS
GAP	TYPE-2	576	PAIRED ENDS
ACCESSION4	?	Scaffold1	MINUS
GAP	TYPE-2	642	PAIRED ENDS
ACCESSION5	?	Scaffold1	MINUS
GAP	TYPE-2	20	PAIRED ENDS
ACCESSION6	?	Scaffold2	MINUS
GAP	TYPE-2	963	PAIRED ENDS
ACCESSION7	?	Scaffold2	MINUS
GAP	TYPE-2	650	PAIRED ENDS
ACCESSION8	?	Scaffold2	MINUS
GAP	TYPE-2	525	PAIRED ENDS
ACCESSION9	?	Scaffold2	MINUS
GAP	TYPE-2	20	PAIRED ENDS
ACCESSION10	?	Scaffold2	MINUS
GAP	TYPE-2	627	PAIRED ENDS
ACCESSION11	?	Scaffold3	MINUS
GAP	TYPE-2	5962	PAIRED ENDS
ACCESSION12	?	Scaffold3	MINUS
GAP	TYPE-2	555	PAIRED ENDS
ACCESSION13	?	Scaffold3	MINUS
GAP	TYPE-2	2517	PAIRED ENDS
ACCESSION14	?	Scaffold3	MINUS
GAP	TYPE-2	999	PAIRED ENDS
ACCESSION15	?	Scaffold3	MINUS
GAP	TYPE-2	1517	PAIRED ENDS
ACCESSION16	?	Scaffold3	MINUS
##=== End of TPF Data ===
);
$tpf = Bio::GenomeUpdate::TPF->new();
$tpf->parse_tpf($tpf_str);

my $scaffold_to_modify="Scaffold2";
my $tpf_copy = Bio::GenomeUpdate::TPF->new();
my $first_line_number = 0;
my $last_line_number = 0;
my $in_scaffold_range = 0;
my $line_number_to_insert_before = 0;
my $scaffold_to_insert_before = "Scaffold3";
my $last_line;

%lines = %{$tpf->get_tpf_lines()};
my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
foreach my $line_key (@sorted_line_numbers) {
  if ($lines{$line_key}->get_line_type() eq "sequence") {
    if ($previous_accession) {
      if ($lines{$line_key}->get_local_contig_identifier() eq $scaffold_to_modify) {
	$in_scaffold_range=1;
	$tpf_copy->add_line_to_end($lines{$line_key});
	if ($first_line_number) {
	  $last_line_number = $line_key;
	} else {
	  $first_line_number = $line_key;
	  $last_line_number = $line_key;
	}
      } else {
	$in_scaffold_range=0;
      }
    } else {
      $in_scaffold_range=0;
    }
    $previous_accession = $lines{$line_key}->get_accession();
  } elsif ($lines{$line_key}->get_line_type() eq "gap") {
    if ($in_scaffold_range==1) {
      $tpf_copy->add_line_to_end($lines{$line_key});
      $last_line_number = $line_key;
    }
  }
}
#remove lines in range
foreach my $line_key (reverse($first_line_number..$last_line_number)){
  $tpf->delete_line($line_key);
}
print STDERR $tpf->get_formatted_tpf();
print STDERR "\n\n\n\n";
print STDERR $tpf_copy->get_formatted_tpf();

%lines = %{$tpf->get_tpf_lines()};
@sorted_line_numbers = sort { $a <=> $b } keys %lines;
foreach my $line_key (@sorted_line_numbers) {
  $last_line = $line_key;
  if ($lines{$line_key}->get_line_type() eq "sequence") {
    if ($line_number_to_insert_before == 0) {
      if ($lines{$line_key}->get_local_contig_identifier() eq $scaffold_to_insert_before) {
	$line_number_to_insert_before = $line_key;
	print STDERR "Insert before: $line_number_to_insert_before\n";
      }
    }
  }
}

my %lines_in_copy = %{$tpf_copy->get_tpf_lines()};
my @sorted_line_numbers_in_copy = sort { $a <=> $b } keys %lines_in_copy;
foreach my $line_key (@sorted_line_numbers_in_copy) {
  $tpf->insert_line_after($line_number_to_insert_before-1,$lines_in_copy{$line_key});
  $line_number_to_insert_before++;
  $last_line++;
}

is($tpf->get_formatted_tpf(),$tpf_str);

print STDERR "\n\n\n\n";
print STDERR $tpf->get_formatted_tpf();

%lines = %{$tpf->get_tpf_lines()};
$last_line=(sort { $a <=> $b } keys %lines)[-1];
@sorted_line_numbers_in_copy = sort { $a <=> $b } keys %lines_in_copy;
my $insert_pos = $last_line;
my $last_line_in_copy;
foreach my $line_key (@sorted_line_numbers_in_copy) {
  #$tpf->insert_line_after($insert_pos+1,$lines_in_copy{$line_key});
  $tpf->add_line_to_end($lines_in_copy{$line_key});
  $insert_pos++;
  $last_line_in_copy = $line_key;
}
$tpf->delete_line($insert_pos+1);
$tpf->insert_line_after($last_line,$lines_in_copy{$last_line_in_copy});

print STDERR "\n\n\n\n";
print STDERR $tpf->get_formatted_tpf();


