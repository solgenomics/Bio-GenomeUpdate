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
ok(my $tpf = Bio::GenomeUpdate::TPF->new());
$tpf->parse_tpf($input_tpf);
my $out_str_from_tpf_parsed = $tpf->get_formatted_tpf();
write_file('/home/jeremy/test_files/processed_tpf.txt',$out_str_from_tpf_parsed);

my $scaffold_to_modify="SL2.40sc04133";


%lines = %{$tpf->get_tpf_lines()};
my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
foreach my $line_key (@sorted_line_numbers) {
  if ($lines{$line_key}->get_line_type() eq "sequence") {
    if ($previous_accession) {
      if ($lines{$line_key}->get_local_contig_identifier() eq $scaffold_to_modify) {
      }
    }
    else {
    }
    $previous_accession = $lines{$line_key}->get_accession();
  }
}
