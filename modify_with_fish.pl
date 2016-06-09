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


my $input_file = "chr1_tpf.txt";
my $input_tpf = read_file($input_file);
my $tpf = TPF->new();
$tpf->parse_tpf($input_tpf);
my $out_str_from_tpf_parsed = $tpf->get_formatted_tpf();
write_file('processed_tpf.txt',$out_str_from_tpf_parsed);


