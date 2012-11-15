#!/usr/bin/perl

=head1 NAME

agpconvert.t
A test for Bio::GenomcdeUpdate::AGP class

=cut

=head1 SYNOPSIS

perl agpconvert.t


=head1 DESCRIPTION



=head2 Author

Jeremy Edwards <jde22@cornell.edu>
=cut

use strict;
use warnings;
use autodie;
use File::Slurp;

use Test::More tests => 16;
BEGIN {use_ok( 'Bio::GenomeUpdate::AGP::AGPConvert' ); }
require_ok( 'Bio::GenomeUpdate::TPF::TPFSequenceLine' );
require_ok( 'Bio::GenomeUpdate::TPF::TPFGapLine' );
require_ok( 'Bio::GenomeUpdate::TPF' );
require_ok( 'Bio::GenomeUpdate::AGP::AGPSequenceLine' );
require_ok( 'Bio::GenomeUpdate::AGP::AGPGapLine' );
require_ok( 'Bio::GenomeUpdate::AGP' );

my %ncbi_to_local;
my %local_to_ncbi;
#my $input_file = "/home/jeremy/Code/component_localID2acc";
#my @lines = read_file($input_file);
#foreach my $line (@lines){
#    chomp($line);
#    my @tab_parsed_line = split(/\t/, $line);
#    $ncbi_to_local{$tab_parsed_line[1]} = $tab_parsed_line[0];
#    $local_to_ncbi{$tab_parsed_line[0]} = $tab_parsed_line[1];
$ncbi_to_local{'ncbi1'} = 'local1';
$ncbi_to_local{'ncbi2'} = 'local2';
$local_to_ncbi{'local1'} = 'ncbi1';
$local_to_ncbi{'local2'} = 'ncbi2';

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
Chromosome 1	1	1000	1	A	local1	1	1000	+
Chromosome 1	1001	1400	2	N	400	scaffold	yes	paired-ends
Chromosome 1	1401	1800	3	A	local2	1	400	+
);

ok(my $agp = Bio::GenomeUpdate::AGP->new());
ok($agp->parse_agp($compare_str));
ok(my $agp_convert = Bio::GenomeUpdate::AGP::AGPConvert->new());
ok($agp_convert->set_agp_to_convert($agp));
ok($agp_convert->set_component_id_is_local(1));
ok($agp_convert->set_ncbi_to_local_conversion(\%ncbi_to_local));
ok($agp_convert->set_local_to_ncbi_conversion(\%local_to_ncbi));
my @converted_tpfs;
ok(@converted_tpfs = $agp_convert->to_tpf());
ok (my $out_str_from_tpf_parsed = $converted_tpfs[0]->get_formatted_tpf());
print STDERR $out_str_from_tpf_parsed;







