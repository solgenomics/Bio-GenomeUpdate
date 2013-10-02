#!/usr/bin/perl

=head1 NAME

gff.t
A test for Bio::GenomeUpdate::GFF and Bio::GenomeUpdate::GFF::GFFRearrange classes

=cut

=head1 SYNOPSIS

perl gff.t


=head1 DESCRIPTION



=head2 Author

Surya Saha <suryasaha@cornell.edu , @SahaSurya> 

=cut

use strict;
use warnings;
use autodie;

use Test::More tests => 15;
BEGIN { use_ok('Bio::GenomeUpdate::GFF'); }
BEGIN { use_ok('Bio::GenomeUpdate::AGP'); }
require_ok('Bio::GenomeUpdate::GFF::GFFRearrange');#uses AGP module
require_ok('Bio::GFF3::LowLevel');

#load orig AGP
my $agp_orig_file = q(##agp-version	2.0
# ORGANISM: Organism
# TAX_ID: 1234
# ASSEMBLY NAME: name of assembly
# ASSEMBLY DATE: 01-January-2012
# GENOME CENTER: genome center
# DESCRIPTION: original AGP
# COMMENTS:
# first comment
# next comment
# last comment
CHR01_FISH2_GAPS	1	200	1	F	SL2.40SC04133	1	200	0
CHR01_FISH2_GAPS	201	300	2	N	99	contig	no	
CHR01_FISH2_GAPS	301	600	3	F	SL2.40SC04191	1	300	+
CHR01_FISH2_GAPS	601	800	4	N	200	contig	no	
CHR01_FISH2_GAPS	801	1000	5	F	SL2.40SC03666	1	200	+
CHR01_FISH2_GAPS	1001	1100	4	N	200	contig	no	
CHR01_FISH2_GAPS	1101	1500	5	F	SL2.40SC03661	1	200	-
);
ok( my $agp_orig = Bio::GenomeUpdate::AGP->new(),'create orig AGP obj');
ok( $agp_orig->parse_agp($agp_orig_file),'parse orig AGP obj');

#load updated AGP, gap size change, scaffold flip 
my $agp_fish_file = q(##agp-version	2.0
# ORGANISM: Organism
# TAX_ID: 1234
# ASSEMBLY NAME: name of assembly
# ASSEMBLY DATE: 01-January-2012
# GENOME CENTER: genome center
# DESCRIPTION: original AGP
# COMMENTS:
# first comment
# next comment
# last comment
CHR01_FISH2_GAPS	1	200	1	F	SL2.40SC04133	1	200	0
CHR01_FISH2_GAPS	201	350	2	N	150	contig	no	
CHR01_FISH2_GAPS	351	650	3	F	SL2.40SC04191	1	300	+
CHR01_FISH2_GAPS	651	800	4	N	150	contig	no	
CHR01_FISH2_GAPS	801	1000	5	F	SL2.40SC03666	1	200	-
CHR01_FISH2_GAPS	1001	1100	4	N	200	contig	no	
CHR01_FISH2_GAPS	1101	1500	5	F	SL2.40SC03661	1	200	+
);
ok(my $agp_fish = Bio::GenomeUpdate::AGP->new(),'create fish AGP obj');
ok( $agp_fish->parse_agp($agp_fish_file),'parse fish AGP obj');
  
# gff 
my $gff_file = q(##gff3
CHR01_FISH2_GAPS	src	CDS	50	150	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	320	400	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	450	500	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	550	600	0	-	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	801	900	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	850	900	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	925	940	0	-	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	950	969	0	-	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	1200	1300	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	1400	1475	0	-	0	ID=attributes
);


ok( my $gff = Bio::GenomeUpdate::GFF->new(),'create GFF obj');
ok( $gff->parse_gff($gff_file),'parse GFF obj');

##update coordinates, calls GFFRearrange object
#ok( my %coords = $gff->get_reordered_coordinates($agp_orig,$agp_fish),'get_reordered_coordinates');
#ok( my %flips = $gff->get_flipped_coordinates($agp_orig,$agp_fish),'get_flipped_coordinates');
#
##print STDERR "flips: ",scalar keys %flips,"\n";
#
##remap GFF
#ok( $gff->remap_coordinates_hash(\%coords,\%flips),'remap_coordinates');

#validate remapping
my $compare_str = q(##gff3
CHR01_FISH2_GAPS	src	CDS	50	150	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	370	450	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	500	550	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	600	650	0	-	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	901	1000	0	-	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	901	951	0	-	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	861	876	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	832	851	0	+	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	1301	1401	0	-	0	ID=attributes
CHR01_FISH2_GAPS	src	CDS	1126	1201	0	+	0	ID=attributes
);
#ok( my $gff_fish = $gff->get_formatted_gff(),'get_formatted_gff');
#is( $gff_fish, $compare_str, 'GFF remapping using hash based method is as expected');

#remap GFF using optimized method
ok( $gff->remap_coordinates($agp_orig,$agp_fish),'remap_coordinates using optimized method');
ok( my $gff_fish = $gff->get_formatted_gff(),'get_formatted_gff');
is( $gff_fish, $compare_str, 'GFF remapping using optimized method is as expected');
