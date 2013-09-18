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

use Test::More tests => 14;
BEGIN { use_ok('Bio::GenomeUpdate::GFF'); }
BEGIN { use_ok('Bio::GenomeUpdate::AGP'); }
require_ok('Bio::GenomeUpdate::GFF::GFFRearrange');#uses AGP module
require_ok('Bio::GFF3::LowLevel::Parser');

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
CHR01_FISH2_GAPS	1	200	1	F	SL2.40SC04133	1	200	+
CHR01_FISH2_GAPS	201	300	2	N	99	contig	no	
CHR01_FISH2_GAPS	301	600	3	F	SL2.40SC04191	1	300	+
CHR01_FISH2_GAPS	601	800	4	N	200	contig	no	
CHR01_FISH2_GAPS	801	1000	5	F	SL2.40SC03666	1	200	+
);
ok(my $agp_orig = Bio::GenomeUpdate::AGP->new());
ok( $agp_orig->parse_agp($agp_orig_file));

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
CHR01_FISH2_GAPS	1	200	1	F	SL2.40SC04133	1	200	+
CHR01_FISH2_GAPS	201	350	2	N	150	contig	no	
CHR01_FISH2_GAPS	351	650	3	F	SL2.40SC04191	1	300	+
CHR01_FISH2_GAPS	651	800	4	N	150	contig	no	
CHR01_FISH2_GAPS	801	1000	5	F	SL2.40SC03666	1	200	-
);
ok(my $agp_fish = Bio::GenomeUpdate::AGP->new() );
ok( $agp_fish->parse_agp($agp_fish_file));
  
# gff 
my $gff_file = q(##gff3
CHR01_FISH2_GAPS	src	CDS	50	150	0	+	0	attributes
CHR01_FISH2_GAPS	src	CDS	320	400	0	+	0	attributes
CHR01_FISH2_GAPS	src	CDS	450	500	0	+	0	attributes
CHR01_FISH2_GAPS	src	CDS	550	600	0	-	0	attributes
CHR01_FISH2_GAPS	src	CDS	801	900	0	+	0	attributes
CHR01_FISH2_GAPS	src	CDS	950	969	0	-	0	attributes
);


ok( my $gff = Bio::GenomeUpdate::GFF->new());
ok( $gff->parse_gff($gff_file));

#update coordinates, calls GFFRearrange object
ok( my %coords = $gff->reorder_coordinates_AGP($agp_orig,$agp_fish));

#remap GFF
ok( $gff->remap_coordinates_AGP(%coords));

#validate remapping
my $compare_str = q(##gff3
CHR01_FISH2_GAPS	src	CDS	50	150	0	+	0	attributes
CHR01_FISH2_GAPS	src	CDS	370	450	0	+	0	attributes
CHR01_FISH2_GAPS	src	CDS	500	550	0	+	0	attributes
CHR01_FISH2_GAPS	src	CDS	600	650	0	-	0	attributes
CHR01_FISH2_GAPS	src	CDS	901	1000	0	-	0	attributes
CHR01_FISH2_GAPS	src	CDS	831	850	0	+	0	attributes
);
ok( my $gff_fish = $gff->get_updated_gff());
is( $gff_fish, $compare_str, 'GFF remapping is as expected');