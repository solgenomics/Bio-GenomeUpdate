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
SL2.40ch01	1	32987597	1	W	SL2.40sc04133	1	32987597	-
SL2.40ch01	32987598	32987697	2	U	100	contig	no
SL2.40ch01	32987698	35621724	3	W	SL2.40sc03666	1	2634027	+
SL2.40ch01	35621725	35621824	4	U	100	contig	no
SL2.40ch01	35621825	37344421	5	W	SL2.40sc04191	1	1722597	+
SL2.40ch01	37344422	37344521	6	U	100	contig	no
SL2.40ch01	37344522	43638707	7	W	SL2.40sc03594	1	6294186	+
SL2.40ch01	43638708	43638807	8	U	100	contig	no
SL2.40ch01	43638808	64098855	9	W	SL2.40sc05010	1	20460048	-
SL2.40ch01	64098856	64098955	10	U	100	contig	no
SL2.40ch01	64098956	70260226	11	W	SL2.40sc05941	1	6161271	-
SL2.40ch01	70260227	70260326	12	U	100	contig	no
SL2.40ch01	70260327	70630432	13	W	SL2.40sc06903	1	370106	+
SL2.40ch01	70630433	70630532	14	U	100	contig	no
SL2.40ch01	70630533	73352533	15	W	SL2.40sc06917	1	2722001	+
SL2.40ch01	73352534	73352633	16	U	100	contig	no
SL2.40ch01	73352634	90304244	17	W	SL2.40sc04323	1	16951611	+
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
SL2.40ch01	1	32987597	1	F	SL2.40sc04133	1	32987597	+
SL2.40ch01	32987598	35267597	2	N	2280000	contig	no	
SL2.40ch01	35267598	36990194	3	F	SL2.40sc04191	1	1722597	+
SL2.40ch01	36990195	39120194	4	N	2130000	contig	no	
SL2.40ch01	39120195	41754221	5	F	SL2.40sc03666	1	2634027	-
SL2.40ch01	41754222	42324221	6	N	570000	contig	no	
SL2.40ch01	42324222	48618407	7	F	SL2.40sc03594	1	6294186	+
SL2.40ch01	48618408	50738407	8	N	2120000	contig	no	
SL2.40ch01	50738408	71198455	9	F	SL2.40sc05010	1	20460048	+
SL2.40ch01	71198456	71708455	10	N	510000	contig	no	
SL2.40ch01	71708456	77869726	11	F	SL2.40sc05941	1	6161271	+
SL2.40ch01	77869727	78119726	12	N	250000	contig	no	
SL2.40ch01	78119727	80841727	13	F	SL2.40sc06917	1	2722001	+
SL2.40ch01	80841728	81011727	14	N	170000	contig	no	
SL2.40ch01	81011728	81381833	15	F	SL2.40sc06903	1	370106	+
SL2.40ch01	81381834	81591833	16	N	210000	contig	no	
SL2.40ch01	81591834	98543444	17	F	SL2.40sc04323	1	16951611	+
);
ok(my $agp_fish = Bio::GenomeUpdate::AGP->new(),'create fish AGP obj');
ok( $agp_fish->parse_agp($agp_fish_file),'parse fish AGP obj');
  
# gff 
my $gff_file = q(##gff3
SL2.40ch01	ITAG_infernal	transcript	598896	599053	81.85	+	.	Alias=U1;ID=Solyc01r005915.1.1;Name=RF00003;e-value=1.817e-16;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	610656	610727	68.71	-	.	Alias=tRNA;ID=Solyc01r005935.1.1;Name=RF00005;e-value=3.034e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	870450	870522	57.91	+	.	Alias=tRNA;ID=Solyc01r006265.1.1;Name=RF00005;e-value=2.98e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1298413	1298603	40.65	+	.	Alias=U3;ID=Solyc01r006715.1.1;Name=RF00012;e-value=4.27e-08;rna_type=snRNA
SL2.40ch01	ITAG_infernal	transcript	1346414	1346485	64.77	-	.	Alias=tRNA;ID=Solyc01r006755.1.1;Name=RF00005;e-value=3.232e-13;rna_type=tRNA
);


ok( my $gff = Bio::GenomeUpdate::GFF->new(),'create GFF obj');
ok( $gff->parse_gff($gff_file),'parse GFF obj');

#update coordinates, calls GFFRearrange object
ok( my %coords = $gff->get_reordered_coordinates($agp_orig,$agp_fish),'get_reordered_coordinates');
#ok( my %flips = $gff->get_flipped_coordinates($agp_orig,$agp_fish),'get_flipped_coordinates');
#
##print STDERR "flips: ",scalar keys %flips,"\n";
#
##remap GFF
#ok( $gff->remap_coordinates(\%coords,\%flips),'remap_coordinates');
#
##validate remapping
#my $compare_str = q(##gff3
#CHR01_FISH2_GAPS	src	CDS	50	150	0	+	0	ID=attributes
#CHR01_FISH2_GAPS	src	CDS	370	450	0	+	0	ID=attributes
#CHR01_FISH2_GAPS	src	CDS	500	550	0	+	0	ID=attributes
#CHR01_FISH2_GAPS	src	CDS	600	650	0	-	0	ID=attributes
#CHR01_FISH2_GAPS	src	CDS	901	1000	0	-	0	ID=attributes
#CHR01_FISH2_GAPS	src	CDS	901	951	0	-	0	ID=attributes
#CHR01_FISH2_GAPS	src	CDS	861	876	0	+	0	ID=attributes
#CHR01_FISH2_GAPS	src	CDS	832	851	0	+	0	ID=attributes
#);
#ok( my $gff_fish = $gff->get_formatted_gff(),'get_formatted_gff');
##is( $gff_fish, $compare_str, 'GFF remapping is as expected');

