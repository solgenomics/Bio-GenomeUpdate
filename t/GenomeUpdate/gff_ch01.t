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
use Proc::ProcessTable;

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
SL2.40ch01	ITAG_infernal	transcript	1747195	1747311	93.77	-	.	Alias=U5;ID=Solyc01r007185.1.1;Name=RF00020;e-value=2.185e-18;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	1748343	1748503	142.78	-	.	Alias=U1;ID=Solyc01r007186.1.1;Name=RF00003;e-value=6.947e-26;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	1749044	1749160	95.05	-	.	Alias=U5;ID=Solyc01r007187.1.1;Name=RF00020;e-value=1.28e-18;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	1750194	1750354	139.74	-	.	Alias=U1;ID=Solyc01r007188.1.1;Name=RF00003;e-value=2.197e-25;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	1750947	1751064	95.19	-	.	Alias=U5;ID=Solyc01r007189.1.1;Name=RF00020;e-value=1.209e-18;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	1751971	1752131	121.52	-	.	Alias=U1;ID=Solyc01r113630.1.1;Name=RF00003;e-value=2.229e-22;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	1753277	1753436	127.52	-	.	Alias=U1;ID=Solyc01r113640.1.1;Name=RF00003;e-value=2.286e-23;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	1756114	1756230	94.54	-	.	Alias=U5;ID=Solyc01r113650.1.1;Name=RF00020;e-value=1.117e-18;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	1757310	1757470	148.73	-	.	Alias=U1;ID=Solyc01r113660.1.1;Name=RF00003;e-value=1.928e-27;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	1758422	1758538	89.46	-	.	Alias=U5;ID=Solyc01r113670.1.1;Name=RF00020;e-value=9.302e-18;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	1762314	1762430	94.54	-	.	Alias=U5;ID=Solyc01r113680.1.1;Name=RF00020;e-value=7.76e-19;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	1763038	1763195	79.7	-	.	Alias=U1;ID=Solyc01r113690.1.1;Name=RF00003;e-value=8.822e-16;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	1764439	1764599	139.59	-	.	Alias=U1;ID=Solyc01r113700.1.1;Name=RF00003;e-value=1.166e-25;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	1765528	1765644	85.63	-	.	Alias=U5;ID=Solyc01r113710.1.1;Name=RF00020;e-value=3.194e-17;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	1839905	1839978	60.61	-	.	Alias=Intron_gpII;ID=Solyc01r007275.1.1;Name=RF00029;e-value=5.648e-13;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1841621	1841707	52.74	+	.	Alias=tRNA;ID=Solyc01r007276.1.1;Name=RF00005;e-value=3.445e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1843009	1843081	64.28	-	.	Alias=tRNA;ID=Solyc01r007285.1.1;Name=RF00005;e-value=2.552e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1844791	1844863	65.55	+	.	Alias=tRNA;ID=Solyc01r007286.1.1;Name=RF00005;e-value=1.482e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1849119	1849191	61.19	+	.	Alias=tRNA;ID=Solyc01r007315.1.1;Name=RF00005;e-value=9.461e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1863651	1863724	65.3	-	.	Alias=tRNA;ID=Solyc01r007435.1.1;Name=RF00005;e-value=1.651e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1863889	1863962	60.86	-	.	Alias=tRNA;ID=Solyc01r007436.1.1;Name=RF00005;e-value=1.089e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1868463	1868523	42.53	-	.	Alias=Intron_gpII;ID=Solyc01r007485.1.1;Name=RF00029;e-value=3.657e-09;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1872968	1873031	50.09	+	.	Alias=Intron_gpII;ID=Solyc01r007535.1.1;Name=RF00029;e-value=9.322e-11;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1874539	1874613	45.81	+	.	Alias=Intron_gpII;ID=Solyc01r007536.1.1;Name=RF00029;e-value=7.447e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1878776	1878852	54.95	-	.	Alias=Intron_gpII;ID=Solyc01r007595.1.1;Name=RF00029;e-value=8.808e-12;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1881903	1881980	49.19	-	.	Alias=Intron_gpII;ID=Solyc01r007625.1.1;Name=RF00029;e-value=1.442e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1883426	1883499	61.58	-	.	Alias=tRNA;ID=Solyc01r007645.1.1;Name=RF00005;e-value=8.025e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1891036	1891116	52.81	-	.	Alias=tRNA;ID=Solyc01r007655.1.1;Name=RF00005;e-value=3.34e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1892412	1892490	56.82	-	.	Alias=Intron_gpII;ID=Solyc01r007656.1.1;Name=RF00029;e-value=3.553e-12;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1894700	1894773	54.01	-	.	Alias=Intron_gpII;ID=Solyc01r007675.1.1;Name=RF00029;e-value=1.393e-11;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1895468	1895544	58.12	-	.	Alias=Intron_gpII;ID=Solyc01r007676.1.1;Name=RF00029;e-value=1.894e-12;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1897101	1897172	58.56	+	.	Alias=tRNA;ID=Solyc01r007695.1.1;Name=RF00005;e-value=2.899e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	1899840	1899958	39.64	+	.	Alias=Intron_gpII;ID=Solyc01r007725.1.1;Name=RF00029;e-value=1.488e-08;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	1900794	1900906	39.36	+	.	Alias=Intron_gpII;ID=Solyc01r007726.1.1;Name=RF00029;e-value=1.703e-08;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	2690856	2690931	42.28	-	.	Alias=Intron_gpII;ID=Solyc01r008605.1.1;Name=RF00029;e-value=3.806e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	2885353	2885504	83.41	-	.	Alias=U4;ID=Solyc01r008855.1.1;Name=RF00015;e-value=1.881e-15;rna_type=U4_snRNA
SL2.40ch01	ITAG_infernal	transcript	3040190	3040262	58.69	+	.	Alias=tRNA;ID=Solyc01r009035.1.1;Name=RF00005;e-value=2.139e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	3274222	3274293	72.46	+	.	Alias=tRNA;ID=Solyc01r009245.1.1;Name=RF00005;e-value=6.155e-15;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	3280560	3280639	56.87	-	.	Alias=tRNA;ID=Solyc01r009255.1.1;Name=RF00005;e-value=4.966e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	3482631	3482704	66.89	+	.	Alias=tRNA;ID=Solyc01r009305.1.1;Name=RF00005;e-value=6.564e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	3629373	3629455	54.11	-	.	Alias=tRNA;ID=Solyc01r009425.1.1;Name=RF00005;e-value=1.6e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	3657477	3657557	57.02	+	.	Alias=tRNA;ID=Solyc01r009435.1.1;Name=RF00005;e-value=4.641e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	3903981	3904052	66.91	-	.	Alias=tRNA;ID=Solyc01r009665.1.1;Name=RF00005;e-value=6.522e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	4410471	4410542	58.35	-	.	Alias=tRNA;ID=Solyc01r009875.1.1;Name=RF00005;e-value=2.473e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	4576869	4576942	58.58	+	.	Alias=tRNA;ID=Solyc01r009975.1.1;Name=RF00005;e-value=2.247e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	4917379	4917451	64.94	-	.	Alias=tRNA;ID=Solyc01r010145.1.1;Name=RF00005;e-value=1.506e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	5279545	5279617	65.92	-	.	Alias=tRNA;ID=Solyc01r010395.1.1;Name=RF00005;e-value=1.981e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	5767999	5768070	45.57	+	.	Alias=tRNA;ID=Solyc01r010735.1.1;Name=RF00005;e-value=5.672e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	5813051	5813153	113.97	-	.	Alias=U6;ID=Solyc01r010755.1.1;Name=RF00026;e-value=2.497e-21;rna_type=U6_snRNA
SL2.40ch01	ITAG_infernal	transcript	6214416	6214487	52.13	-	.	Alias=tRNA;ID=Solyc01r010815.1.1;Name=RF00005;e-value=3.484e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	7029360	7029430	43.39	-	.	Alias=tRNA;ID=Solyc01r011055.1.1;Name=RF00005;e-value=1.43e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	7041221	7041303	52.5	-	.	Alias=tRNA;ID=Solyc01r011065.1.1;Name=RF00005;e-value=5.961e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	7230813	7230910	21.84	+	.	Alias=mir-576;ID=Solyc01r114010.1.1;Name=RF00984;e-value=2.462e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	7391353	7391448	31.96	+	.	Alias=mir-576;ID=Solyc01r114020.1.1;Name=RF00984;e-value=2.108e-07;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	8143769	8143866	20.4	+	.	Alias=mir-576;ID=Solyc01r114030.1.1;Name=RF00984;e-value=1.271e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	9710402	9710473	51.48	-	.	Alias=tRNA;ID=Solyc01r012585.1.1;Name=RF00005;e-value=4.597e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	12603183	12603280	33.28	-	.	Alias=mir-576;ID=Solyc01r114040.1.1;Name=RF00984;e-value=3.624e-08;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	13376038	13376171	69.45	+	.	Alias=MIR1027;ID=Solyc01r114050.1.1;Name=RF00925;e-value=1.639e-15;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	14110320	14110447	83.2	+	.	Alias=MIR1027;ID=Solyc01r114060.1.1;Name=RF00925;e-value=1.203e-18;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	14110327	14110451	57.35	-	.	Alias=MIR1023;ID=Solyc01r114070.1.1;Name=RF01043;e-value=3.017e-11;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	14110343	14110440	46.8	+	.	Alias=mir-576;ID=Solyc01r114080.1.1;Name=RF00984;e-value=1.043e-11;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	14110345	14110437	47.37	-	.	Alias=mir-552;ID=Solyc01r114090.1.1;Name=RF00990;e-value=1.029e-11;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	14110346	14110440	36.56	+	.	Alias=mir-653;ID=Solyc01r114100.1.1;Name=RF00937;e-value=1.645e-09;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	21598666	21598744	56.82	+	.	Alias=Intron_gpII;ID=Solyc01r016715.1.1;Name=RF00029;e-value=3.185e-13;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	23563684	23563779	31.96	-	.	Alias=mir-576;ID=Solyc01r114110.1.1;Name=RF00984;e-value=1.468e-07;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	23563693	23563771	20.7	-	.	Alias=mir-653;ID=Solyc01r114120.1.1;Name=RF00937;e-value=4.564e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	23563694	23563783	42.13	+	.	Alias=mir-450;ID=Solyc01r114130.1.1;Name=RF00708;e-value=1.592e-08;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	23888732	23888806	38.94	-	.	Alias=Intron_gpII;ID=Solyc01r017195.1.1;Name=RF00029;e-value=1.698e-08;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	23893282	23893354	47.01	+	.	Alias=tRNA;ID=Solyc01r017215.1.1;Name=RF00005;e-value=2.213e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	23896536	23896606	42.73	+	.	Alias=tRNA;ID=Solyc01r017245.1.1;Name=RF00005;e-value=1.364e-08;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	23898820	23898898	47.31	+	.	Alias=Intron_gpII;ID=Solyc01r017246.1.1;Name=RF00029;e-value=2.916e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	23900191	23900271	40.96	+	.	Alias=tRNA;ID=Solyc01r017265.1.1;Name=RF00005;e-value=2.893e-08;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	23942759	23942830	35.63	+	.	Alias=Intron_gpII;ID=Solyc01r017375.1.1;Name=RF00029;e-value=8.433e-08;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	23950285	23950360	36.17	+	.	Alias=Intron_gpII;ID=Solyc01r017425.1.1;Name=RF00029;e-value=6.509e-08;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	23954942	23955058	37.76	-	.	Alias=5S_rRNA;ID=Solyc01r017485.1.1;Name=RF00001;e-value=1.715e-08;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	23955919	23956026	82.59	-	.	Alias=PK-G12rRNA;ID=Solyc01r017486.1.1;Name=RF01118;e-value=9.541e-17;rna_type=rRNA
SL2.40ch01	ITAG_infernal	transcript	25452441	25452538	34.66	-	.	Alias=mir-576;ID=Solyc01r114140.1.1;Name=RF00984;e-value=5.146e-08;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	25453629	25453705	41.33	-	.	Alias=Intron_gpII;ID=Solyc01r017805.1.1;Name=RF00029;e-value=5.862e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	27503335	27503432	24.16	-	.	Alias=mir-576;ID=Solyc01r114150.1.1;Name=RF00984;e-value=1.425e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	28500178	28500261	21.86	+	.	Alias=mir-576;ID=Solyc01r114160.1.1;Name=RF00984;e-value=1.56e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	30487744	30487846	112.32	-	.	Alias=U6;ID=Solyc01r020435.1.1;Name=RF00026;e-value=8.872e-21;rna_type=U6_snRNA
SL2.40ch01	ITAG_infernal	transcript	30939310	30939383	47.9	+	.	Alias=tRNA;ID=Solyc01r020455.1.1;Name=RF00005;e-value=4.358e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	30944556	30944630	48.47	+	.	Alias=tRNA;ID=Solyc01r020475.1.1;Name=RF00005;e-value=3.407e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	32649681	32649855	57.01	-	.	Alias=U2;ID=Solyc01r021715.1.1;Name=RF00004;e-value=3.965e-12;rna_type=U2_snRNA
SL2.40ch01	ITAG_infernal	transcript	32709133	32709327	95.55	-	.	Alias=U2;ID=Solyc01r021725.1.1;Name=RF00004;e-value=5.875e-18;rna_type=U2_snRNA
SL2.40ch01	ITAG_infernal	transcript	32713219	32713414	134.05	-	.	Alias=U2;ID=Solyc01r021726.1.1;Name=RF00004;e-value=4.08e-24;rna_type=U2_snRNA
SL2.40ch01	ITAG_infernal	transcript	32714794	32714968	58.39	-	.	Alias=U2;ID=Solyc01r021727.1.1;Name=RF00004;e-value=5.177e-12;rna_type=U2_snRNA
SL2.40ch01	ITAG_infernal	transcript	32729830	32729900	47.22	+	.	Alias=tRNA;ID=Solyc01r021728.1.1;Name=RF00005;e-value=5.806e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	32730748	32730901	83.91	-	.	Alias=U2;ID=Solyc01r021729.1.1;Name=RF00004;e-value=4.282e-16;rna_type=U2_snRNA
SL2.40ch01	ITAG_infernal	transcript	32818742	32818812	43.01	+	.	Alias=tRNA;ID=Solyc01r021745.1.1;Name=RF00005;e-value=1.684e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	32819199	32819361	68.08	-	.	Alias=U2;ID=Solyc01r021746.1.1;Name=RF00004;e-value=1.977e-13;rna_type=U2_snRNA
SL2.40ch01	ITAG_infernal	transcript	32838877	32839070	102.45	-	.	Alias=U2;ID=Solyc01r021747.1.1;Name=RF00004;e-value=6.295e-19;rna_type=U2_snRNA
SL2.40ch01	ITAG_infernal	transcript	32944996	32945114	69.19	-	.	Alias=5S_rRNA;ID=Solyc01r022755.1.1;Name=RF00001;e-value=3.185e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32945345	32945463	81.84	-	.	Alias=5S_rRNA;ID=Solyc01r022756.1.1;Name=RF00001;e-value=1.839e-14;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32945693	32945810	64.39	-	.	Alias=5S_rRNA;ID=Solyc01r022757.1.1;Name=RF00001;e-value=2.251e-11;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32947060	32947178	74.39	-	.	Alias=5S_rRNA;ID=Solyc01r022775.1.1;Name=RF00001;e-value=3.822e-13;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32947405	32947522	68.95	-	.	Alias=5S_rRNA;ID=Solyc01r022776.1.1;Name=RF00001;e-value=3.508e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32947739	32947859	51.82	-	.	Alias=5S_rRNA;ID=Solyc01r022777.1.1;Name=RF00001;e-value=3.776e-09;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32948090	32948207	77.76	-	.	Alias=5S_rRNA;ID=Solyc01r022778.1.1;Name=RF00001;e-value=9.715e-14;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32948436	32948554	88.51	-	.	Alias=5S_rRNA;ID=Solyc01r022779.1.1;Name=RF00001;e-value=1.213e-15;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32948781	32948897	69.51	-	.	Alias=5S_rRNA;ID=Solyc01r113720.1.1;Name=RF00001;e-value=2.793e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32949127	32949236	51.66	-	.	Alias=5S_rRNA;ID=Solyc01r113730.1.1;Name=RF00001;e-value=4.029e-09;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32949465	32949583	95.4	-	.	Alias=5S_rRNA;ID=Solyc01r113740.1.1;Name=RF00001;e-value=7.331e-17;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32949730	32949848	53.1	-	.	Alias=5S_rRNA;ID=Solyc01r113750.1.1;Name=RF00001;e-value=2.244e-09;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32950080	32950198	62.01	-	.	Alias=5S_rRNA;ID=Solyc01r113760.1.1;Name=RF00001;e-value=5.945e-11;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32950424	32950541	70.02	-	.	Alias=5S_rRNA;ID=Solyc01r113770.1.1;Name=RF00001;e-value=2.269e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32950770	32950888	83.67	-	.	Alias=5S_rRNA;ID=Solyc01r113780.1.1;Name=RF00001;e-value=8.718e-15;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32951344	32951465	55.47	-	.	Alias=5S_rRNA;ID=Solyc01r113790.1.1;Name=RF00001;e-value=8.523e-10;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32951694	32951816	76.36	-	.	Alias=5S_rRNA;ID=Solyc01r113800.1.1;Name=RF00001;e-value=1.716e-13;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32952048	32952166	56.03	-	.	Alias=5S_rRNA;ID=Solyc01r113810.1.1;Name=RF00001;e-value=6.773e-10;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32952399	32952520	53.15	-	.	Alias=5S_rRNA;ID=Solyc01r113820.1.1;Name=RF00001;e-value=2.194e-09;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32952750	32952868	77.14	-	.	Alias=5S_rRNA;ID=Solyc01r113830.1.1;Name=RF00001;e-value=1.251e-13;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32953098	32953216	70.77	-	.	Alias=5S_rRNA;ID=Solyc01r113840.1.1;Name=RF00001;e-value=1.673e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32953446	32953564	68.06	-	.	Alias=5S_rRNA;ID=Solyc01r113850.1.1;Name=RF00001;e-value=5.043e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32953794	32953912	73.64	-	.	Alias=5S_rRNA;ID=Solyc01r113860.1.1;Name=RF00001;e-value=5.209e-13;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32954141	32954258	84.18	-	.	Alias=5S_rRNA;ID=Solyc01r113870.1.1;Name=RF00001;e-value=7.083e-15;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32954490	32954608	55.29	-	.	Alias=5S_rRNA;ID=Solyc01r113880.1.1;Name=RF00001;e-value=9.182e-10;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32954838	32954956	82.52	-	.	Alias=5S_rRNA;ID=Solyc01r113890.1.1;Name=RF00001;e-value=1.398e-14;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32955186	32955303	70.77	-	.	Alias=5S_rRNA;ID=Solyc01r113900.1.1;Name=RF00001;e-value=1.675e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32955538	32955656	72.04	-	.	Alias=5S_rRNA;ID=Solyc01r113910.1.1;Name=RF00001;e-value=9.963e-13;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32955886	32956003	71.54	-	.	Alias=5S_rRNA;ID=Solyc01r113920.1.1;Name=RF00001;e-value=1.221e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32965246	32965362	69.09	-	.	Alias=5S_rRNA;ID=Solyc01r113930.1.1;Name=RF00001;e-value=1.684e-13;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32965384	32965502	48.66	-	.	Alias=5S_rRNA;ID=Solyc01r113940.1.1;Name=RF00001;e-value=6.935e-10;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32965732	32965848	62.5	-	.	Alias=5S_rRNA;ID=Solyc01r113950.1.1;Name=RF00001;e-value=2.466e-12;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32966078	32966195	71.02	-	.	Alias=5S_rRNA;ID=Solyc01r113960.1.1;Name=RF00001;e-value=7.657e-14;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32967435	32967553	78.03	-	.	Alias=5S_rRNA;ID=Solyc01r113970.1.1;Name=RF00001;e-value=2.73e-15;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32967784	32967902	70.28	-	.	Alias=5S_rRNA;ID=Solyc01r022785.1.1;Name=RF00001;e-value=6.418e-14;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32968134	32968252	88.03	-	.	Alias=5S_rRNA;ID=Solyc01r022786.1.1;Name=RF00001;e-value=4.638e-17;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	32987420	32987538	76.44	-	.	Alias=5S_rRNA;ID=Solyc01r022805.1.1;Name=RF00001;e-value=5.218e-15;rna_type=rRNA_5S
SL2.40ch01	ITAG_infernal	transcript	33448026	33448118	40.37	+	.	Alias=mir-552;ID=Solyc01r114170.1.1;Name=RF00990;e-value=3.924e-10;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	37089725	37089819	41.58	+	.	Alias=mir-552;ID=Solyc01r114180.1.1;Name=RF00990;e-value=7.136e-10;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	45850689	45850786	27.73	-	.	Alias=mir-576;ID=Solyc01r114190.1.1;Name=RF00984;e-value=1.648e-06;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	46779554	46779633	44.08	-	.	Alias=tRNA;ID=Solyc01r056365.1.1;Name=RF00005;e-value=1.137e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	47249385	47249453	41.14	+	.	Alias=tRNA;ID=Solyc01r056475.1.1;Name=RF00005;e-value=3.725e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	47429616	47429689	54.75	+	.	Alias=tRNA;ID=Solyc01r056515.1.1;Name=RF00005;e-value=1.143e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	50372508	50372588	52.81	-	.	Alias=tRNA;ID=Solyc01r056915.1.1;Name=RF00005;e-value=2.783e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	51718919	51719014	21.37	-	.	Alias=mir-576;ID=Solyc01r114200.1.1;Name=RF00984;e-value=3.134e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	52733783	52733896	61.12	-	.	Alias=U5;ID=Solyc01r057215.1.1;Name=RF00020;e-value=4.209e-13;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	52734774	52734933	74.02	-	.	Alias=U1;ID=Solyc01r057216.1.1;Name=RF00003;e-value=3.556e-15;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	53884347	53884439	25.92	-	.	Alias=mir-576;ID=Solyc01r114210.1.1;Name=RF00984;e-value=3.957e-06;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	54170317	54170429	40.2	+	.	Alias=Intron_gpII;ID=Solyc01r057525.1.1;Name=RF00029;e-value=1.979e-09;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	54348818	54348910	27.88	+	.	Alias=mir-653;ID=Solyc01r114220.1.1;Name=RF00937;e-value=1.86e-06;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	55725447	55725544	20.79	-	.	Alias=mir-576;ID=Solyc01r114230.1.1;Name=RF00984;e-value=4.14e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	58347196	58347283	42.13	+	.	Alias=mir-450;ID=Solyc01r114240.1.1;Name=RF00708;e-value=3.256e-08;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	58347209	58347289	21.02	-	.	Alias=mir-576;ID=Solyc01r114250.1.1;Name=RF00984;e-value=6.174e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	59682590	59682745	118.46	-	.	Alias=5_8S_rRNA;ID=Solyc01r058495.1.1;Name=RF00002;e-value=6.788e-23;rna_type=rRNA_5.8S
SL2.40ch01	ITAG_infernal	transcript	59921481	59921574	20.77	-	.	Alias=mir-576;ID=Solyc01r114260.1.1;Name=RF00984;e-value=7.432e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	62040288	62040360	41.5	+	.	Alias=tRNA;ID=Solyc01r059895.1.1;Name=RF00005;e-value=3.2e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	62148140	62148237	26.31	-	.	Alias=mir-576;ID=Solyc01r114270.1.1;Name=RF00984;e-value=7.327e-06;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	62786462	62786590	22.25	+	.	Alias=MIR475;ID=Solyc01r114280.1.1;Name=RF00721;e-value=1.756e-05;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	62786885	62787009	44.43	+	.	Alias=MIR475;ID=Solyc01r114290.1.1;Name=RF00721;e-value=2.418e-09;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	62840993	62841293	155.75	-	.	Alias=SRP_euk_arch;ID=Solyc01r060285.1.1;Name=RF00017;e-value=3.198e-23;rna_type=SRP_RNA
SL2.40ch01	ITAG_infernal	transcript	62881752	62882052	181.63	-	.	Alias=SRP_euk_arch;ID=Solyc01r060305.1.1;Name=RF00017;e-value=1.365e-26;rna_type=SRP_RNA
SL2.40ch01	ITAG_infernal	transcript	62916004	62916303	120.69	-	.	Alias=SRP_euk_arch;ID=Solyc01r060306.1.1;Name=RF00017;e-value=1.172e-18;rna_type=SRP_RNA
SL2.40ch01	ITAG_infernal	transcript	62965859	62966159	147.42	-	.	Alias=SRP_euk_arch;ID=Solyc01r060315.1.1;Name=RF00017;e-value=3.886e-22;rna_type=SRP_RNA
SL2.40ch01	ITAG_infernal	transcript	63029133	63029421	76.82	-	.	Alias=SRP_euk_arch;ID=Solyc01r060335.1.1;Name=RF00017;e-value=6.033e-13;rna_type=SRP_RNA
SL2.40ch01	ITAG_infernal	transcript	63174491	63174564	34.71	+	.	Alias=Intron_gpII;ID=Solyc01r060365.1.1;Name=RF00029;e-value=5.911e-08;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	63238999	63239071	37.9	+	.	Alias=Intron_gpII;ID=Solyc01r060366.1.1;Name=RF00029;e-value=1.258e-08;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	63276696	63276769	47.68	+	.	Alias=Intron_gpII;ID=Solyc01r060367.1.1;Name=RF00029;e-value=1.091e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	63337698	63337771	35.58	+	.	Alias=Intron_gpII;ID=Solyc01r060368.1.1;Name=RF00029;e-value=3.874e-08;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	63669554	63669627	41.09	+	.	Alias=Intron_gpII;ID=Solyc01r060375.1.1;Name=RF00029;e-value=6.601e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	64342927	64343005	43.57	+	.	Alias=Intron_gpII;ID=Solyc01r065585.1.1;Name=RF00029;e-value=6.159e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	64344983	64345054	41.67	+	.	Alias=tRNA;ID=Solyc01r065586.1.1;Name=RF00005;e-value=6.147e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	64356426	64356497	43.47	-	.	Alias=tRNA;ID=Solyc01r065605.1.1;Name=RF00005;e-value=2.859e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	65146227	65146299	60.64	+	.	Alias=tRNA;ID=Solyc01r065975.1.1;Name=RF00005;e-value=9.366e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	66991087	66991173	34.86	-	.	Alias=mir-653;ID=Solyc01r114300.1.1;Name=RF00937;e-value=7.512e-09;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	66991090	66991174	41.44	+	.	Alias=mir-552;ID=Solyc01r114310.1.1;Name=RF00990;e-value=2.132e-10;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	66991090	66991178	41.28	-	.	Alias=mir-450;ID=Solyc01r114320.1.1;Name=RF00708;e-value=4.657e-09;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	69155802	69155873	50.78	+	.	Alias=tRNA;ID=Solyc01r067785.1.1;Name=RF00005;e-value=6.177e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	69480937	69481009	60.61	+	.	Alias=tRNA;ID=Solyc01r068065.1.1;Name=RF00005;e-value=9.46e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	69492167	69492266	21.36	-	.	Alias=mir-576;ID=Solyc01r114330.1.1;Name=RF00984;e-value=3.181e-06;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	70995588	70995679	45.82	+	.	Alias=mir-576;ID=Solyc01r114340.1.1;Name=RF00984;e-value=1.048e-10;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	71072861	71072982	58.98	+	.	Alias=U6atac;ID=Solyc01r079515.1.1;Name=RF00619;e-value=1.138e-12;rna_type=U6atac_snRNA
SL2.40ch01	ITAG_infernal	transcript	71979651	71979756	85.97	+	.	Alias=snoZ279_R105_R108;ID=Solyc01r080245.1.1;Name=RF00304;e-value=1.181e-17;rna_type=snoRNA
SL2.40ch01	ITAG_infernal	transcript	71981969	71982074	84.02	+	.	Alias=snoZ279_R105_R108;ID=Solyc01r080246.1.1;Name=RF00304;e-value=2.781e-17;rna_type=snoRNA
SL2.40ch01	ITAG_infernal	transcript	72405845	72405928	44.87	+	.	Alias=tRNA;ID=Solyc01r080645.1.1;Name=RF00005;e-value=8.158e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	72472502	72472625	23.07	+	.	Alias=MIR475;ID=Solyc01r114350.1.1;Name=RF00721;e-value=6.314e-06;rna_type=miRNA
SL2.40ch01	ITAG_infernal	transcript	72540083	72540381	104.41	-	.	Alias=SRP_euk_arch;ID=Solyc01r080785.1.1;Name=RF00017;e-value=2.003e-17;rna_type=SRP_RNA
SL2.40ch01	ITAG_infernal	transcript	73820366	73820447	53.18	-	.	Alias=tRNA;ID=Solyc01r087135.1.1;Name=RF00005;e-value=2.382e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	73980964	73981036	65.92	+	.	Alias=tRNA;ID=Solyc01r087285.1.1;Name=RF00005;e-value=9.903e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	74374736	74374807	68.78	-	.	Alias=tRNA;ID=Solyc01r087685.1.1;Name=RF00005;e-value=6.062e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	74384766	74384837	71.04	-	.	Alias=tRNA;ID=Solyc01r087715.1.1;Name=RF00005;e-value=2.327e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	75278726	75278797	55.78	-	.	Alias=tRNA;ID=Solyc01r088755.1.1;Name=RF00005;e-value=7.399e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	75388031	75388176	58.44	-	.	Alias=5_8S_rRNA;ID=Solyc01r089915.1.1;Name=RF00002;e-value=8.498e-13;rna_type=rRNA_5.8S
SL2.40ch01	ITAG_infernal	transcript	76222602	76222674	51.59	+	.	Alias=tRNA;ID=Solyc01r090785.1.1;Name=RF00005;e-value=4.377e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	76721050	76721122	64.05	-	.	Alias=tRNA;ID=Solyc01r091315.1.1;Name=RF00005;e-value=6.874e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	76922204	76922277	59.88	+	.	Alias=tRNA;ID=Solyc01r091545.1.1;Name=RF00005;e-value=4.047e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	77186886	77186966	57.02	+	.	Alias=tRNA;ID=Solyc01r091885.1.1;Name=RF00005;e-value=1.363e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	77653028	77653102	42.36	+	.	Alias=Intron_gpII;ID=Solyc01r094375.1.1;Name=RF00029;e-value=6.95e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	77715764	77715837	60.71	-	.	Alias=tRNA;ID=Solyc01r094505.1.1;Name=RF00005;e-value=6.61e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	77959967	77960052	41.71	+	.	Alias=tRNA;ID=Solyc01r094795.1.1;Name=RF00005;e-value=2.125e-08;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	78014364	78014437	60.08	-	.	Alias=tRNA;ID=Solyc01r094845.1.1;Name=RF00005;e-value=8.638e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	78338139	78338222	56.15	+	.	Alias=tRNA;ID=Solyc01r095235.1.1;Name=RF00005;e-value=4.593e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	78368283	78368354	68.78	-	.	Alias=tRNA;ID=Solyc01r095295.1.1;Name=RF00005;e-value=2.132e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	78899314	78899386	47.38	+	.	Alias=tRNA;ID=Solyc01r096045.1.1;Name=RF00005;e-value=1.906e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	79027616	79027689	59.88	-	.	Alias=tRNA;ID=Solyc01r096165.1.1;Name=RF00005;e-value=9.385e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	79150585	79150664	52.17	-	.	Alias=tRNA;ID=Solyc01r096315.1.1;Name=RF00005;e-value=3.66e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	79836665	79836737	65.92	+	.	Alias=tRNA;ID=Solyc01r097115.1.1;Name=RF00005;e-value=9.903e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	80107961	80108076	94.21	+	.	Alias=U5;ID=Solyc01r097525.1.1;Name=RF00020;e-value=1.285e-18;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	80108963	80109074	76.94	+	.	Alias=U5;ID=Solyc01r097535.1.1;Name=RF00020;e-value=1.718e-15;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	80109516	80109632	86.17	+	.	Alias=U5;ID=Solyc01r097536.1.1;Name=RF00020;e-value=3.665e-17;rna_type=U5_snRNA
SL2.40ch01	ITAG_infernal	transcript	81458498	81458568	59.34	+	.	Alias=tRNA;ID=Solyc01r099445.1.1;Name=RF00005;e-value=1.982e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	81489040	81489110	44.91	+	.	Alias=tRNA;ID=Solyc01r099495.1.1;Name=RF00005;e-value=9.159e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	81494414	81494484	59.34	+	.	Alias=tRNA;ID=Solyc01r099505.1.1;Name=RF00005;e-value=1.982e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	81646849	81646920	62.9	-	.	Alias=tRNA;ID=Solyc01r099725.1.1;Name=RF00005;e-value=4.374e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	82002324	82002397	68.58	+	.	Alias=tRNA;ID=Solyc01r100175.1.1;Name=RF00005;e-value=3.912e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	82155687	82155758	54.54	-	.	Alias=tRNA;ID=Solyc01r100325.1.1;Name=RF00005;e-value=1.527e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	82166851	82166981	45.04	+	.	Alias=Intron_gpII;ID=Solyc01r100335.1.1;Name=RF00029;e-value=1.046e-10;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	82172486	82172558	57.91	-	.	Alias=tRNA;ID=Solyc01r100345.1.1;Name=RF00005;e-value=3.636e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	82174840	82174961	78.21	-	.	Alias=snoR83;ID=Solyc01r100346.1.1;Name=RF01227;e-value=1.632e-17;rna_type=snoRNA
SL2.40ch01	ITAG_infernal	transcript	82519487	82519559	60.61	-	.	Alias=tRNA;ID=Solyc01r100845.1.1;Name=RF00005;e-value=1.154e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	82672487	82672559	60.61	+	.	Alias=tRNA;ID=Solyc01r101055.1.1;Name=RF00005;e-value=1.154e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	82876921	82876993	65.92	-	.	Alias=tRNA;ID=Solyc01r102295.1.1;Name=RF00005;e-value=1.208e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	82896654	82896814	147.46	-	.	Alias=U1;ID=Solyc01r102315.1.1;Name=RF00003;e-value=3.125e-27;rna_type=U1_snRNA
SL2.40ch01	ITAG_infernal	transcript	82898374	82898593	135.3	+	.	Alias=U3;ID=Solyc01r102316.1.1;Name=RF00012;e-value=5.686e-21;rna_type=snRNA
SL2.40ch01	ITAG_infernal	transcript	83912063	83912165	112.5	+	.	Alias=U6;ID=Solyc01r103575.1.1;Name=RF00026;e-value=8.564e-21;rna_type=U6_snRNA
SL2.40ch01	ITAG_infernal	transcript	84159211	84159284	66.11	-	.	Alias=tRNA;ID=Solyc01r103895.1.1;Name=RF00005;e-value=1.889e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	84486358	84486430	57.91	+	.	Alias=tRNA;ID=Solyc01r104305.1.1;Name=RF00005;e-value=2.762e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	84496273	84496344	56.67	-	.	Alias=tRNA;ID=Solyc01r104325.1.1;Name=RF00005;e-value=4.693e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	84517994	84518065	60.36	-	.	Alias=tRNA;ID=Solyc01r104355.1.1;Name=RF00005;e-value=9.755e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	84518345	84518416	62.9	-	.	Alias=tRNA;ID=Solyc01r104356.1.1;Name=RF00005;e-value=3.322e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	84519750	84519821	62.9	-	.	Alias=tRNA;ID=Solyc01r104357.1.1;Name=RF00005;e-value=3.322e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	84521411	84521482	62.9	-	.	Alias=tRNA;ID=Solyc01r104358.1.1;Name=RF00005;e-value=3.322e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	84539974	84540046	64.05	+	.	Alias=tRNA;ID=Solyc01r104405.1.1;Name=RF00005;e-value=2.033e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	85125658	85125730	64.94	+	.	Alias=tRNA;ID=Solyc01r105095.1.1;Name=RF00005;e-value=1.396e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	85519650	85519733	50.48	+	.	Alias=tRNA;ID=Solyc01r105665.1.1;Name=RF00005;e-value=1.209e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86118905	86118978	60.57	+	.	Alias=tRNA;ID=Solyc01r106465.1.1;Name=RF00005;e-value=1.658e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86412401	86412472	53.04	+	.	Alias=tRNA;ID=Solyc01r106955.1.1;Name=RF00005;e-value=4.069e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86413197	86413268	50.68	+	.	Alias=tRNA;ID=Solyc01r106956.1.1;Name=RF00005;e-value=1.111e-09;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86413569	86413640	57.37	+	.	Alias=tRNA;ID=Solyc01r106957.1.1;Name=RF00005;e-value=6.461e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86415417	86415488	54.63	+	.	Alias=tRNA;ID=Solyc01r106958.1.1;Name=RF00005;e-value=2.069e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86415658	86415729	57.04	+	.	Alias=tRNA;ID=Solyc01r106959.1.1;Name=RF00005;e-value=7.445e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86417232	86417303	52.23	+	.	Alias=tRNA;ID=Solyc01r106965.1.1;Name=RF00005;e-value=5.755e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86417587	86417658	51.64	+	.	Alias=tRNA;ID=Solyc01r106966.1.1;Name=RF00005;e-value=7.378e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86418261	86418332	57.04	+	.	Alias=tRNA;ID=Solyc01r106967.1.1;Name=RF00005;e-value=7.445e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86420227	86420298	57.04	+	.	Alias=tRNA;ID=Solyc01r106968.1.1;Name=RF00005;e-value=7.445e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86420599	86420670	57.04	+	.	Alias=tRNA;ID=Solyc01r106969.1.1;Name=RF00005;e-value=7.445e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86421190	86421261	54.73	+	.	Alias=tRNA;ID=Solyc01r113980.1.1;Name=RF00005;e-value=1.983e-10;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86430464	86430535	68.78	-	.	Alias=tRNA;ID=Solyc01r113990.1.1;Name=RF00005;e-value=5.045e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86430968	86431039	57.04	+	.	Alias=tRNA;ID=Solyc01r114000.1.1;Name=RF00005;e-value=7.445e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86437193	86437264	68.78	-	.	Alias=tRNA;ID=Solyc01r106975.1.1;Name=RF00005;e-value=5.045e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	86501800	86501873	66.89	+	.	Alias=tRNA;ID=Solyc01r107115.1.1;Name=RF00005;e-value=6.564e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	87433427	87433500	66.93	+	.	Alias=tRNA;ID=Solyc01r108315.1.1;Name=RF00005;e-value=6.466e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	87900296	87900368	63.19	-	.	Alias=tRNA;ID=Solyc01r109085.1.1;Name=RF00005;e-value=1.287e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	87974211	87974287	37.55	+	.	Alias=Intron_gpII;ID=Solyc01r109165.1.1;Name=RF00029;e-value=3.677e-09;rna_type=group_II_intron
SL2.40ch01	ITAG_infernal	transcript	88123643	88123714	68.78	+	.	Alias=tRNA;ID=Solyc01r109395.1.1;Name=RF00005;e-value=1.193e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	88292984	88293055	66.47	+	.	Alias=tRNA;ID=Solyc01r109645.1.1;Name=RF00005;e-value=3.193e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	88751202	88751273	54.54	-	.	Alias=tRNA;ID=Solyc01r110265.1.1;Name=RF00005;e-value=2.586e-11;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	88815175	88815246	65.76	+	.	Alias=tRNA;ID=Solyc01r110325.1.1;Name=RF00005;e-value=2.189e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	88840972	88841044	63.51	+	.	Alias=tRNA;ID=Solyc01r110345.1.1;Name=RF00005;e-value=2.756e-13;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	88944247	88944330	56.1	+	.	Alias=tRNA;ID=Solyc01r110485.1.1;Name=RF00005;e-value=6.882e-12;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	89125498	89125569	66.91	+	.	Alias=tRNA;ID=Solyc01r110835.1.1;Name=RF00005;e-value=6.522e-14;rna_type=tRNA
SL2.40ch01	ITAG_infernal	transcript	89737229	89737335	35.79	+	.	Alias=snoR103;ID=Solyc01r111805.1.1;Name=RF01213;e-value=1.499e-08;rna_type=snoRNA
SL2.40ch01	ITAG_infernal	transcript	89758031	89758128	22.07	+	.	Alias=mir-576;ID=Solyc01r114360.1.1;Name=RF00984;e-value=2.924e-05;rna_type=miRNA
);


ok( my $gff = Bio::GenomeUpdate::GFF->new(),'create GFF obj');
ok( $gff->parse_gff($gff_file),'parse GFF obj');
print STDERR "Files read..\n";
mem_used();

#update coordinates, calls GFFRearrange object
ok( my %coords = $gff->get_reordered_coordinates($agp_orig,$agp_fish),'get_reordered_coordinates');
ok( my %flips = $gff->get_flipped_coordinates($agp_orig,$agp_fish),'get_flipped_coordinates');
print STDERR "Hashes populated..\n";
mem_used();

print STDERR "\nflips: ",scalar keys %flips,"\n";

#remap GFF
ok( $gff->remap_coordinates(\%coords,\%flips),'remap_coordinates');
print STDERR "Coords remapped..\n";
mem_used();
	
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
);
ok( my $gff_fish = $gff->get_formatted_gff(),'get_formatted_gff');
is( $gff_fish, $compare_str, 'GFF remapping is as expected');
print STDERR "GFF written..\n";
mem_used();

#----------------------------------------------------------------------------
	
sub mem_used{
	my ($i,$t); 
	$t = new Proc::ProcessTable;
	foreach my $got ( @{$t->table} ) {
		next if not $got->pid eq $$; $i=$got->size;
	}
	print STDERR "Process id=",$$,"\n"; print "Memory used(MB)=", $i/1024/1024, "\n";
}