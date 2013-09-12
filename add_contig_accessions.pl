#!/usr/bin/perl

=head1

add_contig_accessions.pl

=head1 SYNOPSIS

    add_contig_accessions.pl -a [AGP file] -n [Allcontig AGP file] -s [Scaffold from contig AGP file]-m [Scaffold names file]

=head1 COMMAND-LINE OPTIONS

 -a  Current (build 2.40) AGP file (mandatory)
 -s  Scaffold from contig AGP from NCBI FTP. Scaffold ID to contig ID (mandatory)
 -n  All Contig AGP file from NCBI FTP. Scaffold acc# to contig acc# (mandatory)
 -m  Scaffold ID to Refseq and Genbank ID mapping file (mandatory) 
 -o  Updated output AGP file 
 -h  Help

=cut

use strict;
use Bio::GenomeUpdate::AGP;
use File::Slurp;
use Getopt::Std;

our ( $opt_a, $opt_s, $opt_n, $opt_m, $opt_o, $opt_h );
getopts('a:s:n:m:o:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_a || !$opt_s || !$opt_n || !$opt_m ) {
	print
	"\nAGP, Scaffold from contig AGP, Allcontig AGP and Mapping file are required. See help below\n\n\n";
	help();
}

#get input files
my $input_curr_agp_file = $opt_a;
my $input_curr_agp      = read_file($input_curr_agp_file)
  or die "Could not open current AGP input file: $input_curr_agp_file \n";
my $input_scaffold_agp_file = $opt_s;
my $input_scaffold_agp      = read_file($input_scaffold_agp_file)
  or die "Could not open NCBI scaffold AGP input file: $input_scaffold_agp_file \n";
my $input_allcontig_agp_file = $opt_n;
my $input_allcontig_agp      = read_file($input_allcontig_agp_file)
  or die "Could not open NCBI allcontig AGP input file: $input_allcontig_agp_file \n";
my $input_mapping_file = $opt_m;
my $output_agp_file;
if (defined($opt_o)){
	$output_agp_file = $opt_o;
}else{
	$output_agp_file = "updated\.$input_curr_agp_file";
}

#objects for AGPs and mapping file
my $curr_agp = Bio::GenomeUpdate::AGP->new();
$curr_agp->parse_agp($input_curr_agp);
my $scaffold_agp = Bio::GenomeUpdate::AGP->new();
$scaffold_agp->parse_agp($input_scaffold_agp);
my $allcontig_agp = Bio::GenomeUpdate::AGP->new();
$allcontig_agp->parse_agp($input_allcontig_agp);
my $output_agp = Bio::GenomeUpdate::AGP->new();

# prepare mapping data
my %scaffoldID_accession = get_scaffold_accession_mapping($input_mapping_file);
my %contigID_scaffoldID = get_component_object_mapping($scaffold_agp);
my %scaffoldACC_contigACC = get_object_components_mapping($allcontig_agp);




#----------------------------------------------------------------------------

=head2 Methods

=over

=item C<get_scaffold_accession_mapping ( $input_mapping_file )>

Gets the hash with key/value pair as Scaffold IDs and Refseq Accession numbers.

=cut
sub get_scaffold_accession_mapping{
	my $input_mapping_file = shift;
	die "not able to open $input_mapping_file\n\n" unless(open(INMAP,"<$input_mapping_file"));
	
	my %scaffoldID_accession;
	while (my $line = <INMAP>){
		my @temp = split ("\t", $line);
		#record only Refseq accessions for Scaffold IDs
		$scaffoldID_accession{$temp[1]} = $temp[2];
	}
	close(INMAP);
	return %scaffoldID_accession;
}

=item C<get_object_components_mapping_hash ( $input_agp )>

Gets the hash with key/value pair as AGP object and array of components in the object.

=cut

sub get_object_components_mapping_hash{
	my $agp = shift;
	my %object_components;
	while ( my $agp_line = $scaffold_agp->get_next_agp_line()){
		#no component info in gap lines
		if ($agp_line->get_line_type() eq "gap"){ next;}
		
		my $agp_line_object_id = $agp_line->get_object_being_assembled();
		my $agp_line_component_id = $agp_line->get_component_id();
		push ( @{$object_components{$agp_line_object_id}}, $agp_line_component_id)
		
		#my @temp; this also works
#		if (!exists $object_components{$agp_line_object_id}){
#			#push (@temp, $agp_line_component_id);
#			push ( @{$object_components{$agp_line_object_id}}, $agp_line_component_id)
#		}
#		else{
#			@temp = $object_components{$agp_line_object_id};
#			push (@temp, $agp_line_component_id);
#		}
#		$object_components{$agp_line_object_id} = \@temp;
	}
	return %object_components;
}

=item C<get_component_object_mapping_hash ( $input_agp )>

Gets the hash with key/value pair as component and parent object.

=cut

sub get_component_object_mapping_hash{
	my $agp = shift;
	my %component_object;
	while ( my $agp_line = $scaffold_agp->get_next_agp_line()){
		#no component info in gap lines
		if ($agp_line->get_line_type() eq "gap"){ next;}
		$component_object{$agp_line->get_component_id()} = $agp_line->get_object_being_assembled();
	}
	return %component_object;
}

=item C<get_contig_accession_mapping ( $curr_agp, %scaffoldID_accession, %contigID_scaffoldID, %scaffoldACC_contigACC )>

Gets the hash with key/value pair as Contig IDs and Accession numbers.

=cut
sub get_contig_accession_mapping{
	my ($curr_agp, %scaffoldID_accession, %contigID_scaffoldID, %scaffoldACC_contigACC ) = @_;
	
	while ( my $agp_line = $curr_agp->get_next_agp_line()){
		#no component info in gap lines
		if ($agp_line->get_line_type() eq "gap"){ next;}
		my $contigID = $agp_line->get_component_id();
		my $scaffoldID = $contigID_scaffoldID{$contigID};
		my $scaffoldACC = $scaffoldID_accession{$scaffoldID};
		my @contigACC_arr = $scaffoldACC_contigACC{$scaffoldACC};
		
		#return $contigACC_arr[] CHK IF LINE NUMBER == COMPONENT/PART NUMBER
		
	}	

}

sub help {
	print STDERR <<EOF;
  $0:

    add_contig_accessions.pl

    Description:

     This script add the Accession numbers of Contigs to an Accessioned Golden Path (AGP) file.

    Usage:
    
     add_contig_accessions.pl -a [AGP file] -n [Allcontig AGP file] -s [Scaffold from contig AGP file]-m [Scaffold names file]

    Flags:

        -a  current (build 2.40) AGP file (mandatory)
        -s  Scaffold from contig AGP from NCBI FTP. Scaffold ID to contig ID (mandatory)
        -n  All Contig AGP file from NCBI FTP. Scaffold acc# to contig acc# (mandatory)
        -m  Scaffold ID to Refseq and Genbank ID mapping file (mandatory) 
        -o  Updated output AGP file
        -h  Help

EOF
	exit(1);
}
=back

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <ss2489@cornell.edu>

=cut

__END__

# Current 2.40 AGP
SL2.40ch00	21801621	21801720	8946	U	100	contig	no
SL2.40ch00	21801721	21803721	8947	W	SL2.40ct09637	1	2001	+
SL2.40ch00	21803722	21803821	8948	U	100	contig	no
SL2.40ch00	21803822	21805821	8949	W	SL2.40ct16352	1	2000	+
SL2.40ch01	1	3610	1	W	SL2.40ct13037	1	3610	-
SL2.40ch01	3611	3863	2	N	253	fragment	yes	
SL2.40ch01	3864	12848	3	W	SL2.40ct13036	1	8985	-
SL2.40ch01	12849	12868	4	N	20	fragment	yes	
SL2.40ch01	12869	15910	5	W	SL2.40ct13035	1	3042	-

# Scaffold AGP from NCBI, all contigs on POS strand
SL2.40sc03923	7523342	7524540	154	N	1199	fragment	yes	
SL2.40sc03923	7524541	7525668	155	W	SL2.40ct09626	1	1128	+
SL2.40sc03923	7525669	7527534	156	N	1866	fragment	yes	
SL2.40sc03923	7527535	7528599	157	W	SL2.40ct09627	1	1065	+
SL2.40sc03771	1	2611	1	W	SL2.40ct06816	1	2611	+
SL2.40sc03771	2612	4591	2	N	1980	fragment	yes	
SL2.40sc03771	4592	18382	3	W	SL2.40ct06817	1	13791	+
SL2.40sc03771	18383	18402	4	N	20	fragment	yes	
SL2.40sc03771	18403	20331	5	W	SL2.40ct06818	1	1929	+
SL2.40sc03771	20332	23613	6	N	3282	fragment	yes	
SL2.40sc03771	23614	45878	7	W	SL2.40ct06819	1	22265	+

# Allcontig AGP from NCBI
NW_004194292.1	32971688	32974729	2359	W	AEKE02013045.1	1	3042	+
NW_004194292.1	32974730	32974749	2360	N	20	scaffold	yes	paired-ends
NW_004194292.1	32974750	32983734	2361	W	AEKE02013046.1	1	8985	+
NW_004194292.1	32983735	32983987	2362	N	253	scaffold	yes	paired-ends
NW_004194292.1	32983988	32987597	2363	W	AEKE02013047.1	1	3610	+
NW_004194293.1	1	1952	1	W	AEKE02001187.1	1	1952	+
NW_004194293.1	1953	2521	2	N	569	scaffold	yes	paired-ends
NW_004194293.1	2522	4629	3	W	AEKE02001188.1	1	2108	+
NW_004194293.1	4630	5831	4	N	1202	scaffold	yes	paired-ends
NW_004194293.1	5832	11217	5	W	AEKE02001189.1	1	5386	+
NW_004194293.1	11218	12963	6	N	1746	scaffold	yes	paired-ends

# Mapping file from NCBI
#Assembly	Genome Center name	RefSeq Accession.version	GenBank Accession.version	NCBI name
SL2.40	SL2.40sc04133	NW_004194292.1	GL758110.1	GPS_001221024.1
SL2.40	SL2.40sc03666	NW_004194293.1	GL758111.1	GPS_001221025.1
SL2.40	SL2.40sc04191	NW_004194294.1	GL758112.1	GPS_001221026.1
SL2.40	SL2.40sc03594	NW_004194295.1	GL758113.1	GPS_001221027.1
SL2.40	SL2.40sc05010	NW_004194296.1	GL758114.1	GPS_001221028.1
SL2.40	SL2.40sc05941	NW_004194297.1	GL758115.1	GPS_001221029.1
