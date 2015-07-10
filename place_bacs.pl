#!/usr/bin/perl

=head1

place_bacs.pl

=head1 SYNOPSIS

    place_bacs.pl -t [TPF file] -s [scaffold AGP file] -c [chromosome AGP file] -n [chromosome number] -l [file with BAC names and coordinates] -a  ACE [file of assembled BACs]

=head1 COMMAND-LINE OPTIONS

 -t  TPF file (required)
 -s  scaffold AGP file (required)
 -c  chromosome AGP file (required)
 -n  chromosome number. Valid values are 1-12 or Un for tomato (required)
 -l  Tab-delimited file containing BAC names and coordinates from align_BACends_group_coords.pl or group_coords.pl (required)
 -a  ACE file of assembled BACs generated by phrap. see scripts/run_phrap.sh (required) 
 -h  Help

=head1 TODO

  See github issues.
  
=cut

use strict;
use warnings;

use Getopt::Std;
use File::Slurp;
use Bio::Assembly::IO;
use Bio::Assembly::Contig;
use Bio::Assembly::Scaffold;
use Bio::SeqFeature::Generic;
use Bio::Assembly::Singlet;
use Bio::Seq;
use Bio::LocatableSeq;

use Bio::GenomeUpdate::TPF;
use Bio::GenomeUpdate::TP;
use Bio::GenomeUpdate::SP;
use Bio::GenomeUpdate::AGP;

use Data::Dumper;#for debugging

our ( $opt_t, $opt_l, $opt_s, $opt_c, $opt_h, $opt_a, $opt_n );
getopts('t:s:c:l:h:a:n:');#h or boolean param should be the last, any params accepting strings should be followed by a : 
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_t || !$opt_l || !$opt_s || !$opt_c || !$opt_a || !$opt_n ) {
	print STDERR
	  "\nTPF, AGP, group_coords BAC name/coordinates, chr number and assembled ACE filenames are required.\n\n\n";
	help();
}

my $tpf_input_file = $opt_t;
my $input_tpf      = read_file($tpf_input_file)
  or die "Could not open TPF input file: $tpf_input_file\n";
my $scaffold_agp_input_file = $opt_s;
my $scaffold_input_agp      = read_file($scaffold_agp_input_file)
  or die "Could not open scaffold AGP input file: $scaffold_agp_input_file\n";
my $chromosome_agp_input_file = $opt_c;
my $chromosome_input_agp      = read_file($chromosome_agp_input_file)
  or die
  "Could not open chromosome AGP input file: $chromosome_agp_input_file\n";
my $group_coords_input_file = $opt_l;
my $input_group_coords     = read_file($group_coords_input_file)
  or die "Could not open BAC name and coordinates file: $group_coords_input_file\n";
my $assembly = Bio::Assembly::IO->new( -file => $opt_a, -format => 'ace'); 

my $chromosome;
if ($opt_n eq '1' || $opt_n eq '2' || $opt_n eq '3' || $opt_n eq '4' || $opt_n eq '5' || $opt_n eq '6' || $opt_n eq '7' || $opt_n eq '8' || $opt_n eq '9' || $opt_n eq '10' || $opt_n eq '11' || $opt_n eq '12' || $opt_n eq "Un"){
	$chromosome = $opt_n;
}
else{
	help();
}


##process assembly
my $scaffold = $assembly->next_assembly();
scaffold_summary($scaffold);
my %scaffold_component_contigs;
my %scaffold_component_contig_directions;
my $contig_ctr = 0;
foreach my $contig ($scaffold->all_contigs()){
##	print Dumper ($contig);
	my ($contig_component_bacs_ref,$contig_component_directions_ref) = contig_component_id_direction($contig);
	$scaffold_component_contigs{$contig->id()}= $contig_component_bacs_ref;  
	$scaffold_component_contig_directions{$contig->id()}  = $contig_component_directions_ref;
	$contig_ctr++;
}
print STDERR "\n$contig_ctr contigs processed from $opt_a\n";

#print STDERR "\%scaffold_component_contigs\n";
#print STDERR Dumper \%scaffold_component_contigs;
#print STDERR "\%scaffold_component_contig_directions\n";
#print STDERR Dumper \%scaffold_component_contig_directions;


my $tpf            = Bio::GenomeUpdate::TPF->new();
my $scaffold_agp   = Bio::GenomeUpdate::AGP->new();
my $chromosome_agp = Bio::GenomeUpdate::AGP->new();
my @group_coords_lines;
my @bacs;

my $ordered_tpf;
my $ordered_tpf_with_bacs_inserted_in_sequences_and_gaps;
my $sp_with_bacs_inserted_in_sequences_and_gaps;
my $tp_with_bacs_inserted_in_sequences_and_gaps;
my %chromosome_agp_offset;
my %scaffold_agp_coords;

$tpf->parse_tpf($input_tpf);
$scaffold_agp->parse_agp($scaffold_input_agp);
$chromosome_agp->parse_agp($chromosome_input_agp);

my %scaffold_agp_lines;
if ( $scaffold_agp->has_agp_lines() ) {
	%scaffold_agp_lines = %{ $scaffold_agp->get_agp_lines() };
}

my %chromosome_agp_lines;
if ( $chromosome_agp->has_agp_lines() ) {
	%chromosome_agp_lines = %{ $chromosome_agp->get_agp_lines() };
}

foreach my $agp_line_key ( keys %chromosome_agp_lines ) {
	if ( $chromosome_agp_lines{$agp_line_key}->get_line_type() eq 'sequence' ) {
		my $offset = $chromosome_agp_lines{$agp_line_key}->get_object_begin();
		my $agp_accession = $chromosome_agp_lines{$agp_line_key}->get_component_id();
		$chromosome_agp_offset{$agp_accession} = $offset;
	}
}

foreach my $agp_line_key ( keys %scaffold_agp_lines ) {
	my %coords;
	if ( $scaffold_agp_lines{$agp_line_key}->get_line_type() eq 'sequence' ) {
		my $scaffold_name = $scaffold_agp_lines{$agp_line_key}->get_object_being_assembled();
		my $offset = $chromosome_agp_offset{$scaffold_name};
		#offset for the start position of the scaffold
		$coords{'start'} = $scaffold_agp_lines{$agp_line_key}->get_object_begin() + $offset - 1;
		$coords{'end'} =  $scaffold_agp_lines{$agp_line_key}->get_object_end() + $offset - 1;
		$coords{'orientation'} = $scaffold_agp_lines{$agp_line_key}->get_orientation();
		my $agp_accession = $scaffold_agp_lines{$agp_line_key}->get_component_id();
		$agp_accession =~ s/\.\d+//;
		$scaffold_agp_coords{$agp_accession} = \%coords;
	}
}

#print STDERR "\%scaffold_agp_coords\n";
#print STDERR Dumper \%scaffold_agp_coords;

my $startline   = 2;
my $currentline = 0;

@group_coords_lines = split( /\n/, $input_group_coords );
foreach my $line (@group_coords_lines) {
	$currentline++;
	if ( $currentline < $startline ) {
		next;
	}
	chomp($line);
	my @bac_end_coordinates = split( /\t/, $line );
	my $bac_name            = $bac_end_coordinates[0]; #query
	if ( is_ncbi_format($bac_name) ){ $bac_name = get_accession($bac_name);}
	my $bac_ref_start       = $bac_end_coordinates[2]; #ref_start
	my $bac_ref_end         = $bac_end_coordinates[3]; #ref_end
	my $bac_query_start     = $bac_end_coordinates[5]; #query_start
	my $bac_query_end       = $bac_end_coordinates[6]; #query_end
	my $bac_query_length    = $bac_end_coordinates[7]; #query_length
	my $direction           = $bac_end_coordinates[9]; #direction (+1 if in same orientation, -1 otherwise)
	my @bac_array           = ( $bac_name, $bac_ref_start, $bac_ref_end, $bac_query_start, $bac_query_end, $bac_query_length, $direction );
	push( @bacs, \@bac_array );
}

#print STDERR "\@bacs\n";
#print STDERR Dumper \@bacs;

#old calls
#$ordered_tpf_with_bacs_inserted_in_gaps = $tpf->get_tpf_with_bacs_inserted_in_gaps( \@bacs, \%scaffold_agp_coords );
#my $out_str_from_tpf_ordered = $ordered_tpf->get_formatted_tpf();
#print $out_str_from_tpf_ordered. "\n";

#$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps = $tpf->get_tpf_with_bacs_inserted_in_sequences_and_gaps( \@bacs, \%scaffold_agp_coords, \%scaffold_component_contigs, \%scaffold_component_contig_directions );
#my $out_str_from_ordered_tpf_with_bacs_inserted_in_sequences_and_gaps = $ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->get_formatted_tpf();
#print $out_str_from_ordered_tpf_with_bacs_inserted_in_sequences_and_gaps. "\n";


#create switch point object
my $sp = Bio::GenomeUpdate::SP->new(
	taxid => '4081',
	assembly_group => 'TGP',
	assembly_unit => 'Primary',
	tpf_type => 'chromosome');

#create trim point object
my $tp = Bio::GenomeUpdate::TP->new(
	taxid => '4081',
	assembly_group => 'TGP',
	assembly_unit => 'Primary',
	tpf_type => 'chromosome');

($ordered_tpf_with_bacs_inserted_in_sequences_and_gaps, $sp_with_bacs_inserted_in_sequences_and_gaps, $tp_with_bacs_inserted_in_sequences_and_gaps) 
	= $tpf->get_tpf_sp_tp_with_bacs_inserted_in_sequences_and_gaps( $chromosome, $sp, $tp, \@bacs, \%scaffold_agp_coords, \%scaffold_component_contigs, \%scaffold_component_contig_directions );

#setting tomato genome values, change for other genomes
$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->set_assembly_version('3.0');
#$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->set_organism('Solanum lycopersicum');
$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->set_assembly_name('SL3.0');
$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->set_chromosome($chromosome);
$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->set_strain_haplotype_cultivar('Heinz 1706');
$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->set_type('Chromosome');
$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->set_comment('Includes HTGS phase 3 BACs');

my $out_str_from_ordered_tpf_with_bacs_inserted_in_sequences_and_gaps = $ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->get_formatted_tpf();
my $out_str_from_sp_with_bacs_inserted_in_sequences_and_gaps = $sp_with_bacs_inserted_in_sequences_and_gaps->get_formatted_sp();


unless (open(MODTPF,">updated.${tpf_input_file}")) { print STDERR "Cannot open output TPF file to write. Exiting..\n"; exit 1;}
print MODTPF $out_str_from_ordered_tpf_with_bacs_inserted_in_sequences_and_gaps;
close (MODTPF);

unless (open(SP,">switchpoints.updated.${tpf_input_file}")) { print STDERR "Cannot open output TPF file to write. Exiting..\n"; exit 1;}
print SP $out_str_from_sp_with_bacs_inserted_in_sequences_and_gaps;
close (SP);

=item C<is_ncbi_format ( $bac_name )>

Checks if BAC name is in NCBI format, i.e. gi|118344469|gb|AC193777.1|

=cut

sub is_ncbi_format{
	my $name = shift;
	my @name_arr = split (/\|/, $name);
	if (( $name_arr[0] eq 'gi' ) && ( $name_arr[2] eq 'gb' )){
		return 1;
	}
	else{
		return 0;
	}
}

=item C<get_accession>

Returns the accession number from a BAC name after trimming off the version number. NCBI will automatically use the latest version.

=cut

sub get_accession{
	my $name = shift;
	my @name_arr = split (/\|/, $name);
	if ($name_arr[3] =~ /\./){
		my @accession_arr = split (/\./, $name_arr[3]);
		return $accession_arr[0]; # return without version number, NCBI GRC will get latest version
	}
	else{
		return $name_arr[3];
	}
}

=item C<scaffold_summary  ( Bio::Assembly::Scaffold  )>

Accepts a scaffold object from an assembly. Prints basis statistics.
=cut

sub scaffold_summary{
	my $scaffold = shift;
	print STDERR "Scaffold or Assembly name: ".$scaffold->id()."\n";
	print STDERR "Number of BACs in scaffold: ".$scaffold->get_nof_seqs()."\n";
	print STDERR "Number of Contigs: ".$scaffold->get_nof_contigs()."\n";
	print STDERR "Number of BACs in Contigs: ".$scaffold->get_nof_contig_seqs()."\n";
	print STDERR "Number of Singlets: ".$scaffold->get_nof_singlets()."\n";
	
}

=item C<contig_component_id_direction ( Bio::Assembly::Contig  )>

Accepts a single contig object from an assembly. Returns array of components and orientation (+1,-1) in order of occurence in contig.

=cut
sub contig_component_id_direction {
	my $contig = shift;
	#ACE file spec http://bozeman.mbt.washington.edu/consed/distributions/README.29.0.txt
	#AF <read name> <C or U> <padded start consensus position>
	
	my $number_of_sequences = $contig->num_sequences();
	my $sequence_counter = 1;
	my $sequences_processed = 0;
	my (@contig_component_sequence_arr, @contig_component_directions_arr); 
	while ( $sequences_processed < $number_of_sequences ){#hacky logic
		if(defined $contig->get_seq_by_pos($sequence_counter)){#does not return value for $sequence_counter = 1, weird BioPerl 
			my $locatableleq = $contig->get_seq_by_pos($sequence_counter); #returns Bio::LocatableSeq
			my $contig_bac_name;
			if ( is_ncbi_format($locatableleq->id()) ){ $contig_bac_name = get_accession($locatableleq->id());}

			push @contig_component_sequence_arr, $contig_bac_name ;
			push @contig_component_directions_arr, $locatableleq->strand(); # returns the value of the strandedness (-1, 0 or 1)
			$sequences_processed++;
		}
		$sequence_counter++;
	}
	
	return (\@contig_component_sequence_arr, \@contig_component_directions_arr);
}


sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script inserts genome-aligned BACs into a Tiling Path File (TPF) for ONE chromosome. It also generates switch point and trim files for submission to NCBI GRC. Use filter_group_coords_output.pl to sort group_coords output by ref start coordinate and remove redundant BACs (when >1 BAC covers tha same region) to avoid problems with NCBI TPF pipeline.

    Usage:
      place_bacs.pl -t [TPF file] -s [scaffold AGP file] -c [chromosome AGP file] -n [chromosome number] -l [file with BAC names and coordinates] -a  ACE [file of assembled BACs]      

    Flags:

      -t <TPF_file>                    Original TPF file (mandatory)
      -s <scaffold AGP_file>           scaffold AGP file (mandatory)
      -c <chromosome AGP_file>         chromosome AGP file (mandatory)
      -n <chromosome number>           chromosome number. Valid values are 1-12 or Un for tomato (mandatory)
      -l <order_and_orientation_file>  Tab-delimited file containing BAC names and coordinates from align_BACends_group_coords.pl or group_coords.pl (mandatory)
      -a <ACE file>                    ACE file of assembled BACs generated by phrap. see scripts/run_phrap.sh (mandatory) 
      -h <help>                        Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>
  Surya Saha <suryasaha at cornell.edu, @SahaSurya>

=cut

