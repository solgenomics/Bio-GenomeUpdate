#!/usr/bin/perl

=head1

place_bacs.pl

=head1 SYNOPSIS

    place_bacs.pl -t [TPF file] -s [scaffold AGP file] -c [chromosome AGP file] -l [file with BAC names and coordinates]

=head1 COMMAND-LINE OPTIONS

 -t  TPF file (required)
 -s  scaffold AGP file (required)
 -c  chromosome AGP file (required)
 -l  Tab-delimited file containing BAC names and start and end coordinates from align_BACends_group_coords.pl or group_coords.pl(required) 
 -h  Help

=head1 TODO

  See github issues.
  
=cut

use strict;
use Bio::GenomeUpdate::TPF::TPFSequenceLine;
use Bio::GenomeUpdate::TPF::TPFGapLine;
use Bio::GenomeUpdate::TPF;
use Bio::GenomeUpdate::AGP;
use File::Slurp;
use Getopt::Std;

use Getopt::Std;
our ( $opt_t, $opt_l, $opt_s, $opt_c, $opt_h );
getopts('t:s:c:l:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_t || !$opt_l || !$opt_s || !$opt_c ) {
	print STDERR
	  "\nTPF, AGP and group_coords BAC name/coordinates filenames are required.\n\n\n";
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
my $tpf            = Bio::GenomeUpdate::TPF->new();
my $scaffold_agp   = Bio::GenomeUpdate::AGP->new();
my $chromosome_agp = Bio::GenomeUpdate::AGP->new();
my $out_tpf;
my @group_coords_lines;
my @bacs;

#my @sorted_bacs;
#my $ordered_tpf_with_bacs_inserted_in_gaps;
my $ordered_tpf;
my $ordered_tpf_with_bacs_inserted_in_sequences_and_gaps;
my %chromosome_agp_offset;
my %scaffold_agp_coords;

$tpf->parse_tpf($input_tpf);
$scaffold_agp->parse_agp($scaffold_input_agp);
$chromosome_agp->parse_agp($chromosome_input_agp);

#create a TPF for output by copying the original and clearing the lines.  This preserves the other TPF info.
#$out_tpf = $tpf;
#$out_tpf->clear_tpf_lines();

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
		$coords{'start'} = $scaffold_agp_lines{$agp_line_key}->get_object_begin() + $offset;
		$coords{'end'} =  $scaffold_agp_lines{$agp_line_key}->get_object_end() + $offset;
		$coords{'orientation'} = $scaffold_agp_lines{$agp_line_key}->get_orientation();
		my $agp_accession = $scaffold_agp_lines{$agp_line_key}->get_component_id();
		$agp_accession =~ s/\.\d+//;
		$scaffold_agp_coords{$agp_accession} = \%coords;
	}
}

my $startline   = 2;
my $currentline = 0;

@group_coords_lines = split( /\n/, $input_group_coords );
foreach my $line (@group_coords_lines) {
	$currentline++;
	if ( $currentline < $startline ) {
		next;
	}
	chomp($line);
	my @bac_and_coordinates = split( /\t/, $line );
	my $bac_name            = $bac_and_coordinates[0]; #query
	if ( is_ncbi_format($bac_name) ){ $bac_name = get_accession($bac_name);}
	my $bac_start           = $bac_and_coordinates[2]; #ref_start
	my $bac_end             = $bac_and_coordinates[3]; #ref_end
	my @bac_array           = ( $bac_name, $bac_start, $bac_end );
	push( @bacs, \@bac_array );
}

#sort @bacs

#$ordered_tpf_with_bacs_inserted_in_gaps = $tpf->get_tpf_with_bacs_inserted_in_gaps( \@bacs, \%scaffold_agp_coords );
#my $out_str_from_tpf_ordered = $ordered_tpf->get_formatted_tpf();
#print $out_str_from_tpf_ordered. "\n";

#$ordered_tpf_with_bacs_inserted_in_gaps_and_sequences = $ordered_tpf_with_bacs_inserted_in_gaps->get_tpf_with_bacs_inserted_in_sequences( \@bacs, \%scaffold_agp_coords );
$ordered_tpf_with_bacs_inserted_in_sequences_and_gaps = $tpf->get_tpf_with_bacs_inserted_in_sequences_and_gaps( \@bacs, \%scaffold_agp_coords );
my $out_str_from_ordered_tpf_with_bacs_inserted_in_sequences_and_gaps = $ordered_tpf_with_bacs_inserted_in_sequences_and_gaps->get_formatted_tpf();
print $out_str_from_ordered_tpf_with_bacs_inserted_in_sequences_and_gaps. "\n";

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

=item C<get_gap_methods>

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

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script inserts genome-aligned BACs into a Tiling Path File (TPF). It also generates switch point and trim files for submission to NCBI GRC.

    Usage:
      place_bacs.pl -t [TPF file] -a [AGP file] -l [file with BAC names and start and end coordinates]

    Flags:

      -t <TPF_file>                    Original TPF file (mandatory)
      -s <scaffold AGP_file>           scaffold AGP file (mandatory)
      -c <chromosome AGP_file>         chromosome AGP file (mandatory)
      -l <order_and_orientation_file>  Tab-delimited file containing tab-delimited lines of BAC names and start and end coordinates.
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

