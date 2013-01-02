#!/usr/bin/perl

=head1

place_bacs.pl

=head1 SYNOPSIS

    place_bacs.pl -t [TPF file] -a [AGP file] -l [file with BAC names and coordinates]

=head1 COMMAND-LINE OPTIONS

 -t  TPF file
 -a  AGP file
 -l  Tab-delimited file containing BAC names and start and end coordinates.
 -h  Help

=cut

use strict;
use Bio::GenomeUpdate::TPF::TPFSequenceLine;
use Bio::GenomeUpdate::TPF::TPFGapLine;
use Bio::GenomeUpdate::TPF;
use Bio::GenomeUpdate::AGP;
use File::Slurp;
use Getopt::Std;

use Getopt::Std;
our ($opt_t, $opt_l, $opt_a, $opt_h);
getopts('t:a:l:h');
if ($opt_h){
  help();
  exit;
}
if (!$opt_t || !$opt_l || !$opt_a) {
    print "\nTPF, AGP and BAC name/coordinates filenames are required.\n\n\n";
    help();
}

my $tpf_input_file = $opt_t;
my $input_tpf = read_file($tpf_input_file) or die "Could not open TPF input file: $tpf_input_file\n";
my $agp_input_file = $opt_a;
my $input_agp = read_file($agp_input_file) or die "Could not open AGP input file: $agp_input_file\n";
my $bac_input_file = $opt_l;
my $input_bacs = read_file($bac_input_file) or die "Could not open BAC name and coordinates file: $bac_input_file\n";
my $tpf = Bio::GenomeUpdate::TPF->new();
my $agp = Bio::GenomeUpdate::AGP->new();
my $out_tpf;
my @bac_lines;
my @bacs;
#my @sorted_bacs;
my $ordered_tpf;
my %agp_coords;

$tpf->parse_tpf($input_tpf);
$agp->parse_agp($input_agp);


#create a TPF for output by copying the original and clearing the lines.  This preserves the other TPF info.
#$out_tpf = $tpf;
#$out_tpf->clear_tpf_lines();


my %agp_lines;
if ($agp->has_agp_lines()) {
  %agp_lines = %{$agp->get_agp_lines()};
}

foreach my $agp_line_key (keys %agp_lines){
  my %coords;
  if ($agp_lines{$agp_line_key}->get_line_type() eq 'sequence') {
    $coords{'start'} = $agp_lines{$agp_line_key}->get_object_begin();
    $coords{'end'} = $agp_lines{$agp_line_key}->get_object_end();
    my $agp_accession = $agp_lines{$agp_line_key}->get_component_id();
    $agp_accession =~ s/\.\d+//;
    $agp_coords{$agp_accession} = \%coords;
  }

}

@bac_lines = split (/\n/, $input_bacs);
foreach my $line (@bac_lines) {
  chomp($line);
  my @bac_and_coordinates = split(/\t/,$line);
  my $bac_name = $bac_and_coordinates[0];
  my $bac_start = $bac_and_coordinates[1];
  my $bac_end = $bac_and_coordinates[2];
  my @bac_array = ($bac_name, $bac_start, $bac_end);
  push(@bacs,\@bac_array);
}

#sort @bacs

$ordered_tpf = $tpf->get_tpf_with_bacs_inserted(\@bacs, \%agp_coords);
my $out_str_from_tpf_ordered = $ordered_tpf->get_formatted_tpf();
print $out_str_from_tpf_ordered."\n";


sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script inserts genome-aligned BACs into a Tiling Path File (TPF).

    Usage:
      place_bacs.pl -t [TPF file] -a [AGP file] -l [file with BAC names and start and end coordinates]

    Flags:

      -t <TPF_file>                    Original TPF file (mandatory)
      -a <AGP_file>                    AGP file (mandatory)
      -l <order_and_orientation_file>  Tab-delimited file containing tab-delimited lines of BAC names and start and end coordinates.
      -h <help>                        Help

EOF
exit (1);
}


=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

=cut
