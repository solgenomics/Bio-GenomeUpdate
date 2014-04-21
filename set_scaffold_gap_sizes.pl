#!/usr/bin/perl

=head1

set_scaffold_gap_sizes.pl

=head1 SYNOPSIS

         set_scaffold_gap_sizes.pl -t [TPF file] -l [file with scaffold pairs, size, and method]

=head1 COMMAND-LINE OPTIONS

 -t  TPF file
 -l  Tab-delimited file containing tab-delimited lines of scaffold pairs,  gap size, and gap method
 -h  Help
=cut

use strict;
use Bio::GenomeUpdate::TPF::TPFSequenceLine;
use Bio::GenomeUpdate::TPF::TPFGapLine;
use Bio::GenomeUpdate::TPF;
use File::Slurp;
use Getopt::Std;

use Getopt::Std;
our ( $opt_t, $opt_l, $opt_h );
getopts('t:l:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_t || !$opt_l ) {
	print "\nTPF and gap info filenames are required.\n\n\n";
	help();
}

my $tpf_input_file = $opt_t;
my $input_tpf      = read_file($tpf_input_file)
  or die "Could not open TPF input file: $tpf_input_file\n";
my $gap_input_file = $opt_l;
my $gap_info       = read_file($gap_input_file)
  or die "Could not open gap info file: $gap_input_file\n";
my $tpf = Bio::GenomeUpdate::TPF->new();
my $out_tpf;
my @gap_info_lines;
my @gaps_sizes_and_methods;
my $output_tpf;

$tpf->parse_tpf($input_tpf);

#create a TPF for output by copying the original and clearing the lines.  This preserves the other TPF info.
#$out_tpf = $tpf;
#$out_tpf->clear_tpf_lines();

@gap_info_lines = split( /\n/, $gap_info );
foreach my $line (@gap_info_lines) {
	chomp($line);
	my @gap_size_and_method = split( /\t/, $line );
	my $first_scaffold      = $gap_size_and_method[0];
	my $second_scaffold     = $gap_size_and_method[1];
	my $gap_size            = $gap_size_and_method[2];
	my $gap_method          = $gap_size_and_method[3];

#  my @gap_size_and_method_array = ($first_scaffold, $second_scaffold, $gap_size, $gap_method);
#  push(@gaps_sizes_and_methods,\@gap_size_and_method_array);
	$tpf->change_gap_size_between_scaffolds( $first_scaffold, $second_scaffold,
		$gap_size, $gap_method );
}

#$output_tpf = $tpf->change_gap_size_between_scaffolds(\@gaps_sizes_and_methods);
my $out_str_from_tpf_altered = $tpf->get_formatted_tpf();
print $out_str_from_tpf_altered. "\n";

###write_file('processed_tpf.txt',$out_str_from_tpf_parsed);

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script changes the sizes and methods of gaps between scaffold pairs in a Tiling Path File (TPF).

    Usage:
      set_scaffold_gap_sizes.pl -t [TPF file] -l [file with scaffold pairs, size, and method]

    Flags:

      -t <TPF_file>                    Original TPF file (mandatory)
      -l <gap_size_and_method>  Tab-delimited file containing tab-delimited lines of scaffold pairs, gap sizes, and gap methods
      -h <help>                        Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

=cut

