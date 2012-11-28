#!/usr/bin/perl

=head1

reorder_scaffolds.pl

=head1 SYNOPSIS

    reorder_scaffolds.pl -t [TPF file] -l [file with order and orientation]

=head1 COMMAND-LINE OPTIONS

 -t  TPF file
 -l  Tab-delimited file containing tab-delimited lines of scaffolds orientation as + or - listed in the desired new order.
 -h  Help

=cut

use strict;
use Bio::GenomeUpdate::TPF::TPFSequenceLine;
use Bio::GenomeUpdate::TPF::TPFGapLine;
use Bio::GenomeUpdate::TPF;
use File::Slurp;
use Getopt::Std;

use Getopt::Std;
our ($opt_t, $opt_l,$opt_h);
getopts('t:l:h');
if ($opt_h){
  help();
  exit;
}
if (!$opt_t || !$opt_l) {
    print "\nTPF and order/orientation filenames are required.\n\n\n";
    help();
}

my $tpf_input_file = $opt_t;
my $input_tpf = read_file($tpf_input_file) or die "Could not open TPF input file: $tpf_input_file\n";
my $order_input_file = $opt_l;
my $input_order = read_file($order_input_file) or die "Could not open scaffold order and orientation file: $order_input_file\n";
my $tpf = Bio::GenomeUpdate::TPF->new();
my $out_tpf;
my @order_lines;
my @ordered_and_oriented_scaffolds;
my $ordered_tpf;

$tpf->parse_tpf($input_tpf);

#create a TPF for output by copying the original and clearing the lines.  This preserves the other TPF info.
#$out_tpf = $tpf;
#$out_tpf->clear_tpf_lines();

@order_lines = split (/\n/, $input_order);
foreach my $line (@order_lines) {
  chomp($line);
  my @scaffold_and_orientation = split(/\t/,$line);
  my $scaffold_name = $scaffold_and_orientation[0];
  my $orientation = $scaffold_and_orientation[1];
  if (($orientation eq "+") || ($orientation eq "-")) {
    my @scaffold_and_orientation_array = ($scaffold_name, $orientation);
    push(@ordered_and_oriented_scaffolds,\@scaffold_and_orientation_array);
   }
  else {
    die "Orientation must be specified as + or -\n";
  }
}

$ordered_tpf = $tpf->get_tpf_in_new_scaffold_order(\@ordered_and_oriented_scaffolds);
my $out_str_from_tpf_ordered = $ordered_tpf->get_formatted_tpf();
print $out_str_from_tpf_ordered."\n";

###write_file('processed_tpf.txt',$out_str_from_tpf_parsed);


sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script chages the order and orientation of scaffolds in a Tiling Path File (TPF).

    Usage:
      reorder_scaffolds.pl -t [TPF file] -l [file with order and orientation]

    Flags:

      -t <TPF_file>                    Original TPF file (mandatory)
      -l <order_and_orientation_file>  Tab-delimited file containing tab-delimited lines of scaffolds orientation as + or - listed in the desired new order.
      -h <help>                        Help

EOF
exit (1);
}


=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

=cut
