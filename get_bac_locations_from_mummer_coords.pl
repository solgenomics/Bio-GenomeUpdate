#!/usr/bin/perl

=head1

get_bac_locatons_from_mummer_coords.pl

=head1 SYNOPSIS

    get_bac_locatons_from_mummer_coords.pl -i [MUMmer show-coords file]

=head1 COMMAND-LINE OPTIONS

 -i  MUMmer show-coords file
 -h  Help

=cut

use strict;
use File::Slurp;
use Getopt::Std;

use Getopt::Std;
our ($opt_i, $opt_h);
getopts('i:h');
if ($opt_h){
  help();
  exit;
}
if (!$opt_i) {
    print "\nMUMmer show-coords filename is required\n\n\n";
    help();
}

my $input_file = $opt_i;
my $input_coords_contents = read_file($input_file) or die "Could not open show-coords input file: $input_file\n";
my @input_lines;
my %query_sequences;


@input_lines = split (/\n/, $input_coords_contents);
foreach my $line (@input_lines) {
  chomp($line);
  my @coord_line_array = split(/\t/,$line);
  if ((!$coord_line_array[14]) || $coord_line_array[0] =~ /\D/) {
    next;
  }
  my $reference_name = $coord_line_array[13];
  my $query_name = $coord_line_array[14];
  my $reference_start = $coord_line_array[0];
  my $reference_end = $coord_line_array[1];
  my $query_start = $coord_line_array[2];
  my $query_end = $coord_line_array[3];
  my $reference_alignment_length = $coord_line_array[4];
  my $query_alignment_length = $coord_line_array[5];
  my $percent_identity = $coord_line_array[6];
  my $reference_length = $coord_line_array[7];
  my $query_length = $coord_line_array[8];
  my $query_orientation = $coord_line_array[12];
  if ($query_sequences{$query_name}) {
    if ($query_orientation eq 1) {
      if ($reference_start  < $query_sequences{$query_name}{'end'}){
	$query_sequences{$query_name}{'alignment_overlaps'} = 1;
      }
    }
    else {
      if ($reference_end  < $query_sequences{$query_name}{'end'}){
	$query_sequences{$query_name}{'alignment_overlaps'} = 1;
      }
    }
    if ($reference_start < $query_sequences{$query_name}{'start'}) {
      $query_sequences{$query_name}{'start'} = $reference_start;
    }
    if ($reference_end > $query_sequences{$query_name}{'end'}) {
      $query_sequences{$query_name}{'end'} = $reference_end;
    }
    if (!$query_orientation == $query_sequences{$query_name}{'orientation'}) {
      $query_sequences{$query_name}{'mixed_orientations'} = 1;
    }
    $query_sequences{$query_name}{'length'} += $reference_alignment_length;
  }
  else {
    my %query_info;
    $query_info{'reference_name'} = $reference_name;
    $query_info{'start'} = $reference_start;
    $query_info{'end'} = $reference_end;
    $query_info{'orientation'} = $query_orientation;
    $query_info{'mixed_orientations'} = 0;
    $query_info{'length'} = $reference_alignment_length;
    $query_info{'aligns_to_beginning'} = 0;
    $query_info{'aligns_to_end'} = 0;
    $query_info{'query_length'} = $query_length;
    $query_info{'alignment_overlaps'} = 0;
    $query_sequences{$query_name} = {%query_info};
  }
  if ($query_start == 1 || $query_end == 1) {
    $query_sequences{$query_name}{'aligns_to_beginning'} = 1;
  }
  if ($query_start == $query_length || $query_end == $query_length) {
    $query_sequences{$query_name}{'aligns_to_end'} = 1;
  }
}

# for my $query_name ( keys %query_sequences ) {
#   print $query_sequences{$query_name}{'reference_name'}."\t";
#   print $query_sequences{$query_name}{'start'}."\t";
#   print $query_sequences{$query_name}{'end'}."\t";
#   print $query_sequences{$query_name}{'length'}."\t";
#   print $query_name."\t";
#   print $query_sequences{$query_name}{'orientation'}."\t";
#   print $query_sequences{$query_name}{'mixed_orientations'}."\t";
#   print $query_sequences{$query_name}{'aligns_to_beginning'}."\t";
#   print $query_sequences{$query_name}{'aligns_to_end'}."\t";
#   print $query_sequences{$query_name}{'alignment_overlaps'}."\n";
# }

#foreach my $key (    #
#    sort { $hash{$a}->{Make} cmp $hash{$b}->{Make} }    #
#    keys %hash
#    )

foreach my $query_name ( sort {$query_sequences{$a}->{'reference_name'} cmp $query_sequences{$b}->{'reference_name'} ||
$query_sequences{$a}->{'start'} <=> $query_sequences{$b}->{'start'}} keys %query_sequences) {
  print $query_sequences{$query_name}{'reference_name'}."\t";
  print $query_sequences{$query_name}{'start'}."\t";
  print $query_sequences{$query_name}{'end'}."\t";
  print $query_sequences{$query_name}{'length'}."\t";
  print $query_name."\t";
  print $query_sequences{$query_name}{'orientation'}."\t";
  print $query_sequences{$query_name}{'mixed_orientations'}."\t";
  print $query_sequences{$query_name}{'aligns_to_beginning'}."\t";
  print $query_sequences{$query_name}{'aligns_to_end'}."\t";
  print $query_sequences{$query_name}{'alignment_overlaps'}."\n";
}

sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script reads coordinates of BACs aligned to a genome from a MUMmer show-coords file and determines start and end positions and does filtering for low quality/aberrant alignments.

    Usage:
      get_bac_locatons_from_mummer_coords.pl -i [MUMmer show-coords file]

    Flags:

      -i <show-coords file>            MUMmer show-coords file (mandatory)
      -h <help>                        Help

EOF
exit (1);
}


=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

=cut
