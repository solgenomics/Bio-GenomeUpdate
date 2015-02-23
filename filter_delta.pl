#!/usr/bin/perl
=head1

filter_delta.pl

=head1 SYNOPSIS

    filter_delta.pl -c [coords file] -d [delta file] -o [delta file]

=head1 COMMAND-LINE OPTIONS

 -c  COORDS file created by show-coords (required)
 -d  DELTA file created by nucmer (required) 
 -o  Output filtered DELTA file 
 -h  Help

=head1 TODO

  See github issues.

=cut

use strict;
use warnings;

use Getopt::Std;
use File::Slurp;
use Bio::GenomeUpdate::AlignmentCoords;
use Bio::GenomeUpdate::AlignmentCoordsGroup;


our ($opt_c, $opt_d, $opt_o, $opt_h);
getopts("c:d:o:h");
if (!$opt_c || !$opt_d ) {
  print STDERR "Required files missing!! exiting..\n";
  help();
}
if ($opt_h) {
  help();
}
unless (-e $opt_c){ print STDERR "COORDS file $opt_c not found. exiting..\n"; exit;}
unless (-e $opt_d){ print STDERR "DELTA file $opt_d not found. exiting..\n"; exit;}

my $input_coords_file = $opt_c || die("-i input_file required\n");
my $input_delta_file = $opt_d || die("-i input_file required\n");

my @lines = read_file($input_coords_file);
my $startline = 5;
my $currentline = 0;
my @alignment_coords_array;
my @query_valid_hits;
my @query_invalid_hits;
my $last_line_query_id;
my $last_query_id;
my $last_query_length;

#parse coords file
#[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[FRM]	[TAGS]
foreach my $line (@lines) {
  $currentline++;
  if ($currentline < $startline) {
    next;
  }
  my @row;
  @row = split('\t',$line);
  my $current_query_id = $row[14];
  my $current_query_length = $row[8];
  if (!defined($last_line_query_id)) {
    $last_line_query_id = $current_query_id;
  }
  #exec if query ID changes, i.e., coords for next assembled or singleton BAC aligned to chr
  if (!($current_query_id eq $last_line_query_id)) { 
    #find alignments that align in mixed orientation or in non co-linear order
    #print DELTA to STDOUT and messages to STDERR
    calc_and_print_info(\@alignment_coords_array, $last_query_id, $last_query_length);
    
    @alignment_coords_array = ();
  }
  my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
  $aln_coords->set_reference_id($row[13]);
  $aln_coords->set_query_id($row[14]);
  $aln_coords->set_reference_start_coord($row[0]);
  $aln_coords->set_reference_end_coord($row[1]);
  $aln_coords->set_query_start_coord($row[2]);
  $aln_coords->set_query_end_coord($row[3]);
  push(@alignment_coords_array, $aln_coords);
  
  #deal with last row since no more alignments for query after this
  if ($currentline==scalar(@lines)) {
    #find alignments that align in mixed orientation or in non co-linear order
    #print DELTA to STDOUT and messages to STDERR
    calc_errors_and_print_delta(\@alignment_coords_array, $last_query_id, $last_query_length);
    @alignment_coords_array = ();
  }
  $last_line_query_id = $current_query_id;
  $last_query_id = $current_query_id;
  $last_query_length = $current_query_length;
}

=item C<calc_errors_and_print_delta (@alignment_coords_array, $last_query_id, $last_query_length)>

If non co-linear alignment is found, it prints error info to STDERR and returns 1. Otherwise it returns a 0.

=cut

sub calc_errors_and_print_delta {
  my ($aref,$q_id,$q_length) = @_;
  my $align_group =  Bio::GenomeUpdate::AlignmentCoordsGroup->new();
  
  #assign all coords for query/assembled or singleton BAC to obj
  $align_group->set_array_of_alignment_coords($aref);
  
  
  
}

=item C<help>

duh..

=cut


sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script filters the DELTA file produced by nucmer by removing BAC alignments that align in non co-linear order. Use the 
     following commands to generate the COORDS and DELTA files where $1 is chr number
     
     nucmer -l 100 -c 500 --noextend -p 500bp.ch0${1}_asm_BACs__SL2.50ch0${1} SL2.50ch0${1}.fa ch0${1}_asm_BACs.fas
     delta-filter -l 500 -u 99 500bp.ch0${1}_asm_BACs__SL2.50ch0${1}.delta > 500bp.ch0${1}_asm_BACs__SL2.50ch0${1}.delta.filtered
     show-coords -c -d -l -q -T -o 500bp.ch0${1}_asm_BACs__SL2.50ch0${1}.delta.filtered > 500bp.ch0${1}_asm_BACs__SL2.50ch0${1}.delta.filtered.coords
     

    Usage:
      filter_delta.pl -c [coords file] -d [delta file] -o [delta file]

    Flags:

      -c  COORDS file created by show-coords (required)
      -d  DELTA file created by nucmer (required) 
      -o  Output filtered DELTA file 
      -h  Help
 
EOF
exit (1);
}


=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha at cornell.edu, @SahaSurya>

=cut

__END__