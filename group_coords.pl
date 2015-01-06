#!/usr/bin/perl
=head1

group_coords.pl

=head1 SYNOPSIS

    group_coords.pl -i [coords file] -g [gap size]

=head1 COMMAND-LINE OPTIONS

 -i  COORDS file created by show-coords (required)
 -u  Sequence ID of chromosome with unmapped contigs (required) 
 -g  Gap size allowed between aligned clusters in the reference sequence, typically the mean/median scaffold gap
 -t  Print header
 -h  Help

=head1 TODO

  Print fasta of query seqs that did not align
  Print # of Ns on ref covered by alignments of query seqs

=cut

use strict;
use warnings;

use Getopt::Std;
use File::Slurp;
use Bio::GenomeUpdate::AlignmentCoords;
use Bio::GenomeUpdate::AlignmentCoordsGroup;

our ($opt_i, $opt_g, $opt_u, $opt_t, $opt_h);
getopts("i:g:u:t:h");
if (!$opt_i || !$opt_g ) {
  help();
}
if ($opt_h) {
  help();
}
my $input_file;
my $gap_size_allowed;
my $unmapped_ID;
my $print_header;
$input_file = $opt_i || die("-i input_file required\n");
if ($opt_g) {
  $gap_size_allowed=$opt_g;
  print STDERR "Gg: $opt_u\n";
}
if ($opt_u) {
  $unmapped_ID = $opt_u;
}else{
	$unmapped_ID = 'NA';#not required
}
if ($opt_t) {
  if ($opt_t eq "T") {
    $print_header="T";
  } elsif ($opt_t eq "F") {
    $print_header="F";
  } else {
    die ("-t must be T or F\n");
  }
}

my $total=0;
my $total_smaller_than_20k=0; #for alignments covering < 20k on ref
my $total_mixed=0;
my $total_over=0;
my $total_alt=0;
my $total_full_length=0;
my $total_to_end=0;
my $total_extend=0;
my $total_ref_covered=0;

print STDERR "G0: $gap_size_allowed\n";
my @lines = read_file($input_file);
my $startline = 5;
my $currentline = 0;
my @alignment_coords_array;
my $last_line_query_id;
my $last_query_id;
my $last_query_length;
print "query\treference\tref_start\tref_end\tlength\tq_start\tq_end\tq_length\tseq_in_clusters\tdirection\tref_count\tincludes_0\tfull_length\tfrom_start\tfrom_end\tinternal_gap\tis_overlapping\tsize_of_alt\talternates\t\n";

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
  	#print info for last query (assembled or singleton BAC ) to STDOUT
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
    #calc_and_print_info(\@alignment_coords_array, $current_query_id, $current_query_length,$gap_size_allowed);
    calc_and_print_info(\@alignment_coords_array, $current_query_id, $current_query_length);
    @alignment_coords_array = ();
  }
  $last_line_query_id = $current_query_id;
  $last_query_id = $current_query_id;
  $last_query_length = $current_query_length;
}

sub calc_and_print_info {
  my ($aref,$q_id,$q_length) = @_;
  my $align_group =  Bio::GenomeUpdate::AlignmentCoordsGroup->new();
  #assign all coords for query/assembled or singleton BAC to obj
  $align_group->set_array_of_alignment_coords($aref);
  my $zero_chromosome_id = $unmapped_ID;
  #Returns IDs, start and end coordinates, total aligned sequence length, and direction for the longest proximity-grouped 
  #alignment clusters sorted by longest to shortest length of non-overlapping sequence covered by alignment clusters.  
  #The proximity grouping is done using the specified length of an allowed gap between aligned clusters in the reference sequence ($gap allowed).
  #$sequence_aligned_in_clusters is the total length of ref covered by alignment grp
  my ($ref_id, $query_id, $ref_start, $ref_end, $query_start, $query_end, $sequence_aligned_in_clusters,$direction,$is_overlapping,
  	$size_of_next_largest_match,$alternates) = $align_group->get_id_coords_and_direction_of_longest_alignment_cluster_group($gap_size_allowed);
  my $is_full_length;
  my $start_gap_length = $query_start - 1;#region of query before aligned part
  my $end_gap_length = $q_length - $query_end;#region of query after aligned part
  my $internal_gap_length = ($q_length - $sequence_aligned_in_clusters) - ($start_gap_length + $end_gap_length);#not clear how its used, see other scripts
  if (($query_start == 1) && ($query_end == $q_length)) {
    $is_full_length = "Contains";#entire query is covered in the alignment grp
  } else {
    $is_full_length = "Partial";
  }
  print $q_id."\t";
  print $ref_id."\t";
  print $ref_start."\t";
  print $ref_end."\t";
  print $ref_end - $ref_start."\t";
  print $query_start."\t";
  print $query_end."\t";
  print $q_length."\t";
  print $sequence_aligned_in_clusters."\t";
  print $direction."\t";#strand
  print $align_group->get_count_of_reference_sequence_ids()."\t";;
  print $align_group-> includes_reference_id($zero_chromosome_id)."\t";
  print $is_full_length."\t";    
  print $start_gap_length."\t";
  print $end_gap_length."\t";
  print $internal_gap_length."\t";
  print $is_overlapping."\t";
  print $size_of_next_largest_match."\t";
  print $alternates."\t";
  #if (defined($second_id)){
  #print $second_id."\t";
  #print $second_size."\t"
  #}
  #else {
  #print "None\tNone\t";
  #}
  print "\n";

  my $flagged = 0;#potential problem

  $total++;
  if ($direction == 0) {#query aligns to bot + and - strand of ref
    $total_mixed++;
    $flagged=1;
  }
  if ($is_overlapping == 1) {#query alignment tiles overlap
    $total_over++;
    $flagged=1;
  }
  if ($size_of_next_largest_match > 10000) {
    $total_alt++;
    $flagged=1;
  }
  if ($ref_end - $ref_start < 20000) {
    $total_smaller_than_20k++;
    #$flagged=1;#short alignment may be for a BAC end so not always a negative
  }
  if ($start_gap_length < 10 && $end_gap_length < 10 && $flagged==0) {#alignments cover query
    $total_full_length++;
  }
  if (($start_gap_length < 10 || $end_gap_length < 10) && $flagged==0) {#alignments cover till one end of query
    $total_to_end++;
  }
  if ($flagged==0) {
    $total_extend += $start_gap_length + $end_gap_length;
  }
  if ($flagged==0) {
    $total_ref_covered += $sequence_aligned_in_clusters;
  }
}

  ##summary info
print STDERR "Total:\t$total\n";
print STDERR "Total smaller than 20,000 on ref:\t$total_smaller_than_20k\n";
print STDERR "Total with mixed orientation:\t$total_mixed\n";
print STDERR "Total with overlapping alignment clusters:\t$total_over\n";
print STDERR "Total with alternate alignments > 10,000:\t$total_alt\n";
print STDERR "Total full length:\t$total_full_length\n";
print STDERR "Total with alignment to at least one end:\t$total_to_end\n";
print STDERR "Total reference extended by BACs:\t$total_extend\n";#new seqs from query
print STDERR "Total reference covered by BACs:\t$total_ref_covered\n";#includes gaps ($gap_size_allowed) between alignment clusters


sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script groups aligned clusters and creates a tab delimited file with BAC alignment details.

    Usage:
      group_coords.pl -i [coords file] -g [gap size]

    Flags:

       -i  <coords file>             COORDS file created by show-coords (required)
       -g  <int>                     Gap size allowed. Recommended 10000 (required)
       -u  <str>                     Sequence ID(seqid) of chromosome with unmapped contigs/scaffolds. Typically chromosome 0.
       -t  <T/F>                     Print header. Must be T or F. Default is T
       -h  <help>
EOF
exit (1);
}

__END__


[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[FRM]	[TAGS]
794498	870595	1	76098	76098	76098	99.99	65875088	76098	0.12	100.00	1	1	SL2.50ch05	gi|108743801|gb|AC187148.1|	[CONTAINS]
632723	750165	117443	1	117443	117443	100.00	65875088	117443	0.18	100.00	1	-1	SL2.50ch05	gi|110431384|gb|AC188778.1|	[CONTAINS]
38351613	38505235	153623	1	153623	153623	100.00	68045021	153623	0.23	100.00	1	-1	SL2.50ch07	gi|117165387|emb|CU024881.7|	[CONTAINS]
84085410	84193187	107778	1	107778	107778	99.99	98543444	107778	0.11	100.00	1	-1	SL2.50ch01	gi|118344475|gb|AC193779.1|	[CONTAINS]
7850690	7969210	1	118521	118521	118521	100.00	56302525	118521	0.21	100.00	1	1	SL2.50ch11	gi|118344479|gb|AC171734.2|	[CONTAINS]
34049026	34050620	117470	115876	1595	1595	99.69	49751636	137379	0.00	1.16	1	-1	SL2.50ch06	gi|119371448|dbj|AP009271.1|	
64063464	64221573	1	158110	158110	158110	100.00	65866657	158110	0.24	100.00	1	1	SL2.50ch08	gi|119371464|dbj|AP009287.1|	[CONTAINS]
50117537	50118617	57136	56056	1081	1081	99.72	68045021	153195	0.00	0.71	1	-1	SL2.50ch07	gi|126153618|emb|CU074307.9|	
