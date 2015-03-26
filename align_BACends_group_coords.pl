#!/usr/bin/perl

=head1

group_coords.pl

=head1 SYNOPSIS

    align_BACends_group_coords.pl -l [BAC end length] -r [fasta] -q [fasta] -c [agp] -s [agp]

=head1 COMMAND-LINE OPTIONS

 -u  Sequence ID of chromosome with unmapped contigs 
 -l  BAC end length (required)
 -r  Fasta file of reference (required)
 -q  Fasta file of query (assembled and singleton BACs, required)
 -c  Contig or component AGP file for reference (includes scaffold gaps)
 -s  Chromosome AGP file for reference (with only scaffolds and gaps) 
 -t  Print header
 -h  Help

=head1 TODO

  See github issues.

=cut

use strict;
use warnings;

use Getopt::Std;
use File::Slurp;
use Scalar::Util qw(looks_like_number);
use Bio::GenomeUpdate::AlignmentCoords;
use Bio::GenomeUpdate::AlignmentCoordsGroup;
use Bio::GenomeUpdate::AGP;
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Seq;


our ( $opt_l, $opt_r, $opt_q, $opt_u, $opt_t, $opt_h, $opt_c, $opt_s );
getopts("l:r:q:u:c:s:t:h");
if ( !$opt_l || !$opt_c || !$opt_s || !$opt_r || !$opt_q ) {
	print STDERR "\nRequired parameters missing!! exiting..\n\n";
	help();
}
if ($opt_h) {
	help();
}
unless ( -e $opt_r ) {
	print STDERR "Fasta file $opt_r not found. exiting..\n";
	exit;
}
unless ( -e $opt_q ) {
	print STDERR "Fasta file $opt_q not found. exiting..\n";
	exit;
}
if ( defined $opt_c ) {
	unless ( -e $opt_c ) {
		print STDERR
		  "Contig or component AGP file $opt_c not found. exiting..\n";
		exit;
	}
}
if ( defined $opt_s ) {
	unless ( -e $opt_s ) {
		print STDERR "Chromosome AGP file $opt_s not found. exiting..\n";
		exit;
	}
}

my $unmapped_ID;
my $print_header;
my $bacend_length;
if ($opt_u) {
	$unmapped_ID = $opt_u;
	print STDERR "Gg: $opt_u\n";
}
else {
	$unmapped_ID = 'NA';    #not required
}
if ($opt_t) {
	if ( $opt_t eq "T" ) {
		$print_header = "T";
	}
	elsif ( $opt_t eq "F" ) {
		$print_header = "F";
	}
	else {
		die("-t must be T or F\n");
	}
}
else{
	$print_header = "T";
}

if (looks_like_number($opt_l)){
	$bacend_length = $opt_l;
}
else{
	die("$opt_l is not a number");
}
  

open( MIXED, ">mixed_end_${bacend_length}_align_BACends_group_coords.out" )
  or die
"Could not create mixed_end_${bacend_length}_align_BACends_group_coords.out for writing out BACs aligned to ref chr in mixed orientation";
open( NONCOLINEAR, ">noncolinear_end_${bacend_length}_align_BACends_group_coords.out" )
  or die
"Could not create noncolinear_end_${bacend_length}_align_BACends_group_coords.out for writing out BACs aligned non co-linearly to ref chr, i.e. different order of aligned tiles on BAC and ref chr. Shows miassembly on ref chr or BAC";

my ($contig_agp,$chr_agp);

if ($opt_c && $opt_s){
	my $contig_agp_input_file = $opt_c;
	my $contig_input_agp      = read_file($contig_agp_input_file)
	  or die "Could not open contig AGP input file: $contig_agp_input_file\n";
	$contig_agp = Bio::GenomeUpdate::AGP->new();
	$contig_agp->parse_agp($contig_input_agp);
	my $chr_agp_input_file = $opt_s;
	my $chr_input_agp      = read_file($chr_agp_input_file)
	  or die "Could not open chr AGP input file: $chr_agp_input_file\n";
	$chr_agp = Bio::GenomeUpdate::AGP->new();
	$chr_agp->parse_agp($chr_input_agp);
}

my $ref_db   = Bio::DB::Fasta->new( $opt_r, '-reindex' => 1 );
die "\nMultiple references will break the BAC end comparison logic. Please run for individual chromosomes. Exiting...\n" if (scalar $ref_db->get_all_ids() > 1);
my $query_db = Bio::DB::Fasta->new( $opt_q, '-reindex' => 1 );
my %query_lengths;

print STDERR "Number of reference sequences: ";
print STDERR scalar $ref_db->get_all_ids();
print STDERR "\nNumber of query sequences: ";
print STDERR scalar $query_db->get_all_ids();
print STDERR "\n\n";


###### create fasta file of BACends ###### 
my $return_value = archive_old('query_bacends.fasta', 'Could not remove old query_bacends.fasta file.');
if (!$return_value) {print STDERR "Old query_bacends.fasta moved to query_bacends.fasta.old.\n";}
my $query_bacends_fasta = Bio::SeqIO->new( -file => ">query_bacends.fasta", -format => 'Fasta');

# Loop through sequence objects
my $stream  = $query_db->get_PrimarySeq_stream();
while (my $query_seq_obj = $stream->next_seq()) { # returns Bio::PrimarySeqI obj
	my ($bacend_source, $bacend_name, $bacend_seq, $bacend_start, $bacend_stop);
	#5' end
	$bacend_source = $query_seq_obj->display_id();
	$bacend_name .=$bacend_source."_left_".$bacend_length;
	$bacend_start = 1;
	$bacend_stop = $bacend_length;
	$bacend_seq = $query_seq_obj->subseq($bacend_start, $bacend_stop);
	my $left_end_seq = Bio::Seq->new( -display_id => $bacend_name,
										-seq => $bacend_seq);
	$query_bacends_fasta->write_seq($left_end_seq);
	
	#3' end
	$bacend_name = $bacend_source."_right_".$bacend_length;
	#$query_seq_length = length  $query_seq_obj->seq();
	$bacend_start = $query_seq_obj->length() - $bacend_length + 1 ;
	$bacend_stop = $bacend_start + $bacend_length - 1;  
	#print STDERR "\n**\nname:",$query_seq_obj->display_id(),"\nlength:",$query_seq_obj->length(),"\nstart: $bacend_start \nend: $bacend_stop\n";
	$bacend_seq = $query_seq_obj->subseq( $bacend_start, $bacend_stop );
	my $right_end_seq = Bio::Seq->new( -display_id => $bacend_name,
										-seq => $bacend_seq);
	$query_bacends_fasta->write_seq($right_end_seq);
	
	#remember the name, length of query
	$query_lengths{$bacend_source} = $query_seq_obj->length(); 
}


###### run mummer with optimized paramaters ######

print STDERR "\nStarting mummer run\n";
archive_old('nucmer.coords','Could not remove old nucmer.coords');
#system('nucmer', '-l 100', "-c $bacend_length", '-p nucmer.coords', $opt_r, $opt_q );
#anchors matches that are unique in in the reference but not necessarily unique in the query
#so only 1 match reported for both BAC ends, it one exists 
system("nucmer --noextend -l 100 -c $bacend_length -p nucmer.coords $opt_r query_bacends.fasta");
die("\nCould not run nucmer. $!\nExiting...\n\n") if ($? == -1);

#is this required any more?? cluster size in nucmer should only report alignment >= $bacend_length 
archive_old('nucmer.coords.delta.filtered','Could not remove old nucmer.coords.delta.filtered');
system("delta-filter -l $bacend_length nucmer.coords.delta > nucmer.coords.delta.filtered");
die("\nCould not run delta-filter. $!\nExiting...\n\n") if ($? == -1);

archive_old('nucmer.coords.delta.filtered.coords','Could not remove old nucmer.coords.delta.filtered.coords');
system("show-coords -c -d -l -q -T -o nucmer.coords.delta.filtered > nucmer.coords.delta.filtered.coords");
die("\nCould not run show-coords. $!\nExiting...\n\n") if ($? == -1);

print STDERR "\n";

###### compute alignment and coverage statistics ######
my $total                  = 0;
my $total_mixed_end_orientation = 0;
my $total_unhandled        = 0;
my $total_smaller_than_20k = 0;    #for alignments covering < 20k on ref
my $total_mixed            = 0;
my $total_noncolinear      = 0;
my $total_over             = 0;
my $total_alt              = 0;
my $total_full_length      = 0;
my $total_to_end           = 0;
my $total_extend           = 0;
my $total_ref_covered      = 0;
my $total_ref_Ns_covered   = 0;
my $total_query_Ns_covered   = 0;
my $total_complete_contig_gaps_covered       = 0;
my $total_complete_contig_gap_length_covered = 0;
my $total_partial_contig_gaps_covered        = 0;
my $total_partial_contig_gap_length_covered  = 0;
my $total_complete_chr_gaps_covered          = 0;
my $total_complete_chr_gap_length_covered    = 0;
my $total_partial_chr_gaps_covered           = 0;
my $total_partial_chr_gap_length_covered     = 0;

my @lines       = read_file('nucmer.coords.delta.filtered.coords');
my $startline   = 5;
my $currentline = 0;
my @alignment_coords_array;
my @query_valid_hits;
my @query_invalid_hits;
my $last_line_query_id;
my $last_query_id;
my $last_query_length;

if ($print_header eq 'T'){
	print
"query\treference\tref_start\tref_end\tlength\tq_start\tq_end\tq_length\tseq_in_clusters\tdirection\tref_count\tincludes_0\tfull_length\tfrom_start\tfrom_end\tinternal_gap\tis_overlapping\tsize_of_alt\talternates\t\n";
	print MIXED
	"query\treference\tref_start\tref_end\tlength\tq_start\tq_end\tq_length\tseq_in_clusters\tdirection\tref_count\tincludes_0\tfull_length\tfrom_start\tfrom_end\tinternal_gap\tis_overlapping\tsize_of_alt\talternates\t\n";
	print NONCOLINEAR
	"query\treference\tref_start\tref_end\tlength\tq_start\tq_end\tq_length\tseq_in_clusters\tdirection\tref_count\tincludes_0\tfull_length\tfrom_start\tfrom_end\tinternal_gap\tis_overlapping\tsize_of_alt\talternates\t\n";
		
}

#parse nucmer.coords.delta.filtered.coords, convert BAC right end coordinates to BAC   
#[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[FRM]	[TAGS]
my ($left_reference_start_coord, $left_reference_end_coord, $left_aligned, $left_aligned_reference_direction, $left_aligned_query_direction);
$left_aligned = 0;
foreach my $line (@lines) {
	$currentline++;
	if ( $currentline < $startline ) {
		next;
	}
	my @row;
	@row = split( '\t', $line );
	
	my $query_name =  $row[14];
	$query_name =~ s/_\w+_\d+$//; #trim BAC end suffix 
	
	my $current_query_id     = $query_name;
	#my $current_query_length = $row[8];
	my $current_query_length = $query_lengths{$query_name};
	
	if ( !defined($last_line_query_id) ) {
		$last_line_query_id = $current_query_id;
	}

	#exec if query ID changes, i.e., coords for next assembled or singleton BAC aligned to chr
	if ( !( $current_query_id eq $last_line_query_id ) ) {

		#print info for prev query (assembled or singleton BAC ) to STDOUT
		if (scalar @alignment_coords_array > 0 ){
			calc_and_print_info( \@alignment_coords_array, $last_query_id,
							 $last_query_length, $bacend_length );
		}
		@alignment_coords_array = ();
	}
	
	if ( $row[14] =~ /_left_\d+/){
		$left_reference_start_coord = $row[0];
		$left_reference_end_coord = $row[1];
				
		#package BAC end coords as BAC coords
		my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
		$aln_coords->set_reference_id( $row[13] );
		$aln_coords->set_query_id( $query_name );	
		$aln_coords->set_reference_start_coord( $row[0] );
		$aln_coords->set_reference_end_coord( $row[1] );
		$aln_coords->set_reference_strand( $row[11] );
		$aln_coords->set_query_start_coord( $row[2] );
		$aln_coords->set_query_end_coord( $row[3] );
		$aln_coords->set_query_strand( $row[12] );
		push( @alignment_coords_array, $aln_coords );
		
		$left_aligned = 1;#found left BAC end alignment
		$left_aligned_reference_direction = $row[11];
		$left_aligned_query_direction =  $row[12];
	}
	#convert right BAC end coordinates into BAC coordinate space
	#right BAC end is always after the left BAC end in query_bacends.fasta
	elsif( $row[14] =~ /_right_\d+/){
		my ( $right_reference_start_coord, $right_reference_end_coord, $right_aligned_reference_direction, $right_aligned_query_direction, $query_start_coord, $query_end_coord);
		$right_reference_start_coord = $row[0];
		$right_reference_end_coord = $row[1];
		$right_aligned_reference_direction = $row[11]; 
		$right_aligned_query_direction = $row[12];

		if ($right_aligned_query_direction == 1){
			$query_start_coord = $query_lengths{$query_name} - $bacend_length + 1 ;
			$query_end_coord = $query_lengths{$query_name};	
		}
		elsif($right_aligned_query_direction == -1){
			$query_start_coord = $query_lengths{$query_name};
			$query_end_coord = $query_lengths{$query_name} - $bacend_length + 1 ;
		}
		else{
			die "Invalid direction for query. Exiting...\n";
		}
		

		if ($left_aligned){
			#if BAC ends from +ive strand
			#alignments on the reference will be on +ive strand
#			if (($left_aligned_reference_direction == 1 )
#				&& ($left_aligned_query_direction == 1 )
#				&& ($right_aligned_reference_direction == 1 )
#				&& ($right_aligned_query_direction == 1 )){
			#if BAC ends align on the same strand
			#alignments on the reference on same strand
			if (($left_aligned_reference_direction == $right_aligned_reference_direction )
				&& ($left_aligned_query_direction == $right_aligned_query_direction )){
				
				#check if BAC ends align within range +- 5% of BAC length
#				my ($min_reference_aligned_length, $max_reference_aligned_length, $reference_aligned_length);
#				$reference_aligned_length =  $right_reference_end_coord - $right_reference_start_coord + 1 ;
#				
#				$min_reference_aligned_length = $query_lengths{$query_name} - (0.05 * $query_lengths{$query_name});
#				$max_reference_aligned_length = $query_lengths{$query_name} + (0.05 * $query_lengths{$query_name});
#				
#				print STDERR "\n******** ",join("\t",$min_reference_aligned_length, $max_reference_aligned_length, $reference_aligned_length),"\n";
#				
#				#if valid query alignment
#				if (($reference_aligned_length >= $min_reference_aligned_length)
#				&& ($reference_aligned_length <= $max_reference_aligned_length)){
#					my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
#					$aln_coords->set_reference_id( $row[13] );
#					$aln_coords->set_query_id( $query_name );
#					$aln_coords->set_reference_start_coord( $right_reference_start_coord );
#					$aln_coords->set_reference_end_coord( $right_reference_end_coord );
#					$aln_coords->set_query_start_coord( $query_start_coord );
#					$aln_coords->set_query_end_coord( $query_end_coord );
#					push( @alignment_coords_array, $aln_coords );					
#				}
				my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
				$aln_coords->set_reference_id( $row[13] );
				$aln_coords->set_query_id( $query_name );
				$aln_coords->set_reference_start_coord( $right_reference_start_coord );
				$aln_coords->set_reference_end_coord( $right_reference_end_coord );
				$aln_coords->set_reference_strand( $right_aligned_reference_direction );
				$aln_coords->set_query_start_coord( $query_start_coord );
				$aln_coords->set_query_end_coord( $query_end_coord );
				$aln_coords->set_query_strand( $right_aligned_query_direction );
				push( @alignment_coords_array, $aln_coords );	
				
				if(($left_aligned_query_direction == -1 )
				|| ($right_aligned_query_direction == -1 )){
					print STDERR "\nComplimentary strand of BAC end from $query_name aligns to ref chr.\n";
				}

			}
			#if BAC was from -ive strand
			#alignments on the reference will be on -ive strand
#			elsif(($left_aligned_reference_direction == -1 )
#				&& ($left_aligned_query_direction == 1 )
#				&& ($right_aligned_reference_direction == -1 )
#				&& ($right_aligned_query_direction == 1 )){
#				
				#check if BAC ends align within range +- 5% of BAC length
#				my ($min_reference_aligned_length, $max_reference_aligned_length, $reference_aligned_length);
#				$reference_aligned_length =  $right_reference_end_coord - $right_reference_start_coord + 1 ;
#				
#				$min_reference_aligned_length = $query_lengths{$query_name} - (0.05 * $query_lengths{$query_name});
#				$max_reference_aligned_length = $query_lengths{$query_name} + (0.05 * $query_lengths{$query_name});
#				
#				#if valid query alignment
#				if (($reference_aligned_length >= $min_reference_aligned_length)
#				&& ($reference_aligned_length <= $max_reference_aligned_length)){
#					my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
#					$aln_coords->set_reference_id( $row[13] );
#					$aln_coords->set_query_id( $query_name );
#					$aln_coords->set_reference_start_coord( $right_reference_start_coord );
#					$aln_coords->set_reference_end_coord( $right_reference_end_coord );
#					$aln_coords->set_query_start_coord( $query_start_coord );
#					$aln_coords->set_query_end_coord( $query_end_coord );
#					push( @alignment_coords_array, $aln_coords );					
#				}
#				
#				my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
#				$aln_coords->set_reference_id( $row[13] );
#				$aln_coords->set_query_id( $query_name );
#				$aln_coords->set_reference_start_coord( $right_reference_start_coord );
#				$aln_coords->set_reference_end_coord( $right_reference_end_coord );
#				$aln_coords->set_query_start_coord( $query_start_coord );
#				$aln_coords->set_query_end_coord( $query_end_coord );
#				push( @alignment_coords_array, $aln_coords );
#			}
#			#if BAC ends align -ive strand
#			#should be fine as long as ref dir is same for both ends			
#			elsif(($left_aligned_reference_direction == $right_aligned_reference_direction )
#				&& ($left_aligned_query_direction == -1 )
#				&& ($right_aligned_query_direction == -1 )){
#				print STDERR "\nComplimentary strand of BAC end from $query_name aligns to ref chr. See mummer output files.\n";
#				my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
#				$aln_coords->set_reference_id( $row[13] );
#				$aln_coords->set_query_id( $query_name );
#				$aln_coords->set_reference_start_coord( $right_reference_start_coord );
#				$aln_coords->set_reference_end_coord( $right_reference_end_coord );
#				$aln_coords->set_query_start_coord( $query_start_coord );
#				$aln_coords->set_query_end_coord( $query_end_coord );
#				push( @alignment_coords_array, $aln_coords );				
#				
#				## TODO: Combine the 3 conditions as no diff in logic??
#			}
			elsif(($left_aligned_query_direction != $right_aligned_query_direction)
				|| ($left_aligned_reference_direction != $right_aligned_reference_direction)){
				print STDERR "\nBAC ends or ref orientation is opposite for $query_name. Ignoring. Can be misassembly in ref chr.";
				$total_mixed_end_orientation++;
			}
			else{
				print STDERR "\nUnhandled case from $query_name";
				$total_unhandled++;
			}
				
			#prep for next BAC end pair
			$left_aligned = 0;
			undef $left_reference_start_coord;
			undef $left_reference_start_coord;
			undef $left_reference_end_coord; 
			undef $left_aligned;
			undef $left_aligned_reference_direction;
			undef $left_aligned_query_direction;
		}
		else{# solo right BAC end alignment
			my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
			$aln_coords->set_reference_id( $row[13] );
			$aln_coords->set_query_id( $query_name );
			$aln_coords->set_reference_start_coord( $right_reference_start_coord );
			$aln_coords->set_reference_end_coord( $right_reference_end_coord );
			$aln_coords->set_reference_strand( $right_aligned_reference_direction );
			$aln_coords->set_query_start_coord( $query_start_coord );
			$aln_coords->set_query_end_coord( $query_end_coord );
			$aln_coords->set_query_strand( $right_aligned_query_direction );
			push( @alignment_coords_array, $aln_coords );				
		}
	}
		
	#deal with last row since no more alignments for query after this
	if ( $currentline == scalar(@lines) ) {
		if (scalar @alignment_coords_array > 0 ){
			calc_and_print_info( \@alignment_coords_array, $current_query_id,
							 $current_query_length, $bacend_length );
		}
		@alignment_coords_array = ();
	}
	$last_line_query_id = $current_query_id;
	$last_query_id      = $current_query_id;
	$last_query_length  = $current_query_length;
}

#cleanup
unlink "${opt_r}.index";
unlink "${opt_q}.index";
close(MIXED);
if ($total_mixed == 0) { unlink "mixed_qry_${opt_q}_ref_${opt_r}_group_coords.out";}
close(NONCOLINEAR);
if ($total_noncolinear == 0) { unlink "noncolinear_qry_${opt_q}_ref_${opt_r}_group_coords.out";}

##summary info
print STDERR "\n\nTotal queries:\t\t\t\t\t\t\t$total\n";
print STDERR "Total ignored ends/ref dir is opposite:\t\t\t\t$total_mixed_end_orientation\n";
print STDERR "Total unhandled and ignored queries:\t\t\t\t$total_unhandled\n";
print STDERR 
"Total queries with alignments smaller than 20,000 on ref:\t$total_smaller_than_20k\n";
print STDERR "Total queries with mixed orientation:\t\t\t\t$total_mixed\n";
print STDERR "Total queries with non co-linear alignments:\t\t\t$total_noncolinear\n";
print STDERR
  "Total queries with overlapping alignment clusters:\t\t$total_over\n";
print STDERR
  "Total queries with alternate alignments > 10,000:\t\t$total_alt\n";
print STDERR "Total queries aligned full length:\t\t\t\t$total_full_length\n";
print STDERR
  "Total queries with alignment to at least one end:\t\t$total_to_end\n";

print STDERR "Total reference extended by valid BAC hits:\t\t\t$total_extend\n"
  ;    #new seqs from query
print STDERR
  "Total reference covered by valid BAC hits:\t\t\t$total_ref_covered\n"
  ;    #includes gaps ($gap_size_allowed) between alignment clusters so BAC ends + seq in middle
print STDERR
"Total N's within reference covered by valid BAC hits:\t\t$total_ref_Ns_covered\n"
  ;    #includes gaps ($gap_size_allowed) between alignment clusters so BAC ends + seq in middle
print STDERR
"Total N's within queries covered by valid reference hits:\t$total_query_Ns_covered\n"
  ;    #includes gaps ($gap_size_allowed) between alignment clusters
  
#if ( $total_extend > $total_ref_Ns_covered ){ #new sequence beyond ends of chromosome
#	print STDERR "Total novel sequence beyond chr ends from valid BAC hits:\t",$total_extend - $total_ref_Ns_covered,"\n";
#}

if ($opt_c && $opt_s){
	print STDERR "\nStatistics from AGPs\n";
	print STDERR "Contig or component AGP (contigs and contig gaps)\n";
	print STDERR
	"\tTotal gaps completely covered from contig AGP:\t\t\t$total_complete_contig_gaps_covered\n";
	print STDERR
	"\tTotal length of gaps completely covered from contig AGP:\t$total_complete_contig_gap_length_covered\n";
	if ( $total_complete_contig_gaps_covered > 0){
		print STDERR "\tAvg length of gaps completely covered from contig AGP:\t\t"
	  . $total_complete_contig_gap_length_covered /
	  $total_complete_contig_gaps_covered . "\n";
	}
	print STDERR
	"\tTotal gaps partially covered from contig AGP:\t\t\t$total_partial_contig_gaps_covered\n";
	print STDERR
	"\tTotal length of gaps partially covered from contig AGP:\t\t$total_partial_contig_gap_length_covered\n";
	
	print STDERR "Chromosome AGP (scaffolds and scaffold gaps)\n";
	print STDERR
	"\tTotal gaps completely covered from chr AGP:\t\t\t$total_complete_chr_gaps_covered\n";
	print STDERR
	"\tTotal length of gaps completely covered from chr AGP:\t\t$total_complete_chr_gap_length_covered\n";
	if ($total_complete_chr_gaps_covered > 0){
		print STDERR "Avg length of gaps completely covered from chr AGP:\t\t\t"
	  . $total_complete_chr_gap_length_covered / $total_complete_chr_gaps_covered
	  . "\n";
	}
	print STDERR
	"\tTotal gaps partially covered from chr AGP:\t\t\t$total_partial_chr_gaps_covered\n";
	print STDERR
	"\tTotal length of gaps partial covered from chr AGP:\t\t$total_partial_chr_gap_length_covered\n";
}


=item C<calc_and_print_info (@alignment_coords_array, $last_query_id, $last_query_length, $bacend_length)>

Prints info to STDOUT. No return value.

=cut

sub calc_and_print_info {
	my ( $aref, $q_id, $q_length, $bacend_length ) = @_;
	my $align_group = Bio::GenomeUpdate::AlignmentCoordsGroup->new();

	#assign all coords for query/assembled or singleton BAC to obj
	$align_group->set_array_of_alignment_coords($aref);
	my $zero_chromosome_id = $unmapped_ID;
	#only checking for max alignment length, no check for BAC ends that align too close
	my $gap_size_allowed = $q_length + ( 0.05 * $q_length) - (2 * $bacend_length);#allow gap of up to query length + 5% between BAC ends
	
	#print STDERR "\n******* gap_size_allowed $gap_size_allowed\n";

	#Returns IDs, start and end coordinates, total aligned sequence length, and direction for the longest proximity-grouped
	#alignment clusters sorted by longest to shortest length of non-overlapping sequence covered by alignment clusters.
	#The proximity grouping is done using the specified length of an allowed gap between aligned clusters in the reference sequence ($gap allowed).
	#$sequence_aligned_in_clusters is the total length of ref covered by alignment grp
	my (
		 $ref_id,                       $query_id,
		 $ref_start,                    $ref_end,
		 $query_start,                  $query_end,
		 $sequence_aligned_in_clusters, $direction,
		 $colinear_order_check,		
		 $is_overlapping,               $size_of_next_largest_match,
		 $alternates
	  )
	  = $align_group
	  ->get_id_coords_and_direction_of_longest_alignment_cluster_group(
															 $gap_size_allowed);
	my $is_full_length;
	my $start_gap_length =
	  $query_start - 1;    #region of query before aligned part
	my $end_gap_length =
	  $q_length - $query_end;    #region of query after aligned part
	 #for calculating any gaps within the BAC alignment. It doesn’t count gaps at the beginning or end of the BAC
	 #(those may be from the BAC extending beyond the contig. Large gaps within the BAC could indicate a problem
	my $internal_gap_length =
	  ( $q_length - $sequence_aligned_in_clusters ) -
	  ( $start_gap_length + $end_gap_length );

	if ( ( $query_start == 1 ) && ( $query_end == $q_length ) ) {
		$is_full_length =
		  "Contains";    #entire query is covered in the alignment grp
	}
	else {
		$is_full_length = "Partial";
	}

	my $flagged = 0;    #flag 1 for potential problem

	$total++;

	#problem alignments
	if ( $direction == 0 ) {    #query aligns to both + and - strand of ref
		$total_mixed++;
		$flagged = 1;

		print MIXED $q_id . "\t";
		print MIXED $ref_id . "\t";
		print MIXED $ref_start . "\t";
		print MIXED $ref_end . "\t";
		print MIXED $ref_end - $ref_start . "\t";
		print MIXED $query_start . "\t";
		print MIXED $query_end . "\t";
		print MIXED $q_length . "\t";
		print MIXED $sequence_aligned_in_clusters . "\t";
		print MIXED $direction . "\t";    #strand
		print MIXED $align_group->get_count_of_reference_sequence_ids() . "\t";
		print MIXED $align_group->includes_reference_id($zero_chromosome_id)
		  . "\t";
		print MIXED $is_full_length . "\t";
		print MIXED $start_gap_length . "\t";
		print MIXED $end_gap_length . "\t";
		print MIXED $internal_gap_length . "\t";
		print MIXED $is_overlapping . "\t";
		print MIXED $size_of_next_largest_match . "\t";
		print MIXED $alternates . "\t";
		print MIXED "\n";

	}
	if ( $colinear_order_check == 1 ) {    #query BACs aligned non co-linearly to ref chr, i.e. different order of aligned tiles on BAC and ref chr. Shows miassembly on ref chr or BAC
		$total_noncolinear++;
		$flagged = 1;

		print NONCOLINEAR $q_id . "\t";
		print NONCOLINEAR $ref_id . "\t";
		print NONCOLINEAR $ref_start . "\t";
		print NONCOLINEAR $ref_end . "\t";
		print NONCOLINEAR $ref_end - $ref_start . "\t";
		print NONCOLINEAR $query_start . "\t";
		print NONCOLINEAR $query_end . "\t";
		print NONCOLINEAR $q_length . "\t";
		print NONCOLINEAR $sequence_aligned_in_clusters . "\t";
		print NONCOLINEAR $direction . "\t";    #strand
		print NONCOLINEAR $align_group->get_count_of_reference_sequence_ids() . "\t";
		print NONCOLINEAR $align_group->includes_reference_id($zero_chromosome_id)
		  . "\t";
		print NONCOLINEAR $is_full_length . "\t";
		print NONCOLINEAR $start_gap_length . "\t";
		print NONCOLINEAR $end_gap_length . "\t";
		print NONCOLINEAR $internal_gap_length . "\t";
		print NONCOLINEAR $is_overlapping . "\t";
		print NONCOLINEAR $size_of_next_largest_match . "\t";
		print NONCOLINEAR $alternates . "\t";
		print NONCOLINEAR "\n";

	}
	if ( $is_overlapping == 1 ) {    #query alignment tiles overlap
		$total_over++;
		$flagged = 1;
	}
	if ( $size_of_next_largest_match > 10000 ) {
		$total_alt++;
		$flagged = 1;
	}
	if ( $ref_end - $ref_start < 20000 ) {
		$total_smaller_than_20k++;

		#$flagged=1;#short alignment may be for a BAC end so not always a negative
	}

	#good alignments
	if ( $start_gap_length < 10 && $end_gap_length < 10 && $flagged == 0 )
	{    #alignments cover query
		$total_full_length++;
	}
	if ( ( $start_gap_length < 10 || $end_gap_length < 10 ) && $flagged == 0 )
	{    #alignments cover till one end of query
		$total_to_end++;
	}
	if ( $flagged == 0 ) {

		#print STDERR "**",join(' ',$query_id, $ref_start, $ref_end, $query_start, $query_end, $sequence_aligned_in_clusters, $start_gap_length,
		#	$end_gap_length, $internal_gap_length, $start_gap_length + $end_gap_length + $internal_gap_length),"\n\n";

		$total_extend +=
		  $start_gap_length + $end_gap_length + $internal_gap_length;

		$total_ref_covered += $sequence_aligned_in_clusters;

		my $ref_aligned_seq = $ref_db->seq( $ref_id, $ref_start, $ref_end );
		$total_ref_Ns_covered += ( $ref_aligned_seq =~ tr/N// );
		$total_ref_Ns_covered += ( $ref_aligned_seq =~ tr/n// );
		
		my $query_aligned_seq = $query_db->seq( $query_id, $query_start, $query_end );
		$total_query_Ns_covered += ( $query_aligned_seq =~ tr/N// );
		$total_query_Ns_covered += ( $query_aligned_seq =~ tr/n// );

		if ($opt_c && $opt_s){
		my ($cov_gap_count,     $cov_gap_length,
			 $par_cov_gap_count, $par_cov_gap_length
			  ) = $contig_agp->get_gap_overlap( $ref_start, $ref_end );
			$total_complete_contig_gaps_covered       += $cov_gap_count;
			$total_complete_contig_gap_length_covered += $cov_gap_length;
			$total_partial_contig_gaps_covered        += $par_cov_gap_count;
			$total_partial_contig_gap_length_covered  += $par_cov_gap_length;
	
			($cov_gap_count,     $cov_gap_length,
			   $par_cov_gap_count, $par_cov_gap_length
			  ) = $chr_agp->get_gap_overlap( $ref_start, $ref_end );
			$total_complete_chr_gaps_covered       += $cov_gap_count;
			$total_complete_chr_gap_length_covered += $cov_gap_length;
			$total_partial_chr_gaps_covered        += $par_cov_gap_count;
			$total_partial_chr_gap_length_covered  += $par_cov_gap_length;
		}
		
		print $q_id. "\t";
		print $ref_id. "\t";
		print $ref_start. "\t";
		print $ref_end. "\t";
		print $ref_end - $ref_start . "\t";
		print $query_start. "\t";
		print $query_end. "\t";
		print $q_length. "\t";
		print $sequence_aligned_in_clusters. "\t";
		print $direction. "\t";    #strand
		print $align_group->get_count_of_reference_sequence_ids() . "\t";
		print $align_group->includes_reference_id($zero_chromosome_id) . "\t";
		print $is_full_length. "\t";
		print $start_gap_length. "\t";
		print $end_gap_length. "\t";
		print $internal_gap_length. "\t";
		print $is_overlapping. "\t";
		print $size_of_next_largest_match. "\t";
		print $alternates. "\t";
		print "\n";
			
	}
}

=item C<archive_old ($file, $message)>

Checks if file exists. If yes, moves to file.old. Returns 1 if no file found and 0 if found and moved.

=cut
sub archive_old{
	my ($file, $message)= @_;
	if ( -e $file){ 
		system("mv $file ${file}.old");
		die("\n$message. $!\nExiting...\n\n") if ($? == -1);
		return 0;
	}
	else{
		return 1;
	}
}

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script creates a fasta of BAC ends, aligns to reference with mummer and groups aligned clusters to create a tab delimited file with BAC alignment details. Mixed and out of order alignments are written to separate files. Only removing BACs whose ends align beyond range on ref, no check for BAC ends that align too close.

    Usage:
      align_BACends_group_coords.pl -r [fasta] -q [fasta] -c [agp] -s [agp]

    Flags:

    -u  Sequence ID of chromosome with unmapped contigs 
    -l  BAC end length (required)
    -r  Fasta file of reference (required)
    -q  Fasta file of query (assembled and singleton BACs)
    -c  Contig or component AGP file for reference (includes scaffold gaps)
    -s  Chromosome AGP file for reference (with only scaffolds and gaps) 
    -t  Print header
    -h  Help
 
EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha at cornell.edu, @SahaSurya>

=cut

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

