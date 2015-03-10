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
if (looks_like_number($opt_l)){
	$bacend_length = $opt_l;
}
else{
	die("$opt_l is not a number");
}
  

open( MIXED, ">mixed_qry_${opt_q}_ref_${opt_r}_group_coords.out" )
  or die
"Could not create mixed_qry_${opt_q}_ref_${opt_r}_group_coords.out for writing out BACs aligned to ref chr in mixed orientation";
open( NONCOLINEAR, ">noncolinear_qry_${opt_q}_ref_${opt_r}_group_coords.out" )
  or die
"Could not create noncolinear_qry_${opt_q}_ref_${opt_r}_group_coords.out for writing out BACs aligned non co-linearly to ref chr, i.e. different order of aligned tiles on BAC and ref chr. Shows miassembly on ref chr or BAC";

my $contig_agp_input_file = $opt_c;
my $contig_input_agp      = read_file($contig_agp_input_file)
  or die "Could not open contig AGP input file: $contig_agp_input_file\n";
my $contig_agp = Bio::GenomeUpdate::AGP->new();
$contig_agp->parse_agp($contig_input_agp);
my $chr_agp_input_file = $opt_s;
my $chr_input_agp      = read_file($chr_agp_input_file)
  or die "Could not open chr AGP input file: $chr_agp_input_file\n";
my $chr_agp = Bio::GenomeUpdate::AGP->new();
$chr_agp->parse_agp($chr_input_agp);

my $ref_db   = Bio::DB::Fasta->new( $opt_r, '-reindex' => 1 );
my $query_db = Bio::DB::Fasta->new( $opt_q, '-reindex' => 1 );
my %query_lengths;


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

archive_old('nucmer.coords','Could not remove old nucmer.coords');
#system('nucmer', '-l 100', "-c $bacend_length", '-p nucmer.coords', $opt_r, $opt_q );
system("nucmer --noextend -l 100 -c $bacend_length -p nucmer.coords $opt_r query_bacends.fasta");
die("\nCould not run nucmer. $!\nExiting...\n\n") if ($? == -1);

#is this required any more?? cluster size in nucmer should only report alignment >= $bacend_length 
archive_old('nucmer.coords.delta.filtered','Could not remove old nucmer.coords.delta.filtered');
system("delta-filter -l $bacend_length nucmer.coords.delta > nucmer.coords.delta.filtered");
die("\nCould not run delta-filter. $!\nExiting...\n\n") if ($? == -1);

archive_old('nucmer.coords.delta.filtered.coords','Could not remove old nucmer.coords.delta.filtered.coords');
system("show-coords -c -d -l -q -T -o nucmer.coords.delta.filtered > nucmer.coords.delta.filtered.coords");
die("\nCould not run show-coords. $!\nExiting...\n\n") if ($? == -1);




###### compute alignment and coverage statistics ######


#read nucmer.coords.delta.filtered.coords
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

#parse coords file, convert BAC right end coordinates to BAC   
#[S1]	[E1]	[S2]	[E2]	[LEN 1]	[LEN 2]	[% IDY]	[LEN R]	[LEN Q]	[COV R]	[COV Q]	[FRM]	[TAGS]
foreach my $line (@lines) {
	$currentline++;
	if ( $currentline < $startline ) {
		next;
	}
	my @row;
	@row = split( '\t', $line );
	my $current_query_id     = $row[14];
	my $current_query_length = $row[8];
	if ( !defined($last_line_query_id) ) {
		$last_line_query_id = $current_query_id;
	}

#exec if query ID changes, i.e., coords for next assembled or singleton BAC aligned to chr
	if ( !( $current_query_id eq $last_line_query_id ) ) {

		#print info for prev query (assembled or singleton BAC ) to STDOUT
		calc_and_print_info( \@alignment_coords_array, $last_query_id,
							 $last_query_length );
		@alignment_coords_array = ();
	}
	my $aln_coords = Bio::GenomeUpdate::AlignmentCoords->new();
	$aln_coords->set_reference_id( $row[13] );
	$aln_coords->set_query_id( $row[14] );
	$aln_coords->set_reference_start_coord( $row[0] );
	$aln_coords->set_reference_end_coord( $row[1] );
	$aln_coords->set_query_start_coord( $row[2] );
	$aln_coords->set_query_end_coord( $row[3] );
	push( @alignment_coords_array, $aln_coords );

	#deal with last row since no more alignments for query after this
	if ( $currentline == scalar(@lines) ) {

#calc_and_print_info(\@alignment_coords_array, $current_query_id, $current_query_length,$gap_size_allowed);
		calc_and_print_info( \@alignment_coords_array, $current_query_id,
							 $current_query_length );
		@alignment_coords_array = ();
	}
	$last_line_query_id = $current_query_id;
	$last_query_id      = $current_query_id;
	$last_query_length  = $current_query_length;
}









my $total                  = 0;
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
my $total_complete_contig_gaps_covered       = 0;
my $total_complete_contig_gap_length_covered = 0;
my $total_partial_contig_gaps_covered        = 0;
my $total_partial_contig_gap_length_covered  = 0;
my $total_complete_chr_gaps_covered          = 0;
my $total_complete_chr_gap_length_covered    = 0;
my $total_partial_chr_gaps_covered           = 0;
my $total_partial_chr_gap_length_covered     = 0;



print STDERR "Number of reference sequences: ";
print STDERR scalar $ref_db->get_all_ids();
print STDERR "\nNumber of query sequences: ";
print STDERR scalar $query_db->get_all_ids();
print STDERR "\n\n";

#cleanup
unlink "${opt_r}.index";
unlink "${opt_q}.index";
close(MIXED);
if ($total_mixed == 0) { unlink "mixed_qry_${opt_q}_ref_${opt_r}_group_coords.out";}
close(NONCOLINEAR);
if ($total_noncolinear == 0) { unlink "noncolinear_qry_${opt_q}_ref_${opt_r}_group_coords.out";}

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

     This script creates a fasta of BAC ends, aligns to reference with mummer and groups aligned clusters to create a tab delimited file with BAC alignment details. Mixed and out of order alignments are written to separate files.

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

