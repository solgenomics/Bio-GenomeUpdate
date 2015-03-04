#!/usr/bin/perl

=head1

group_coords.pl

=head1 SYNOPSIS

    align_BACends_group_coords.pl -g [gap size] -r [fasta] -q [fasta] -c [agp] -s [agp]

=head1 COMMAND-LINE OPTIONS

 -u  Sequence ID of chromosome with unmapped contigs 
 -r  Fasta file of reference (required)
 -q  Fasta file of query (assembled and singleton BACs, required)
 -g  Gap size allowed between aligned clusters in the reference sequence, typically the mean/median scaffold gap (required)
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
use Bio::GenomeUpdate::AlignmentCoords;
use Bio::GenomeUpdate::AlignmentCoordsGroup;
use Bio::GenomeUpdate::AGP;
use Bio::DB::Fasta;

our ( $opt_r, $opt_q, $opt_u, $opt_t, $opt_h, $opt_c, $opt_s );
getopts("r:q:u:c:s:t:h");
if ( !$opt_c || !$opt_s || !$opt_r || !$opt_q ) {
	print STDERR "Required parameters missing!! exiting..\n";
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
my $ref_db   = Bio::DB::Fasta->new( $opt_r, '-reindex' => 1 );
my $query_db = Bio::DB::Fasta->new( $opt_q, '-reindex' => 1 );


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


sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script groups aligned clusters and creates a tab delimited file with BAC alignment details.

    Usage:
      group_coords.pl -i [coords file] -g [gap size] -r [fasta] -q [fasta] -c [agp] -s [agp]

    Flags:

    -i  COORDS file created by show-coords (required)
    -u  Sequence ID of chromosome with unmapped contigs 
    -r  Fasta file of reference (required)
    -q  Fasta file of query (assembled and singleton BACs)
    -g  Gap size allowed between aligned clusters in the reference sequence, typically the mean/median scaffold gap (required)
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

