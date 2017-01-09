#!/usr/bin/perl

=head1 NAME

update_maker_names_gff.pl

=head1 SYNOPSIS

update_maker_names_gff.pl -i [old GFF file] -o [new GFF file] -n [New Solyc ids] -d 0

=head1

Gets Solyc id from the mRNA record and assigns to gene and mRNA records. Generates a new Solyc id if no old id was passed through and assigns to gene and mRNA records. Generates a list of new Solyc ids. Need to check if any new id overlaps with a deprecated gene.

=head1 COMMAND-LINE OPTIONS

 -i  old GFF file (required)
 -o  new GFF file with updated names for genes and mRNA (required)
 -n  List of new Solyc ids
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;
use Utilities;
use File::Slurp;
use Getopt::Std;
use Bio::GFF3::LowLevel
  qw (gff3_parse_feature  gff3_format_feature gff3_parse_attributes);


our ( $opt_i, $opt_o, $opt_n, $opt_d, $opt_h );
getopts('i:o:n:d:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_i || !$opt_o ) {
	print
"\nOld and new GFF3 are required. 
See help below\n\n\n";
	help();
}

#get input files
my $old_gff_input_file = $opt_i;
my $input_old_gff      = read_file($old_gff_input_file)
  or die "Could not open old GFF input file: $old_gff_input_file\n";
my $gff_output_file = $opt_o;
if (defined $opt_n){ 
	my $new_id_output_file   = $opt_n;
}
else{
	my $new_id_output_file   = $opt_o_{new_ids.names};
}


my $util = Utilities->new();
if ($opt_d) {
	print STDERR "Params parsed..\n";
	$util->mem_used();
	$util->run_time();
}

#object for gff

my @lines        = split( /\n/, $input_old_gff );
my $line_count   = scalar(@lines);
my $line_counter = 0;
my $gene_flag    = 0;
my @gene_gff_line_hashes;

foreach my $line (@lines) {
	$line_counter++;
	chomp($line);

	if ( $line =~ m/^#/ ) {
		next;
	}
	else {
		my $gff_line_hash = gff3_parse_feature($line);
	}
	print STDERR "\rParsing GFF3 line "
	  . $line_counter . " of "
	  . $line_count;
	
	## if first gene
	if ( ($gff_line_hash{'type'} eq 'gene') && ( $gene_flag == 0) ){
		$gene_flag = 1;
		push (@gene_gff_line_hashes, $gff_line_hash);
	}
	## if next gene
	elsif( ($gff_line_hash{'type'} eq 'gene') && ( $gene_flag == 1) ){
		#process last gene
		
		@gene_gff_line_hashes = (); # reset
	}
	
}

## last gene


if ($opt_d) {
	print STDERR "Input GFF file read..\n";
	$util->mem_used();
	$util->run_time();
}




###TODO

#print the GFF ($gff_output_file)
if ( !defined($new_gff) ) {
	print STDERR "No valid GFF records found..\n";
	exit 1;
} else {
	unless ( open( OGFF, ">$gff_output_file" ) ) {
		print STDERR "Cannot open $gff_output_file\n";
		exit 1;
	}
	print OGFF $new_gff;
	close(OGFF);
	if ($opt_d) {
		print STDERR "GFF written..\n";
		$util->mem_used();
		$util->run_time();
	}
}

#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     Gets Solyc id from the mRNA record and assigns to gene and mRNA records. Generates a new Solyc id if no old id was passed through and assigns to gene 
     and mRNA records. Generates a list of new Solyc ids. Need to check if any new id overlaps with a deprecated gene.

     
    Usage:
      update_maker_names_gff.pl -i [old GFF file] -o [new GFF file] -n [New Solyc ids] -d 0
      
    Flags:

     -i  old GFF file (required)
     -o  new GFF file with updated names for genes and mRNA (required)
     -n  List of new Solyc ids
     -d  debugging messages (1 or 0)
     -h  Help


EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

__END__


# MAKER GFF
##gff-version 3
##sequence-region   SL3.0ch00 16480 20797619
..
##sequence-region   SL3.0ch12 1065 68124245
#Gap=M130
SL3.0ch00	maker	gene	16480	17940	.	+	.	ID=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0;Name=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0
SL3.0ch00	maker	mRNA	16480	17940	.	+	.	ID=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1;Parent=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0;Name=Solyc00g005000.2.1;_AED=0.08;_eAED=0.08;_QI=0|0|0|1|1|1|2|0|458;Alias=1_t,Solyc00g005000.2.1
SL3.0ch00	maker	exon	16480	16794	.	+	.	ID=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1:1;Parent=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1
SL3.0ch00	maker	CDS	16480	16794	.	+	0	ID=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1:cds;Parent=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1
SL3.0ch00	maker	exon	16879	17940	.	+	.	ID=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1:2;Parent=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1
SL3.0ch00	maker	CDS	16879	17940	.	+	0	ID=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1:cds;Parent=maker-SL3.0ch00-pred_gff_SL3.0_RMmasked-gene-0.0-mRNA-1
###
SL3.0ch00	maker	gene	328352	334459	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4;Name=maker-SL3.0ch00-snap-gene-3.4
SL3.0ch00	maker	mRNA	328352	334459	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1;Parent=maker-SL3.0ch00-snap-gene-3.4;Name=maker-SL3.0ch00-snap-gene-3.4-mRNA-1;_AED=0.56;_eAED=0.64;_QI=0|0|0|0.57|0.16|0.57|7|0|266
SL3.0ch00	maker	exon	328352	328372	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:1;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	CDS	328352	328372	.	+	0	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:cds;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	exon	328439	328507	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:2;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	CDS	328439	328507	.	+	0	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:cds;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	exon	328538	328702	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:3;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	CDS	328538	328702	.	+	0	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:cds;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	exon	328940	329026	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:4;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	CDS	328940	329026	.	+	0	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:cds;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	exon	329196	329318	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:5;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	CDS	329196	329318	.	+	0	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:cds;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	exon	333732	333782	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:6;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	CDS	333732	333782	.	+	0	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:cds;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	exon	334175	334459	.	+	.	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:7;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
SL3.0ch00	maker	CDS	334175	334459	.	+	0	ID=maker-SL3.0ch00-snap-gene-3.4-mRNA-1:cds;Parent=maker-SL3.0ch00-snap-gene-3.4-mRNA-1
###
SL3.0ch00	maker	gene	526682	527116	.	-	.	ID=augustus-SL3.0ch00-processed-gene-5.3;Name=augustus-SL3.0ch00-processed-gene-5.3
SL3.0ch00	maker	mRNA	526682	527116	.	-	.	ID=augustus-SL3.0ch00-processed-gene-5.3-mRNA-1;Parent=augustus-SL3.0ch00-processed-gene-5.3;Name=augustus-SL3.0ch00-processed-gene-5.3-mRNA-1;_AED=0.19;_eAED=0.19;_QI=0|-1|0|1|-1|1|1|0|144
SL3.0ch00	maker	exon	526682	527116	.	-	.	ID=augustus-SL3.0ch00-processed-gene-5.3-mRNA-1:1;Parent=augustus-SL3.0ch00-processed-gene-5.3-mRNA-1
SL3.0ch00	maker	CDS	526682	527116	.	-	0	ID=augustus-SL3.0ch00-processed-gene-5.3-mRNA-1:cds;Parent=augustus-SL3.0ch00-processed-gene-5.3-mRNA-1
###

# ITAG GFF3
##gff-version 3
##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.93
##sequence-region SL2.50ch00 1 21805821
SL2.50ch00	ITAG_eugene	gene	16437	18189	.	+	.	Alias=Solyc00g005000;ID=gene:Solyc00g005000.2;Name=Solyc00g005000.2;from_BOGAS=1;length=1753
SL2.50ch00	ITAG_eugene	mRNA	16437	18189	.	+	.	ID=mRNA:Solyc00g005000.2.1;Name=Solyc00g005000.2.1;Note=Aspartic proteinase nepenthesin I (AHRD V1 **-- A9ZMF9_NEPAL)%3B contains Interpro domain(s)  IPR001461  Peptidase A1 ;Ontology_term=GO:0006508;Parent=gene:Solyc00g005000.2;from_BOGAS=1;interpro2go_term=GO:0006508;length=1753;nb_exon=2
SL2.50ch00	ITAG_eugene	exon	16437	17275	.	+	.	ID=exon:Solyc00g005000.2.1.1;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.50ch00	ITAG_eugene	five_prime_UTR	16437	16479	.	+	.	ID=five_prime_UTR:Solyc00g005000.2.1.0;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.50ch00	ITAG_eugene	CDS	16480	17275	.	+	0	ID=CDS:Solyc00g005000.2.1.1;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.50ch00	ITAG_eugene	intron	17276	17335	.	+	.	ID=intron:Solyc00g005000.2.1.1;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.50ch00	ITAG_eugene	exon	17336	18189	.	+	0	ID=exon:Solyc00g005000.2.1.2;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.50ch00	ITAG_eugene	CDS	17336	17940	.	+	2	ID=CDS:Solyc00g005000.2.1.2;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.50ch00	ITAG_eugene	three_prime_UTR	17941	18189	.	+	.	ID=three_prime_UTR:Solyc00g005000.2.1.0;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
