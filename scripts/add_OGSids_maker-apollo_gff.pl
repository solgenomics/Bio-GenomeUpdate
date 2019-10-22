#!/usr/bin/perl

=head1 NAME

add_OGSids_maker-apollo_gff.pl

=head1 SYNOPSIS

add_OGSids_maker-apollo_GFF.pl -g [old Merged Maker and Apollo GFF file] -a [AHRD file for mRNA with Maker and Apollo ids] -p [species acronym or prefix] -c [Prefix for chromosome in GFF]  -s [starting value for naming] -o [output formatted gff file with OGS ids]

=head1 COMMAND-LINE OPTIONS

 -g  Merged Maker and Apollo GFF file (required)
 -a  AHRD tab separated file for mRNA with Maker and Apollo ids, e.g. maker-DC3.0sc00-snap-gene-91.44-mRNA-1       Glycerol-3-phosphate dehydrogenase [NAD(P)+] (AHRD V3.11 *** tr|A0A0A9WH09|A0A0A9WH09_LYGHE) (required)
 -f  Curated descriptions from Apollo, e.g. 18712462-6e8c-478a-83c8-a40c9d0be977  Dihydrolipoyl transacetylase (required)
 -p  Prefix for name, e.g DcitrP (required)
 -c  Prefix for chromosome in GFF e.g. Dc3.0sc (required)
 -s  Starting seed for gene number for each chromosome, e.g. 1000 (required)
 -o  output GFF file with OGS ids
 -h  Help

=cut

use strict;
use warnings;

use File::Slurp;
use Getopt::Std;
use Bio::GFF3::LowLevel qw (gff3_parse_feature gff3_format_feature gff3_parse_attributes);
use Switch;

our ( $opt_g, $opt_a, $opt_f, $opt_p, $opt_s, $opt_c, $opt_o, $opt_h );
getopts('g:a:f:p:s:c:o:h');
if ($opt_h) {
  help();
  exit;
}
# if ( !$opt_g || !$opt_a || !$opt_f || !$opt_p || !$opt_c || !$opt_s || !$opt_o) {
#   print
# "\nOld GFF file, AHRD file, Apollo descriptions, name prefix, chr prefix, starting seed, output GFF file is required.
# See help below\n\n\n";
#   help();  
# }

# get input files
my $gff_input_file = $opt_g;
my $gff_input   = read_file($gff_input_file)
  or die "Could not open old gff input file: $gff_input_file\n";
my $ahrd_input_file = $opt_a;
my $ahrd_input      = read_file($ahrd_input_file)
  or die "Could not open AHRD input file: $ahrd_input_file\n";
my $apollo_desc_input_file = $opt_f;
my $apollo_desc_input      = read_file($apollo_desc_input_file)
  or die "Could not open Apollo description input file: $apollo_desc_input_file\n";
chomp $opt_p; my $prefix = $opt_p;
chomp $opt_c; my $chrprefix = $opt_c;
# my $seed = $opt_s;
# if ($seed !~ /^[0-9]+$/){ die "$seed should be a number\n"; }
my $seed = 990; #hack to start from 990 to match names in /export/species2/Diaphorina_citri/annotation/OGSv3.0/ahrd_final_OGS3_IDS.txt

# hash of Apollo descriptions
my %apollo_curated_function;
my @lines = split( /\n/, $apollo_desc_input );
foreach my $line (@lines) {
	chomp($line);
	my @line_arr = split ("\t", $line);
	$apollo_curated_function{$line_arr[0]}=$line_arr[1];            #apollo id = curated function
}

# hash of AHRD descriptions
my (%ahrd_function, %ahrd_domain);
@lines = split( /\n/, $ahrd_input );
foreach my $line (@lines) {
	chomp($line);
	my @line_arr = split ("\t", $line);
	#$ahrd_function{$line_arr[0]}=$line_arr[1];                      # maker id = AHRD function string
	$ahrd_function{$line_arr[0]}=$line_arr[1];                      # temp hack new OGSv3 id = AHRD function string
	if ( defined $line_arr[2] ){
		#$ahrd_domain{$line_arr[0]}=$line_arr[2];					# maker id = AHRD comma separated domain string
		$ahrd_domain{$line_arr[0]}=$line_arr[2];					# temp hack new OGSv3 id = AHRD comma separated domain string
	}

}

# output variables
my ($gff_output, $index_output, $desc_output);
# tracking variables
my (%scaffold_last_gene_id, %gene_old_id_mrna_last_rank, %gene_old_new_index, %mrna_old_new_index);
my (%mrna_old_id_cds_last_rank, %mrna_old_id_exon_last_rank, %mrna_old_id_threeprimeutr_last_rank, %mrna_old_id_fiveprimeutr_last_rank);

$gff_output   = '';                                                 # initialize
$index_output = '';
$desc_output  = '';



@lines = split( /\n/, $gff_input );
foreach my $line (@lines){
	# comments line
	if ( $line=~ m/^#/ ){
		$gff_output .= $line."\n";

		if ( $line=~ m/^##sequence-region/ ){						#add seqid to hash assuming that sequence-region is present in a properly formatted GFF file
			##sequence-region   DC3.0sc00 1 34736251
			my @seq_region = split ( ' ', $line);
			$seq_region[1] =~ s/$chrprefix//;						#remove genome version and chr, so only chr number
			$scaffold_last_gene_id{$seq_region[1]} = $seed; 		#init counter from 990 to match with logic in update_maker_names_fasta_chrlocs.pl???
		}

		next;
	}

	my $gff_features = gff3_parse_feature ($line);

	if ( $gff_features->{'type'} eq 'gene' ){						# if gene, create the OGS id Dcitr00g00990.1.1
		my $scaffold = $gff_features->{'seq_id'};
		$scaffold    =~ s/$chrprefix//;								# remove genome version and chr, so only chr number

		my $gene_id;
		if ( $scaffold_last_gene_id{$scaffold} == $seed ) {			# first gene on scaffold
			$gene_id = $seed;
		}
		else{
			$scaffold_last_gene_id{$scaffold} += 10;				# namespace for 9 new genes now
			$gene_id = $scaffold_last_gene_id{$scaffold};
		}

		my $gene_new_id;

		my @chars = split //,$gene_id;
		my $padding_count = 5 - scalar @chars;						#presuming a gene space of <1,00,000 per scaffold so 6 characters, e.g. 00001 - 99999, was 1 mill earlier
		
		$gene_new_id = $prefix.$scaffold . 'g';						# prefix with chr number and 0's 
		foreach (1..$padding_count){
			$gene_new_id = $gene_new_id . '0';
		}
		$gene_new_id = $gene_new_id . $gene_id . '.1';				# assign version 1

		my $old_id = $gff_features->{'attributes'}->{'ID'}->[0];	# add gene_new_id to gene index
		$gene_old_new_index{$old_id} = $gene_new_id;

		# create the new GFF record
		my $gene_attributes_hashref = gff3_parse_attributes ("ID=$gene_new_id;Name=$gene_new_id");
		$gff_features->{'attributes'} = $gene_attributes_hashref;

		$gff_output = $gff_output . gff3_format_feature ($gff_features);
	}
	elsif ( $gff_features->{'type'} eq 'mRNA' ){					# if mRNA, increment the mRNA isoform count if this is the not the first mRNA for this gene
		
		#create mrna id
		my $mrna_parent_new_id = $gene_old_new_index{$gff_features->{'attributes'}->{'Parent'}->[0]};
		my $mrna_rank;
		if ( exists $gene_old_id_mrna_last_rank{$gff_features->{'attributes'}->{'Parent'}->[0]} ){
			$mrna_rank = $gene_old_id_mrna_last_rank{$gff_features->{'attributes'}->{'Parent'}->[0]}++;
		}
		else{
			$gene_old_id_mrna_last_rank{$gff_features->{'attributes'}->{'Parent'}->[0]} = $mrna_rank = 1;
		}
		my $mrna_new_id = $mrna_parent_new_id . '.' . $mrna_rank;


		# # create the func description string
		# # Only Apollo or AHRD desc, adding domains after | separator without IPR,GO,Pfam prefixes
		# # can have special characters, not adding -RA -RB automatically for curated genes, that should be done by curators
		my $mrna_desc;
		if ( exists $apollo_curated_function{$gff_features->{'attributes'}->{'ID'}->[0]} ) {		# get the Apollo description if it exists
			$mrna_desc = $apollo_curated_function{$gff_features->{'attributes'}->{'ID'}->[0]};
		}
		else{
			#$mrna_desc = $ahrd_function{$gff_features->{'attributes'}->{'ID'}->[0]};				# get AHRD description
			$mrna_desc = $ahrd_function{$mrna_new_id};										# temp hack to get AHRD description using new OGSv3 id
		}

		my $mrna_domain;
		if ( exists $ahrd_domain{$gff_features->{'attributes'}->{'ID'}->[0]} ){					# add domains
			$mrna_domain = $ahrd_domain{$gff_features->{'attributes'}->{'ID'}->[0]};
			$mrna_domain =~ s/^,//;															# remove extra , e.g. ,PF12698
			chomp $mrna_domain;
		}
		if ( defined $mrna_domain ){
			$mrna_desc = $mrna_desc. ' | ' .$mrna_domain;
		}

		my $mrna_old_id = $gff_features->{'attributes'}->{'ID'}->[0];								# add mrna new id to mrna index
		$mrna_old_new_index{$mrna_old_id} = $mrna_new_id;

		if ( $gff_features->{'source'} eq 'maker'){													# maker mRNA
			my $mrna_aed  = $gff_features->{'attributes'}->{'_AED'}->[0];
			my $mrna_eaed = $gff_features->{'attributes'}->{'_eAED'}->[0];
			my $mrna_qi   = $gff_features->{'attributes'}->{'_QI'}->[0];

			# write the new mRNA record
			my $mrna_attributes_hashref = gff3_parse_attributes ("ID=$mrna_new_id;Name=$mrna_new_id;Note=$mrna_desc;Parent=$mrna_parent_new_id;_AED=$mrna_aed;_eAED=$mrna_eaed;;_QI=$mrna_qi");
			$gff_features->{'attributes'} = $mrna_attributes_hashref;
		}
		else{																						# apollo mRNA
			# write the new mRNA record
			my $mrna_attributes_hashref = gff3_parse_attributes ("ID=$mrna_new_id;Name=$mrna_new_id;Note=$mrna_desc;Parent=$mrna_parent_new_id");
			$gff_features->{'attributes'} = $mrna_attributes_hashref;
		}

		$gff_output = $gff_output . gff3_format_feature ($gff_features);
	}
	# if any other child record (exon, UTRs, CDS)
	# use the mRNA parent to create the gff record, can have multiple mRNA parents
	elsif ( ($gff_features->{'type'} eq 'CDS') || ($gff_features->{'type'} eq 'exon')
			|| ($gff_features->{'type'} eq 'three_prime_UTR') || ($gff_features->{'type'} eq 'five_prime_UTR')){

		my $child_parent_new_id = '';
		my $child_single_mrna_parent_old_id;												# get first or only parent
		if ( scalar @{$gff_features->{'attributes'}->{'Parent'}} > 1 ){						# multiple parents in arrayref
			my $child_parents_arref = $gff_features->{'attributes'}->{'Parent'};

			foreach my $child_parent (@{$child_parents_arref}){
				if ( length $child_parent_new_id > 0 ){										# if not first parent
					$child_parent_new_id = $child_parent_new_id. ',' .$mrna_old_new_index{$child_parent};
				}
				else{																		# if first parent
					$child_parent_new_id = $mrna_old_new_index{$child_parent};
					$child_single_mrna_parent_old_id  = $child_parent;
				}
			}
		}
		else{																				# single parent in arrayref
			$child_parent_new_id = $mrna_old_new_index{$gff_features->{'attributes'}->{'Parent'}->[0]};
			$child_single_mrna_parent_old_id  = $gff_features->{'attributes'}->{'Parent'}->[0] ;
		}

		# create child new id and record
		my $child_rank;
		my $child_attributes_hashref;
		switch ( $gff_features->{'type'} ){
			case 'CDS'{
				if ( exists $mrna_old_id_cds_last_rank {$child_single_mrna_parent_old_id} ){
					$child_rank = $mrna_old_id_cds_last_rank{$child_single_mrna_parent_old_id}++;
				}
				else{
					$mrna_old_id_cds_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				}

				# create id
				my $cds_new_id = 'CDS:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$cds_new_id;Name=$cds_new_id;Parent=$child_parent_new_id");
			}
			case 'exon'{
				if ( exists $mrna_old_id_exon_last_rank {$child_single_mrna_parent_old_id} ){
					$child_rank = $mrna_old_id_exon_last_rank{$child_single_mrna_parent_old_id}++;
				}
				else{
					$mrna_old_id_exon_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				}

				# create id
				my $exon_new_id = 'exon:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$exon_new_id;Name=$exon_new_id;Parent=$child_parent_new_id");
			}
			case 'three_prime_UTR'{
				if ( exists $mrna_old_id_threeprimeutr_last_rank {$child_single_mrna_parent_old_id} ){
					$child_rank = $mrna_old_id_threeprimeutr_last_rank{$child_single_mrna_parent_old_id}++;
				}
				else{
					$mrna_old_id_threeprimeutr_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				}

				# create id
				my $threeprimeutr_new_id = 'three_prime_UTR:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$threeprimeutr_new_id;Name=$threeprimeutr_new_id;Parent=$child_parent_new_id");
			}
			case 'five_prime_UTR'{
				if ( exists $mrna_old_id_fiveprimeutr_last_rank {$child_single_mrna_parent_old_id} ){
					$child_rank = $mrna_old_id_fiveprimeutr_last_rank{$child_single_mrna_parent_old_id}++;
				}
				else{
					$mrna_old_id_fiveprimeutr_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				}

				# create id
				my $fiveprimeutr_new_id = 'five_prime_UTR:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$fiveprimeutr_new_id;Name=$fiveprimeutr_new_id;Parent=$child_parent_new_id");
			}
		}

		# write the new child record
		$gff_features->{'attributes'} = $child_attributes_hashref;
		$gff_output = $gff_output . gff3_format_feature ($gff_features);
	}
}





# write output files
chomp $opt_o;
open my $OGFF, '>', "$opt_o" or die "Cannot open $opt_o\n";
open my $OINDEX, '>', "index_mRNA.$opt_o" or die  "Cannot open index_mRNA.$opt_o\n";
open my $ODESC, '>', "desc_mRNA.$opt_o" or die "Cannot open desc_mRNA.$opt_o\n";

print $OGFF $gff_output;
close $OGFF;
print $OINDEX $index_output;
close $OINDEX;
print $ODESC $desc_output;
close $ODESC;


#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

	Description:

	 Renames maker and Apollo assigned gene names to gene name with version numbers, e.g.  maker-ScVcwli_1-pred_gff_maker-gene-0.0-mRNA-1 becomes DcitrP00001.1.1. Creates an index file with old and new mRNA ids. This is hard coded for <99,999 mRNAs per scaffold. Counter skips over 10 gene models so manually curated genes can be added later.
	 
	Output:
	 Index: maker/Apollo ids -> OGSv3 ids
	 Functional desc: OGSv3 id -> Func desc string
	 GFF file with formatted functional description and OGS ids


	Usage:
	  add_OGSids_maker-apollo_GFF.pl -g [old Merged Maker and Apollo GFF file] -a [AHRD file for mRNA with Maker and Apollo ids] -p [species acronym or prefix] -c [Prefix for chromosome in GFF]  -s [starting value for naming] -o [output formatted gff file with OGS ids]

	  
	Flags:

	  -g  Merged Maker and Apollo GFF file (required)
	  -a  AHRD tab separated file for mRNA with Maker and Apollo ids, e.g. maker-DC3.0sc00-snap-gene-91.44-mRNA-1       Glycerol-3-phosphate dehydrogenase [NAD(P)+] (AHRD V3.11 *** tr|A0A0A9WH09|A0A0A9WH09_LYGHE) (required)
	  -f  Curated descriptions from Apollo, e.g. 18712462-6e8c-478a-83c8-a40c9d0be977  Dihydrolipoyl transacetylase (required)
	  -p  Prefix for name, e.g DcitrP (required)
	  -c  Prefix for chromosome in GFF e.g. Dc3.0sc (required)
	  -s  Starting seed for gene number for each chromosome, e.g. 1000 (required)
	  -o  output GFF file with OGS ids
	  -h  Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

