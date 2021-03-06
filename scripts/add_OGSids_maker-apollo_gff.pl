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
if ( !$opt_g || !$opt_a || !$opt_f || !$opt_p || !$opt_c || !$opt_s || !$opt_o) {
  print
"\nOld GFF file, AHRD file, Apollo descriptions, name prefix, chr prefix, starting seed, output GFF file is required.
See help below\n\n\n";
  help();  
}

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
my $seed = $opt_s;
if ($seed !~ /^[0-9]+$/){ die "$seed should be a number\n"; }
#my $seed = 990; #hack to start from 990 to match names in /export/species2/Diaphorina_citri/annotation/OGSv3.0/ahrd_final_OGS3_IDS.txt

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
	$ahrd_function{$line_arr[0]}=$line_arr[1];                      # maker id = AHRD function string
	if ( defined $line_arr[2] ){
		$ahrd_domain{$line_arr[0]}=$line_arr[2];					# maker id = AHRD comma separated domain string
	}

}

# output variables
my ($gff_output, $index_output, $curated_desc_output);
# tracking variables
my (%scaffold_last_gene_id, %gene_old_id_mrna_last_rank, %gene_old_new_index, %mrna_old_new_index, $gene_counter, $mrna_counter);
my (%mrna_old_id_exon_last_rank);
my (%mrna_old_id_noncanonical_threeprimesplicesite_last_rank, %mrna_old_id_noncanonical_fiveprimesplicesite_last_rank);

$gff_output   = '';                                                 # initialize
$index_output = '';
$gene_counter = 0;
$mrna_counter = 0;
$curated_desc_output  = '';


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

		die "$scaffold not defined using ##sequence-region for $line\n" if ( !exists $scaffold_last_gene_id{$scaffold} );

		my $gene_id;
		if ( $scaffold_last_gene_id{$scaffold} == $seed ) {			# first gene on scaffold
			$gene_id = $seed;
		}
		else{
			$gene_id = $scaffold_last_gene_id{$scaffold};
		}
		$scaffold_last_gene_id{$scaffold} += 10;					# increment for 2nd gene, namespace for 9 new genes now

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
		my $gene_attributes_hashref;
		if ( !defined $gff_features->{'source'}){					#Apollo as source as no source in apollo exported GFF
			$gff_features->{'source'} = 'Apollo';
			$gene_attributes_hashref = gff3_parse_attributes ("ID=$gene_new_id;Name=$gene_new_id;method=ManualCuration");
		}
		else{
			$gene_attributes_hashref = gff3_parse_attributes ("ID=$gene_new_id;Name=$gene_new_id");
		}
		
		$gff_features->{'attributes'} = $gene_attributes_hashref;
		$gff_output = $gff_output . gff3_format_feature ($gff_features);

		$gene_counter++;

		# if ( $gene_counter % 10 == 0 ){
		# 	print STDERR "\rParsing gene number $gene_counter";
		# }
		print STDERR "\rParsing gene number $gene_counter";
	}
	elsif ( $gff_features->{'type'} eq 'mRNA' ){					# if mRNA, increment the mRNA isoform count if this is the not the first mRNA for this gene
		
		die "Multiple parents for mRNA in $line\n" if ( scalar @{$gff_features->{'attributes'}->{'Parent'}} > 1 );

		#create mrna id
		my $mrna_parent_new_id = $gene_old_new_index{$gff_features->{'attributes'}->{'Parent'}->[0]};
		my $mrna_rank;
		if ( exists $gene_old_id_mrna_last_rank{$gff_features->{'attributes'}->{'Parent'}->[0]} ){
			$mrna_rank = $gene_old_id_mrna_last_rank{$gff_features->{'attributes'}->{'Parent'}->[0]};
		}
		else{
			$gene_old_id_mrna_last_rank{$gff_features->{'attributes'}->{'Parent'}->[0]} = $mrna_rank = 1;
		}
		$gene_old_id_mrna_last_rank{$gff_features->{'attributes'}->{'Parent'}->[0]}++;				# increment for next isoform
		my $mrna_new_id = $mrna_parent_new_id . '.' . $mrna_rank;

		$index_output = $index_output . $gff_features->{'attributes'}->{'ID'}->[0] . "\t" . $mrna_new_id . "\n";


		# create the func description string
		# Only Apollo or AHRD desc, adding domains after | separator without IPR,GO,Pfam prefixes
		# can have special characters, not adding -RA -RB automatically for curated genes, that should be done by curators
		my $mrna_desc;
		if ( exists $apollo_curated_function{$gff_features->{'attributes'}->{'ID'}->[0]} ) {		# get the Apollo description if it exists
			$mrna_desc = $apollo_curated_function{$gff_features->{'attributes'}->{'ID'}->[0]};
		}
		else{
			# check if mNRA id exists in the AHRD file
			die "No desc for ".$gff_features->{'attributes'}->{'ID'}->[0]." in AHRD file" if ( !exists $ahrd_function{$gff_features->{'attributes'}->{'ID'}->[0]}); 
			$mrna_desc = $ahrd_function{$gff_features->{'attributes'}->{'ID'}->[0]};				# get AHRD description
			#die "No desc for $mrna_new_id in Mirella's OGSv3 AHRD file" if ( !exists $ahrd_function{$mrna_new_id}); # temp hack to get AHRD description using new OGSv3 id
			#$mrna_desc = $ahrd_function{$mrna_new_id};
		}

		my $mrna_domain;
		if ( exists $ahrd_domain{$gff_features->{'attributes'}->{'ID'}->[0]} ){						# add domains
			$mrna_domain = $ahrd_domain{$gff_features->{'attributes'}->{'ID'}->[0]};
			$mrna_domain =~ s/^,//;																	# remove extra , e.g. ,PF12698
			$mrna_domain =~ s/\|/,/g;																# replace | between GO terms with ,
			chomp $mrna_domain;
		}
		if ( defined $mrna_domain ){
			$mrna_desc = $mrna_desc. ' | ' .$mrna_domain;
		}

		my $mrna_old_id = $gff_features->{'attributes'}->{'ID'}->[0];								# add mrna new id to mrna index
		$mrna_old_new_index{$mrna_old_id} = $mrna_new_id;

		if ( (defined $gff_features->{'source'}) && ($gff_features->{'source'} eq 'maker') ){		# maker mRNA
			my $mrna_aed  = $gff_features->{'attributes'}->{'_AED'}->[0];
			my $mrna_eaed = $gff_features->{'attributes'}->{'_eAED'}->[0];
			my $mrna_qi   = $gff_features->{'attributes'}->{'_QI'}->[0];

			# write the new mRNA record
			my $mrna_attributes_hashref = gff3_parse_attributes ("ID=$mrna_new_id;Name=$mrna_new_id;Note=$mrna_desc;Parent=$mrna_parent_new_id;_AED=$mrna_aed;_eAED=$mrna_eaed;;_QI=$mrna_qi");
			$gff_features->{'attributes'} = $mrna_attributes_hashref;
		}
		else{																						# apollo mRNA as no source in apollo exported GFF
			# write curated_desc_output
			$curated_desc_output = $curated_desc_output . $mrna_new_id . "\t" . $mrna_desc . "\t" . 
			$gff_features->{'attributes'}->{'owner'}->[0] . "\n";

			# write the new mRNA record
			$gff_features->{'source'} = 'Apollo';
			my $mrna_attributes_hashref = gff3_parse_attributes ("ID=$mrna_new_id;Name=$mrna_new_id;Note=$mrna_desc;Parent=$mrna_parent_new_id;method=ManualCuration");
			$gff_features->{'attributes'} = $mrna_attributes_hashref;
		}

		$gff_output = $gff_output . gff3_format_feature ($gff_features);
		$mrna_counter++;
	}
	# if any other child record (exon, UTRs, CDS)
	# use the mRNA parent to create the gff record, can have multiple mRNA parents
	elsif ( ($gff_features->{'type'} eq 'CDS') || ($gff_features->{'type'} eq 'exon')
			|| ($gff_features->{'type'} eq 'three_prime_UTR') || ($gff_features->{'type'} eq 'five_prime_UTR')
			|| ($gff_features->{'type'} eq 'non_canonical_three_prime_splice_site')
			|| ($gff_features->{'type'} eq 'non_canonical_five_prime_splice_site') ){

		my $child_parent_new_id = '';
		my $child_single_mrna_parent_old_id;												

		if ( $gff_features->{'type'} ne 'exon'){											# get only parent
			$child_parent_new_id = $mrna_old_new_index{$gff_features->
				{'attributes'}->{'Parent'}->[0]};											# single parent in arrayref for everything except exons
			$child_single_mrna_parent_old_id  = $gff_features->{'attributes'}->{'Parent'}->[0] ;
		}

		if ( (scalar @{$gff_features->{'attributes'}->{'Parent'}} > 1)
			&& ($gff_features->{'type'} ne 'exon') ){										# multiple parents in arrayref of non exon feature!!!
			die "Multiple parents found for non-exon child in $line\n";
		}

		# create child new id and record
		my $child_rank;
		my $child_attributes_hashref;
		switch ( $gff_features->{'type'} ){
			case 'CDS'{
				# if ( exists $mrna_old_id_cds_last_rank {$child_single_mrna_parent_old_id} ){
				# 	$child_rank = $mrna_old_id_cds_last_rank{$child_single_mrna_parent_old_id}++;
				# }
				# else{
				# 	$mrna_old_id_cds_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				# }

				# create id
				# my $cds_new_id = 'CDS:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;
				my $cds_new_id = 'CDS:' . $mrna_old_new_index{$child_single_mrna_parent_old_id};	# no need for rank as multi line
																									# CDS feature with same name for each mRNA

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$cds_new_id;Name=$cds_new_id;Parent=$child_parent_new_id");
			}
			case 'exon'{
				if ( scalar @{$gff_features->{'attributes'}->{'Parent'}} > 1 ){						# multiple parents in arrayref !!!
					
					# print STDERR "Multiple parents found for exon child in $line\n";

					my $child_parents_arref = $gff_features->{'attributes'}->{'Parent'};

					foreach my $child_parent (@{$child_parents_arref}){
						if ( length $child_parent_new_id > 0 ){										# if not first parent
							$child_parent_new_id = $child_parent_new_id. ',' .$mrna_old_new_index{$child_parent};
						}
						else{																		# if first parent
							$child_parent_new_id = $mrna_old_new_index{$child_parent};				# using the first parent even if its the second parent that is the 
							$child_single_mrna_parent_old_id  = $child_parent;						# source of the name since exon will be shared whatever the name
						}
					}
				}
				else{																				# single parent in arrayref
					$child_parent_new_id = $mrna_old_new_index{$gff_features->{'attributes'}->{'Parent'}->[0]};
					$child_single_mrna_parent_old_id  = $gff_features->{'attributes'}->{'Parent'}->[0] ;
				}

				if ( exists $mrna_old_id_exon_last_rank {$child_single_mrna_parent_old_id} ){
					$child_rank = ++$mrna_old_id_exon_last_rank{$child_single_mrna_parent_old_id};	# pre assignment increment
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
				# if ( exists $mrna_old_id_threeprimeutr_last_rank {$child_single_mrna_parent_old_id} ){
				# 	$child_rank = $mrna_old_id_threeprimeutr_last_rank{$child_single_mrna_parent_old_id}++;
				# }
				# else{
				# 	$mrna_old_id_threeprimeutr_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				# }

				# create id
				# my $threeprimeutr_new_id = 'three_prime_UTR:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;
				my $threeprimeutr_new_id = 'three_prime_UTR:' . $mrna_old_new_index{$child_single_mrna_parent_old_id};	# no need for rank as multi line
																														# UTR feature with same name for each mRNA

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$threeprimeutr_new_id;Name=$threeprimeutr_new_id;Parent=$child_parent_new_id");
			}
			case 'non_canonical_three_prime_splice_site'{
				if ( exists $mrna_old_id_noncanonical_threeprimesplicesite_last_rank {$child_single_mrna_parent_old_id} ){
					$child_rank = ++$mrna_old_id_noncanonical_threeprimesplicesite_last_rank{$child_single_mrna_parent_old_id};	# pre assignment increment
				}
				else{
					$mrna_old_id_noncanonical_threeprimesplicesite_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				}

				# create id
				my $noncanonical_threeprimesplicesite_new_id = 'non_canonical_three_prime_splice_site:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$noncanonical_threeprimesplicesite_new_id;Name=$noncanonical_threeprimesplicesite_new_id;Parent=$child_parent_new_id");
			}
			case 'five_prime_UTR'{
				# if ( exists $mrna_old_id_fiveprimeutr_last_rank {$child_single_mrna_parent_old_id} ){
				# 	$child_rank = $mrna_old_id_fiveprimeutr_last_rank{$child_single_mrna_parent_old_id}++;
				# }
				# else{
				# 	$mrna_old_id_fiveprimeutr_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				# }

				# create id
				# my $fiveprimeutr_new_id = 'five_prime_UTR:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;
				my $fiveprimeutr_new_id = 'five_prime_UTR:' . $mrna_old_new_index{$child_single_mrna_parent_old_id};	# no need for rank as multi line
																														# UTR feature with same name for each mRNA

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$fiveprimeutr_new_id;Name=$fiveprimeutr_new_id;Parent=$child_parent_new_id");
			}
			case 'non_canonical_five_prime_splice_site'{
				if ( exists $mrna_old_id_noncanonical_fiveprimesplicesite_last_rank {$child_single_mrna_parent_old_id} ){
					$child_rank = ++$mrna_old_id_noncanonical_fiveprimesplicesite_last_rank{$child_single_mrna_parent_old_id};	# pre assignment increment
				}
				else{
					$mrna_old_id_noncanonical_fiveprimesplicesite_last_rank{$child_single_mrna_parent_old_id} = $child_rank = 1;
				}

				# create id
				my $noncanonical_fiveprimesplicesite_new_id = 'non_canonical_five_prime_splice_site:' . $mrna_old_new_index{$child_single_mrna_parent_old_id} . '.' . $child_rank;

				# create attribute hashref
				$child_attributes_hashref = gff3_parse_attributes ("ID=$noncanonical_fiveprimesplicesite_new_id;Name=$noncanonical_fiveprimesplicesite_new_id;Parent=$child_parent_new_id");
			}
		}

		# write the new child record
		$gff_features->{'attributes'} = $child_attributes_hashref;
		$gff_output = $gff_output . gff3_format_feature ($gff_features);
	}
	else{
		print STDERR "\nUnhandled record in GFF\n$line\n";										# for genome edit records: insertion_artifact, deletion_artifact, stop_codon_read_through, substitution_artifact
	}
}

print STDERR "\nDone!!\n";

# write output files
chomp $opt_o;
open my $OGFF, '>', "$opt_o" or die "Cannot open $opt_o\n";
open my $OINDEX, '>', "index_mRNA.$opt_o" or die  "Cannot open index_mRNA.$opt_o\n";
open my $OCURDESC, '>', "curated_ID_desc_owner.$opt_o" or die  "Cannot open index_mRNA.$opt_o\n";

print $OGFF $gff_output;
close $OGFF;
print $OINDEX $index_output;
close $OINDEX;
print $OCURDESC $curated_desc_output;
close $OCURDESC;


#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

	Description:

	 Renames maker and Apollo assigned gene names to gene name with version numbers, e.g.  maker-ScVcwli_1-pred_gff_maker-gene-0.0-mRNA-1 becomes DcitrP00001.1.1. Creates an index file with old and new mRNA ids. This is hard coded for <99,999 mRNAs per scaffold. Counter skips over 10 gene models so manually curated genes can be added later.
	 
	Output:
	 Index: maker/Apollo ids -> OGS ids
	 Curated genes: OGS id, description, owner
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

__DATA__


Maker and Apollo file
----------------------
##sequence-region   DC3.0sc00 1 34736251
DC3.0sc00       maker   gene    91478   95809   .       -       .       ID=maker-DC3.0sc00-pred_gff_Mikado_loci-gene-1.7;Name=maker-DC3.0sc00-pre
DC3.0sc00       maker   mRNA    91478   95809   .       -       .       ID=maker-DC3.0sc00-pred_gff_Mikado_loci-gene-1.7-mRNA-1;Parent=maker-DC3.
DC3.0sc00       maker   three_prime_UTR 91478   92623   .       -       .       ID=maker-DC3.0sc00-pred_gff_Mikado_loci-gene-1.7-mRNA-1:three_pri
DC3.0sc00       maker   exon    91478   92627   .       -       .       ID=maker-DC3.0sc00-pred_gff_Mikado_loci-gene-1.7-mRNA-1:1;Parent=maker-DC
DC3.0sc00       maker   CDS     92624   92627   .       -       1       ID=maker-DC3.0sc00-pred_gff_Mikado_loci-gene-1.7-mRNA-1:cds;Parent=maker-
DC3.0sc00       maker   CDS     95379   95767   .       -       0       ID=maker-DC3.0sc00-pred_gff_Mikado_loci-gene-1.7-mRNA-1:cds;Parent=maker-
DC3.0sc00       maker   exon    95379   95809   .       -       .       ID=maker-DC3.0sc00-pred_gff_Mikado_loci-gene-1.7-mRNA-1:2;Parent=maker-DC
DC3.0sc00       maker   five_prime_UTR  95768   95809   .       -       .       ID=maker-DC3.0sc00-pred_gff_Mikado_loci-gene-1.7-mRNA-1:five_prim
###
DC3.0sc00       .       gene    3440778 3454726 .       -       .       ID=804d5524-b17a-4c54-afff-fa7625d97f13;owner=kkercher2017@my.fit.edu;dat
DC3.0sc00       .       mRNA    3440778 3454726 .       -       .       ID=8988b2bf-3893-4932-a4cc-31606a2a300c;Parent=804d5524-b17a-4c54-afff-fa
DC3.0sc00       .       CDS     3440778 3440885 .       -       0       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3440778 3440885 .       -       .       ID=3c2dfd5f-7bd0-46fa-aa7e-a81a28cf740a;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3443362 3443580 .       -       .       ID=866eb680-539f-475e-8087-ce7ebb736293;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3443362 3443580 .       -       0       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3444910 3444988 .       -       1       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3444910 3444988 .       -       .       ID=5052b312-14c3-4e07-a7d5-d6db579076df;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3445185 3445354 .       -       .       ID=f4077db7-7e8f-47cb-8a0d-e168718de592;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3445185 3445354 .       -       0       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3446911 3447073 .       -       1       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3446911 3447073 .       -       .       ID=65505803-156d-4716-a79e-9e78caf8f524;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3447452 3447653 .       -       2       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3447452 3447653 .       -       .       ID=1dc6c007-6550-4834-8f5a-c710b0e6ca48;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3447919 3447995 .       -       .       ID=5894da40-ed11-4c1f-a7f1-522bf837a386;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3447919 3447995 .       -       1       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3449837 3449953 .       -       1       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3449837 3449953 .       -       .       ID=2cf185bf-1496-44bc-8349-c1450ef31fb9;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3450976 3451222 .       -       2       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3450976 3451222 .       -       .       ID=0a26640c-d064-4bcf-8153-5c27b5f99412;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       CDS     3451597 3451879 .       -       0       ID=8048c322-47c4-49ed-a7a5-7d75b9040dcc;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3451597 3451900 .       -       .       ID=33c9c3a2-c28b-430d-9829-73e220d04b30;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3453648 3454101 .       -       .       ID=229a7df0-8b92-44c6-bee5-45788ff12974;Parent=8988b2bf-3893-4932-a4cc-31
DC3.0sc00       .       exon    3454130 3454726 .       -       .       ID=ade1d2e1-199a-4ad2-90dd-5bb2d137c20d;Parent=8988b2bf-3893-4932-a4cc-31
###


AHRD description
----------------
augustus-DC3.0sc00-processed-gene-318.13-mRNA-1	Heat shock protein Hsp20 (AHRD V3.11 *** tr|A0A2N3UBE4|A0A2N3UBE4_9BACT)	IPR002068,PF00011
augustus-DC3.0sc00-processed-gene-31.91-mRNA-1	Protein tweety homolog (AHRD V3.11 *** tr|A0A0K8THM1|A0A0K8THM1_LYGHE)	IPR006990,PF04906
augustus-DC3.0sc00-processed-gene-31.93-mRNA-1	Protein tweety homolog (AHRD V3.11 *** tr|A0A0K8THM1|A0A0K8THM1_LYGHE)	IPR006990,PF04906
augustus-DC3.0sc00-processed-gene-325.9-mRNA-1	DNA replication complex GINS protein PSF1 (AHRD V3.11 *** tr|A0A146M5X8|A0A146M5X8_LYGHE)	I
PR021151,PF05916


Apollo curated description
--------------------------
312f7c3a-8a38-4bc9-8e12-9ef46f2203ad    S-phase kinase associated protein 1
22300a7a-3dfa-4b4e-8f84-9e9e899ad2c8    Alcohol Dehydrogenase Class 3
fae5f81b-9d66-46dd-accb-a20f10342ae5    Interference hedgehog
43292adb-c368-4b74-94af-77a897bbeeb2    Timeless 2 part A
