#!/usr/bin/perl

=head1 NAME

update_maker_names_gff.pl

=head1 SYNOPSIS

update_maker_names_gff.pl -i [old GFF file]

=head1 COMMAND-LINE OPTIONS

 -i  Maker GFF file for 1 chr (required)
 -h  Help

=cut

use strict;
use warnings;
use File::Slurp;
use Getopt::Std;

our ( $opt_i, $opt_h );
getopts('i:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_i ) {
	print
"\nOld GFF3 are required.
See help below\n\n\n";
	help();
}

#get input files
my $old_gff_input_file = $opt_i;
my $input_old_gff      = read_file($old_gff_input_file)
  or die "Could not open old GFF input file: $old_gff_input_file\n";

my $new_id_output_file;
$new_id_output_file   = $old_gff_input_file.'_new-ids.names';

my @lines        = split( /\n/, $input_old_gff );
my $line_count   = scalar(@lines);
my $line_counter = 0;
my $gene_flag    = 0;
my @gene_gff_lines_arr;
my $current_mRNA_Solycid;
my $prev_mRNA_Solycid;			#can be old ITAG2.4 Solyc id or new ITAG Solyc id or ID_OUT_OF_RANGE
my $prev_mRNA_ITAG24_Solycid;		#can be only ITAG2.4 Solyc id ending wt 0
my $new_id_output = '';
my $outofrange_gene_counter = 0;
my $new_gene_counter = 0;
my %mRNA_Solycid_hash;
my %mRNA_Solycid_new_gene_count_hash;	#hash of nof new genes after $prev_mRNA_ITAG24_Solycid
my @solycid_new_gene_block_arr;		#arr of Solyc ids for new genes after $prev_mRNA_ITAG24_Solycid and before the next mRNA_ITAG24_Solycid
my @solycid_increment_arr;		#arr of increments to use to create solyc ids for new genes in a block

# PARSING FILE TO FIGURE OUT NUMBER OF NOVEL GENES BETWEEN OLD SOLYC IDS
foreach my $line (@lines) {
	$line_counter++;
	chomp($line);

	if ( $line =~ m/^#/ ) {
		next;
	}

	print STDERR "\rParsing GFF3 line ". $line_counter . " of ". $line_count . ' to count novel genes';
	
	#get Solyc id, if any, from mRNA record
	if ( $line =~ m/\tmRNA\t/ ){
		my @line_arr = split ("\t", $line);
		my @line_attr_arr = split (/\;/, $line_arr[8]);
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if (($key eq 'Name') && ($value =~ /^Solyc/)){
				chomp $value; $current_mRNA_Solycid = $value; last;
			}
			elsif (($key eq 'Alias') && ($value =~ /^Solyc/) && ($value !~ /\,/)){ #only Solyc in Alias
				chomp $value; $current_mRNA_Solycid = $value; last;
			}
			elsif (($key eq 'Alias') && ($value =~ /Solyc/)){ #multiple aliases incl Solyc id
				my @value_arr = split (/,/,$value);
				foreach my $alias (@value_arr){
					if ($alias =~ /^Solyc/) {
						chomp $alias; $current_mRNA_Solycid = $alias;
						last;
					}
				}
			}
		}
		
		#reset if Solyc id already exists as Maker assigns same id to multiple genes during pass through
		if ( defined $current_mRNA_Solycid){
			if (exists ($mRNA_Solycid_hash{$current_mRNA_Solycid})){
				undef $current_mRNA_Solycid;
			}
			else{
				$mRNA_Solycid_hash{$current_mRNA_Solycid} = '';#add to hash
			}
		}
	}
	
	## if first gene
	if (( $line =~ /\tgene\t/ ) && ( $gene_flag == 0) ){
		$gene_flag = 1;
	}
	## if next gene
	elsif (( $line =~ /\tgene\t/ ) && ( $gene_flag == 1) ){
		#IF NO SOLYC ID IN PREV GENE AND PREV SOLYC ID EXISTS, ADD TO COUNT OF NEW GENES AFTER PREV SOLYC ID
		if (( ! defined $current_mRNA_Solycid ) && (defined $prev_mRNA_ITAG24_Solycid)){
			if ( exists $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} ){
				#$mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} = $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid}++;
				$mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} = $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} + 1;
			}
			else{
				$mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} = 1;
			}
			$new_gene_counter++;
		}
		
		if ( (defined $current_mRNA_Solycid) && ($current_mRNA_Solycid =~ /^Solyc/) ){
			$prev_mRNA_Solycid = $current_mRNA_Solycid;
		}
		
		#only remember ITAG2.4 mRNA Solyc ids (ending with 0)
		if ((defined $current_mRNA_Solycid) && ($current_mRNA_Solycid =~ /^Solyc/) && ($current_mRNA_Solycid =~ /0\.\d\.\d$/)){
			$prev_mRNA_ITAG24_Solycid = $current_mRNA_Solycid;
		}

		undef $current_mRNA_Solycid;
	}	
}

#last gene
if ( ! defined $current_mRNA_Solycid ){
	if ( exists $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid}){
		#$mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} = $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid}++;
		$mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} = $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} + 1;
	}
	else{
		$mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid} = 1;
	}
	$new_gene_counter++;
}



# RESET
print STDERR "\nIdentified $new_gene_counter new genes. Now writing modified GFF3\n";
%mRNA_Solycid_hash = ();
$line_counter = 0;
$gene_flag    = 0;
@gene_gff_lines_arr = ();
undef $current_mRNA_Solycid;
undef $prev_mRNA_Solycid;
undef $prev_mRNA_ITAG24_Solycid;

#HASH OF ARRAYS FOR UP TO 9 NEW GENES, > 9 SHOULD HAVE ID_OUT_OF_RANGE_ NAMES
#hard coded for interval of 9 between 2 ITAG2.4 Solyc ids but it can be larger. Elastic ranges are better
my %new_solyc_id_increments = (
	1 => [5],
	2 => [3,7],
	3 => [3,5,7],
	4 => [2,4,6,8],
	5 => [1,2,4,6,8],
	6 => [1,2,3,4,6,8],
	7 => [1,2,3,4,5,6,8],
	8 => [1,2,3,4,5,6,7,8],
	9 => [1,2,3,4,5,6,7,8,9]
);

# WRITING OUT MODIFIED GFF3
foreach my $line (@lines) {
	$line_counter++;
	chomp($line);

	if ( $line =~ m/^#/ ) {
		next;
	}

	print STDERR "\rParsing GFF3 line ". $line_counter . " of ". $line_count . ' to write modified GFF3';

	#get Solyc id ,if any, from mRNA record
	if ( $line =~ m/\tmRNA\t/ ){
		my @line_arr = split ("\t", $line);
		my @line_attr_arr = split (/\;/, $line_arr[8]);
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if (($key eq 'Name') && ($value =~ /^Solyc/)){
				chomp $value; $current_mRNA_Solycid = $value; last;
			}
			elsif (($key eq 'Alias') && ($value =~ /^Solyc/) && ($value !~ /\,/)){ #only Solyc in Alias
				chomp $value; $current_mRNA_Solycid = $value; last;
			}
			elsif (($key eq 'Alias') && ($value =~ /Solyc/)){ #multiple aliases incl Solyc id
				my @value_arr = split (/,/,$value);
				foreach my $alias (@value_arr){
					if ($alias =~ /^Solyc/) {
						chomp $alias; $current_mRNA_Solycid = $alias; last;
					}
				}
			}
		}
		
		#reset if Solyc id already exists
		if ( defined $current_mRNA_Solycid){
			if (exists ($mRNA_Solycid_hash{$current_mRNA_Solycid})){
				undef $current_mRNA_Solycid;
			}
			else{
				$mRNA_Solycid_hash{$current_mRNA_Solycid} = ''; #add to hash
			}		
		}
	}

	## if first gene
	if (( $line =~ /\tgene\t/ ) && ( $gene_flag == 0) ){
		$gene_flag = 1;
	}
	## if next gene
	elsif (( $line =~ /\tgene\t/ ) && ( $gene_flag == 1) ){
		#IF NO SOLYC ID AND PREV GENE WITH SOLYC ID FROM ITAG2.4 EXISTS IN HASH
		#GENERATE A NEW UNIQUE ID BASED UPON PREVIOUS ID
		if ((!defined $current_mRNA_Solycid) 
			&& (defined $prev_mRNA_ITAG24_Solycid) 
			&& (exists $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid})){

			if (scalar @solycid_new_gene_block_arr == 0) {#this is a new gene interval
				my $new_gene_count = $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid};

				#Get Solycid prefix
				#Solyc02g094750.1.1				
				my @prev_mRNA_ITAG24_Solycid_arr = split (/\./, $prev_mRNA_ITAG24_Solycid);
				my $old_count = substr ($prev_mRNA_ITAG24_Solycid_arr[0], -1, 1);#will always be 0
				die "\nprev_mRNA_ITAG24_Solycid does not end with 0 in $prev_mRNA_ITAG24_Solycid_arr[0]\n" if $old_count != 0;
				(my $prev_mRNA_ITAG24_Solycid_prefix = $prev_mRNA_ITAG24_Solycid) =~ s/\d\.\d\.\d$//; #removing gene and mRNA ver #
				
				#create new solyc ids for new genes in block before next ITAG2.4 gene
				#@solycid_increment_arr = ();#reset
				if ($new_gene_count <= 9){
					#@solycid_increment_arr = $new_solyc_id_increments{$new_gene_count};
					for (0 .. ($new_gene_count - 1)){
						#my $new_count   = $old_count + $solycid_increment_arr[$_];
						#my $new_count   = $old_count + @{$new_solyc_id_increments{$new_gene_count}}[$_];
						#my $new_count   = $old_count + @{$new_solyc_id_increments{$new_gene_count}}[$_];
						my $new_suffix   = @{$new_solyc_id_increments{$new_gene_count}}[$_];
						#my $new_Solycid = $prev_mRNA_ITAG24_Solycid_prefix.$new_count.'.1.1';
						my $new_Solycid = $prev_mRNA_ITAG24_Solycid_prefix.$new_suffix.'.1.1';
						push @solycid_new_gene_block_arr,$new_Solycid;#assigns to bottom of arr
					}
					
					
				}
				elsif ($new_gene_count > 9){
					for (0..8){
						#my $new_count   = $old_count + $solycid_increment_arr[$_];
						#my $new_count   = $old_count + @{$new_solyc_id_increments{$new_gene_count}}[$_];
						#my $new_count   = $old_count + @{$new_solyc_id_increments{$new_gene_count}}[$_];
						#print STDERR "For $prev_mRNA_ITAG24_Solycid with $new_gene_count value ".@{$new_solyc_id_increments{$new_gene_count}}[$_]."\n";
						my $new_suffix   = @{$new_solyc_id_increments{9}}[$_];
						#my $new_Solycid = $prev_mRNA_ITAG24_Solycid_prefix.$new_count.'.1.1';
						my $new_Solycid = $prev_mRNA_ITAG24_Solycid_prefix.$new_suffix.'.1.1';
						push @solycid_new_gene_block_arr,$new_Solycid;#assigns to bottom of arr
					}
					for (9 .. ($new_gene_count - 1)){
						$outofrange_gene_counter++;
						my $new_Solycid = 'ID_OUT_OF_RANGE_'.$outofrange_gene_counter;
						push @solycid_new_gene_block_arr,$new_Solycid;#assigns to bottom of arr
					}
				}
			}
			
			#shift and assign id in @solycid_new_gene_block_arr to current_mRNA_Solycid
			$current_mRNA_Solycid = shift @solycid_new_gene_block_arr;# to current_mRNA_Solycid
			$new_id_output = $new_id_output.$current_mRNA_Solycid."\n";
		}
		
		if ((!defined $current_mRNA_Solycid) 
			&& (!defined $prev_mRNA_ITAG24_Solycid)){
			$outofrange_gene_counter++;
			$current_mRNA_Solycid = 'ID_OUT_OF_RANGE_'.$outofrange_gene_counter;
			
		}


		my $exon_count = 1;
		my $cds_count = 1;
		my $three_prime_UTR_count = 0;
		my $five_prime_UTR_count = 0;

		#process prev gene
		foreach my $prev_gff_line ( @gene_gff_lines_arr ){
			my @line_arr = split ("\t", $prev_gff_line);
			$line_arr[1] = 'maker_ITAG'; #using source to reflect ITAG/eugene fed into maker
			my $new_attr;

			if ( $line_arr[2] eq 'gene' ){
				my $length = $line_arr[4] - $line_arr[3];
				my $current_gene_Solycid = $current_mRNA_Solycid;
				$current_gene_Solycid =~ s/\.\d$//;
				my $current_alias_Solycid = $current_gene_Solycid;
				$current_alias_Solycid =~ s/\.\d$//;
				$new_attr = 'Alias='.$current_alias_Solycid.';ID=gene:'.$current_gene_Solycid.';Name='.$current_gene_Solycid.';length='.$length."\n";

				for (0..7){
					print STDOUT $line_arr[$_]."\t";
				}
				print STDOUT $new_attr;
			}
			elsif( $line_arr[2] eq 'mRNA' ){ #add AED
				my $current_gene_Solycid = $current_mRNA_Solycid;
				$current_gene_Solycid =~ s/\.\d$//;
				#$new_attr = 'ID=mRNA:'.$current_mRNA_Solycid.';Name='.$current_mRNA_Solycid.';';
				$new_attr = 'ID=mRNA:'.$current_mRNA_Solycid.';Parent=gene:'.$current_gene_Solycid.';Name='.$current_mRNA_Solycid.';';

				my @line_attr_arr = split (/\;/, $line_arr[8]);
				foreach my $attr (@line_attr_arr){
					my ($key,$value) = split (/=/, $attr);
					if ($key eq '_AED'){
						$new_attr = $new_attr.'_AED='.$value;
					}
				}
				for (0..7){
					print STDOUT $line_arr[$_]."\t";
				}
				print STDOUT $new_attr."\n";
			}
			elsif( $line_arr[2] eq 'exon' ){
				my $current_exon_Solycid = $current_mRNA_Solycid.'.'.$exon_count;
				$exon_count++;
				$new_attr = 'ID=exon:'.$current_exon_Solycid.';Parent=mRNA:'.$current_mRNA_Solycid."\n";
				for (0..7){
					print STDOUT $line_arr[$_]."\t";
				}
				print STDOUT $new_attr;

			}
			elsif( $line_arr[2] eq 'CDS' ){
				my $current_cds_Solycid = $current_mRNA_Solycid.'.'.$cds_count;
				$cds_count++;
				$new_attr = 'ID=CDS:'.$current_cds_Solycid.';Parent=mRNA:'.$current_mRNA_Solycid."\n";
				for (0..7){
					print STDOUT $line_arr[$_]."\t";
				}
				print STDOUT $new_attr;
			}
			elsif( $line_arr[2] eq 'five_prime_UTR' ){
				my $current_fiveprime_Solycid = $current_mRNA_Solycid.'.'.$five_prime_UTR_count;
				$five_prime_UTR_count++;
				$new_attr = 'ID=five_prime_UTR:'.$current_fiveprime_Solycid.';Parent=mRNA:'.$current_mRNA_Solycid."\n";
				for (0..7){
					print STDOUT $line_arr[$_]."\t";
				}
				print STDOUT $new_attr;
			}
			elsif( $line_arr[2] eq 'three_prime_UTR' ){
				my $current_threeprime_Solycid = $current_mRNA_Solycid.'.'.$three_prime_UTR_count;
				$three_prime_UTR_count++;
				$new_attr = 'ID=three_prime_UTR:'.$current_threeprime_Solycid.';Parent=mRNA:'.$current_mRNA_Solycid."\n";
				for (0..7){
					print STDOUT $line_arr[$_]."\t";
				}
				print STDOUT $new_attr;
			}
		}

		$prev_mRNA_Solycid = $current_mRNA_Solycid;
		#only remember ITAG2.4 mRNA Solyc ids (ending with 0)
		if (($current_mRNA_Solycid =~ /^Solyc/) && ($current_mRNA_Solycid =~ /0\.\d\.\d$/)){
			$prev_mRNA_ITAG24_Solycid = $current_mRNA_Solycid;
		}
		undef $current_mRNA_Solycid;
		@gene_gff_lines_arr = (); # reset
	}

	#load for next round
	push (@gene_gff_lines_arr, $line);
}

## last gene
#IF NO SOLYC ID, GENERATE A NEW UNIQUE ID BASED UPON PREVIOUS ID
if ( ! defined $current_mRNA_Solycid ){

	#IF NO SOLYC ID AND PREV GENE WITH SOLYC ID FROM ITAG2.4 EXISTS IN HASH
	#GENERATE A NEW UNIQUE ID BASED UPON PREVIOUS ID
	if ((!defined $current_mRNA_Solycid) 
		&& (defined $prev_mRNA_ITAG24_Solycid) 
		&& (exists $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid})){

		if (scalar @solycid_new_gene_block_arr == 0) {#this is a new gene interval
			my $new_gene_count = $mRNA_Solycid_new_gene_count_hash{$prev_mRNA_ITAG24_Solycid};

			#Get Solycid prefix
			#Solyc02g094750.1.1				
			my @prev_mRNA_ITAG24_Solycid_arr = split (/\./, $prev_mRNA_ITAG24_Solycid);
			my $old_count = substr ($prev_mRNA_ITAG24_Solycid_arr[0], -1, 1);#will always be 0
			die "\nprev_mRNA_ITAG24_Solycid does not end with 0 in $prev_mRNA_ITAG24_Solycid_arr[0]\n" if $old_count!=0;
			(my $prev_mRNA_ITAG24_Solycid_prefix = $prev_mRNA_ITAG24_Solycid) =~ s/\d\.\d\.\d$//; #removing gene and mRNA ver #
			
			#create new solyc ids for new genes in block before next ITAG2.4 gene
			#@solycid_increment_arr = ();#reset <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< FIX
			if ($new_gene_count <= 9){
				#@solycid_increment_arr = $new_solyc_id_increments{$new_gene_count};
				for (0 .. ($new_gene_count - 1)){
					#my $new_count   = $old_count + $solycid_increment_arr[$_];
					#my $new_count   = $old_count + @{$new_solyc_id_increments{$new_gene_count}}[$_];
					#my $new_count   = $old_count + @{$new_solyc_id_increments{$new_gene_count}}[$_];
					my $new_suffix   = @{$new_solyc_id_increments{$new_gene_count}}[$_];
					#my $new_Solycid = $prev_mRNA_ITAG24_Solycid_prefix.$new_count.'.1.1';
					my $new_Solycid = $prev_mRNA_ITAG24_Solycid_prefix.$new_suffix.'.1.1';
					push @solycid_new_gene_block_arr,$new_Solycid;#assigns to bottom of arr
				}
				
				
			}
			elsif ($new_gene_count > 9){
				for (0 .. 8){
					#my $new_count   = $old_count + $solycid_increment_arr[$_];
					#my $new_count   = $old_count + @{$new_solyc_id_increments{$new_gene_count}}[$_];
					#my $new_count   = $old_count + @{$new_solyc_id_increments{$new_gene_count}}[$_];
					#my $new_suffix   = @{$new_solyc_id_increments{$new_gene_count}}[$_];
					my $new_suffix   = @{$new_solyc_id_increments{9}}[$_];
					#my $new_Solycid = $prev_mRNA_ITAG24_Solycid_prefix.$new_count.'.1.1';
					my $new_Solycid = $prev_mRNA_ITAG24_Solycid_prefix.$new_suffix.'.1.1';
					push @solycid_new_gene_block_arr,$new_Solycid;#assigns to bottom of arr
				}
				for (9 .. ($new_gene_count - 1)){
					$outofrange_gene_counter++;
					my $new_Solycid = 'ID_OUT_OF_RANGE_'.$outofrange_gene_counter;
					push @solycid_new_gene_block_arr,$new_Solycid;#assigns to bottom of arr
				}
			}
		}
		
		#shift and assign id in @solycid_new_gene_block_arr to current_mRNA_Solycid
		$current_mRNA_Solycid = shift @solycid_new_gene_block_arr;# to current_mRNA_Solycid
		$new_id_output = $new_id_output.$current_mRNA_Solycid."\n";
	}

	$new_id_output = $new_id_output.$current_mRNA_Solycid."\n";
}

my $exon_count = 1;
my $cds_count = 1;
my $three_prime_UTR_count = 0;
my $five_prime_UTR_count = 0;

#process last gene
foreach my $prev_gff_line ( @gene_gff_lines_arr ){
	my @line_arr = split ("\t", $prev_gff_line);
	$line_arr[1] = 'maker_ITAG'; #using source to reflect ITAG/eugene fed into maker
	my $new_attr;

	if ( $line_arr[2] eq 'gene' ){
		my $length = $line_arr[4] - $line_arr[3];
		my $current_gene_Solycid = $current_mRNA_Solycid;
		$current_gene_Solycid =~ s/\.\d$//;
		my $current_alias_Solycid = $current_gene_Solycid;
		$current_alias_Solycid =~ s/\.\d$//;
		$new_attr = 'Alias='.$current_alias_Solycid.';ID=gene:'.$current_gene_Solycid.';Name='.$current_gene_Solycid.';length='.$length."\n";

		for (0..7){
			print STDOUT $line_arr[$_]."\t";
		}
		print STDOUT $new_attr;
	}
	elsif( $line_arr[2] eq 'mRNA' ){ #add AED
		my $current_gene_Solycid = $current_mRNA_Solycid;
		$current_gene_Solycid =~ s/\.\d$//;
		$new_attr = 'ID=mRNA:'.$current_mRNA_Solycid.';Parent=gene:'.$current_gene_Solycid.';Name='.$current_mRNA_Solycid.';';
		#$new_attr = 'ID=mRNA:'.$current_mRNA_Solycid.';Name='.$current_mRNA_Solycid.';';

		my @line_attr_arr = split (/\;/, $line_arr[8]);
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if ($key eq '_AED'){
				$new_attr = $new_attr.'AED='.$value;
			}
		}
		for (0..7){
			print STDOUT $line_arr[$_]."\t";
		}
		print STDOUT $new_attr."\n";
	}
	elsif( $line_arr[2] eq 'exon' ){
		my $current_exon_Solycid = $current_mRNA_Solycid.'.'.$exon_count;
		$exon_count++;
		$new_attr = 'ID=exon:'.$current_exon_Solycid.';Parent=mRNA:'.$current_mRNA_Solycid."\n";
		for (0..7){
			print STDOUT $line_arr[$_]."\t";
		}
		print STDOUT $new_attr;

	}
	elsif( $line_arr[2] eq 'CDS' ){
		my $current_cds_Solycid = $current_mRNA_Solycid.'.'.$cds_count;
		$cds_count++;
		$new_attr = 'ID=CDS:'.$current_cds_Solycid.';Parent=mRNA:'.$current_mRNA_Solycid."\n";
		for (0..7){
			print STDOUT $line_arr[$_]."\t";
		}
		print STDOUT $new_attr;
	}
	elsif( $line_arr[2] eq 'five_prime_UTR' ){
		my $current_fiveprime_Solycid = $current_mRNA_Solycid.'.'.$five_prime_UTR_count;
		$five_prime_UTR_count++;
		$new_attr = 'ID=five_prime_UTR:'.$current_fiveprime_Solycid.';Parent=mRNA:'.$current_mRNA_Solycid."\n";
		for (0..7){
			print STDOUT $line_arr[$_]."\t";
		}
		print STDOUT $new_attr;
	}
	elsif( $line_arr[2] eq 'three_prime_UTR' ){
		my $current_threeprime_Solycid = $current_mRNA_Solycid.'.'.$three_prime_UTR_count;
		$three_prime_UTR_count++;
		$new_attr = 'ID=three_prime_UTR:'.$current_threeprime_Solycid.';Parent=mRNA:'.$current_mRNA_Solycid."\n";
		for (0..7){
			print STDOUT $line_arr[$_]."\t";
		}
		print STDOUT $new_attr;
	}
}

unless ( open( OID, ">$new_id_output_file" ) ) {
	print STDERR "Cannot open $new_id_output_file\n";
	exit 1;
}
print OID $new_id_output;
close(OID);

print STDERR "\nNumber of genes without Solyc id assignments (ID_OUT_OF_RANGE): $outofrange_gene_counter\n";

#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     Gets Solyc id from the mRNA record and assigns to gene and mRNA records. Generates a new Solyc id (skipping over multiples of 10 to avoid old ids) if no old id was passed through and assigns to gene and mRNA records. Generates a file with list of new Solyc ids. No need to check if any new id overlaps with a deprecated gene. Output GFF contains feature names in ITAG convention


    Usage:
      update_maker_names_gff.pl

    Flags:

     -i  old GFF file (required)
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
