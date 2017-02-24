#!/usr/bin/perl

=head1 NAME

update_maker_names_gff.pl

=head1 SYNOPSIS

update_maker_names_gff.pl -i [Apollo GFF file]

=head1 COMMAND-LINE OPTIONS

 -i  Apollo GFF dump with names for mRNA record specified
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
"\nApollo GFF3 are required.
See help below\n\n\n";
	help();
}

#get input files
my $old_gff_input_file = $opt_i;
my $input_old_gff      = read_file($old_gff_input_file)
  or die "Could not open Apollo GFF input file: $old_gff_input_file\n";

my $new_id_output_file;
$new_id_output_file   = $old_gff_input_file.'_new-ids.names';

my @lines        = split( /\n/, $input_old_gff );
my $line_count   = scalar(@lines);
my $line_counter = 0;
my $gene_flag    = 0;
my @gene_gff_line_arr;
my $current_mRNA_Solycid;
my $prev_mRNA_Solycid;
my $new_id_output = '';
my $outofrange_gene_counter = 0;
my %mRNA_Solycid_hash;

foreach my $line (@lines) {
	$line_counter++;
	chomp($line);

	if ( $line =~ m/^#/ ) {
		next;
	}

	print STDERR "\rParsing GFF3 line ". $line_counter . " of ". $line_count;

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
						chomp $alias; $current_mRNA_Solycid = $alias;
						last;
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
				#add to hash
				$mRNA_Solycid_hash{$current_mRNA_Solycid} = '';
			}
		}
	}

	## if first gene
	if (( $line =~ /\tgene\t/ ) && ( $gene_flag == 0) ){
		$gene_flag = 1;
	}
	## if next gene
	elsif (( $line =~ /\tgene\t/ ) && ( $gene_flag == 1) ){
		#IF NO SOLYC ID, GENERATE A NEW UNIQUE ID BASED UPON PREVIOUS ID
		if ( ! defined $current_mRNA_Solycid ){
			#$current_mRNA_Solycid = 'TODO'; Solyc02g094750.1.1
			my $old_count;

			#if gene has Solyc id
			if ( (defined $prev_mRNA_Solycid) && ($prev_mRNA_Solycid =~ /^Solyc/) ){
				my @prev_mRNA_Solycid_arr = split (/\./, $prev_mRNA_Solycid);
				#print STDERR "\t\t".$prev_mRNA_Solycid_arr[0]."\n\n";
				#$prev_mRNA_Solycid_arr[0] =~ m/\d\d$/ or die "Invalid data in $prev_mRNA_Solycid_arr[0]\n";
				$old_count = substr ($prev_mRNA_Solycid_arr[0], -1, 1);
				#print STDERR "\t\t$old_count\n\n";
				$prev_mRNA_Solycid =~ s/\d\.\d\.\d$//;
			}
			else{#if first gene does not have Solyc id
				$prev_mRNA_Solycid = 'ID_OUT_OF_RANGE_';
			}

			# no multiples of 10
			#print STDERR "old count $old_count prev_mRNA_Solycid $prev_mRNA_Solycid))\n";
			if ( ( $prev_mRNA_Solycid !~ /ID_OUT_OF_RANGE/) && (($old_count + 1) % 10 != 0) ){ #if prev gene was out of range
				my $new_count = $old_count + 1;
				$current_mRNA_Solycid = $prev_mRNA_Solycid.$new_count.'.1.1' ;
			}
			else{
				$outofrange_gene_counter++;
				$current_mRNA_Solycid = 'ID_OUT_OF_RANGE_'.$outofrange_gene_counter; #placeholder for cases where there are > 9 new genes between 2 old Solyc ids
			}

			$new_id_output = $new_id_output.$current_mRNA_Solycid."\n";
		}

		my $exon_count = 1;
		my $cds_count = 1;
		my $three_prime_UTR_count = 0;
		my $five_prime_UTR_count = 0;

		#process prev gene
		foreach my $prev_gff_line ( @gene_gff_line_arr ){
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
			elsif( $line_arr[2] eq 'mRNA' ){ #apollo attribute
				my $current_gene_Solycid = $current_mRNA_Solycid;
				$current_gene_Solycid =~ s/\.\d$//;
				$new_attr = 'ID=mRNA:'.$current_mRNA_Solycid.';Parent=gene:'.$current_gene_Solycid.';Name='.$current_mRNA_Solycid.';from_apollo=1;';

				my @line_attr_arr = split (/\;/, $line_arr[8]);
				foreach my $attr (@line_attr_arr){
					my ($key,$value) = split (/=/, $attr);
					if ($key eq 'description'){
						$new_attr = $new_attr.'Note='.$value;
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
		undef $current_mRNA_Solycid;
		@gene_gff_line_arr = (); # reset
	}

	#load for next round
	push (@gene_gff_line_arr, $line);
}

## last gene
#IF NO SOLYC ID, GENERATE A NEW UNIQUE ID BASED UPON PREVIOUS ID
if ( ! defined $current_mRNA_Solycid ){
	my @prev_mRNA_Solycid_arr = split (/\./, $prev_mRNA_Solycid);
	#print STDERR "\t\t".$prev_mRNA_Solycid_arr[0]."\n\n";
	my $old_count = substr ($prev_mRNA_Solycid_arr[0], -1, 1);
	#print STDERR "\t\t$old_count\n\n";
	$prev_mRNA_Solycid =~ s/\d\.\d\.\d$//;

	# no multiples of 10
	if ( ( $prev_mRNA_Solycid !~ /ID_OUT_OF_RANGE/) && (($old_count + 1) % 10 != 0) ){ #if prev gene was out of range
		my $new_count = $old_count + 1;
		$current_mRNA_Solycid = $prev_mRNA_Solycid.$new_count.'.1.1' ;
	}
	else{
		$outofrange_gene_counter++;
		$current_mRNA_Solycid = 'ID_OUT_OF_RANGE_'.$outofrange_gene_counter; #placeholder for cases where there are > 9 new genes between 2 old Solyc ids
	}


	$new_id_output = $new_id_output.$current_mRNA_Solycid."\n";
}

my $exon_count = 1;
my $cds_count = 1;
my $three_prime_UTR_count = 0;
my $five_prime_UTR_count = 0;

#process last gene
foreach my $prev_gff_line ( @gene_gff_line_arr ){
	my @line_arr = split ("\t", $prev_gff_line);
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
      update_apollo_names_description_gff.pl

    Flags:

     -i  Apollo GFF file (required)
     -h  Help


EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

Prashant Hosmani <psh65@cornell.edu> and Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

__END__


SL3.0ch06	.	mRNA	38915340	38924935	.	-	.	owner=psh65@cornell.edu;Parent=2cdea178-4a56-44d5-acbb-08082b035853;description=phosphatidylserine decarboxylase;ID=19665396-1a2e-48b7-8ca6-7f3a5a09825f;date_last_modified=2017-02-17;Name=Solyc06g060780.4.1;date_creation=2017-02-17
