#!/usr/bin/perl

=head1 NAME

update_maker_ITAG_version_annotation_gff.pl

=head1 SYNOPSIS

update_maker_ITAG_version_annotation_gff.pl -i [old GFF file] -a [version and annotation]

=head1 COMMAND-LINE OPTIONS

 -i  Maker GFF file for 1 chr (required)
 -a  Tab separated file with version and annotation (required)
 -h  Help

=cut

use strict;
use warnings;
use File::Slurp;
use Getopt::Std;

our ( $opt_i, $opt_a, $opt_h );
getopts('i:a:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_i || !$opt_a) {
	print
"\nOld GFF3 and version_annotation files are required.
See help below\n\n\n";
	help();
}

#get input files
my $old_gff_input_file = $opt_i;
my $input_old_gff      = read_file($old_gff_input_file)
  or die "Could not open old GFF input file: $old_gff_input_file\n";

my $new_id_output_file;
$new_id_output_file   = $old_gff_input_file.'_new-ids.names';

my $version_annotation_input_file = $opt_a;
my $input_version_annotation      = read_file($version_annotation_input_file)
  or die "Could not open version and annotation input file: $version_annotation_input_file\n";

my @lines        = split( /\n/, $input_old_gff );
my $line_count   = scalar(@lines);
my $line_counter = 0;
my $current_Solycid;

my @version_annotations        = split( /\n/, $input_version_annotation );
my $annotation_line_count   = scalar(@version_annotations);
my %solycid_new_version_hash;
my %solycid_annotation_hash;

# CREATE HASH OF NEW VERSION AND ANNOTATION
foreach my $line (@version_annotations) {
	$line_counter++;
	print STDERR "\rParsing Annotation line ". $line_counter . " of ". $annotation_line_count . ' to populate hashes';

	chomp($line);

	my @line_arr = split ("\t", $line);
	$solycid_new_version_hash{$line_arr[0]} = $line_arr[1];
	chomp $line_arr[2];
	$solycid_annotation_hash{$line_arr[0]} = $line_arr[2];
}

$line_counter = 0;
# WRITING OUT MODIFIED GFF3
foreach my $line (@lines) {
	$line_counter++;
	chomp($line);

	if ( $line =~ m/^#/ ) {
		next;
	}

	print STDERR "\rParsing GFF3 line ". $line_counter . " of ". $line_count . ' to write modified GFF3';

	my @line_arr = split ("\t", $line);
	my @line_attr_arr = split (/\;/, $line_arr[8]);
	my $new_attr;

	if ( $line_arr[2] eq 'gene' ){
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if (($key eq 'Name') && ($value =~ /^Solyc/)){
				chomp $value; $current_Solycid = $value; last;
			}
		}

		if ( $current_Solycid ne $solycid_new_version_hash{$current_Solycid} ){#if there is an updated version
			my $current_Solycid_new = $solycid_new_version_hash{$current_Solycid};
			my $length = $line_arr[4] - $line_arr[3];
			my $current_alias_Solycid_new = $current_Solycid_new;
			$current_alias_Solycid_new =~ s/\.\d$//;
			$new_attr = 'Alias='.$current_alias_Solycid_new.';ID=gene:'.$current_Solycid_new.';Name='.$current_Solycid_new.';length='.$length."\n";

			for (0..7){
				print STDOUT $line_arr[$_]."\t";
			}
			print STDOUT $new_attr;
		}
		else{
			print STDOUT $line."\n"; #no changes
		}
	}
	elsif ( $line_arr[2] eq 'mRNA' ){
		my $AED;
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if (($key eq 'Name') && ($value =~ /^Solyc/)){
				chomp $value; $current_Solycid = $value;
			}
			elsif ( $key eq '_AED' ){
				chomp $value; $AED = $value;
			}
		}

		$current_Solycid =~ s/\.\d$//;

		if ( $current_Solycid ne $solycid_new_version_hash{$current_Solycid} ){#if there is an updated version
			die "$current_Solycid not present in version hash" if ! exists $solycid_new_version_hash{$current_Solycid};
			my $current_Solycid_new = $solycid_new_version_hash{$current_Solycid};

			my $current_Solycid_gene = $current_Solycid_new;
			$current_Solycid_new = $current_Solycid_new.'.1';
			$new_attr = 'ID=mRNA:'.$current_Solycid_new.';Name='.$current_Solycid_new.';Parent=gene:'.$current_Solycid_gene.';_AED='.$AED;

			for (0..7){
				print STDOUT $line_arr[$_]."\t";
			}
			print STDOUT $new_attr;
		}
		else{
			chomp $line;
			print STDOUT $line;
		}

		die "$current_Solycid not present in annotation hash" if ! exists $solycid_annotation_hash{$current_Solycid};
		my $annotation = $solycid_annotation_hash{$current_Solycid};
		print STDOUT ';'.$annotation."\n"; #only adding annotation


	}
	elsif ( $line_arr[2] eq 'exon' ){
		my $exon_count;
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if ( $key eq 'ID' ){
				chomp $value;
				$value =~ s/^exon://;
				#$exon_count = $value =~ m/[\d]+$/; #get count for any exon count, not working in some cases
				my @val_arr = split (/\./, $value); $exon_count = $val_arr[3];
				$value =~ s/\.\d\.[\d]+$//; # remove mRNA and exon count
				$current_Solycid = $value;
			}
		}

		if ( $current_Solycid ne $solycid_new_version_hash{$current_Solycid} ){#if there is an updated version
			die "$current_Solycid not present in version hash" if ! exists $solycid_new_version_hash{$current_Solycid};
			my $current_Solycid_new = $solycid_new_version_hash{$current_Solycid};
			$new_attr = 'ID=exon:'.$current_Solycid_new.'.1.'.$exon_count.';Parent=mRNA:'.$current_Solycid_new.'.1'."\n";

			for (0..7){
				print STDOUT $line_arr[$_]."\t";
			}
			print STDOUT $new_attr;
		}
		else{
			print STDOUT $line."\n";
		}
	}
	elsif ( $line_arr[2] eq 'CDS' ){
		my $CDS_count;
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if ( $key eq 'ID' ){
				chomp $value;
				$value =~ s/^CDS://;
				#$CDS_count = $value =~ m/[\d]+$/; #get count, not working in some cases
				my @val_arr = split (/\./, $value); $CDS_count = $val_arr[3];
				$value =~ s/\.\d\.[\d]+$//;
				$current_Solycid = $value;
			}
		}

		if ( $current_Solycid ne $solycid_new_version_hash{$current_Solycid} ){#if there is an updated version
			die "$current_Solycid not present in version hash" if ! exists $solycid_new_version_hash{$current_Solycid};
			my $current_Solycid_new = $solycid_new_version_hash{$current_Solycid};
			$new_attr = 'ID=CDS:'.$current_Solycid_new.'.1.'.$CDS_count.';Parent=mRNA:'.$current_Solycid_new.'.1'."\n";

			for (0..7){
				print STDOUT $line_arr[$_]."\t";
			}
			print STDOUT $new_attr;
		}
		else{
			print STDOUT $line."\n";
		}
	}
	elsif ( $line_arr[2] eq 'five_prime_UTR' ){
		my $five_prime_UTR_count;
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if ( $key eq 'ID' ){
				chomp $value;
				$value =~ s/^five_prime_UTR://;
				#$five_prime_UTR_count = $value =~ m/[\d]+$/; #get count, not working in some cases
				my @val_arr = split (/\./, $value); $five_prime_UTR_count = $val_arr[3];
				$value =~ s/\.\d\.[\d]+$//;
				$current_Solycid = $value;
			}
		}

		if ( $current_Solycid ne $solycid_new_version_hash{$current_Solycid} ){#if there is an updated version
			die "$current_Solycid not present in version hash" if ! exists $solycid_new_version_hash{$current_Solycid};
			my $current_Solycid_new = $solycid_new_version_hash{$current_Solycid};
			$new_attr = 'ID=five_prime_UTR:'.$current_Solycid_new.'.1.'.$five_prime_UTR_count.';Parent=mRNA:'.$current_Solycid_new.'.1'."\n";

			for (0..7){
				print STDOUT $line_arr[$_]."\t";
			}
			print STDOUT $new_attr;
		}
		else{
			print STDOUT $line."\n";
		}
	}
	elsif ( $line_arr[2] eq 'three_prime_UTR' ){
		my $three_prime_UTR_count;
		foreach my $attr (@line_attr_arr){
			my ($key,$value) = split (/=/, $attr);
			if ( $key eq 'ID' ){
				chomp $value;
				$value =~ s/^three_prime_UTR://;
				#$three_prime_UTR_count = $value =~ m/[\d]+$/; #get count, not working in some cases
				my @val_arr = split (/\./, $value); $three_prime_UTR_count = $val_arr[3];
				$value =~ s/\.\d\.[\d]+$//;
				$current_Solycid = $value;
			}
		}

		if ( $current_Solycid ne $solycid_new_version_hash{$current_Solycid} ){#if there is an updated version
			die "$current_Solycid not present in version hash" if ! exists $solycid_new_version_hash{$current_Solycid};
			my $current_Solycid_new = $solycid_new_version_hash{$current_Solycid};
			$new_attr = 'ID=three_prime_UTR:'.$current_Solycid_new.'.1.'.$three_prime_UTR_count.';Parent=mRNA:'.$current_Solycid_new.'.1'."\n";

			for (0..7){
				print STDOUT $line_arr[$_]."\t";
			}
			print STDOUT $new_attr;
		}
		else{
			print STDOUT $line."\n";
		}
	}
}

print STDERR "\n";



#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     For a Maker_ITAG GFF file produced by update_maker_names_gff.pl, this script updates the version number of the Solyc id if the gene model has been updated and adds in anotation fields to the attribute column Output GFF contains feature names in ITAG convention


    Usage:
      update_maker_ITAG_version_annotation_gff.pl

    Flags:

     -i  old GFF file (required)
     -a  Tab separated file with version and annotation (required)
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


# Maker ITAG GFF3
SL3.0ch07	maker_ITAG	gene	16549	35133	.	+	.	Alias=Solyc07g005000;ID=gene:Solyc07g005000.2;Name=Solyc07g005000.2;length=18584
SL3.0ch07	maker_ITAG	mRNA	16549	35133	.	+	.	ID=mRNA:Solyc07g005000.2.1;Name=Solyc07g005000.2.1;_AED=0.10
SL3.0ch07	maker_ITAG	five_prime_UTR	16549	16768	.	+	.	ID=five_prime_UTR:Solyc07g005000.2.1.0;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	16549	17005	.	+	.	ID=exon:Solyc07g005000.2.1.1;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	16769	17005	.	+	0	ID=CDS:Solyc07g005000.2.1.1;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	19253	19327	.	+	.	ID=exon:Solyc07g005000.2.1.2;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	19253	19327	.	+	0	ID=CDS:Solyc07g005000.2.1.2;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	19480	19533	.	+	.	ID=exon:Solyc07g005000.2.1.3;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	19480	19533	.	+	0	ID=CDS:Solyc07g005000.2.1.3;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	20043	20108	.	+	.	ID=exon:Solyc07g005000.2.1.4;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	20043	20108	.	+	0	ID=CDS:Solyc07g005000.2.1.4;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	21090	21148	.	+	.	ID=exon:Solyc07g005000.2.1.5;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	21090	21148	.	+	0	ID=CDS:Solyc07g005000.2.1.5;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	21608	21684	.	+	.	ID=exon:Solyc07g005000.2.1.6;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	21608	21684	.	+	1	ID=CDS:Solyc07g005000.2.1.6;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	22176	22301	.	+	.	ID=exon:Solyc07g005000.2.1.7;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	22176	22301	.	+	2	ID=CDS:Solyc07g005000.2.1.7;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	23257	23291	.	+	.	ID=exon:Solyc07g005000.2.1.8;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	23257	23291	.	+	2	ID=CDS:Solyc07g005000.2.1.8;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	24248	24348	.	+	.	ID=exon:Solyc07g005000.2.1.9;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	24248	24348	.	+	0	ID=CDS:Solyc07g005000.2.1.9;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	24419	24514	.	+	.	ID=exon:Solyc07g005000.2.1.10;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	24419	24514	.	+	1	ID=CDS:Solyc07g005000.2.1.10;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	25782	25846	.	+	.	ID=exon:Solyc07g005000.2.1.11;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	25782	25846	.	+	1	ID=CDS:Solyc07g005000.2.1.11;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	26427	26498	.	+	.	ID=exon:Solyc07g005000.2.1.12;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	26427	26498	.	+	2	ID=CDS:Solyc07g005000.2.1.12;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	28705	28875	.	+	.	ID=exon:Solyc07g005000.2.1.13;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	28705	28875	.	+	2	ID=CDS:Solyc07g005000.2.1.13;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	29459	29658	.	+	.	ID=exon:Solyc07g005000.2.1.14;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	29459	29658	.	+	2	ID=CDS:Solyc07g005000.2.1.14;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	29779	29845	.	+	.	ID=exon:Solyc07g005000.2.1.15;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	29779	29845	.	+	0	ID=CDS:Solyc07g005000.2.1.15;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	30863	31038	.	+	.	ID=exon:Solyc07g005000.2.1.16;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	30863	31038	.	+	2	ID=CDS:Solyc07g005000.2.1.16;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	31356	31457	.	+	.	ID=exon:Solyc07g005000.2.1.17;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	31356	31457	.	+	0	ID=CDS:Solyc07g005000.2.1.17;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	32716	32763	.	+	.	ID=exon:Solyc07g005000.2.1.18;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	32716	32763	.	+	0	ID=CDS:Solyc07g005000.2.1.18;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	33107	33250	.	+	.	ID=exon:Solyc07g005000.2.1.19;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	33107	33250	.	+	0	ID=CDS:Solyc07g005000.2.1.19;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	33320	33415	.	+	.	ID=exon:Solyc07g005000.2.1.20;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	33320	33415	.	+	0	ID=CDS:Solyc07g005000.2.1.20;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	34031	34132	.	+	.	ID=exon:Solyc07g005000.2.1.21;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	34031	34132	.	+	0	ID=CDS:Solyc07g005000.2.1.21;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	CDS	34805	34897	.	+	0	ID=CDS:Solyc07g005000.2.1.22;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	exon	34805	35133	.	+	.	ID=exon:Solyc07g005000.2.1.22;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	three_prime_UTR	34898	35133	.	+	.	ID=three_prime_UTR:Solyc07g005000.2.1.0;Parent=mRNA:Solyc07g005000.2.1
SL3.0ch07	maker_ITAG	gene	41762	45808	.	+	.	Alias=Solyc07g005010;ID=gene:Solyc07g005010.2;Name=Solyc07g005010.2;length=4046
SL3.0ch07	maker_ITAG	mRNA	41762	45808	.	+	.	ID=mRNA:Solyc07g005010.2.1;Name=Solyc07g005010.2.1;_AED=0.03
SL3.0ch07	maker_ITAG	exon	41762	42298	.	+	.	ID=exon:Solyc07g005010.2.1.1;Parent=mRNA:Solyc07g005010.2.1
SL3.0ch07	maker_ITAG	CDS	41762	42298	.	+	0	ID=CDS:Solyc07g005010.2.1.1;Parent=mRNA:Solyc07g005010.2.1
SL3.0ch07	maker_ITAG	exon	42918	44907	.	+	.	ID=exon:Solyc07g005010.2.1.2;Parent=mRNA:Solyc07g005010.2.1
SL3.0ch07	maker_ITAG	CDS	42918	44907	.	+	0	ID=CDS:Solyc07g005010.2.1.2;Parent=mRNA:Solyc07g005010.2.1
SL3.0ch07	maker_ITAG	CDS	45174	45478	.	+	2	ID=CDS:Solyc07g005010.2.1.3;Parent=mRNA:Solyc07g005010.2.1
SL3.0ch07	maker_ITAG	exon	45174	45808	.	+	.	ID=exon:Solyc07g005010.2.1.3;Parent=mRNA:Solyc07g005010.2.1
SL3.0ch07	maker_ITAG	three_prime_UTR	45479	45808	.	+	.	ID=three_prime_UTR:Solyc07g005010.2.1.0;Parent=mRNA:Solyc07g005010.2.1
