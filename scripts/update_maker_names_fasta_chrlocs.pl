#!/usr/bin/perl

=head1 NAME

update_maker_names_fasta_chrlocs.pl

=head1 SYNOPSIS

update_maker_names_fasta_chrlocs.pl -i [old fasta file] -p [species acronym or prefix] -c [Prefix for chromosome in GFF]  -s [starting value for naming] -g [gff file with chromosomal location]

=head1 COMMAND-LINE OPTIONS

 -i  Maker Fasta file (required)
 -p  Prefix for name, e.g DcitrP (required)
 -c  Prefix for chromosome in GFF e.g. Dc3.0sc (required)
 -s  Starting seed, e.g. 1 (required)
 -g  GFF file with chromosomal location
 -h  Help

=cut

use strict;
use warnings;
use File::Slurp;
use Getopt::Std;

our ( $opt_i, $opt_p, $opt_s, $opt_c, $opt_g, $opt_h );
getopts('i:p:s:c:g:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_i || !$opt_p || !$opt_c || !$opt_s || !$opt_g) {
	print
"\nOld Fasta file, name prefix, chr prefix, GFF file with location and starting seed number is required.
See help below\n\n\n";
	help();
}

#get input files
my $old_fasta_input_file = $opt_i;
my $input_old_fasta      = read_file($old_fasta_input_file)
	or die "Could not open old fasta input file: $old_fasta_input_file\n";
chomp $opt_p; my $prefix = $opt_p;
chomp $opt_c; my $chrprefix = $opt_c;
my $seed = $opt_s;
if ($seed !~ /^[0-9]+$/){
	die "$seed should be a number\n";
}
my $old_gff_input_file = $opt_g;
my $input_old_gff      = read_file($old_gff_input_file)
	or die "Could not open old fasta input file: $old_gff_input_file\n";

my $new_id_fasta_output_file   = 'renamed.'.${old_fasta_input_file};
my $new_id_index_output_file   = 'index.'.${old_fasta_input_file};
my @lines        = split( /\n/, $input_old_fasta );
my @gff_lines    = split( /\n/, $input_old_gff );
my $last_id      = 0;
my $last_maker_id= '';
my $new_id_fasta_output;
my $new_id_index_output = '';
my $seq_counter  = 0;

#hash of locations from mRNA for simplicity
my %chr_location;
foreach my $line (@gff_lines) {
	if ($line =~ m/^#/){ next;}
	my @gff_line_arr = split ( "\t", $line);
	if ($gff_line_arr[2] eq 'mRNA'){
		my $chr = $gff_line_arr[0];
		$chr =~ s/\$chrprefix//; #remove genome version and sc, so only sc number
		my @gff_line_attr_arr = split (/;/, $gff_line_arr[8]);
		my $name;
		foreach my $attr (@gff_line_attr_arr){
				if ( $attr =~ m/^Name\=/ ){
					$name = $attr;
					$name =~ s/^Name\=//;
					$chr_location{$name} = $chr;
					last;
				}
		}
	}
}

#print number of mRNAs
print STDERR "\nNumber of mRNAs in GFF: ", scalar keys %chr_location, "\n\n";

foreach my $line (@lines) {
	chomp($line);
	my $new_id;

	if ( $line =~ m/^>/ ) {

		print $line."\n";
		$line =~ s/^>//;
		$new_id_index_output = $new_id_index_output.$line;

		my $current_maker_id = $line;
		$current_maker_id =~ s/[0-9]+$//;

		if(($last_maker_id ne $current_maker_id) && ($seq_counter > 0)){
			$last_id += 10; #namespace for 9 new genes now
		}

		$line =~ s/^[\S]+-mRNA-//;
		my $mRNA_rank = $line;
		#print "$mRNA_rank\n";

		#presuming a gene space of <1,00,000 so 6 characters, e.g. 00001 - 99999, was 1 mill earlier
		my @chars = split //,$last_id;
		my $padding_count = 5 - scalar @chars;

		#prefix with chr number
		$new_id = $prefix.$chr_location{$line};
		foreach (1..$padding_count){
			$new_id = $new_id.'0';
		}

		if ($mRNA_rank == 1){
			$new_id = $new_id.$last_id.'.1.'.$mRNA_rank;
			$last_maker_id = $current_maker_id;
		}
		elsif ($mRNA_rank > 1){
			#checks if all the mRNAs of a gene are grouped together
			if ($last_maker_id ne $current_maker_id){
				die "Id of current seq $current_maker_id with mRNA rank > 1 is not equal to previous id $last_maker_id\n";
			}
			$new_id = $new_id.$last_id.'.1.'.$mRNA_rank;
		}

		$new_id_fasta_output = $new_id_fasta_output.">$new_id\n";
		$new_id_index_output = $new_id_index_output."\t$new_id\n";

		$seq_counter++;
	}
	else{
		$new_id_fasta_output = $new_id_fasta_output."$line\n";
	}

}
unless ( open( OID, '>', "$new_id_fasta_output_file" ) ) {
	print STDERR "Cannot open $new_id_fasta_output_file\n";
	exit 1;
}
print OID $new_id_fasta_output;
close(OID);

unless ( open( OIX, '>', "$new_id_index_output_file" ) ) {
	print STDERR "Cannot open $new_id_index_output_file\n";
	exit 1;
}
print OIX $new_id_index_output;
close(OIX);

#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     Renames maker assigned gene names to gene name with version numbers. maker-ScVcwli_1-pred_gff_maker-gene-0.0-mRNA-1 becomes DcitrP00001.1.1. This is hard coded for <1,000,000 mRNAs. Counter skips over 10 gene models so manually curated genes can be added. It checks if all the mRNAs of a gene are grouped together but it does not check if all mRNAs are in order 1,2,3....


    Usage:
      update_maker_names_gff_chrlocs.pl

    Flags:

        -i  Maker Fasta file (required)
        -p  Prefix for name, e.g DcitrP (required)
        -s  Starting seed, e.g. 1 (required)
        -g  GFF file with chromosomal location (required)
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
