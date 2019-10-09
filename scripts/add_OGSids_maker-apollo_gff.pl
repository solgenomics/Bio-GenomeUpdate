#!/usr/bin/perl

=head1 NAME

add_OGSids_maker-apollo_GFF.pl

=head1 SYNOPSIS

add_OGSids_maker-apollo_GFF.pl -g [old Merged Maker and Apollo GFF file] -a [AHRD file for mRNA with Maker and Apollo ids] -p [species acronym or prefix] -c [Prefix for chromosome in GFF]  -s [starting value for naming] -o [output formatted gff file with OGS ids]

=head1 COMMAND-LINE OPTIONS

 -g  Merged Maker and Apollo GFF file (required)
 -a  AHRD tab separated file for mRNA with Maker and Apollo ids, e.g. maker-DC3.0sc00-snap-gene-91.44-mRNA-1       Glycerol-3-phosphate dehydrogenase [NAD(P)+] (AHRD V3.11 *** tr|A0A0A9WH09|A0A0A9WH09_LYGHE) (required)
 -f  Curated descriptions from Apollo, e.g. 18712462-6e8c-478a-83c8-a40c9d0be977  Dihydrolipoyl transacetylase (required)
 -p  Prefix for name, e.g DcitrP (required)
 -c  Prefix for chromosome in GFF e.g. Dc3.0sc (required)
 -s  Starting seed, e.g. 1 (required)
 -o  output GFF file with OGS ids
 -h  Help

=cut

use strict;
use warnings;

use File::Slurp;
use Getopt::Std;
use Bio::GFF3::LowLevel qw (gff3_parse_feature  gff3_format_feature gff3_parse_attributes);

our ( $opt_g, $opt_a, $opt_p, $opt_s, $opt_c, $opt_o, $opt_h );
getopts('g:a:p:s:c:g:h');
if ($opt_h) {
  help();
  exit;
}
if ( !$opt_g || !$opt_a || !$opt_p || !$opt_c || !$opt_s || !$opt_o) {
  print
"\nOld GFF file, AHRD file, name prefix, chr prefix, starting seed, output GFF file is required.
See help below\n\n\n";
  help();
}

#get input files
my $gff_input_file = $opt_g;
my $gff_input   = read_file($gff_input_file)
  or die "Could not open old fasta input file: $gff_input_file\n";
my $ahrd_input_file = $opt_a;
my $ahrd_input      = read_file($ahrd_input_file)
  or die "Could not open AHRD input file: $ahrd_input_file\n";
chomp $opt_p; my $prefix = $opt_p;
chomp $opt_c; my $chrprefix = $opt_c;
my $seed = $opt_s;
if ($seed !~ /^[0-9]+$/){ die "$seed should be a number\n"; }

#output variables
my ($gff_output, $index_output, $desc_output);
















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

     Renames maker and Apollo assigned gene names to gene name with version numbers, e.g.  maker-ScVcwli_1-pred_gff_maker-gene-0.0-mRNA-1 becomes DcitrP00001.1.1. Creates an index file with old and new mRNA ids. This is hard coded for <1,000,000 mRNAs. Counter skips over 10 gene models so manually curated genes can be added.
     
    NOTE:


    Usage:
      cmdline_perldoc.pl -g [?? file] <other params>
      
    Flags:

		  -g  Merged Maker and Apollo GFF file (required)
      -a  AHRD tab separated file for mRNA with Maker and Apollo ids, e.g. maker-DC3.0sc00-snap-gene-91.44-mRNA-1       Glycerol-3-phosphate dehydrogenase [NAD(P)+] (AHRD V3.11 *** tr|A0A0A9WH09|A0A0A9WH09_LYGHE) (required)
      -f  Curated descriptions from Apollo, e.g. 18712462-6e8c-478a-83c8-a40c9d0be977  Dihydrolipoyl transacetylase (required)
      -p  Prefix for name, e.g DcitrP (required)
      -c  Prefix for chromosome in GFF e.g. Dc3.0sc (required)
      -s  Starting seed, e.g. 1 (required)
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

