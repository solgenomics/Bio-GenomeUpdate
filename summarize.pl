#!/usr/bin/perl

=head1 NAME

summarize.pl

=head1 SYNOPSIS

summarize.pl -f [TPF/AGP file] -t [type tpf/agp]

=head1 COMMAND-LINE OPTIONS

 -f  TPF/AGP file (required)
 -t  type (tpf or agp) (required)
 -h  Help

=cut

use strict;
use warnings;

use Bio::GenomeUpdate::AGP;
use Bio::GenomeUpdate::TPF;
use File::Slurp;
use Getopt::Std;
use PDL;

our ( $opt_f, $opt_t, $opt_h );
getopts('f:t:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_f || !$opt_t ) {
	print "\nTPF or AGP file and type are required. See help below\n\n\n";
	help();
}

#prep input data
my $input_file = $opt_f;
my $input      = read_file($input_file)
  or die "Could not open file: $input_file\n";

my $type = $opt_t;
if (($type ne lc 'tpf') && ($type ne lc 'agp')){
	print STDERR "Invalid type. Only tpf or agp are accepted.\n";
	help();
}

if ($type eq lc 'tpf'){
	my $tpf = Bio::GenomeUpdate::TPF->new();
	$tpf->parse_tpf($input);
	print STDERR 'Number of gaps: '.$tpf->get_number_of_gap_lines()."\n";
	print STDERR 'Number of sequences: '.$tpf->get_number_of_sequence_lines()."\n";
	
	my @gap_lengths_from_tpf=$tpf->get_gap_lengths();
	my $piddle = pdl @gap_lengths_from_tpf;
	my ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($piddle);
	print STDERR "\nGap statistics:\n";
	print STDERR "Mean: $mean\nMedian: $median\nMin: $min\nMax: $max\nAvg abs dev: $adev\nStd dev: $rms\n";
	#Std dev is fine here as we know the population of gaps
}
elsif($type eq lc 'agp'){
	my $agp = Bio::GenomeUpdate::AGP->new();
	$agp->parse_agp($input);
	print STDERR 'Number of gaps: '.$agp->get_number_of_gap_lines()."\n";
	print STDERR 'Number of sequences: '.$agp->get_number_of_sequence_lines()."\n";
	
	my @gap_lengths_from_agp=$agp->get_gap_lengths();
	my $piddle = pdl @gap_lengths_from_agp;
	my ($mean,$prms,$median,$min,$max,$adev,$rms) = stats($piddle);
	print STDERR "\nGap statistics:\n";
	print STDERR "Mean: $mean\nMedian: $median\nMin: $min\nMax: $max\nAvg abs dev: $adev\nStd dev: $rms\n";
	
	my @sequence_lengths_from_agp=$agp->get_sequence_lengths();
	$piddle = pdl @sequence_lengths_from_agp;
	($mean,$prms,$median,$min,$max,$adev,$rms) = stats($piddle);
	print STDERR "\nComponent statistics:\n";
	print STDERR "Mean: $mean\nPop std dev: $prms\nMedian: $median\nMin: $min\nMax: $max\nAvg abs dev: $adev\nStd dev: $rms\n";
	#Std dev is fine here as we know the population of gaps
}
else{
	#this should not happen
}




#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script prints a summary of TPF or AGP file.

    Usage:
      summarize.pl -f [TPF/AGP file] -t [type]
      
    Flags:

         -f  TPF/AGP file (required)
         -t  type (tpf or agp) (required)
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