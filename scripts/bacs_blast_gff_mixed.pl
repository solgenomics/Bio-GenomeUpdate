#!/usr/bin/perl

=head1 NAME

bacs_blast_gff_mixed.pl

=head1 SYNOPSIS

bacs_blast_gff_mixed.pl -g [ GFF file] -l [ 1000 ]

=head1 COMMAND-LINE OPTIONS

 -g  GFF file from blast output (required, use Blast2HitGFF3.pl)
 -l  minimum length of feature or HSP (default 1000bp)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;

use Getopt::Std;

our ( $opt_g, $opt_l, $opt_d, $opt_h );
getopts('g:l:d:h');
if ($opt_h) {
	help();
	exit;
}

if ( !$opt_g ) {
	print "\n GFF3 of BAC alignments to reference chromosome is required. Use blastn and Blast2HitGFF3.pl. See help below\n\n\n";
	help();
}

my $minimum_length = defined $opt_l ? $opt_l : 1000;
unless (open (GFF, "<${opt_g}")) {die "Not able to open ${opt_g}\n\n";}
my ($bacs_mixed, $bacs_mixed_gff, $prev_rec_strand, $prev_rec_bac, $prev_rec_ref, $prev_rec_start, $prev_rec_stop);
$bacs_mixed = $bacs_mixed_gff = '';

while (my $rec = <GFF>){
	if( $rec =~ /^#/ ){ next;}
	
	my @rec_arr = split ("\t", $rec );
	
	if ( $rec_arr[4] - $rec_arr[3] + 1 > $minimum_length ){
		my @rec_attributes_arr = split (/\;/, $rec_arr[8]);
		if (!defined($prev_rec_strand) && !defined($prev_rec_bac) && !defined($prev_rec_ref)){
			$prev_rec_strand = $rec_arr[6];
			$prev_rec_bac    = $rec_attributes_arr[1];
			$prev_rec_bac =~ s/Name=//;
			$prev_rec_ref = $rec_arr[0];
		}
		
		my ($rec_strand, $rec_bac, $rec_ref);
		$rec_strand = $rec_arr[6];
		$rec_bac    = $rec_attributes_arr[1];
		$rec_bac =~ s/Name=//;
		$rec_ref = $rec_arr[0];
		
		if ( ( $rec_ref eq $prev_rec_ref ) && ( $rec_bac eq $prev_rec_bac ) && ( $rec_strand eq $prev_rec_strand ) ){
			next;
		}
		elsif ( ( $rec_ref eq $prev_rec_ref ) && ( $rec_bac eq $prev_rec_bac ) && ( $rec_strand ne $prev_rec_strand ) ){
			$bacs_mixed     = $bacs_mixed."\n".$rec_bac;
			
			if ($opt_d){ print STDERR "$rec_bac in mixed orientation at \n$rec";}
			
			$bacs_mixed_gff = $bacs_mixed_gff.$rec;
		}
#		#if bac changes or ref changes
#		elsif ( ( $rec_ref ne $prev_rec_ref ) || ( $rec_bac ne $prev_rec_bac ) ){
#			
#		}
		
		$prev_rec_ref    = $rec_ref;
		$prev_rec_bac    = $rec_bac;
		$prev_rec_strand = $rec_strand;
	}
}
close(GFF);

unless (open (OGFF, ">mixed.${opt_g}")) {die "Not able to open mixed.${opt_g}\n\n";}
print OGFF $bacs_mixed_gff;
close (OGFF);
unless (open (OBACS, ">${opt_g}.mixed_bacs.full_names")) {die "Not able to open ${opt_g}.mixed_bacs.full_names\n\n";}
print OBACS $bacs_mixed;
close (OBACS);




#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

    This script will print a list of BACs that align in mixed orientation to the reference chromosome. These regions could be errors in the reference assembly. This does NOT check if the alignments to the reference are OUT OF ORDER on the BAC. It does NOT check if the HSPs are in the same region in the reference sequence.
     
    Usage:
      bacs_blast_gff_mixed.pl -g [ GFF file] -l [ 1000 ]
      
    Flags:

          -g  GFF file from blast output (required, use Blast2HitGFF3.pl)
          -l  minimum length of feature or HSP (default 1000bp)
          -d  debugging messages (1 or 0)
          -h  Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha at cornell.edu , @SahaSurya>

=cut

