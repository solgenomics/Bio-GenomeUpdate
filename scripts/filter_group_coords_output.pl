#!/usr/bin/perl

=head1 NAME

filter_group_coords_output.pl

=head1 SYNOPSIS

filter_group_coords_output.pl -f [ align_BACends_group_coords.pl or group_coords.pl output file]

=head1 COMMAND-LINE OPTIONS

 -f  align_BACends_group_coords.pl or group_coords.pl output file (required)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;

use Getopt::Std;

our ( $opt_f, $opt_d, $opt_h );
getopts('f:d:h');
if ($opt_h) {
	help();
	exit;
}

if ( !$opt_f ) {
	print "\n align_BACends_group_coords.pl or group_coords.pl output file is required. See help below\n\n\n";
	help();
}

unless (open (IN, "<$opt_f")) { die "Cannot open $opt_f";	}
my @input_group_coords;
my @sorted_input_group_coords;
my $line;

my $header = <IN>; #skip header
while ($line = <IN>){
	my @temp_arr = split ("\t",$line);
	push @input_group_coords, \@temp_arr;
}
close(IN);

#sort by start position on ref
#query	reference	ref_start	ref_end	length .....
@sorted_input_group_coords = sort { $a->[2] <=> $b->[2] } @input_group_coords;

if ($opt_d) { print STDERR "Finished sorting $opt_f by ref start coordinate\n\n";}

#remove multiple coverage
unless (open (OR, ">redundant.${opt_f}")) {die "Not able to open redundant.${opt_f}\n\n";}
print OR $header;
unless (open (OF, ">filtered.${opt_f}")) {die "Not able to open filtered.${opt_f}\n\n";}
print OF $header;

my $line_query;
my $prev_line_query;
my $line_ref_start;
my $prev_line_ref_start;
my $line_ref_end;
my $prev_line_ref_end;
my $line_counter;
my $redundant_line_counter;
my $filtered_line_counter;
$line_counter = $redundant_line_counter = $filtered_line_counter = 1;

foreach $line (@sorted_input_group_coords){
	if ((!defined $prev_line_ref_start) && (!defined $prev_line_ref_end)){ #first line
		$prev_line_query = $line->[0];
		if ( $line->[3] > $line->[2]){ #pos strand
			$prev_line_ref_start = $line->[2];
			$prev_line_ref_end = $line->[3];
		}
		elsif ( $line->[2] > $line->[3]){ #neg strand
			$prev_line_ref_start = $line->[3];
			$prev_line_ref_end = $line->[2];
		}
		else{
			die "Unexpected orientation in $opt_f in line $line_counter\n";
		}
		
		print OF join ("\t",@$line);
		$line_counter++;
		$filtered_line_counter++;
		next;
	}
	
	$line_query = $line->[0];
	if ( $line->[3] > $line->[2]){ #pos strand
		$line_ref_start = $line->[2];
		$line_ref_end = $line->[3];
	}
	elsif ( $line->[2] > $line->[3]){ #neg strand
		$line_ref_start = $line->[3];
		$line_ref_end = $line->[2];
	}
	else{
		die "Unexpected orientation in $opt_f in line $line_counter\n";
	}
	
	# <-----------------------------> ref
	# <--------------->               bac1
	#          <-------------->       bac2
	#     <-------------------------> bac3
	if ( ($line_ref_start >= $prev_line_ref_start) && ($line_ref_end <= $prev_line_ref_end) ){
		if ($opt_d) { 
			print STDERR "Line $line_counter for $line_query BAC is contained in $prev_line_query\n";
			print OR "\nLine $line_counter for $line_query BAC is contained in $prev_line_query\n";
		}
		print OR join ("\t",@$line);
		$redundant_line_counter++;
	}
	else{
		print OF join ("\t",@$line);
		$filtered_line_counter++;
		$prev_line_query      = $line_query;
		$prev_line_ref_start  = $line_ref_start;
		$prev_line_ref_end    = $line_ref_end;
	}
	
	$line_counter++;
}

print STDERR "\n\nLines in $opt_f\t\t:\t$line_counter\n";
print STDERR "Lines in redundant.${opt_f}\t:\t$redundant_line_counter\n";
print STDERR "Lines in filtered.${opt_f}\t:\t$filtered_line_counter\n\n";

#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

    Some BACs and BAC contigs cover the same region on the reference. This script will keep only one BAC or Contig BAC for a reference region and filters out the redundant BAC hits. It also sorts by ref location as group_coords does not do that.
     
    Usage:
      filter_group_coords_output.pl -f [ align_BACends_group_coords.pl or group_coords.pl output file]
      
    Flags:

		 -f  align_BACends_group_coords.pl or group_coords.pl output file (required)
		 -d  debugging messages (1 or 0)
		 -h  Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl. Change??????????

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

