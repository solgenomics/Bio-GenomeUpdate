#!/usr/bin/perl -w

Does not work for mutiple regions in the same seq. Array parsing is broken

=head1 NAME

replace_with_Ns_fasta.pl

=head1 SYNOPSIS

replace_with_Ns_fasta -f [input fasta file] -i [file with seqname, start and end base]

=head1 COMMAND-LINE OPTIONS

 -f  Old Fasta file (required)
 -i  File with seqname, start and end base (required)
 -h  help

=cut

use strict;
use warnings;

use Getopt::Std;
use File::Slurp;
use Bio::Seq;
use Bio::SeqIO;


our ( $opt_f, $opt_i, $opt_h );
getopts('f:i:h');
if ($opt_h) {
	help();
	exit;
}

if ( !$opt_f || !$opt_i ) {
	print "\n fasta file, file with seqname, start and end base is required. See help below\n\n\n";
	help();
}
 
my $in_fasta  = Bio::SeqIO->new(-file => "<$opt_f", -format => 'Fasta');
chomp $opt_i;
my $out_fasta  = Bio::SeqIO->new(-file => ">$opt_i.softmasked.$opt_f", -format => 'Fasta');

# hash of starts and stops
my $region_input_file = $opt_i;
my $region_input   = read_file($region_input_file);

my %region_hash;
my @lines = split( /\n/, $region_input );
foreach my $line (@lines) {
  chomp($line);
  my @line_arr = split ("\t", $line);
  my $region = "$line_arr[1] $line_arr[2]";
  
  if ( !exists $region_hash{$line_arr[0]} ){
  	my @arr;
  	push @arr, $region; 
  	$region_hash{$line_arr[0]} = @arr;				#seq id = startbase stopbase
  }
  else{
  	my @arr = $region_hash{$line_arr[0]};
  	push @arr, $region; 
  	$region_hash{$line_arr[0]} = @arr;				#seq id = startbase stopbase
  }
}

while ( my $seq_object = $in_fasta->next_seq ) {

	if ( exists $region_hash{$seq_object->id} ) {
		my $sequence = $seq_object->seq();

		my @arr = $region_hash{$seq_object->id};
		foreach my $region (@arr) {
			my @region_arr = split (' ', $region);
			my $start      = $region_arr[0];
			my $stop       = $region_arr[1];
			for my $position ( $start-1 .. $stop-1 ){
				substr ($sequence,$position,1) = 'n';
			}
		}

		my $out_seq_object = Bio::Seq->new( -display_id => $seq_object->id, -seq => $sequence);
		$out_fasta->write_seq($out_seq_object);
	}
	else{
		$out_fasta->write_seq($seq_object);
	}
}






#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     Replaces the region in the fasta file with small Ns (n)
    


    Usage:
      replace_with_Ns_fasta -f [input fasta file] -i [file with seqname, start and end base]
      
    Flags:

		-f  Old Fasta file (required)
 		-i  File with seqname, start and end base (required)
 		-h  help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

