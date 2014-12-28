#!/usr/bin/perl

=head1 NAME

process_bac_assembly.pl

=head1 SYNOPSIS

process_bac_assembly.pl -f [ACE file] -m [mismatch %] -o [output directory]

=head1 COMMAND-LINE OPTIONS

 -f  ACE file from Phrap assembly (required)
 -m  Mismatch percentage (recommended 0.5, required)
 -t  Do a test run (no ACE file required)
 -o  Output directory
 -h  Help

=TODO

Add test param and function (see Aure testing prez)

=cut

use strict;
use warnings;
use Getopt::Std;
use Bio::Assembly::IO;
use Bio::Assembly::Contig;
use Bio::Assembly::Scaffold;
use Bio::SeqFeature::Generic;
use Bio::Assembly::Singlet;
use Bio::Seq;
use Bio::LocatableSeq;

use Data::Dumper;

=item C<scaffold_summary  ( Bio::Assembly::Scaffold  )>

Accepts a scaffold object from an assembly. Prints basis statistics.
=cut

sub scaffold_summary{
	my $scaffold = shift;
	print STDERR "Scaffold or Assembly name: ".$scaffold->id()."\n";
	print STDERR "Number of BACs in scaffold: ".$scaffold->get_nof_seqs()."\n";
	print STDERR "Number of Contigs: ".$scaffold->get_nof_contigs()."\n";
	print STDERR "Number of BACs in Contigs: ".$scaffold->get_nof_contig_seqs()."\n";
	print STDERR "Number of Singlets: ".$scaffold->get_nof_singlets()."\n";
	
}


=item C<contig_to_fasta ( Bio::Assembly::Contig  )>

Accepts a single contig object from an assembly. Returns Bio::SeqIO::fasta object and the corresponding meta-data in a string variable. 

=cut

sub contig_to_fasta {
	my $contig = shift;
	
	my @seqs = $contig->get_seq_ids();
	
	if (scalar @seqs == 1){
		my $seq = $contig->get_seq_by_name($seqs[0]);
		print STDERR $seq->seq,"\n";
		#return ($seq[0] , )	
	}
	else{
		
	}
	
	
	
	
}


=item C<contig_to_ace  ( Bio::Assembly::Contig  )>

Accepts a single contig object from an assembly. Returns a Bio::Assembly::IO::ace object for the contig.If the [mismatch %] is more than threshold then a new ACE file error_contigs.ace will be created with only those erroneous 
     contigs.
=cut

=item C<run_tests ()>

Runs a test of all functions with dummy data. No ACE file required in this case.
=cut

sub run_tests{
	#create scaffold
	my (@contigs,@singlets,$scaffold);
	
	#create contig
	my $c1 =  Bio::Assembly::Contig->new(-id => 'contig1');
	my $ls1 = Bio::LocatableSeq->new(-seq=>"ACCG-T", -id=>"bac1", -alphabet=>'dna');
    my $ls2 = Bio::LocatableSeq->new(-seq=>"ACA-CG-T", -id=>"bac2", -alphabet=>'dna');
    my $ls1_coord = Bio::SeqFeature::Generic->new(-start=>3, -end=>8, -strand=>1);
    my $ls2_coord = Bio::SeqFeature::Generic->new(-start=>1, -end=>8, -strand=>1);
	$c1->add_seq($ls1);
	$c1->add_seq($ls2);
	$c1->set_seq_coord($ls1_coord,$ls1);
	$c1->set_seq_coord($ls2_coord,$ls2);

	my $con1 = Bio::LocatableSeq->new(-seq=>"ACACCG-T", -alphabet=>'dna');
	$c1->set_consensus_sequence($con1);

	#create singlet
	my $seq = Bio::Seq->new(-id=>'bac3', -seq=>'ATGGGGGTGGTGGTACCCT');
	my $s1 = Bio::Assembly::Singlet->new(-id=>'singlet1', -seqref=>$seq);
	
	$scaffold = Bio::Assembly::Scaffold->new (-id => 'assembly1',
					 -source => 'test_program',
#					 -contigs => \@contigs, these do not work 
#					 -singlets => \@singlets
					);
	#had to add contig and scaffold manually
	$scaffold->add_contig($c1);
	$scaffold->add_singlet($s1);

	#print summary
	scaffold_summary($scaffold);

	#get seqs and compare
	
	
}

our ( $opt_f, $opt_m, $opt_t, $opt_o, $opt_h );
getopts('f:m:t:o:h');
if ($opt_h) {
	help();
	exit;
}
if ($opt_t){
	run_tests();
	exit;
}

if ( !$opt_f ) {
	print "\nACE file is required. See help below\n\n\n";
	help();
}

#prep input data
my $assembly = Bio::Assembly::IO->new( -file => $opt_f, -format => 'ace'); 


#process contigs
while (my $contig = $assembly->next_contig()){
	#if mismatch % > threshold, add to new ACE file error_contigs.ace 

	#if only 1 BAC, write out original BAC to singleton_BACs.fas
	#If >1 BAC write to contigs_BACs.fas and corresponding meta-data to contigs_BACs.txt. 
	contig_to_fasta($contig);
	
	
	
}





#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script analyzes a ACE file from a Phrap (http://www.phrap.org/phredphrapconsed.html) assembly of BACs. Contigs with 
     1 read (BAC) are written to singleton_BACs.fas. Contigs with multiple reads (BACs) are written to contigs_BACs.fas and the 
     corresponding meta-data to contigs_BACs.txt. 
     
     If the [mismatch %] is more than threshold then a new ACE file error_contigs.ace will be created with only those erroneous 
     contigs. Manually explore the erroneous contigs in a ACE viewer (Tablet http://ics.hutton.ac.uk/tablet/) and remove the 
     misfit BAC(s) from the contigs and treat as singletons. Reassemble the rest of the BACs in the contig and add to final 
     assembled BAC set.
     
     Phrap parameters (recommended) to generate assembly and ACE file
      phrap -new_ace -shatter_greedy -penalty -4 -minmatch 20 FILE.fas

    Usage:
      process_bac_assembly.pl -f [ACE file] -m [mismatch %] -o [output directory]
      
    Flags:

         -f  ACE file from Phrap assembly (required)
         -m  Mismatch percentage (recommended 0.5, required)
         -t  Do a test run (no ACE file required)
         -o  Output directory 
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