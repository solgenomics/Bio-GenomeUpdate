#!/usr/bin/perl

=head1 NAME

process_bac_assembly.pl

=head1 SYNOPSIS

process_bac_assembly.pl -f [ACE file] -m [mismatch %] -o [output directory]

=head1 COMMAND-LINE OPTIONS

 -f  ACE file from Phrap assembly (required)
 -m  Mismatch percentage (recommended 0.5, required)
 -t  Do a test run (e.g. -t 1, no ACE file required in this case)
 -d  Print extra status messages (e.g. -d 1)
 -o  Output directory
 -h  Help

=NOTES

If the [mismatch %] is more than threshold then a new ACE file error_contigs.ace can be created with only those erroneous 
     contigs. Manually explore the erroneous contigs in a ACE viewer (Tablet http://ics.hutton.ac.uk/tablet/) and remove the 
     misfit BAC(s) from the contigs and treat as singletons. Reassemble the rest of the BACs in the contig and add to final 
     assembled BAC set.
     
Many routines mentioned in Bio::Assembly::Contig documentation ARE NOT IMPLEMENTED. 

TODO: 
 - Fix contig_mismatch() to work with ACE files from phrap.
 - Fix add_contig_to_scaffold to work with test ACE object (wrong aligned coordinates reported) and ACE files from phrap.  

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

use Cwd;
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

Accepts a single contig object from an assembly. Returns fasta sequences for contig and assembled BACs.

=cut
sub contig_to_fasta {
	my $contig = shift;
#	print Dumper ($contig);
	my @seqs = $contig->get_seq_ids();
	my $BAC_fasta = '';
	foreach my $seqname (@seqs){
		my $seq = $contig->get_seq_by_name($seqname); #returns Bio::LocatableSeq
#		print STDERR $seq->id(),"\n";
#		print STDERR $seq->seq(),"\n";
		my $cleaned_seq = $seq->seq();
		$cleaned_seq =~ s/-//g;
#		print STDERR $cleaned_seq."\n";
		$BAC_fasta = $BAC_fasta.'>'.$seq->id()."\n".$cleaned_seq."\n";
	}
	my $contig_fasta = '>'.$contig->id()."\n".$contig->get_consensus_sequence()->seq()."\n";
	return ($contig_fasta,$BAC_fasta)
}

=item C<singlet_to_fasta ( Bio::Assembly::Singlet  )>

Accepts a single singlet object from an assembly. Returns string with Fasta sequence. 

=cut

sub singlet_to_fasta {
	my $singlet = shift;
	my $seq = $singlet->seqref();
#	print STDERR $seq->id(),"\n";
#	print STDERR $seq->seq(),"\n";
	
	if ($seq->seq() =~ /-/){
		print STDERR "Singlet ",$seq->id()," has a gapped sequence which is not expected. Please investigate the assembly ACE file.  Exiting...\n";
		exit 0;
	}
	my $fasta = '>'.$seq->id()."\n".$seq->seq()."\n";
	return ($fasta);
}

=item C<contig_mismatch ( Bio::Assembly::Contig  )>

Accepts a single contig object from an assembly. Returns floats containing the mismatch percentage and the gap percentage for the contig.

=cut
sub contig_mismatch {
	my ($contig,$debug) = (@_);
#	print STDERR ref($contig);
	my $consensus_sequence = $contig->get_consensus_sequence()->seq(); #calling Bio::Seq->seq()
	if ($consensus_sequence =~ /-/){
		print STDERR "Contig ",$contig->id()," consensus has a gapped sequence which is not expected. Please investigate the assembly ACE file.  Exiting...\n";
	}
	my @consensus_sequence_arr = split '',$consensus_sequence; 
	
	#get reads (BACs), positions and consensus. Compare reads to consensus positions to get mismatch.
	my $featureDB = $contig->get_features_collection(); # returns Bio::DB::SeqFeature::Store::memory
#	print STDERR ref($featureDB)."\n";
	my $mismatches = 0;
	my $gaps = 0;
	foreach my $feature ($featureDB->get_all_features()){#returns Bio::SeqFeature::Generic for each read or BAC in the contig
		print STDERR ref($feature)."\n";
		my $aligned_start = $feature->start();
		my $aligned_end = $feature->end();
		
		#Workaround for exception when end > length. Happens because start and end in parent consensus seq coordinate space
		#Reset the start and end to valid values to get the sequence out
		#------------- EXCEPTION: Bio::Root::Exception ------------- MSG: Bad end parameter. End must be less than the total length of sequence
		$feature->end($feature->length());
		$feature->start(1);
		if ($debug){
			#print STDERR "Read ".$feature->seq()->id()." start $aligned_start end $aligned_end\n";
			print STDERR "Read ".$feature->seq_id()." start $aligned_start end $aligned_end\n";
		}

		my $aligned_seq = $feature->seq()->seq();
		my @aligned_seq_arr = split '', $aligned_seq;
		my $aligned_seq_pos = 0;
		for (my $con_pos = $aligned_start - 1 ; $con_pos <= $aligned_end - 1; $con_pos++){
			if ($consensus_sequence_arr[$con_pos] ne $aligned_seq_arr[$aligned_seq_pos]){
				$mismatches++;
			}
			if ($aligned_seq_arr[$aligned_seq_pos] eq '-'){
				$gaps++;	
			}
			$aligned_seq_pos++;
		}
		if ($debug){
			print STDERR "cumulative gaps $gaps mismatches $mismatches\n";
		}
	}
	my $mismatch_percent = $mismatches / scalar @consensus_sequence_arr;
	my $gap_percent = $gaps / scalar @consensus_sequence_arr;
	return ($mismatch_percent,$gap_percent);
}

=item C<add_contig_to_scaffold  ( Bio::Assembly::Scaffold, Bio::Assembly::Contig  )>

Accepts a single contig object from an assembly and a test scaffold object. Returns updated scaffold object. ABANDONED for now as new ACE file 
has wrong coordinates for test ACE assembly. Can just write out the fasta for the BACs in the contig for re-assembly with phrap.

=cut
sub add_contig_to_scaffold {
#	my ($ace,$contig,$debug) = shift;
#	my ($scaffold,$contig,$debug) = @_;#need to use for passing >1 objects or it fails
#	print STDERR ref($scaffold)."\n";
#	print STDERR ref($contig)."\n";
#	my $status = $scaffold->add_contig($contig);
#	
#	if ($debug && $status){
#		print STDERR $contig->id()." successfully written\n";
#	}elsif(!$status){
#		print STDERR $contig->id()." was not added to scaffold. Exiting..\n";
#		exit 0;
#	} 
#	return $scaffold;
}

=item C<run_tests ()>

Runs a test of all functions with dummy data. No ACE file required in this case.
=cut

sub run_tests{
	my $debug = shift;
	#create scaffold
	my (@contigs,@singlets,$scaffold);
	print STDERR "Running tests..\ncreating test data..\n";
	#create contig
	my $c1 =  Bio::Assembly::Contig->new(-id => 'contig1');
	my $ls1 = Bio::LocatableSeq->new(-seq=>"ACCG-T", -id=>"bac1", -alphabet=>'dna');
	my $ls2 = Bio::LocatableSeq->new(-seq=>"ACA-CG-T", -id=>"bac2", -alphabet=>'dna');
	my $ls1_coord = Bio::SeqFeature::Generic->new(-start=>3, -end=>8, -strand=>1);
	my $ls2_coord = Bio::SeqFeature::Generic->new(-start=>1, -end=>8, -strand=>1);
	$c1->add_seq($ls1); $c1->set_seq_coord($ls1_coord,$ls1);
	$c1->add_seq($ls2);	$c1->set_seq_coord($ls2_coord,$ls2);

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
	#had to add contig and singlet manually
	$scaffold->add_contig($c1);
	$scaffold->add_singlet($s1);

	#print summary
	scaffold_summary($scaffold);

	#get seqs and compare
	my $ctr = 1;
	foreach my $contig  ($scaffold->all_contigs()){
		if ($debug){
			print STDERR "read contig $ctr\n";	
		}
		my($mismatch_percent,$gap_percent) = contig_mismatch($contig,$debug);
		if ($debug){
			print STDERR "mismatch % $mismatch_percent \ngap % $gap_percent\n";
		}
		if (($mismatch_percent == 0.125) && ($gap_percent == 0.375)){
		    print STDERR "contig_mismatch() test successful\n";
		}
		my ($contig_fasta,$BAC_fasta) = contig_to_fasta($contig);
		if ($debug){
			print STDERR $contig_fasta.$BAC_fasta;
		}
		my $contig_fasta_test = ">contig1\nACACCG-T\n";
		my $BAC_fasta_test = ">bac1\nACCGT\n>bac2\nACACGT\n";
		if (($contig_fasta eq $contig_fasta_test) && ($BAC_fasta eq $BAC_fasta_test)){
		    print STDERR "contig_fasta() test successful\n";
		}
		$ctr++;
	}
	
	$ctr = 1;
	foreach my $singlet  ($scaffold->all_singlets()){
		if ($debug){
			print STDERR "read singlet $ctr\n";
		}
		my $singlet_fasta = singlet_to_fasta($singlet);
		if ($debug){
			print STDERR $singlet_fasta;
		}
		my $singlet_fasta_test = ">bac3\nATGGGGGTGGTGGTACCCT\n";
		if ($singlet_fasta eq $singlet_fasta_test){
			print STDERR "singlet_fasta() test successful\n";
		}
		$ctr++
	}
	
	#create new ACE (for error contigs)
#	my $scaffold_out = Bio::Assembly::Scaffold->new (-id => 'erroneous_contigs',
#					 -source => 'test_program',
#					);
##	$scaffold_out = contig_to_scaffold($scaffold_out,$c1,$debug);
#	$scaffold_out->add_contig($c1);
#	my $ace_out = Bio::Assembly::IO->new(	-file   => ">test.ace",
#											-format => 'ace' );
##	$ace_out->write_assembly( -scaffold => $scaffold_out);
##	print STDERR Dumper($c1)."\n";
#	foreach my $contig ($scaffold_out->all_contigs()){
#		$ace_out->write_contig($contig);
#	}
#	
#	$ace_out->write_header();
#	$ace_out->write_footer();
#	
#	my $ace_test_in = Bio::Assembly::IO->new(-file   => "test.ace",
#										-format => 'ace' );
#	scaffold_summary($ace_test_in->next_assembly());
#	unlink 'test.ace';	
	print STDERR "ALL TESTS PASSED FOR MANUALLY CREATED ACE OBJECT. MAY NOT WORK FOR ACE FILE FROM PHRAP\n\n";
}

our ( $opt_f, $opt_m, $opt_t, $opt_d, $opt_o, $opt_h );
getopts('f:m:t:d:o:h');

if (defined $opt_t) {$opt_t = !defined $opt_t ? 0 : 1 unless $opt_t == 0;} #assign 0 if no value, else assign 1 unless value == 0
if (defined $opt_d) {$opt_d = !defined $opt_d ? 0 : 1 unless $opt_d == 0;}
$opt_o ||= 'process_bac_assembly_out';
if ($opt_h) {
	help();
	exit;
}
if ($opt_t){
	run_tests($opt_d);
	exit;
}

if ( (!$opt_f) || (!$opt_m) ) {
	print "\nACE file and mismatch % are required. See help below\n\n\n";
	help();
}

#prep input data
my $assembly = Bio::Assembly::IO->new( -file => $opt_f, -format => 'ace'); 
my $scaffold = $assembly->next_assembly();

#process assembly
scaffold_summary($scaffold);


#prep dirs
my $cwd = getcwd();
my $path_contigs = "$cwd/${opt_o}/contigs/";
my $path_singletons = "$cwd/${opt_o}/singletons/";
if (!( -d "$cwd/${opt_o}")) { 
	mkdir "$cwd/${opt_o}"; 
	mkdir "$cwd/${opt_o}/contigs/";
	mkdir "$cwd/${opt_o}/singletons/";
}

my $ctr = 0;
foreach my $contig ($scaffold->all_contigs()){
	if ($opt_d){
		print STDERR "read contig $ctr\n";	
	}
	#if mismatch % > threshold, write contig fasta to error_contigs directory 
#	my($mismatch_percent,$gap_percent) = contig_mismatch($contig,$opt_d);
#	if ($opt_d){
#		print STDERR "mismatch % $mismatch_percent \ngap % $gap_percent\n";
#	}

	my ($contig_fasta,$BAC_fasta) = contig_to_fasta($contig);
	if ($opt_d){
		print STDERR $contig_fasta.$BAC_fasta;
	}
	
	#Not using mismatch as Bio::Seq objects not created properly in contig_mismatch() for test ACE file. 
	#Works fine for hand created ACE object.
#	if ($mismatch_percent >= $opt_m){
#		open (BF , ">${opt_o}\/poor\/".$contig->id()."_BAC.fas");
#		print BF $BAC_fasta;
#		close BF;
#	}else{
#		open (CF , ">${opt_o}\/good\/".$contig->id().".fas");
#		print CF $contig_fasta;
#		close CF;
#		
#		open (BF , ">${opt_o}\/good\/".$contig->id()."_BAC.fas");
#		print BF $BAC_fasta;
#		close BF;
#	}

	#writing out all contigs and BACs
	my $contig_file = $path_contigs.$contig->id();
#	print STDERR  "$contig_file.fas\n";
	open (CF , '>',"$contig_file.fas");
	print CF $contig_fasta;
	close CF;
	
	open (BF ,'>', "${contig_file}_BAC.fas");
	print BF $BAC_fasta;
	close BF;
	
	$ctr++;
}
print STDERR "$ctr contigs written to $path_contigs\n";

$ctr = 0;
foreach my $singlet  ($scaffold->all_singlets()){
	if ($opt_d){
		print STDERR "read singlet $ctr\n";
	}
	my $singlet_fasta = singlet_to_fasta($singlet);
	if ($opt_d){
		print STDERR $singlet_fasta;
	}
	my $singleton_file = $path_singletons.$singlet->seqref()->id();
	#open (SF , ">${opt_o}\/good\/singleton_".$singlet->seqref()->id().".fas");
	open (SF , '>',"${singleton_file}.fas");
	print SF $singlet_fasta;
	close SF;

	$ctr++
}
print STDERR "$ctr singletons written to $path_singletons\n";


#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script analyzes a ACE file from a Phrap (http://www.phrap.org/phredphrapconsed.html) assembly of BACs. Contigs with 
     1 read (BAC) are written to singleton_BACs.fas. Contigs with multiple reads (BACs) are written to contigs_BACs.fas and the 
     corresponding meta-data to contigs_BACs.txt. 
     
     If the [mismatch %] is more than threshold then Fasta of BACs in contig will be written to <Contig name> directory for 
     re-assembly with phrap. Not creating new ACE file as Bio::Assembly::IO::ace writes wrong coordinates for BACs.

     Phrap parameters (recommended) to generate assembly and ACE file
      phrap -new_ace -shatter_greedy -penalty -4 -minmatch 20 FILE.fas

    Usage:
      process_bac_assembly.pl -f [ACE file] -m [mismatch %] -o [output directory]
      
    Flags:

         -f  ACE file from Phrap assembly (required)
         -m  Mismatch percentage (recommended 0.5, required)
         -t  Do a test run (e.g. -t 1, no ACE file required)
         -d  Print extra status messages (e.g. -d 1)
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
