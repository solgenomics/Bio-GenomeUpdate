package Bio::GenomeUpdate::AGPConvert;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::TPF;
use Bio::GenomeUpdate::AGP;

=head1 NAME

    AGPConvert - Converts an Accessioned Golden Path (AGP) to other formats.

=head1 SYNOPSIS

    my $variable = TPF->new();

=head1 DESCRIPTION

This class converts Accessioned Golden Path (AGP) information including headers, sequence lines and gap lines, into other formats (e.g., Tiling Path Files (TPF)). 

=head2 Methods

=over 

=item C<set_agp_to_convert ( $agp_object )>

Sets the Accessioned Golden Path (AGP) object to convert.

=cut

has 'agp_to_convert' => (isa => 'AGP', is => 'rw', predicate => 'has_agp_to_convert');

=item C<set_component_id_is_local ( $int )>

Set to 0 if the component IDs in the AGP object are NCBI accessions and set to 1 if they are local IDs (default 0).

=cut

has 'component_id_is_local' => (isa => 'Int', is => 'rw', default => 0);

=item C<set_ncbi_to_local_conversion ( $hash_ref )>

Sets a hash reference to use for converting NCBI accessions to local IDs.

=cut

has 'ncbi_to_local_conversion' => (isa => 'HashRef[Str]', is => 'rw', predicate => 'has_ncbi_to_local_conversion');

=item C<set_local_to_ncbi_conversion ( $hash_ref )>

Sets a hash reference to use for converting NCBI accessions to local IDs.

=cut

has 'local_to_ncbi_conversion' => (isa => 'HashRef[Str]', is => 'rw', predicate => 'has_local_to_ncbi_conversion');

=item C<set_strain_haplotype_cultivar ( $str )>

Sets the strain/haplotype/cultivar to list in converted files.

=cut

has 'strain_haplotype_cultivar' => (isa => 'Str', is => 'rw', default => '?');

sub to_tpf {
  my $self = shift;
  my $agp = $self->get_agp_to_convert();
  

  my %agp_lines = %{$agp->get_agp_lines()};
  my @sorted_agp_line_numbers = sort { $a <=> $b } keys %agp_lines;
  my %seen_agp_objects = ();
  my @uniq_agp_objects = (); 
  foreach my $agp_line_key (@sorted_agp_line_numbers) {
    push(@uniq_agp_objects, $agp_lines{$agp_line_key}->get_object_being_assembled()) unless $seen_agp_objects{$agp_lines{$agp_line_key}->get_object_being_assembled()}++;	
  }
    
  my @tpfs_to_return;
    
  foreach my $unique_agp_object (@uniq_agp_objects) {
    my $tpf = TPF->new();
    $tpf->set_organism($agp->get_organism());
    $tpf->set_assembly_name($agp->get_assembly_name());
    $tpf->set_strain_haplotype_cultivar($self->get_strain_haplotype_cultivar());
    $tpf->set_comment("@{$agp->get_comment_lines()}");
    $tpf->set_chromosome($unique_agp_object);
    foreach my $agp_line_key (@sorted_agp_line_numbers) {
      if ($agp_lines{$agp_line_key}->get_object_being_assembled() eq $unique_agp_object) {
	if ($agp_lines{$agp_line_key}->get_line_type() eq "sequence") {
	  my $tpf_sequence_line = TPFSequenceLine->new();
	  if ($self->get_component_id_is_local() == 0) {
	    print 
	      $tpf_sequence_line->set_accession($agp_lines{$agp_line_key}->get_component_id());
	    my %clone_lookup = %{$self->get_ncbi_to_local_conversion()};
	    if (defined($clone_lookup{$agp_lines{$agp_line_key}->get_component_id()})) {
	      $tpf_sequence_line->set_clone_name($clone_lookup{$agp_lines{$agp_line_key}->get_component_id()});
	      print "local: ".$clone_lookup{$agp_lines{$agp_line_key}->get_component_id()}.": ",$agp_lines{$agp_line_key}->get_component_id()."\n";
	    } else {
	      $tpf_sequence_line->set_clone_name('?');
	      print "ncbi accession for local clone ID not found\n";
	    }
	  } else {
	    $tpf_sequence_line->set_clone_name($agp_lines{$agp_line_key}->get_component_id());
	    my %ncbi_lookup = %{$self->get_local_to_ncbi_conversion()};
	    if (defined($ncbi_lookup{$agp_lines{$agp_line_key}->get_component_id()})) {
	      my $ncbi_id = $ncbi_lookup{$agp_lines{$agp_line_key}->get_component_id()};
	      $ncbi_id =~ s/(.*?)\.\d+/$1/; #remove version number
	      $tpf_sequence_line->set_accession($ncbi_id);
	      print "ncbi: $ncbi_id\n";
	      print "accesson from sequence line: ";
	      print $tpf_sequence_line->get_accession();
	      print "\n";
	    } else {
	      $tpf_sequence_line->set_accession("?"); 
	      print "ncbi accession for local clone ID not found\n";
	    }
			
	  }
	  $tpf_sequence_line->set_local_contig_identifier(""); #add ability to lookup later
	  if ($agp_lines{$agp_line_key}->has_orientation()) {
	    if ($agp_lines{$agp_line_key}->get_orientation() eq '+') {
	      $tpf_sequence_line->set_local_contig_identifier("PLUS");
	    }
	    if ($agp_lines{$agp_line_key}->get_orientation() eq '-') {
	      $tpf_sequence_line->set_local_contig_identifier("MINUS");
	    }
	  }
	  $tpf->add_line_to_end($tpf_sequence_line);
	}
	if ($agp_lines{$agp_line_key}->get_line_type() eq 'gap') {
	  my $tpf_gap_line = TPFGapLine->new();
	  #set to TYPE-2 (clone) gap for gaps of unknown length
	  if ($agp_lines{$agp_line_key}->get_component_type() eq 'U') {  
	    $tpf_gap_line->set_gap_type('TYPE-2'); 
	  }
	  #set to TYPE-3 gap for contig or scaffold gaps
	  elsif (($agp_lines{$agp_line_key}->get_gap_type() eq 'contig') || ($agp_lines{$agp_line_key}->get_gap_type() eq 'scaffold')) {  
	    $tpf_gap_line->set_gap_type('TYPE-3');
	  }
	  #set biological type gaps
	  elsif ($agp_lines{$agp_line_key}->get_gap_type() eq 'centromere') {
	    $tpf_gap_line->set_gap_type('CENTROMERE');
	  } elsif ($agp_lines{$agp_line_key}->get_gap_type() eq 'telomere') {
	    $tpf_gap_line->set_gap_type('TELOMERE');
	  } elsif ($agp_lines{$agp_line_key}->get_gap_type() eq 'heterochromatin') {
	    $tpf_gap_line->set_gap_type('HETEROCHROMATIN');
	  } elsif ($agp_lines{$agp_line_key}->get_gap_type() eq 'short_arm') {
	    $tpf_gap_line->set_gap_type('SHORT-ARM');
	  } else {
	    #add warning missing acceptable gap type
	  }
	  if ($agp_lines{$agp_line_key}->has_gap_length()) {
	    $tpf_gap_line->set_gap_size($agp_lines{$agp_line_key}->get_gap_length());
	  } elsif ($agp_lines{$agp_line_key}->get_gap_type() eq ('centromere' || 'telomere' || 'heterochromatin' || 'short_arm')) {
	    #add warning: gap size required for biological gap
	  }
	  #gap method required if gap size field is populated and the gap type is not biological
	  if ($tpf_gap_line->has_gap_size() && !($agp_lines{$agp_line_key}->get_gap_type() eq 'centromere' ||
						 $agp_lines{$agp_line_key}->get_gap_type() eq 'telomere' || 
						 $agp_lines{$agp_line_key}->get_gap_type() eq 'heterochromatin' || 
						 $agp_lines{$agp_line_key}->get_gap_type() eq 'short_arm')) {
	    if ($agp_lines{$agp_line_key}->has_linkage_evidence()) {
	      my @linkage_evidence = @{$agp_lines{$agp_line_key}->get_linkage_evidence()};
	      foreach my $evidence_item (@linkage_evidence) {
		if ($evidence_item eq 'paired-ends') {
		  $tpf_gap_line->add_gap_method('PAIRED ENDS');
		} elsif ($evidence_item eq 'align_genus') {
		  $tpf_gap_line->add_gap_method('ALIGN GENUS');
		} elsif ($evidence_item eq 'align_xgenus') {
		  $tpf_gap_line->add_gap_method('ALIGN XGENUS');
		} elsif ($evidence_item eq 'align_trnscpt') {
		  $tpf_gap_line->add_gap_method('ALIGN TRNSCPT');
		} else {
		  #add warning: gap method could be translated to a valid TPF method type
		}
	      }
	    } else {
	      #add warning: missing required gap method
	    }
	  }
	  $tpf->add_line_to_end($tpf_gap_line);
	}
      }
    }

    push (@tpfs_to_return, $tpf);
  }

  return @tpfs_to_return;
}


###
1;				#do not remove
###

=back

=head1 LICENSE

    Same as Perl.

=head1 AUTHORS

    Jeremy D. Edwards <jde22@cornell.edu>   

=cut
