package  Bio::GenomeUpdate::TPF;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::TPF::TPFSequenceLine;
use Bio::GenomeUpdate::TPF::TPFGapLine;
use File::Slurp;

=head1 NAME

    TPF - Tiling path information used to generate a Tiling Path File (TPF)

=head1 SYNOPSIS

    my $variable = TPF->new();

=head1 DESCRIPTION

This class stores tiling path information including headers, sequence lines and gap lines, and generates a Tiling Path File (TPF).  The tiling path specifies the set of clones that will provide the best possible sequence coverage for a particular chromosome, the order of the clones along the chromosome, and the location of any gaps. 

=head2 Methods

=over 

=item C<set_version ( $version_number )>

Sets the tiling path file (TPF) specification version number.  Currently the default is version 1.7 (released April 9, 2012) and no other versions are supported at this time.  This is here to support any future changes in the TPF specifications.

=item C<get_version>

Gets the tiling path file (TPF) specification version number.

=cut

subtype 'TPFVersionNumber',
  as 'Str',
  where { $_ eq "1.7"},
  message {"The string, $_, was not a valid version number"};
has 'tpf_version' => (isa => 'TPFVersionNumber', is => 'rw', default => '1.7', clearer => 'clear_tpf_version');


=item C<set_assembly_version ( $organism_string )>

Sets the assembly_version (required).

=cut

has 'assembly_version' => (isa => 'Str', is => 'rw', clearer => 'clear_assembly_version');


=item C<set_organism ( $organism_string )>

Sets the organism (required).

=cut

has 'organism' => (isa => 'Str', is => 'rw', clearer => 'clear_organism');

=item C<set_assembly_name ( $assembly_name_string )>

Sets the assembly name (required).

=cut

has 'assembly_name' => (isa => 'Str', is => 'rw', clearer => 'clear_assembly_name');

=item C<set_chromosome ( $chromosome_number )>

Sets the chromosome (required).

=cut

has 'chromosome' => (isa => 'Str', is => 'rw', clearer => 'clear_chromosome');

=item C<set_strain_haplotype_cultivar ( $strain_haplotype_or_cultivar_string )>

Sets the strain, haplotype, or cultivar (optional).

=cut

has 'strain_haplotype_cultivar' => (isa => 'Str', is => 'rw', clearer => 'clear_strain_haplotype_cultivar');

=item C<set_type ( $type_string )>

Sets the TPF type (required).
The default type is "Complete Chromosome".

=cut

has 'type' => (isa => 'Str', is => 'rw', default => 'Complete Chromosome', clearer => 'clear_type');

=item C<set_comment ( $comment_string )>

Sets a comment for the TPF (optional).

=cut

has 'comment' => (isa => 'Str', is => 'rw', predicate => 'has_comment', clearer => 'clear_comment');

subtype 'TPFLine',
  as 'Bio::GenomeUpdate::TPF::TPFSequenceLine | Bio::GenomeUpdate::TPF::TPFGapLine',
  message {"The object was not a TPF sequence or gap line"};

has 'tpf_lines' => (isa => 'HashRef[TPFLine]',is => 'rw', predicate => 'has_tpf_lines', clearer => 'clear_tpf_lines');
 
sub add_line_to_end {
  my $self = shift;
  my $line_to_add = shift;
  my %lines;
  if ($self->has_tpf_lines()) {
    %lines = %{$self->get_tpf_lines()};
  }
  my $last_line = $self->get_number_of_lines();
  $lines{$last_line+1} = $line_to_add;
  $self->set_tpf_lines({%lines});
}

sub add_line_to_beginning {
  my $self = shift;
  my $line_to_add = shift;
  my %lines;
  if ($self->has_tpf_lines()) {
    %lines = %{$self->get_tpf_lines()};
    my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
    foreach my $line_key (@reverse_sorted_line_numbers) {
      $lines{$line_key+1} = $lines{$line_key};
    }
    $lines{$reverse_sorted_line_numbers[-1]} = $line_to_add;
  } else {
    $lines{1} = $line_to_add;
  }
  $self->set_tpf_lines({%lines});
}

sub insert_line_before {
  my $self = shift;
  my $line_number_to_insert_before = shift;
  my $line_to_add = shift;
  my %lines;
  if ($self->has_tpf_lines()) {
    %lines = %{$self->get_tpf_lines()};
    my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
    foreach my $line_key (@reverse_sorted_line_numbers) {
      if ($line_key >= $line_number_to_insert_before) {
	$lines{$line_key+1} = $lines{$line_key};
      }
    }
    $lines{$line_number_to_insert_before} = $line_to_add;
  } else {
    $lines{1} = $line_to_add;
  }
  $self->set_tpf_lines({%lines});
}

sub insert_line_after {
  my $self = shift;
  my $line_number_to_insert_after = shift;
  my $line_to_add = shift;
  my %lines;
  if ($self->has_tpf_lines()) {
    %lines = %{$self->get_tpf_lines()};
    my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
    foreach my $line_key (@reverse_sorted_line_numbers) {
      if ($line_key > $line_number_to_insert_after) {
	$lines{$line_key+1} = $lines{$line_key};
      }
    }
    $lines{$line_number_to_insert_after+1} = $line_to_add;
  } else {
    $lines{1} = $line_to_add;
  }
  $self->set_tpf_lines({%lines});
}

sub delete_line {
  my $self = shift;
  my $line_number_to_delete = shift;
  my %lines;
  if ($self->has_tpf_lines()) {
    %lines = %{$self->get_tpf_lines()};
    my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
    foreach my $line_key (@sorted_line_numbers) {
      if ($line_key > $line_number_to_delete) {
	$lines{$line_key-1} = $lines{$line_key};
      }
    }
    delete $lines{$sorted_line_numbers[-1]}; 
  } else {
    #add error warning;
  }
  $self->set_tpf_lines({%lines});
}

sub get_number_of_lines {
  my $self = shift;
  my %lines;
  if ($self->has_tpf_lines()) {
    %lines = %{$self->get_tpf_lines()};
    my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
    return $sorted_line_numbers[-1];
  } else {
    return 0;
  }
}

sub parse_tpf {
  my $self=shift;
  my $input_string = shift;
  $self->clear_organism();
  $self->clear_assembly_name();
  $self->clear_chromosome();
  $self->clear_strain_haplotype_cultivar();
  $self->clear_type();
  $self->clear_assembly_version();
  $self->clear_comment();
  $self->clear_tpf_lines();
  my @lines = split (/\n/, $input_string);
  my $tpf_data_has_begun = 0;
  my $tpf_data_has_ended = 0;
  foreach my $line (@lines) {
    chomp($line);
    if ($line =~ m/^\s*$/) {	#skip blank lines
      next;
    }
    if ($line =~m/^##/) {	#identify comment lines
      if ($line =~ m/ORGANISM: /i) {
	my @organism_line = split(/ORGANISM: /i,$line);
	my $organism = $organism_line[-1];
	$organism =~ s/^\s+|\s+$//g;
	$self->set_organism($organism);
	next;
      }
      if ($line =~ m/ASSEMBLY NAME: /i) {
	my @assembly_name_line = split(/ASSEMBLY NAME: /i,$line);
	my $assembly_name = $assembly_name_line[-1];
	$assembly_name =~ s/^\s+|\s+$//g;
	$self->set_assembly_name($assembly_name);
	next;
      }
      if ($line =~ m/CHROMOSOME: /i) {
	my @chromosome_line = split(/CHROMOSOME: /i,$line);
	my $chromosome = $chromosome_line[-1];
	$chromosome =~ s/^\s+|\s+$//g;
	$self->set_chromosome($chromosome);
	next;
      }
      if ($line =~ m/STRAIN\/HAPLOTYPE\/CULTIVAR: /i) {
	my @strain_haplotype_cultivar_line = split(/STRAIN\/HAPLOTYPE\/CULTIVAR: /i,$line);
	my $strain_haplotype_cultivar = $strain_haplotype_cultivar_line[-1];
	$strain_haplotype_cultivar =~ s/^\s+|\s+$//g;
	$self->set_strain_haplotype_cultivar($strain_haplotype_cultivar);
	next;
      }
      if ($line =~ m/TYPE: /i) {
	my @type_line = split(/TYPE: /i,$line);
	my $type = $type_line[-1];
	$type =~ s/^\s+|\s+$//g;
	$self->set_type($type);
	next;
      }
      if ($line =~ m/VERSION: /i) {
	my @assembly_version_line = split(/VERSION: /i,$line);
	my $assembly_version = $assembly_version_line[-1];
	$assembly_version =~ s/^\s+|\s+$//g;
	$self->set_assembly_version($assembly_version);
	next;
      }
      if ($line =~ m/COMMENT: /i) {
	my @comment_line = split(/COMMENT: /i,$line);
	my $comment = $comment_line[-1];
	$comment =~ s/^\s+|\s+$//g;
	$self->set_comment($comment);
	next;
      }
      if ($line =~ m/=== Beginning of TPF Data ===/) {
	$tpf_data_has_begun = 1;
	next;
      }
      if ($line =~ m/=== End of TPF Data ===/) {
	$tpf_data_has_ended = 1;
	next;
      }
    }

    if (($tpf_data_has_begun == 1) && ($tpf_data_has_ended == 0)) {
      my @tab_parsed_line = split(/\t/, $line);

      if (!defined($tab_parsed_line[0])) {
	#die with error line information
      }

      $tab_parsed_line[0] =~ s/^\s+|\s+$//g;

      if ($tab_parsed_line[0] eq "GAP") {
	my $tpf_gap_line = Bio::GenomeUpdate::TPF::TPFGapLine->new();
	if (!defined($tab_parsed_line[1])) {
	  #die with error missing gap line information
	  print STDERR "error in tpf\n";
	  print $line."\n";
	}
	$tab_parsed_line[1] =~ s/^\s+|\s+$//g;
	$tpf_gap_line->set_gap_type($tab_parsed_line[1]);
	if (defined($tab_parsed_line[2])) {
	  $tab_parsed_line[2] =~ s/^\s+|\s+$//g;
	  $tpf_gap_line->set_gap_size($tab_parsed_line[2]);
	}
	if (defined($tab_parsed_line[3])) {	
	  my @methods = split(/;/,$tab_parsed_line[3]);
	  foreach my $method (@methods) {
	    $method  =~ s/^\s+|\s+$//g;
	    $tpf_gap_line->add_gap_method($method);
	  }
	}
	$self->add_line_to_end($tpf_gap_line);
      } else {
	my $tpf_sequence_line = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
	if (defined($tab_parsed_line[0])) {
	  $tab_parsed_line[0] =~ s/^\s+|\s+$//g;
	  $tpf_sequence_line->set_accession($tab_parsed_line[0]);
	}
	if (defined($tab_parsed_line[1])) {
	  $tab_parsed_line[1] =~ s/^\s+|\s+$//g;
	  $tpf_sequence_line->set_clone_name($tab_parsed_line[1]);
	}
	if (defined($tab_parsed_line[2])) {
	  $tab_parsed_line[2] =~ s/^\s+|\s+$//g;
	  $tpf_sequence_line->set_local_contig_identifier($tab_parsed_line[2]);
	}
	if (defined($tab_parsed_line[3])) {
	  $tab_parsed_line[3] =~ s/^\s+|\s+$//g;
	  if ($tab_parsed_line[3] eq 'PLUS' || $tab_parsed_line[3] eq 'MINUS') {
	    $tpf_sequence_line->set_orientation($tab_parsed_line[3]);
	  } else {
	    $tpf_sequence_line->set_contains($tab_parsed_line[3]);
	    if (defined($tab_parsed_line[4])) {
	      $tab_parsed_line[4] =~ s/^\s+|\s+$//g;
	      $tpf_sequence_line->set_containing_accesion($tab_parsed_line[4]);
	    }
	    if (defined($tab_parsed_line[5])) {
	      $tab_parsed_line[5] =~ s/^\s+|\s+$//g;
	      $tpf_sequence_line->set_containing_clone_name($tab_parsed_line[4]);
	    }
	  }
	}
	$self->add_line_to_end($tpf_sequence_line);
      }
      next;
    }
  }
  return $self;
}

sub get_formatted_tpf {
  my $self = shift;
  my %lines;
  my $out_str;
  #Print header info
  $out_str .= "##ORGANISM: ".$self->get_organism()."\n";
  $out_str .= "##ASSEMBLY NAME: ".$self->get_assembly_name()."\n";
  $out_str .= "##CHROMOSOME: ".$self->get_chromosome()."\n";
  $out_str .= "##STRAIN/HAPLOTYPE/CULTIVAR: ".$self->get_strain_haplotype_cultivar()."\n";
  $out_str .= "##TYPE: ".$self->get_type()."\n";
  if ($self->has_comment()) {
    $out_str .= "##COMMENT: ".$self->get_comment()."\n";
  }
  $out_str .= "##=== Beginning of TPF Data ===\n";
  if ($self->has_tpf_lines()) {
    %lines = %{$self->get_tpf_lines()};
    my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
    foreach my $line_key (@sorted_line_numbers) {
      if ($lines{$line_key}->get_line_type() eq "sequence") {
	if ($lines{$line_key}->has_accession()) {
	  $out_str .= $lines{$line_key}->get_accession()."\t";
	} else {
	  $out_str .= "??\t";
	  print "accession not found\n";
	}
	if ($lines{$line_key}->has_clone_name()) {
	  $out_str .= $lines{$line_key}->get_clone_name()."\t";
	} else {
	  $out_str .= "?\t";
	}

	$out_str .= $lines{$line_key}->get_local_contig_identifier()."\t";
	if ($lines{$line_key}->has_contains()) {
	  $out_str .= $lines{$line_key}->get_contains()."\t";
	  if ($lines{$line_key}->has_containing_accession()) {
	    $out_str .= $lines{$line_key}->get_containing_accession()."\t";
	  } else {
	    $out_str .= "?\t";
	  }
	  if ($lines{$line_key}->has_containing_clone_name()) {
	    $out_str .= $lines{$line_key}->get_containing_clone_name();
	  } else {
	    $out_str .= "?";
	  }
	} else {
	  if ($lines{$line_key}->has_orientation()) {
	    $out_str .= $lines{$line_key}->get_orientation();
	  }
	}
      } elsif ($lines{$line_key}->get_line_type() eq "gap") {
	$out_str .= $lines{$line_key}->get_gap_identifier()."\t";
	$out_str .= $lines{$line_key}->get_gap_type();
	if ($lines{$line_key}->has_gap_size()) {
	  $out_str .= "\t".$lines{$line_key}->get_gap_size();
		    
	  if (!(($lines{$line_key}->get_gap_type() eq "CENTROMERE"|| 
		 $lines{$line_key}->get_gap_type() eq "TELOMERE" || 
		 $lines{$line_key}->get_gap_type() eq "HETEROCHROMATIN" ||
		 $lines{$line_key}->get_gap_type() eq "SHORT-ARM"))) {
	    if ( $lines{$line_key}->has_gap_methods()) {
	      $out_str .= "\t";
	      my $first_gmethod = 1;
	      foreach my $gmethod (@{$lines{$line_key}->get_gap_methods()}) {
		if ($first_gmethod == 0) {
		  $out_str .= ";";
		}
		$first_gmethod = 0;
		$out_str .= $gmethod; 
	      } 
	    } else {
	      #add error warning method required
	    }
	  }
	} else {
	  if ($lines{$line_key}->get_gap_type() eq "CENTROMERE"||
	      $lines{$line_key}->get_gap_type() eq "TELOMERE" ||
	      $lines{$line_key}->get_gap_type() eq "HETEROCHROMATIN"|| 
	      $lines{$line_key}->get_gap_type() eq "SHORT-ARM") {
	    #add error warning biological gap must have size
	  }
	}
		
		
      } else {
	#add error warning;
      }
      $out_str .= "\n";
    }
  } else {
    #add error warning;
  }
  #$out_str .= footer
  $out_str .= "##=== End of TPF Data ===\n";
  return $out_str;
}

###
  1;				#do not remove
###

=pod

  =back

  =head1 LICENSE

  Same as Perl.

  =head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

  =cut
