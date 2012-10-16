package Bio::GenomeUpdate::AGP;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::AGP::AGPSequenceLine;
use Bio::GenomeUpdate::AGP::AGPGapLine;
use Bio::GenomeUpdate::AGP::AGPConvert;
use File::Slurp;

=head1 NAME

    AGP - Accessioned Golden Path (AGP) describes the assembly of a larger sequence object from smaller objects.

=head1 SYNOPSIS

    my $variable = TPF->new();

=head1 DESCRIPTION

    This class stores Accessioned Golden Path (AGP) information including headers, sequence lines and gap lines, and generates a formatted AGP file.

=head2 Methods

=over

=item C<set_version ( $version_number )>

Sets the Accessioned Golden Path (AGP) specification version number.  The default is version 2.0.  Version 1.1 is currently the only other supported version.

=item C<get_version>

Gets the tiling path file (TPF) specification version number.

=cut

subtype 'AGPVersionNumber',
  as 'Str',
  where { $_ eq "2.0" ||  $_ eq "1.1" },
  message {"The string, $_, was not a valid version number"};
has 'version' => (isa => 'AGPVersionNumber', is => 'rw', default => '2.0', predicate => 'has_version', clearer => 'clear_version');


=item C<set_organism ( $organism_string )>

Sets the organism (required).

=cut

has 'organism' => (isa => 'Str', is => 'rw', predicate => 'has_organism', clearer => 'clear_organism');

=item C<set_tax_id ( $tax_id_string )>

Sets the NCBI taxonomic id (required).

=cut

has 'tax_id' => (isa => 'Str', is => 'rw', predicate => 'has_tax_id', clearer => 'clear_tax_id');

=item C<set_assembly_name ( $assembly_name_string )>

Sets the assembly name (required).

=cut

has 'assembly_name' => (isa => 'Str', is => 'rw', predicate => 'has_assembly_name', clearer => 'clear_assembly_name');

=item C<set_assembly_date ( $assembly_date_string )>

Sets the assembly date (required).

=cut

has 'assembly_date' => (isa => 'Str', is => 'rw', predicate => 'has_assembly_date', clearer => 'clear_assembly_date');

=item C<set_genome_center ( $genome_center_string )>

Sets the genome center (required).

=cut

has 'genome_center' => (isa => 'Str', is => 'rw', predicate => 'has_genome_center', clearer => 'clear_genome_center');

=item C<set_description ( $type_string )>

Sets the description (optional).

=cut

has 'description' => (isa => 'Str', is => 'rw', predicate => 'has_description', clearer => 'clear_description');

=item C<set_comment_lines ( @array_of_strings )>

Sets comment lines from an array of strings (optional).

=cut

has 'comment_lines' => (isa => 'ArrayRef[Str]', is => 'rw', predicate => 'has_comment_lines', clearer => 'clear_comment_lines');

=item C<add_comment_line ( @comment_string )>

Adds a comment line (optional).

=cut

sub add_comment_line {
  my $self = shift;
  my $comment_to_add = shift;
  my @lines;
  if ($self->has_comment_lines()) {
    @lines = @{$self->get_comment_lines()};
  }
  push (@lines, $comment_to_add);
  $self->set_comment_lines([@lines]);
}

subtype 'AGPLine',
  as 'AGPSequenceLine | AGPGapLine',
  message {"The object was not an AGP sequence or gap line"};

has 'agp_lines' => (isa => 'HashRef[AGPLine]',is => 'rw', predicate => 'has_agp_lines', clearer => 'clear_agp_lines');

sub add_line_to_end {
  my $self = shift;
  my $line_to_add = shift;
  my $line_number = shift;
  my %lines;
  if ($self->has_agp_lines()) {
    %lines = %{$self->get_agp_lines()};
  }
  if (defined($line_number)) {
    $lines{$line_number} = $line_to_add;
  } else {
    my $last_line = $self->get_number_of_lines();
    $lines{$last_line+1} = $line_to_add;
    $lines{$last_line+1}->set_line_number($last_line+1);
  }
  $self->set_agp_lines({%lines});
}


sub add_line_to_beginning {
  my $self = shift;
  my $line_to_add = shift;
  $line_to_add->set_line_number(1);
  my %lines;
  if ($self->has_agp_lines()) {
    %lines = %{$self->get_agp_lines()};
    my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
    foreach my $line_key (@reverse_sorted_line_numbers) {
      $lines{$line_key+1} = $lines{$line_key};
      $lines{$line_key+1}->set_line_number($line_key+1);
    }
    $lines{$reverse_sorted_line_numbers[-1]} = $line_to_add;
  } else {
    $lines{1} = $line_to_add;
  }
  # $lines{1}->set_line_number(1);
  $self->set_agp_lines({%lines});
}

sub insert_line_before {
  my $self = shift;
  my $line_number_to_insert_before = shift;
  my $line_to_add = shift;
  $line_to_add->set_line_number($line_number_to_insert_before);
  my %lines;
  if ($self->has_agp_lines()) {
    %lines = %{$self->get_agp_lines()};
    my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
    foreach my $line_key (@reverse_sorted_line_numbers) {
      if ($line_key >= $line_number_to_insert_before) {
	$lines{$line_key+1} = $lines{$line_key};
	$lines{$line_key+1}->set_line_number($line_key+1);
      }
    }
    $lines{$line_number_to_insert_before} = $line_to_add;
    $lines{$line_number_to_insert_before}->set_line_number($line_number_to_insert_before);
  } else {
    $lines{1} = $line_to_add;
    $lines{1}->set_line_number(1);
  }
  $self->set_agp_lines({%lines});
}

sub insert_line_after {
  my $self = shift;
  my $line_number_to_insert_after = shift;
  my $line_to_add = shift;
  $line_to_add->set_line_number($line_number_to_insert_after+1);
  my %lines;
  if ($self->has_agp_lines()) {
    %lines = %{$self->get_agp_lines()};
    my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
    foreach my $line_key (@reverse_sorted_line_numbers) {
      if ($line_key > $line_number_to_insert_after) {
	$lines{$line_key+1} = $lines{$line_key};
	$lines{$line_key+1}->set_line_number($line_key+1);
      }
    }
    $lines{$line_number_to_insert_after+1} = $line_to_add;
    $lines{$line_number_to_insert_after}->set_line_number($line_number_to_insert_after);
  } else {
    $lines{1} = $line_to_add;
    $lines{1}->set_line_number(1);
  }
  $self->set_agp_lines({%lines});
}

sub delete_line {
  my $self = shift;
  my $line_number_to_delete = shift;
  my %lines;
  if ($self->has_agp_lines()) {
    %lines = %{$self->get_agp_lines()};
    my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
    foreach my $line_key (@sorted_line_numbers) {
      if ($line_key > $line_number_to_delete) {
	$lines{$line_key-1} = $lines{$line_key};
	$lines{$line_key-1}->set_line_number($line_key-1);
      }
    }
    delete $lines{$sorted_line_numbers[-1]}; 
  } else {
    #add error warning;
  }
  $self->set_agp_lines({%lines});
}

sub get_number_of_lines {
  my $self = shift;
  my %lines;
  if ($self->has_agp_lines()) {
    %lines = %{$self->get_agp_lines()};
    return scalar(keys %lines);
  } else {
    return 0;
  }
}

sub get_largest_line_number {
  my $self = shift;
  my %lines;
  if ($self->has_agp_lines()) {
    %lines = %{$self->get_agp_lines()};
    my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
    return $sorted_line_numbers[-1];
  } else {
    return 0;
  }
}

sub print_formatted_agp {
  my $self = shift;
  my $out_file = shift;
  $out_file = '>'.$out_file;
  open(OUTFILE, $out_file) || die ("Unable to open output file $out_file\n");

  my %lines;
  my $default_unknown_gap_size = 100;
    
  #Print header info
  print OUTFILE "##agp-version\t".$self->get_version()."\n";
  print OUTFILE "# ORGANISM: ".$self->get_organism()."\n";
  print OUTFILE "# TAX_ID: ".$self->get_tax_id()."\n";
  print OUTFILE "# ASSEMBLY NAME: ".$self->get_assembly_name()."\n";
  #print OUTFILE "# ASSEMBLY DATE: ".$self->get_assembly_date()->format_cldr("DD-MMMM-YYYY")."\n";
  print OUTFILE "# ASSEMBLY DATE: ".$self->get_assembly_date()."\n";
  print OUTFILE "# GENOME CENTER: ".$self->get_genome_center()."\n";
  print OUTFILE "# DESCRIPTION: ".$self->get_description()."\n";
  if ($self->has_comment_lines()) {
    print OUTFILE "# COMMENTS:\n";
    foreach my $comment (@{$self->get_comment_lines()}) {
      print OUTFILE "# ".$comment."\n";
    }
  }
  if ($self->has_agp_lines()) {
    %lines = %{$self->get_agp_lines()};
    my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
    foreach my $line_key (@sorted_line_numbers) {
      print OUTFILE $lines{$line_key}->get_object_being_assembled()."\t";
      print OUTFILE $lines{$line_key}->get_object_begin()."\t";
      print OUTFILE $lines{$line_key}->get_object_end()."\t";
      #print OUTFILE $lines{$line_key}->get_line_number()."\t";
      print OUTFILE $line_key."\t";
      print OUTFILE $lines{$line_key}->get_component_type()."\t";
      if ($lines{$line_key}->get_line_type() eq "sequence") {
	print OUTFILE $lines{$line_key}->get_component_id()."\t";
	print OUTFILE $lines{$line_key}->get_component_begin()."\t";
	print OUTFILE $lines{$line_key}->get_component_end()."\t";
	print OUTFILE $lines{$line_key}->get_orientation();

      } elsif ($lines{$line_key}->get_line_type() eq "gap") {
	if ($lines{$line_key}->get_component_type() eq "N") {
	  print OUTFILE $lines{$line_key}->get_gap_length()."\t";
	} else {
		    
	  print OUTFILE $default_unknown_gap_size."\t";
	}
	print OUTFILE $lines{$line_key}->get_gap_type()."\t";
	print OUTFILE $lines{$line_key}->get_linkage()."\t";
	my $first_evidence_line=1;
	if ($lines{$line_key}->get_linkage() eq "yes") {
	  foreach my $evidence (@{$lines{$line_key}->get_linkage_evidence()}) {
	    if ($first_evidence_line != 1) {
	      print OUTFILE ";";
	    }
	    $first_evidence_line=0;
	    print OUTFILE $evidence;
	  }
	} else {
	  print OUTFILE "na";
	}
      } else {
	#add warning that type is not sequence or gap
      }
      print OUTFILE "\n";
    }
  }

}

  sub parse_agp {
    my $self = shift;
    my $input_file = shift;
    $self->clear_version();
    $self->clear_organism();
    $self->clear_tax_id();
    $self->clear_assembly_name();
    $self->clear_assembly_date();
    $self->clear_genome_center();
    $self->clear_description();
    $self->clear_comment_lines();
    $self->clear_agp_lines();
    
    my @lines = read_file($input_file);
    my $load_comment_lines = 0;
    my $line_count = scalar(@lines);
    my $line_counter = 0;
    foreach my $line (@lines) {
      chomp($line);
      if ($line =~m/^#/) {
	if ($load_comment_lines == 1) {
	  my @comment_line = split(/#/,$line);
	  my $comment = $comment_line[-1];
	  $comment =~ s/^\s+|\s+$//g;
	  $self->add_comment_line($comment);
	  next;
	}
	if (($line =~ m/^##/) && ($line =~ m/agp-version/i)) {
	  my @version_line = split(/\s+/,$line);
	  my $version = $version_line[-1];
	  $version =~ s/^\s+|\s+$//g;
	  $self->set_version($version);
	  next;
	}
	if ($line =~ m/ORGANISM:/i) {
	  my @organism_line = split(/ORGANISM:/i,$line);
	  my $organism = $organism_line[-1];
	  $organism =~ s/^\s+|\s+$//g;
	  $self->set_organism($organism);
	  next;
	}
	if ($line =~ m/TAX_ID:/i) {
	  my @tax_id_line = split(/TAX_ID:/i,$line);
	  my $tax_id = $tax_id_line[-1];
	  $tax_id =~ s/^\s+|\s+$//g;
	  $self->set_tax_id($tax_id);
	  next;
	}
	if ($line =~ m/ASSEMBLY NAME:/i) {
	  my @assembly_name_line = split(/ASSEMBLY NAME:/i,$line);
	  my $assembly_name = $assembly_name_line[-1];
	  $assembly_name =~ s/^\s+|\s+$//g;
	  $self->set_assembly_name($assembly_name);
	  next;
	}
	if ($line =~ m/ASSEMBLY DATE:/i) {
	  my @assembly_date_line = split(/ASSEMBLY DATE:/i,$line);
	  my $assembly_date = $assembly_date_line[-1];
	  $assembly_date =~ s/^\s+|\s+$//g;
	  $self->set_assembly_date($assembly_date);
	  next;
	}
	if ($line =~ m/GENOME CENTER:/) {
	  my @genome_center_line = split(/GENOME CENTER:/,$line);
	  my $genome_center = $genome_center_line[-1];
	  $genome_center =~ s/^\s+|\s+$//g;
	  $self->set_genome_center($genome_center);
	  next;
	}
	if ($line =~ m/DESCRIPTION:/i) {
	  my @description_line = split(/DESCRIPTION:/i,$line);
	  my $description = $description_line[-1];
	  $description =~ s/^\s+|\s+$//g;
	  $self->set_description($description);
	  next;
	}
	if ($line =~ m/COMMENTS:/i) {
	  $load_comment_lines = 1;
	  next;
	}
      } else {
	$line_counter++;
	my @tab_parsed_line = split(/\t/, $line);
	if (!defined($tab_parsed_line[0])) {
	  next;
	}
	print "Parsing line ".$line_counter." of ".$line_count."\t".$tab_parsed_line[0]."\t".$tab_parsed_line[1]."\n";
	if (($tab_parsed_line[4] eq "N") | ($tab_parsed_line[4] eq "U")) {
	  my $agp_gap_line = AGPGapLine->new();
	  $agp_gap_line->set_object_being_assembled($tab_parsed_line[0]);
	  $agp_gap_line->set_object_begin($tab_parsed_line[1]);
	  $agp_gap_line->set_object_end($tab_parsed_line[2]);
	  $agp_gap_line->set_line_number($tab_parsed_line[3]);
	  $agp_gap_line->set_component_type($tab_parsed_line[4]);
	  $agp_gap_line->set_gap_length($tab_parsed_line[5]);
	  #convert from 1.1 to 2.0 AGP specification for 'fragment' gap type
	  if ($tab_parsed_line[6] eq 'fragment') {
	    $agp_gap_line->set_gap_type('scaffold');
	    $agp_gap_line->set_linkage('yes');
	    if ($tab_parsed_line[7] eq 'yes') {
	      $agp_gap_line->add_linkage_evidence('paired-ends');
	    } elsif ($tab_parsed_line[7] eq 'no') {
	      $agp_gap_line->add_linkage_evidence('within_clone');
	    } else {
	      # add warning linkage should be yes or no
		    
	    }
	  }
	  #convert from 1.1 to 2.0 AGP specification for 'clone' gap type
	  elsif ($tab_parsed_line[6] eq 'clone') {
	    if ($tab_parsed_line[7] eq 'yes') {
	      $agp_gap_line->set_gap_type('scaffold');
	      $agp_gap_line->set_linkage('yes');
	      $agp_gap_line->add_linkage_evidence('clone_contig');
	    } elsif ($tab_parsed_line[7] eq 'no') {
	      $agp_gap_line->set_gap_type('contig');
	      $agp_gap_line->set_linkage('no');
	    } else {
	      #add warning 
	    }
	  } else {
	    $agp_gap_line->set_gap_type($tab_parsed_line[6]);
	    $agp_gap_line->set_linkage($tab_parsed_line[7]);
	    if ($tab_parsed_line[7] eq "yes") {
	      my @linkage_evidence = split(/;/,$tab_parsed_line[8]);
	      foreach my $evidence_item (@linkage_evidence) {
		$agp_gap_line->add_linkage_evidence($evidence_item);
	      }
	    }
	  }
	  $self->add_line_to_end($agp_gap_line, $line_counter);
	} else {
	  my $agp_sequence_line = AGPSequenceLine->new();
	  $agp_sequence_line->set_object_being_assembled($tab_parsed_line[0]);
	  $agp_sequence_line->set_object_begin($tab_parsed_line[1]);
	  $agp_sequence_line->set_object_end($tab_parsed_line[2]);
	  $agp_sequence_line->set_line_number($tab_parsed_line[3]);
	  $agp_sequence_line->set_component_type($tab_parsed_line[4]);
	  $agp_sequence_line->set_component_id($tab_parsed_line[5]);
	  $agp_sequence_line->set_component_begin($tab_parsed_line[6]);
	  $agp_sequence_line->set_component_end($tab_parsed_line[7]);
	  $agp_sequence_line->set_orientation($tab_parsed_line[8]);
	  $self->add_line_to_end($agp_sequence_line, $line_counter);
	}
      }
    }
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
