package  Bio::GenomeUpdate::TPF;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::TPF::TPFSequenceLine;
use Bio::GenomeUpdate::TPF::TPFGapLine;
use File::Slurp;

use Data::Dumper;#for debugging

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

subtype 'TPFVersionNumber', as 'Str', where { $_ eq "1.7" },
  message { "The string, $_, was not a valid version number" };
has 'tpf_version' => (
	isa     => 'TPFVersionNumber',
	is      => 'rw',
	default => '1.7',
	clearer => 'clear_tpf_version'
);

=item C<set_assembly_version ( $organism_string )>

Sets the assembly_version (required).

=cut

has 'assembly_version' =>
  ( isa => 'Str', is => 'rw', clearer => 'clear_assembly_version' );

=item C<set_organism ( $organism_string )>

Sets the organism (required).

=cut

has 'organism' => ( isa => 'Str', is => 'rw', clearer => 'clear_organism' );

=item C<set_assembly_name ( $assembly_name_string )>

Sets the assembly name (required).

=cut

has 'assembly_name' =>
  ( isa => 'Str', is => 'rw', clearer => 'clear_assembly_name' );

=item C<set_chromosome ( $chromosome_number )>

Sets the chromosome (required).

=cut

has 'chromosome' => ( isa => 'Str', is => 'rw', clearer => 'clear_chromosome' );

=item C<set_strain_haplotype_cultivar ( $strain_haplotype_or_cultivar_string )>

Sets the strain, haplotype, or cultivar (optional).

=cut

has 'strain_haplotype_cultivar' =>
  ( isa => 'Str', is => 'rw', clearer => 'clear_strain_haplotype_cultivar' );

=item C<set_type ( $type_string )>

Sets the TPF type (required).
The default type is "Complete Chromosome".

=cut

has 'type' => (
	isa     => 'Str',
	is      => 'rw',
	default => 'Complete Chromosome',
	clearer => 'clear_type'
);

=item C<set_comment ( $comment_string )>

Sets a comment for the TPF (optional).

=cut

has 'comment' => (
	isa       => 'Str',
	is        => 'rw',
	predicate => 'has_comment',
	clearer   => 'clear_comment'
);

subtype 'TPFLine',
  as 'Bio::GenomeUpdate::TPF::TPFSequenceLine | Bio::GenomeUpdate::TPF::TPFGapLine',
  message { "The object was not a TPF sequence or gap line" };

has 'tpf_lines' => (
	isa       => 'HashRef[TPFLine]',
	is        => 'rw',
	predicate => 'has_tpf_lines',
	clearer   => 'clear_tpf_lines'
);

sub add_line_to_end {
	my $self        = shift;
	my $line_to_add = shift;
	my %lines;
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
	}
	my $last_line = $self->get_number_of_lines();
	$lines{ $last_line + 1 } = $line_to_add;#key is just the index or line number
	$self->set_tpf_lines( {%lines} );
}

sub add_line_to_beginning {
	my $self        = shift;
	my $line_to_add = shift;
	my %lines;
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
		foreach my $line_key (@reverse_sorted_line_numbers) {
			$lines{ $line_key + 1 } = $lines{$line_key};
		}
		$lines{ $reverse_sorted_line_numbers[-1] } = $line_to_add;
	}
	else {
		$lines{1} = $line_to_add;
	}
	$self->set_tpf_lines( {%lines} );
}

sub insert_line_before {
	my $self                         = shift;
	my $line_number_to_insert_before = shift;
	my $line_to_add                  = shift;
	my %lines;
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
		foreach my $line_key (@reverse_sorted_line_numbers) {
			if ( $line_key >= $line_number_to_insert_before ) {
				$lines{ $line_key + 1 } = $lines{$line_key};
			}
		}
		$lines{$line_number_to_insert_before} = $line_to_add;
	}
	else {
		$lines{1} = $line_to_add;
	}
	$self->set_tpf_lines( {%lines} );
}

sub insert_line_after {
	my $self                        = shift;
	my $line_number_to_insert_after = shift;
	my $line_to_add                 = shift;
	my %lines;
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
		foreach my $line_key (@reverse_sorted_line_numbers) {
			if ( $line_key > $line_number_to_insert_after ) {
				$lines{ $line_key + 1 } = $lines{$line_key};
			}
		}
		$lines{ $line_number_to_insert_after + 1 } = $line_to_add;
	}
	else {
		$lines{1} = $line_to_add;
	}
	$self->set_tpf_lines( {%lines} );
}

sub delete_line {
	my $self                  = shift;
	my $line_number_to_delete = shift;
	my %lines;
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
		foreach my $line_key (@sorted_line_numbers) {
			if ( $line_key > $line_number_to_delete ) {
				$lines{ $line_key - 1 } = $lines{$line_key};
			}
		}
		delete $lines{ $sorted_line_numbers[-1] };
	}
	else {

		#add error warning;
	}
	$self->set_tpf_lines( {%lines} );
}

sub get_number_of_lines {
	my $self = shift;
	my %lines;
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
		return $sorted_line_numbers[-1];
	}
	else {
		return 0;
	}
}

sub parse_tpf {
	my $self         = shift;
	my $input_string = shift;
	$self->clear_organism();
	$self->clear_assembly_name();
	$self->clear_chromosome();
	$self->clear_strain_haplotype_cultivar();
	$self->clear_type();
	$self->clear_assembly_version();
	$self->clear_comment();
	$self->clear_tpf_lines();

	#reading in input TPF into array, element=line
	my @lines               = split( /\n/, $input_string );
	my $tpf_data_has_begun  = 0;
	my $tpf_data_has_ended  = 0;
	my $current_line_number = 1;
	my $total_line_number   = scalar(@lines);
	foreach my $line (@lines) {
		print STDERR
		  "\rParsing TPF line $current_line_number of $total_line_number";
		$current_line_number++;
		chomp($line);
		if ( $line =~ m/^\s*$/ ) {    #skip blank lines
			next;
		}
		#now parsing TPF meta data
		if ( $line =~ m/^##/ )
		{    #identify comment lines and get value, split on key
			if ( $line =~ m/ORGANISM: /i ) {
				my @organism_line = split( /ORGANISM: /i, $line );
				my $organism =
				  $organism_line[-1]
				  ;  #get last element of array, could have used $#organism_line
				$organism =~ s/^\s+|\s+$//g;
				$self->set_organism($organism);
				next;
			}
			if ( $line =~ m/ASSEMBLY NAME: /i ) {
				my @assembly_name_line = split( /ASSEMBLY NAME: /i, $line );
				my $assembly_name = $assembly_name_line[-1];
				$assembly_name =~ s/^\s+|\s+$//g;
				$self->set_assembly_name($assembly_name);
				next;
			}
			if ( $line =~ m/CHROMOSOME: /i ) {
				my @chromosome_line = split( /CHROMOSOME: /i, $line );
				my $chromosome = $chromosome_line[-1];
				$chromosome =~ s/^\s+|\s+$//g;
				$self->set_chromosome($chromosome);
				next;
			}
			if ( $line =~ m/STRAIN\/HAPLOTYPE\/CULTIVAR: /i ) {
				my @strain_haplotype_cultivar_line =
				  split( /STRAIN\/HAPLOTYPE\/CULTIVAR: /i, $line );
				my $strain_haplotype_cultivar =
				  $strain_haplotype_cultivar_line[-1];
				$strain_haplotype_cultivar =~ s/^\s+|\s+$//g;
				$self->set_strain_haplotype_cultivar(
					$strain_haplotype_cultivar);
				next;
			}
			if ( $line =~ m/TYPE: /i ) {
				my @type_line = split( /TYPE: /i, $line );
				my $type = $type_line[-1];
				$type =~ s/^\s+|\s+$//g;
				$self->set_type($type);
				next;
			}
			if ( $line =~ m/VERSION: /i ) {
				my @assembly_version_line = split( /VERSION: /i, $line );
				my $assembly_version = $assembly_version_line[-1];
				$assembly_version =~ s/^\s+|\s+$//g;
				$self->set_assembly_version($assembly_version);
				next;
			}
			if ( $line =~ m/COMMENT: /i ) {
				my @comment_line = split( /COMMENT: /i, $line );
				my $comment = $comment_line[-1];
				$comment =~ s/^\s+|\s+$//g;
				$self->set_comment($comment);
				next;
			}
			if ( $line =~ m/=== Beginning of TPF Data ===/ ) {
				$tpf_data_has_begun = 1;
				next;
			}
			if ( $line =~ m/=== End of TPF Data ===/ ) {
				$tpf_data_has_ended = 1;
				next;
			}
		}

		#now parsing TPF data
		if ( ( $tpf_data_has_begun == 1 ) && ( $tpf_data_has_ended == 0 ) ) {
			my @tab_parsed_line = split( /\t/, $line );

			if ( !defined( $tab_parsed_line[0] ) ) {
				die "error in TPF line information\n";
			}

			$tab_parsed_line[0] =~ s/^\s+|\s+$//g;

			if ( $tab_parsed_line[0] eq "GAP" ) {
				my $tpf_gap_line = Bio::GenomeUpdate::TPF::TPFGapLine->new();#GAP line object
				if ( !defined( $tab_parsed_line[1] ) ) {

					#die with error missing gap line information
					die "Error in TPF line $line\n";
				}
				
				#get GAP data, get/set methods implemented by Moose/PBP
				$tab_parsed_line[1] =~ s/^\s+|\s+$//g;
				$tpf_gap_line->set_gap_type( $tab_parsed_line[1] );
				if ( defined( $tab_parsed_line[2] ) ) {
					$tab_parsed_line[2] =~ s/^\s+|\s+$//g;
					$tpf_gap_line->set_gap_size( $tab_parsed_line[2] );
				}
				if ( defined( $tab_parsed_line[3] ) ) {
					my @methods = split( /;/, $tab_parsed_line[3] );
					foreach my $method (@methods) {#can have multiple methods
						$method =~ s/^\s+|\s+$//g;
						$tpf_gap_line->add_gap_method($method);
					}
				}
				$self->add_line_to_end($tpf_gap_line);#add to TPF object
			}
			else {#For sequence line
				my $tpf_sequence_line =
				  Bio::GenomeUpdate::TPF::TPFSequenceLine->new();#Scaffold line object
				#Get scaffold data
				if ( defined( $tab_parsed_line[0] ) ) {
					$tab_parsed_line[0] =~ s/^\s+|\s+$//g;
					$tpf_sequence_line->set_accession( $tab_parsed_line[0] );
				}
				if ( defined( $tab_parsed_line[1] ) ) {
					$tab_parsed_line[1] =~ s/^\s+|\s+$//g;
					$tpf_sequence_line->set_clone_name( $tab_parsed_line[1] );
				}
				if ( defined( $tab_parsed_line[2] ) ) {
					$tab_parsed_line[2] =~ s/^\s+|\s+$//g;
					$tpf_sequence_line->set_local_contig_identifier(
						$tab_parsed_line[2] );
				}
				if ( defined( $tab_parsed_line[3] ) ) {
					$tab_parsed_line[3] =~ s/^\s+|\s+$//g;
					if (   $tab_parsed_line[3] eq 'PLUS'
						|| $tab_parsed_line[3] eq 'MINUS' )
					{
						$tpf_sequence_line->set_orientation(
							$tab_parsed_line[3] );
					}
					#if contig is contained within another accession, 
					#CONTAINED or CONTAINED_TURNOUT
					else { 
						$tpf_sequence_line->set_contains( $tab_parsed_line[3] );
						if ( defined( $tab_parsed_line[4] ) ) {
							$tab_parsed_line[4] =~ s/^\s+|\s+$//g;
							$tpf_sequence_line->set_containing_accession(
								$tab_parsed_line[4] );
						}
						if ( defined( $tab_parsed_line[5] ) ) {
							$tab_parsed_line[5] =~ s/^\s+|\s+$//g;
							$tpf_sequence_line->set_containing_clone_name(
								$tab_parsed_line[4] );
						}
					}
				}
				$self->add_line_to_end($tpf_sequence_line);
			}
			next;
		}
	}
	print STDERR "\nParsed $total_line_number TPF lines\n";
	return $self;
}

=item C<get_number_of_gap_lines ()>

Return nof gaps in TPF.

=cut

sub get_number_of_gap_lines {
	my $self = shift;
	my %lines;
	my $count=0;
	
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
		foreach my $line_key (@sorted_line_numbers) {
			if ( $lines{$line_key}->get_line_type() eq "gap" ) {
				$count++;
			}
		}
	}
	return $count;
}


=item C<get_number_of_sequence_lines ()>

Return nof sequences in TPF.

=cut

sub get_number_of_sequence_lines {
	my $self = shift;
	my %lines;
	my $count=0;
	
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
		foreach my $line_key (@sorted_line_numbers) {
			if ( $lines{$line_key}->get_line_type() eq "sequence" ) {
				$count++;
			}
		}
	}
	return $count;
}


=item C<get_gap_lengths ()>

Return array with all gap lengths in TPF.

=cut

sub get_gap_lengths{
	my $self = shift;
	my (@gap_lengths,%lines);
	
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
		foreach my $line_key (@sorted_line_numbers) {
			if ( $lines{$line_key}->get_line_type() eq "gap" ) {
				push @gap_lengths,$lines{$line_key}->get_gap_size();
			}
		}
	}
	
	return @gap_lengths;
}

sub get_formatted_tpf {
	my $self = shift;
	my %lines;
	my $out_str;

	#Print header info
	if ( $self->get_organism() ) {
		$out_str .= "##Organism: " . $self->get_organism() . "\n";
	}
	if ( $self->get_chromosome() ) {
		$out_str .= "##Chromosome: " . $self->get_chromosome() . "\n";
	}
	if ( $self->get_assembly_name() ) {
		$out_str .= "##Assembly Name: " . $self->get_assembly_name() . "\n";
	}
	if ( $self->get_strain_haplotype_cultivar() ) {
		$out_str .=
		  "##Strain/Haplotype/Cultivar: "
		  . $self->get_strain_haplotype_cultivar() . "\n";
	}
	if ( $self->get_type() ) {
		$out_str .= "##Type: " . $self->get_type() . "\n";
	}
	if ( $self->has_comment() ) {
		$out_str .= "##Comment: " . $self->get_comment() . "\n\n";
	}
	$out_str .= "\n##=== Beginning of TPF Data ===\n\n";
	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
		my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
		foreach my $line_key (@sorted_line_numbers) {
			if ( $lines{$line_key}->get_line_type() eq "sequence" ) {
				if ( $lines{$line_key}->has_accession() ) {
					$out_str .= $lines{$line_key}->get_accession() . "\t";
				}
				else {
					$out_str .= "??\t";
					print STDERR "accession not found\n";
				}
				if ( $lines{$line_key}->has_clone_name() ) {
					$out_str .= $lines{$line_key}->get_clone_name() . "\t";
				}
				else {
					$out_str .= "?\t";
				}
				if ( $lines{$line_key}->has_local_contig_identifier() ) {
					$out_str .=
					  $lines{$line_key}->get_local_contig_identifier() . "\t";
				}
				else {
					$out_str .= "?\t";
				}
				if ( $lines{$line_key}->has_contains() ) {
					$out_str .= $lines{$line_key}->get_contains() . "\t";
					if ( $lines{$line_key}->has_containing_accession() ) {
						$out_str .=
						  $lines{$line_key}->get_containing_accession() . "\t";
					}
					else {
						$out_str .= "?\t";
					}
					if ( $lines{$line_key}->has_containing_clone_name() ) {
						$out_str .=
						  $lines{$line_key}->get_containing_clone_name();
					}
					else {
						$out_str .= "?";
					}
				}
				else {
					if ( $lines{$line_key}->has_orientation() ) {
						$out_str .= $lines{$line_key}->get_orientation();
					}
				}
			}
			elsif ( $lines{$line_key}->get_line_type() eq "gap" ) {
				$out_str .= $lines{$line_key}->get_gap_identifier() . "\t";
				$out_str .= $lines{$line_key}->get_gap_type();
				if ( $lines{$line_key}->has_gap_size() ) {
					$out_str .= "\t" . $lines{$line_key}->get_gap_size();

					if (
						!(
							(
								$lines{$line_key}->get_gap_type() eq "CENTROMERE"
								|| $lines{$line_key}->get_gap_type() eq "TELOMERE"
								|| $lines{$line_key}->get_gap_type() eq "HETEROCHROMATIN"
								|| $lines{$line_key}->get_gap_type() eq "SHORT-ARM"
							)
						)
					  )
					{
						if ( $lines{$line_key}->has_gap_methods() ) {
							$out_str .= "\t";
							my $first_gmethod = 1;
							foreach my $gmethod (@{$lines{$line_key}->get_gap_methods()})
							{
								if ( $first_gmethod == 0 ) {
									$out_str .= ";";
								}
								$first_gmethod = 0;
								$out_str .= $gmethod;
							}
						}
						else {
							#add error warning method required
						}
					}
				}
				else {
					if (   $lines{$line_key}->get_gap_type() eq "CENTROMERE"
						|| $lines{$line_key}->get_gap_type() eq "TELOMERE"
						|| $lines{$line_key}->get_gap_type() eq "HETEROCHROMATIN"
						|| $lines{$line_key}->get_gap_type() eq "SHORT-ARM" )
					{

						#add error warning biological gap must have size
					}
				}

			}
			else {

				#add error warning;
			}
			$out_str .= "\n";
		}
	}
	else {

		#add error warning;
	}

	#$out_str .= footer
	$out_str .= "##=== End of TPF Data ===\n";
	return $out_str;
}

sub change_gap_size_between_scaffolds {
	my $self            = shift;
	my $first_scaffold  = shift;
	my $second_scaffold = shift;
	my $gap_size        = shift;
	my $gap_method      = shift;
	my %lines;

	#  if ($self->has_tpf_lines()) {
	%lines = %{ $self->get_tpf_lines() };

	#  }
	my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
	foreach my $line_key (@sorted_line_numbers) {

		#skip if on first or last line of TPF
		if (   ( $line_key == 1 )
			|| ( $line_key >= scalar(@sorted_line_numbers) ) )
		{
			next;
		}
		if ( $lines{$line_key}->get_line_type() eq "gap" ) {
			if (
				(
					$lines{ $line_key - 1 }->get_local_contig_identifier() eq
					$first_scaffold
				)
				&& ( $lines{ $line_key + 1 }->get_local_contig_identifier() eq
					$second_scaffold )
			  )
			{

				#print STDERR "FOUND  \n";
				if ($gap_size) {
					$lines{$line_key}->set_gap_size($gap_size);
				}
				if ($gap_method) {
					$lines{$line_key}->clear_gap_methods();
					$lines{$line_key}->add_gap_method($gap_method);
				}
			}
		}
	}

	#$self->clear_tpf_lines(%lines);
	$self->set_tpf_lines( {%lines} );
	return $self;
}

sub get_tpf_in_new_scaffold_order {
	my $self                           = shift;#getting calling TPF object
	my $a_ref                          = shift;#getting the order array 
	my @ordered_and_oriented_scaffolds = @$a_ref;
	my %lines;
	my $last_line_key;
	my $in_scaffold_range             = 0;
	my $out_tpf                       = Bio::GenomeUpdate::TPF->new();
	my $last_gap_line_number          = 0;
	my $last_sequence_line_number     = 0;
	my $is_first_scaffold             = 1;
	my $last_gap_line_number_in_range = "none";
	#get old TPF lines read in from the input file
	if ( $self->has_tpf_lines() ) {
		# returns hash of 
		# Bio::GenomeUpdate::TPF::TPFSequenceLine | Bio::GenomeUpdate::TPF::TPFGapLine
		%lines = %{ $self->get_tpf_lines() };
	}
	# Sorts on line numbers
	# key is just the index or line number
	my @forward_sorted_line_numbers = sort { $a <=> $b } keys %lines;
	my @reverse_sorted_line_numbers = sort { $b <=> $a } keys %lines;
	my @sorted_line_numbers;
	foreach my $scaffold_and_order_ref (@ordered_and_oriented_scaffolds) {

		#insert a TYPE 3 contig gap before the scaffold unless it is the first one.
		if ( $is_first_scaffold == 0 ) {
			my $tpf_gap_line = Bio::GenomeUpdate::TPF::TPFGapLine->new();
			$tpf_gap_line->set_gap_type("TYPE-3");
			$out_tpf->add_line_to_end($tpf_gap_line);
		}
		else {
			$is_first_scaffold = 0;# was init to 1
		}
		
		my @scaffold_and_order   = @$scaffold_and_order_ref;
		my $scaffold_to_move     = $scaffold_and_order[0];
		my $scaffold_orientation = $scaffold_and_order[1];
		
		if ( $scaffold_orientation eq "+" ) {
			@sorted_line_numbers = @forward_sorted_line_numbers;
		}
		elsif ( $scaffold_orientation eq "-" ) {
			@sorted_line_numbers = @reverse_sorted_line_numbers;
		}
		else {
			die "Orientation not specified\n";
		}
		
		foreach my $line_key (@sorted_line_numbers) {
			if ( $lines{$line_key}->get_line_type() eq "sequence" ) {
				if ( $lines{$line_key}->get_local_contig_identifier() eq
					$scaffold_to_move ){
					$in_scaffold_range = 1;
					unless ( $last_gap_line_number_in_range eq "none" ) {## ???????
						$out_tpf->add_line_to_end($lines{$last_gap_line_number_in_range});
					}
					$last_gap_line_number_in_range = "none";
					my $line_orientation = $lines{$line_key}->get_orientation();
					# -,PLUS = MINUS and -,MINUS = PLUS
					# if -, flip the orientation, nothing if +
					if ( $scaffold_orientation eq "-" ) {
						if ( $line_orientation eq "PLUS" ) {
							$lines{$line_key}->set_orientation("MINUS");
						}
						elsif ( $line_orientation eq "MINUS" ) {
							$lines{$line_key}->set_orientation("PLUS");
						}
						else {
							die "TPF sequence line missing PLUS or MINUS information\n";
						}
					}
					$out_tpf->add_line_to_end( $lines{$line_key} );

					$last_sequence_line_number = $line_key;
				}
				else {
					$in_scaffold_range             = 0;
					$last_gap_line_number_in_range = "none";
				}
			}
			elsif ( $lines{$line_key}->get_line_type() eq "gap" ) {
				if ( $in_scaffold_range == 1 ) {

					#$out_tpf->add_line_to_end($lines{$line_key});
					$last_gap_line_number          = $line_key;
					$last_gap_line_number_in_range = $line_key;
				}
			}
		}

#    if ((($scaffold_and_order[1] eq "+") && ($last_gap_line_number > $last_sequence_line_number)) ||
#(($scaffold_and_order[1] eq "-") && ($last_gap_line_number < $last_sequence_line_number)))  {
#      $out_tpf->delete_line($last_gap_line_number);
#    }
	}
	$self->set_tpf_lines( $out_tpf->get_tpf_lines() );
	return $self;
}

=item C<get_tpf_with_bacs_inserted_in_gaps ( @bacs, %scaffold_agp_coords )>

DEPRECATED. Does not account for BAC regions that start or end within gaps when resizing gaps. See issue #32. Returns a full TPF with the BAC accessions inserted in order that replace gaps. Components are now CONTAINED in BACs that encompass them. Using get_tpf_sp_tp_with_bacs_inserted_in_sequences_and_gaps() but retaining for historical reasons.

=cut

sub get_tpf_with_bacs_inserted_in_gaps {
#	my $self           = shift;
#	my $bacs_ref       = shift;
#	my $agp_coords_ref = shift;
#	my @bacs           = @$bacs_ref; # ref to array of arrays with bac names and coordinates
#	my %agp_coords     = %$agp_coords_ref;
#	
#	print STDERR "Please note this method is deprecated now. Please use get_tpf_sp_tp_with_bacs_inserted_in_sequences_and_gaps().\n";

#	#make sure BACs are sorted by position
#	foreach my $bac_ref (@bacs) {
#		my @bac      = @$bac_ref;
#		my $bac_name = $bac[0];
#		my $bac_start;
#		my $bac_end;
#		my $bac_to_insert = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
#		my %tpf_lines;
#		if ( $self->has_tpf_lines() ) {
#			%tpf_lines = %{ $self->get_tpf_lines() };
#		}
#		my @sorted_tpf_line_numbers = sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive
#		$bac_to_insert->set_accession($bac_name);
#		if ( $bac[1] < $bac[2] ) {
#			$bac_to_insert->set_orientation('PLUS'); #records the orientation of ref region that aligned to bac
#			$bac_start = $bac[1];
#			$bac_end   = $bac[2];
#		}
#		elsif ( $bac[1] > $bac[2] ) {
#			$bac_to_insert->set_orientation('MINUS'); #records the orientation of ref region that aligned to bac
#			$bac_start = $bac[2];
#			$bac_end   = $bac[1];
#		}
#		else {
#			die	"Error in BAC coordinates for BAC $bac_start Start: $bac_start End: $bac_end\n";
#		}
#		my $prev_agp_start = 0;
#		my $prev_agp_end   = 0;
#		my $prev_accession = 'none';
#		my $prev_line_key;
#		my $bac_is_contained = 0;
#		my $agp_start;
#		my $agp_end;
#		my $bac_is_inserted = 0;
#		my %gaps_to_resize;  #key will be line number and value will be new size
#		my @sorted_gaps_to_resize;
#		my %gaps_to_remove;    #key is line number value is undef
#		my @rev_sorted_gaps_to_remove;
#		my $insert_before_or_after = undef;
#		my $insert_line_number     = undef;
#		my %contained_contigs;    #key will be line number and value will be the contig accession
#		my $past_bac = 0;
#		my $line_key = 1;

#		#add BAC coordinates to AGP info (not saved/ output to STDOUT or file)
#		my %add_agp_coords;
#		$add_agp_coords{'start'} = $bac_start;
#		$add_agp_coords{'end'}   = $bac_end;
#		if ( $bac_to_insert->get_orientation() eq 'PLUS' ) {
#			$add_agp_coords{'orientation'} = '+';
#		}
#		elsif ( $bac_to_insert->get_orientation() eq 'MINUS' ) {
#			$add_agp_coords{'orientation'} = '-';
#		}
#		else {
#			die "No orientation specified for BAC: $bac_name\n";
#		}
#		$agp_coords{$bac_name} = \%add_agp_coords;
#		
#		#print STDERR "* sorted_tpf_line_numbers: ",@sorted_tpf_line_numbers + 1,"\n";

#		# the +1 breaks the code
#		#while ( $past_bac == 0 && $line_key <= @sorted_tpf_line_numbers + 1 ) {
#		while ( $past_bac == 0 && $line_key <= @sorted_tpf_line_numbers ) {
#			#print STDERR "** processing line $line_key\n";
#			if (!exists $tpf_lines{$line_key}){print STDERR "No TPF line for $line_key\n";}
#			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
#				my $accession = $tpf_lines{$line_key}->get_accession();
#				my $agp_line_coords_ref = $agp_coords{$accession};
#				my %line_coords         = %$agp_line_coords_ref;
#				$agp_start = $line_coords{'start'};
#				$agp_end   = $line_coords{'end'};

#				#check if past the BAC
#				if ( $bac_end < $prev_agp_start ) {
#					$past_bac = 1;
#				}

#				#check if current contig is contained in the BAC
#				if ( $agp_start >= $bac_start && $agp_end <= $bac_end ) {
#					$tpf_lines{$line_key}->set_contains('CONTAINED');
#					$tpf_lines{$line_key}->set_containing_accession($bac_name);
#					$self->set_tpf_lines( \%tpf_lines );

#					#$contained_contigs{$line_key}=$bac_name;
#				}

#				#check if current BAC is contained in the contig
#				if ( $bac_start >= $agp_start && $bac_end <= $agp_end ) {
#					$bac_to_insert->set_contains('CONTAINED');
#					$bac_to_insert->set_containing_accession($accession);
#				}

#				#check if gap is spanned by the BAC
#				if (   $prev_line_key
#					&& $bac_start <= $prev_agp_end
#					&& $bac_end >= $agp_start
#					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
#				{
#					$gaps_to_remove{ $line_key - 1 } = 'delete';
#					my $gap_location = $line_key - 1;
#					print STDERR "Removing gap at line $gap_location between $accession and $prev_accession\n";
#				}

#				#shrink gaps when partially spanned by a BAC
#				if (   $prev_line_key
#					&& $bac_start < $agp_start
#					&& $bac_start > $prev_agp_end
#					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
#				{
#					$gaps_to_resize{ $line_key - 1 } =
#					  $bac_start - $prev_agp_end;
#				}
#				if (   $prev_line_key
#					&& $bac_end < $agp_start
#					&& $bac_end > $prev_agp_end
#					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
#				{
#					$gaps_to_resize{ $line_key - 1 } = $agp_start - $bac_end;
#				}
#				$prev_line_key  = $line_key;
#				$prev_accession = $accession;
#				$prev_agp_start = $agp_start;
#				$prev_agp_end   = $agp_end;
#			}
#			$line_key++;
#		}

#		@sorted_gaps_to_resize     = sort { $a <=> $b } keys %gaps_to_resize;
#		@rev_sorted_gaps_to_remove = sort { $b <=> $a } keys %gaps_to_remove;

#		foreach my $line_number (@sorted_gaps_to_resize) {
#			$tpf_lines{$line_number}->set_gap_size( $gaps_to_resize{$line_number} );
#		}
#		$self->set_tpf_lines( \%tpf_lines );
#		foreach my $line_number (@rev_sorted_gaps_to_remove) {
#			$self->delete_line($line_number);
#		}
#		%tpf_lines               = %{ $self->get_tpf_lines() };
#		@sorted_tpf_line_numbers =
#		  sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive

#		$line_key = 1;
#		while ($bac_is_inserted == 0 && $line_key <= @sorted_tpf_line_numbers + 1 )
#		{
#			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
#				my $accession = $tpf_lines{$line_key}->get_accession();
#				my $agp_line_coords_ref = $agp_coords{$accession};
#				my %line_coords         = %$agp_line_coords_ref;
#				$agp_start = $line_coords{'start'};
#				$agp_end   = $line_coords{'end'};
#				if ( $line_key == 1 ) {          #deal with first one
#					if ( $bac_start <= 0 ) {
#						$insert_before_or_after = 'before';
#						$insert_line_number     = $line_key;
#						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier()
#						);
#						$bac_is_inserted = 1;
#					}
#				}
#				elsif ( $line_key == @sorted_tpf_line_numbers + 1 )
#				{                                #deal with last one
#					if ( $bac_start >= $agp_start ) {
#						$insert_before_or_after = 'after';
#						$insert_line_number     = $line_key;
#						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier()
#						);
#						$bac_is_inserted = 1;
#					}
#				}
#				elsif ($bac_start >= $prev_agp_start && $bac_start < $agp_start )
#				{
#					if ( $bac_start <= $prev_agp_end ) {
#						$insert_before_or_after = 'after';
#						$insert_line_number     = $prev_line_key;
#						$bac_to_insert->set_local_contig_identifier($tpf_lines{$prev_line_key}->get_local_contig_identifier() );
#						$bac_is_inserted = 1;
#					}
#					elsif ( $bac_start > $prev_agp_end ) {
#						$insert_before_or_after = 'before';
#						$insert_line_number     = $line_key;
#						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier()
#						);
#						$bac_is_inserted = 1;
#					}
#				}
#				$prev_line_key  = $line_key;
#				$prev_accession = $accession;
#				$prev_agp_start = $agp_start;
#				$prev_agp_end   = $agp_end;
#			}
#			$line_key++;
#		}

#		if ( $insert_before_or_after eq 'before' ) {
#			$self->insert_line_before( $insert_line_number, $bac_to_insert );
#		}
#		elsif ( $insert_before_or_after eq 'after' ) {
#			$self->insert_line_after( $insert_line_number, $bac_to_insert );
#		}
#		else {
#			die "BAC $bac_name not inserted\n";
#		}
#		%tpf_lines = %{ $self->get_tpf_lines() };
#	}
#	return $self;
}

=item C<get_tpf_with_bacs_inserted_in_sequences_and_gaps ( @bacs, %scaffold_agp_coords, %scaffold_component_contigs, %scaffold_component_contig_directions )>

DEPRECATED. Returns a full TPF with the BAC accessions inserted in order that replace gaps AND sequences. The sequence and gap components that are encompassed by a BAC are now deleted from the TPF. The assembled BACs start with ContigX in the group_coords.out file. These are substituted with the BACs as they are ordered in the ACE file. Each member BAC will have its own TPF line. Using get_tpf_sp_tp_with_bacs_inserted_in_sequences_and_gaps() but retaining for historical reasons.

=cut

sub get_tpf_with_bacs_inserted_in_sequences_and_gaps {
#	my $self           = shift;
#	my $bacs_ref       = shift;
#	my $agp_coords_ref = shift;
#	my $scaffold_component_contigs_ref = shift;
#	my $scaffold_component_contig_directions_ref = shift;
#	my @bacs           = @$bacs_ref; # ref to array of arrays with bac names and coordinates
#	my %agp_coords     = %$agp_coords_ref;
#	my %scaffold_component_contigs = %$scaffold_component_contigs_ref;
#	my %scaffold_component_contig_directions = %$scaffold_component_contig_directions_ref;
#	my %bac_inserted_accessions; 
#	my %sequence_accessions_to_remove; #key is accession and value is delete, can be undef
#	
#	print STDERR "Please note this method is deprecated now. Please use get_tpf_sp_tp_with_bacs_inserted_in_sequences_and_gaps().\n";

#	#make sure BACs are sorted by position
#	foreach my $bac_ref (@bacs) {
#		my @bac      = @$bac_ref;
#		my $bac_name = $bac[0];
#		my $bac_ref_start;
#		my $bac_ref_end;
#		my $bac_query_start;
#		my $bac_query_end;
#		my $bac_query_length;
#		
#		my $bac_to_insert = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
#		my %tpf_lines;
#		if ( $self->has_tpf_lines() ) {
#			%tpf_lines = %{ $self->get_tpf_lines() };
#		}
#		my @sorted_tpf_line_numbers = sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive

#		#set BAC variables
#		$bac_to_insert->set_accession($bac_name);
#		if ( $bac[1] < $bac[2] ) {
#			$bac_to_insert->set_orientation('PLUS'); #records the orientation of ref region that aligned to bac
#			$bac_ref_start = $bac[1];
#			$bac_ref_end   = $bac[2];
#		}
#		elsif ( $bac[1] > $bac[2] ) {#as mummer flips coords for alignments on MINUS strand
#			$bac_to_insert->set_orientation('MINUS'); #records the orientation of ref region that aligned to bac
#			$bac_ref_start = $bac[2];
#			$bac_ref_end   = $bac[1];
#		}
#		else {
#			die	"Error in BAC ref coordinates for BAC $bac_name Start: $bac_ref_start End: $bac_ref_end\n";
#		}
#		if ( $bac[3] < $bac[4] ) {
#			$bac_query_start = $bac[3];
#			$bac_query_end   = $bac[4];
#		}
#		elsif ( $bac[3] > $bac[4] ) {
#			$bac_query_start = $bac[4];
#			$bac_query_end   = $bac[3];
#		}
#		else {
#			die	"Error in BAC query coordinates for BAC $bac_name Start: $bac_query_start End: $bac_query_end\n";
#		}
#		$bac_query_length = $bac[5];
#		
#		my $prev_agp_sequence_start = 0;
#		my $prev_agp_sequence_end   = 0;
#		my $prev_accession = 'none';
#		my $prev_line_key;
#		my $bac_is_contained = 0;
#		my $agp_sequence_start;
#		my $agp_sequence_end;
#		my $bac_is_inserted = 0;
#		my %gaps_to_resize;  #key will be line number and value will be new size
#		my @sorted_gaps_to_resize;
#		my %sequences_and_gaps_to_remove; #key is line number and value is delete, can be undef
#		my @rev_sorted_sequences_and_gaps_to_remove;
#		my $insert_before_or_after = undef;
#		my $insert_line_number     = undef;
#		my %contained_contigs;    #key will be line number and value will be the contig accession
#		my $past_bac = 0;
#		my $line_key = 1; # starting from line 1 of TPF

#		#add BAC coordinates to AGP info (NOT saved or output to STDOUT or file, maybe write to inserted.scaffolds.agp??)
#		my %add_agp_coords;
#		$add_agp_coords{'start'} = $bac_ref_start;
#		$add_agp_coords{'end'}   = $bac_ref_end;
#		if ( $bac_to_insert->get_orientation() eq 'PLUS' ) {
#			$add_agp_coords{'orientation'} = '+';
#		}
#		elsif ( $bac_to_insert->get_orientation() eq 'MINUS' ) {
#			$add_agp_coords{'orientation'} = '-';
#		}
#		else {
#			die "No orientation specified for BAC: $bac_name\n";
#		}
#		$agp_coords{$bac_name} = \%add_agp_coords;
#		
#		#print STDERR "* sorted_tpf_line_numbers: ",@sorted_tpf_line_numbers + 1,"\n";
#		
#		#Parse through the TPF line by line and modify as required 
#		#while ( $past_bac == 0 && $line_key <= @sorted_tpf_line_numbers + 1 ) { # the +1 breaks the code
#		while ( $past_bac == 0 && $line_key <= @sorted_tpf_line_numbers ) {
#			#print STDERR "** processing line $line_key\n";
#			if ( !exists $tpf_lines{$line_key} ){print STDERR "No TPF line for $line_key\n";}
#			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
#				my $accession = $tpf_lines{$line_key}->get_accession();
#				
#				#skip BACs/accessions already inserted in TPF with no scaffold AGP records
#				if ( exists $bac_inserted_accessions{$accession} ){ $line_key++; next; }  
#				
#				#print STDERR "** processing accession $accession\n";
#				my $agp_line_coords_ref = $agp_coords{$accession};
#				my %line_coords         = %$agp_line_coords_ref;
#				$agp_sequence_start = $line_coords{'start'};
#				$agp_sequence_end   = $line_coords{'end'};

#				#check if past the BAC
#				if ( $bac_ref_end < $prev_agp_sequence_start ) {
#					$past_bac = 1;
#				}

#				#check if current contig is contained in the BAC
#				if ( $agp_sequence_start >= $bac_ref_start && $agp_sequence_end <= $bac_ref_end ) {
#					#$tpf_lines{$line_key}->set_contains('CONTAINED');
#					#$tpf_lines{$line_key}->set_containing_accession($bac_name);
#					#$self->set_tpf_lines( \%tpf_lines );
#					#$contained_contigs{$line_key}=$bac_name;
#					
#					$sequences_and_gaps_to_remove{ $line_key } = 'delete';
#					$sequence_accessions_to_remove{ $accession } = 'delete' ; # remember so that accession is not used in contained clause
#					#print STDERR "Added $accession to \%sequence_accessions_to_remove\n";
#					my $sequence_location = $line_key;
#					
#					print STDERR "Removing sequence at line $sequence_location for $accession\n"; 
#				}

#				#check if current BAC is contained in the contig
#				#happens a lot in SL2.50->SL3.0, final sequence will come from BAC but leads to conflicts in GRC alignment
#				if ( ($bac_ref_start >= $agp_sequence_start) && ($bac_ref_end <= $agp_sequence_end) ){
#					$bac_to_insert->set_contains('CONTAINED');
#					$bac_to_insert->set_containing_accession($accession);
#					print STDERR "Setting $bac_name contained in $accession\n";
#				}

#				#check if gap is spanned by the BAC
#				if (   $prev_line_key
#					&& $bac_ref_start <= $prev_agp_sequence_end
#					&& $bac_ref_end >= $agp_sequence_start
#					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
#				{
#					#$gaps_to_remove{ $line_key - 1 } = 'delete'; # subtract one as we want to remove the prev TPF line
#					$sequences_and_gaps_to_remove{ $line_key - 1 } = 'delete'; # subtract one as we want to remove the prev TPF line  
#					my $gap_location = $line_key - 1;
#					print STDERR "Removing gap at line $gap_location between $accession and $prev_accession w.r.t. original TPF\n";
#				}

#				#shrink gaps when partially spanned by a BAC
#				#no need to shrink sequence line as GRC aligner will figure out the switch over points
#				
#				#if BAC start is in the middle of a gap
#				# $bac_ref_start < $agp_sequence_start works as there is a offset of 1 even if mummer does not align at N 
#				if (   $prev_line_key
#					&& $bac_ref_start < $agp_sequence_start
#					&& $bac_ref_start > $prev_agp_sequence_end
#					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
#				{
#					if ( ($bac_query_end - $bac_query_start + 1) == $bac_query_length ){ #whole BAC aligned to ref
#						$gaps_to_resize{ $line_key - 1 } = $bac_ref_start - $prev_agp_sequence_end;
#						my $gap_location = $line_key - 1;
#						print STDERR "BAC $bac_name starts within a gap. Resizing gap at line $gap_location between $accession and $prev_accession to ".
#							$gaps_to_resize{$gap_location}." bp w.r.t. original TPF\n";
#					}
#					elsif( ($bac_query_end - $bac_query_start + 1) < $bac_query_length ){ #partial BAC aligned to ref
#						print STDERR "Partial alignment for $bac_name query start: $bac_query_start query end: $bac_query_end\n";
#						if ( $bac_query_start > 1 ){#fix overhang at begininng
#							my $bac_start_overhang =  $bac_query_start - 1;
#							
#							$gaps_to_resize{ $line_key - 1 } = $bac_ref_start - $bac_start_overhang - $prev_agp_sequence_end;
#							my $gap_location = $line_key - 1;
#							print STDERR "BAC $bac_name starts within a gap. Resizing gap at line $gap_location between $accession and $prev_accession to ".
#								$gaps_to_resize{$gap_location}." bp w.r.t. original TPF\n";
#						}
#					}
#					else{
#						print STDERR "Error for aligned region calculation for $bac_name query start: $bac_query_start query end: $bac_query_end\n";
#					}
#				}
#				
#				#if BAC end is in the middle of a gap
#				#$bac_ref_end < $agp_sequence_end works as there is a offset of 1 even if mummer does not align at N
#				if (   $prev_line_key
#					&& $bac_ref_end < $agp_sequence_start
#					&& $bac_ref_end > $prev_agp_sequence_end
#					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
#				{
#					if ( ($bac_query_end - $bac_query_start + 1) == $bac_query_length ){ #whole BAC aligned to ref
#						$gaps_to_resize{ $line_key - 1 } = $agp_sequence_start - $bac_ref_end;
#						my $gap_location = $line_key - 1;
#						print STDERR "BAC $bac_name ends within a gap. Resizing gap at line $gap_location between $accession and $prev_accession to ".
#							$gaps_to_resize{$gap_location}." bp w.r.t. original TPF\n";
#					}
#					elsif( ($bac_query_end - $bac_query_start + 1) < $bac_query_length ){ #partial BAC aligned to ref
#						print STDERR "Partial alignment for $bac_name query start: $bac_query_start query end: $bac_query_end\n";
#						my $bac_end_overhang = $bac_query_length - $bac_query_end + 1;
#						
#						$gaps_to_resize{ $line_key - 1 } = $agp_sequence_start - $bac_ref_end - $bac_end_overhang ;
#						my $gap_location = $line_key - 1;
#						print STDERR "BAC $bac_name ends within a gap. Resizing gap at line $gap_location between $accession and $prev_accession to ".
#							$gaps_to_resize{$gap_location}." bp w.r.t. original TPF\n";
#					}
#					else{
#						print STDERR "Error for aligned region calculation for $bac_name query start: $bac_query_start query end: $bac_query_end\n";
#					}
#				}
#				$prev_line_key  = $line_key;
#				$prev_accession = $accession;
#				$prev_agp_sequence_start = $agp_sequence_start;
#				$prev_agp_sequence_end   = $agp_sequence_end;
#			}
#			$line_key++;
#		}
#		
#		#editing lines already in TPF
#		#resizing gap lines first before order is distrupted
#		@sorted_gaps_to_resize     = sort { $a <=> $b } keys %gaps_to_resize;
#		foreach my $line_number (@sorted_gaps_to_resize) {
#			$tpf_lines{$line_number}->set_gap_size( $gaps_to_resize{$line_number} );
#		}
#		$self->set_tpf_lines( \%tpf_lines );

#		#DELETING sequence lines AND gaps
#		@rev_sorted_sequences_and_gaps_to_remove = sort { $b <=> $a } keys %sequences_and_gaps_to_remove;
#		foreach my $line_number (@rev_sorted_sequences_and_gaps_to_remove) {
#			$self->delete_line($line_number);
#		}
#		%tpf_lines = %{ $self->get_tpf_lines() };

#		@sorted_tpf_line_numbers = sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive
#		
#		#set the containing scaffold for the BAC, i.e. the local contig identifier
#		$line_key = 1;
#		while ($bac_is_inserted == 0 && $line_key <= @sorted_tpf_line_numbers + 1 )
#		{
#			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
#				my $accession = $tpf_lines{$line_key}->get_accession();
#				
#				#skip BACs/accessions already inserted in TPF with no scaffold AGP records
#				if ( exists $bac_inserted_accessions{$accession} ){ $line_key++; next; }
#				 
#				my $agp_line_coords_ref = $agp_coords{$accession};
#				my %line_coords         = %$agp_line_coords_ref;
#				$agp_sequence_start = $line_coords{'start'};
#				$agp_sequence_end   = $line_coords{'end'};
#				if ( $line_key == 1 ) {          #deal with first one
#					if ( $bac_ref_start <= 0 ) {
#						$insert_before_or_after = 'before';
#						$insert_line_number     = $line_key;
#						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier());
#						$bac_is_inserted = 1;
#					}
#				}
#				elsif ( $line_key == @sorted_tpf_line_numbers + 1 ) #deal with last one
#				{
#					if ( $bac_ref_start >= $agp_sequence_start ) {
#						$insert_before_or_after = 'after';
#						$insert_line_number     = $line_key;
#						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier()
#						);
#						$bac_is_inserted = 1;
#					}
#				}
#				elsif ($bac_ref_start >= $prev_agp_sequence_start
#					&& $bac_ref_start < $agp_sequence_start )
#				{
#					if ( $bac_ref_start <= $prev_agp_sequence_end ) {
#						$insert_before_or_after = 'after';
#						$insert_line_number     = $prev_line_key;
#						$bac_to_insert->set_local_contig_identifier($tpf_lines{$prev_line_key}->get_local_contig_identifier() );
#						$bac_is_inserted = 1;
#					}
#					elsif ( $bac_ref_start > $prev_agp_sequence_end ) {
#						$insert_before_or_after = 'before';
#						$insert_line_number     = $line_key;
#						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier()
#						);
#						$bac_is_inserted = 1;
#					}
#				}
#				$prev_line_key  = $line_key;
#				$prev_accession = $accession;
#				$prev_agp_sequence_start = $agp_sequence_start;
#				$prev_agp_sequence_end   = $agp_sequence_end;
#			}
#			$line_key++;
#		}
#		
#		#finally INSERTING the BAC TPF line
#		#print STDERR Dumper \%sequence_accessions_to_remove;
#		
#		if ($bac_name =~ /^Contig/ ){
#			print STDERR "Substituting in BACs for assembled contig $bac_name\n";
#			my $component_accessions_ref = $scaffold_component_contigs{$bac_name};
#			my @component_accessions_arr = @$component_accessions_ref;
#			my $component_accession_directions_ref = $scaffold_component_contig_directions{$bac_name};# orientation (+1,-1)
#			my @component_accession_directions_arr = @$component_accession_directions_ref; 
#			if (!(exists $scaffold_component_contigs{$bac_name}) || !(exists $scaffold_component_contig_directions{$bac_name})){
#				print STDERR "$bac_name not found in user supplied ACE file. Exiting ....\n\n"; exit 1;
#			}
#			
#			my $component_accessions_count = scalar @component_accessions_arr;
#			my $contig_bac_loop_counter;
#			
#			if ( $bac_to_insert->get_orientation() eq 'PLUS' ){
#				#simple order, same orientation
#				$contig_bac_loop_counter = 0;
#				while ($contig_bac_loop_counter < $component_accessions_count){
#					my $contig_bac_to_insert = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
#					print STDERR "******** inserting ";
#					print STDERR $component_accessions_arr[$contig_bac_loop_counter];
#					print STDERR "\n";
#					#print Dumper $component_accessions_arr[$contig_bac_loop_counter];
#					$contig_bac_to_insert->set_accession($component_accessions_arr[$contig_bac_loop_counter]);
#					$bac_inserted_accessions{$component_accessions_arr[$contig_bac_loop_counter]} = 'inserted'; #recording accession name
#					$contig_bac_to_insert->set_local_contig_identifier($bac_to_insert->get_local_contig_identifier() );
#					
#					#set contained, do NOT set contained if the accession was already deleted from TPF
#					if (($bac_to_insert->has_contains()) &&
#						(!exists $sequence_accessions_to_remove{$bac_to_insert->get_containing_accession()})){
#						$contig_bac_to_insert->set_contains('CONTAINED');
#						$contig_bac_to_insert->set_containing_accession($bac_to_insert->get_containing_accession());
#						print STDERR "Added containing clause for BAC ";
#						print STDERR $contig_bac_to_insert->get_accession();
#						print STDERR " with containing accession ";
#						print STDERR $contig_bac_to_insert->get_containing_accession();
#						print STDERR "\n";
#					}
#					
#					if ($component_accession_directions_arr[$contig_bac_loop_counter] == 1){
#						$contig_bac_to_insert->set_orientation('PLUS');
#					}
#					elsif ($component_accession_directions_arr[$contig_bac_loop_counter] == -1){
#						$contig_bac_to_insert->set_orientation('MINUS');	
#					}
#					
#					if ( $insert_before_or_after eq 'before' ){
#						$self->insert_line_before( $insert_line_number, $contig_bac_to_insert );	
#					}
#					elsif( $insert_before_or_after eq 'after' ){
#						$self->insert_line_after( $insert_line_number, $contig_bac_to_insert );	
#					}
#					print STDERR "Inserted BAC: ";
#					print STDERR $component_accessions_arr[$contig_bac_loop_counter];
#					print STDERR " for assembled contig $bac_name\n";
#					$contig_bac_loop_counter++;
#				}
#			}
#			elsif ( $bac_to_insert->get_orientation() eq 'MINUS' ){
#				#reverse order, flip orientation
#				$contig_bac_loop_counter = $component_accessions_count - 1 ;
#				while ($contig_bac_loop_counter >= 0 ){
#					my $contig_bac_to_insert = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
#					$contig_bac_to_insert->set_accession($component_accessions_arr[$contig_bac_loop_counter]);
#					$bac_inserted_accessions{$component_accessions_arr[$contig_bac_loop_counter]} = 'inserted'; #recording accession name
#					$contig_bac_to_insert->set_local_contig_identifier($bac_to_insert->get_local_contig_identifier() );
#					
#					#set contained, do NOT set contained if the accession was already deleted from TPF
#					if (($bac_to_insert->has_contains()) &&
#						(!exists $sequence_accessions_to_remove{$bac_to_insert->get_containing_accession()})){
#						$contig_bac_to_insert->set_contains('CONTAINED');
#						$contig_bac_to_insert->set_containing_accession($bac_to_insert->get_containing_accession());
#						
#						print STDERR "Added containing clause for BAC ";
#						print STDERR $contig_bac_to_insert->get_accession();
#						print STDERR " with containing accession ";
#						print STDERR $contig_bac_to_insert->get_containing_accession();
#						print STDERR "\n";
#					}
#					
#					if ($component_accession_directions_arr[$contig_bac_loop_counter] == 1){
#						$contig_bac_to_insert->set_orientation('MINUS');#flip
#					}
#					elsif ($component_accession_directions_arr[$contig_bac_loop_counter] == -1){
#						$contig_bac_to_insert->set_orientation('PLUS');#flip	
#					}
#					
#					if ( $insert_before_or_after eq 'before' ){
#						$self->insert_line_before( $insert_line_number, $contig_bac_to_insert );	
#					}
#					elsif( $insert_before_or_after eq 'after' ){
#						$self->insert_line_after( $insert_line_number, $contig_bac_to_insert );	
#					} 
#					print STDERR "Inserted BAC: ";
#					print STDERR $component_accessions_arr[$contig_bac_loop_counter];
#					print STDERR " for assembled contig $bac_name\n";
#					$contig_bac_loop_counter--;
#				}				
#			}
#			else{
#				die "BAC $bac_name not inserted. Exiting ....\n";
#			}
#		}
#		else{#if its a singleton
#			$bac_inserted_accessions{$bac_to_insert->get_accession()} = 'inserted'; #recording accession name
#			if ( $insert_before_or_after eq 'before' ) {
#				$self->insert_line_before( $insert_line_number, $bac_to_insert );
#			}
#			elsif ( $insert_before_or_after eq 'after' ) {
#				$self->insert_line_after( $insert_line_number, $bac_to_insert );
#			}
#			else {
#				die "BAC $bac_name not inserted\n";
#			}
#			print STDERR "Inserted BAC: $bac_name\n";
#		}
#		%tpf_lines = %{ $self->get_tpf_lines() };
#		#print STDERR $self->get_formatted_tpf();
#	}
#	
#	#remove contained if the accession (typically WGS from prev assembly) was deleted from TPF
#	my %lines = %{ $self->get_tpf_lines() };
#	my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
#	foreach my $line_key (@sorted_line_numbers) {
#		if ( $lines{$line_key}->get_line_type() eq "sequence" ) {
#			if (($lines{$line_key}->has_contains()) &&
#				(exists $sequence_accessions_to_remove{$lines{$line_key}->get_containing_accession()})){
#				my $accession_to_delete = $lines{$line_key}->get_containing_accession();
#				$lines{$line_key}->clear_contains();
#				$lines{$line_key}->clear_containing_accession();
#				print STDERR "Removed containing accession $accession_to_delete for BAC ";
#				print STDERR $lines{$line_key}->get_accession();
#				print STDERR "\n";
#			}
#		}
#	}
#	return $self;
}

=item C<get_tpf_sp_tp_with_bacs_inserted_in_sequences_and_gaps ( $chr, $switch_points, $trim_points, @bacs, %scaffold_agp_coords, %scaffold_component_contigs, %scaffold_component_contig_directions )>

Returns SP, TP and full TPF objects with BAC accessions inserted in order that replace gaps AND sequences. The sequence and gap components that are encompassed by a BAC are now deleted from the TPF. The assembled BACs start with ContigX in the group_coords.out file. These are substituted with the BACs as they are ordered in the ACE file. Each member BAC will have its own TPF line. You may need to remove redundant BACs as they confuse the GRC end-to-end aligner. Switch point lines are created and added to switch point object but no trim point lines are added as we don't use them right now.

=cut

sub get_tpf_sp_tp_with_bacs_inserted_in_sequences_and_gaps {
	my $self           = shift;
	my $chromosome     = shift;
	my $switch_points  = shift;
	my $trim_points    = shift;
	my $bacs_ref       = shift;
	my $scaffold_agp_coords_ref = shift;
	my $scaffold_component_contigs_ref = shift;
	my $scaffold_component_contig_directions_ref = shift;
	my @bacs           = @$bacs_ref;   # ref to array of arrays with bac names and coordinates
	my %scaffold_agp_coords     = %$scaffold_agp_coords_ref;
	my %scaffold_component_contigs = %$scaffold_component_contigs_ref;
	my %scaffold_component_contig_directions = %$scaffold_component_contig_directions_ref;
	my %bac_inserted_accessions; 
	my %sequence_accessions_to_remove; #key is accession and value is delete, can be undef
	my %before_insertion_counter;      #track insertions before a TPF line to calculate offsets for subsequent insertions
	my %after_insertion_counter;       #track insertions after a TPF line to calculate offsets for subsequent insertions

	
=item C<_get_accession_coordinate ($accession, $chromosome_coordinate)> 
	
Returns the local coordinate of the accession. Added for use in switchover and trim files
	
=cut
	
	local *_get_accession_coordinate = sub {
		my $accession = shift;
		my $chromosome_coordinate = shift;
		
		my $accession_start = $scaffold_agp_coords{$accession}->{start};
		
		return $chromosome_coordinate - $accession_start;
	};

=item C<_get_accession_length ($accession)> 
	
Returns the length of the accession. Added for use in switchover and trim files
	
=cut
	
	local *_get_accession_length = sub {
		my $accession = shift;
		
		return ( $scaffold_agp_coords{$accession}->{end} - $scaffold_agp_coords{$accession}->{start} + 1 );
	};
	
	#make sure BACs are sorted by position
	foreach my $bac_ref (@bacs) {
		my @bac      = @$bac_ref;
		my $bac_name = $bac[0];
		my $bac_ref_start;
		my $bac_ref_end;
		my $bac_query_start;
		my $bac_query_end;
		my $bac_query_length;
		my $ref_orientation_groupcoords;
		my $qry_orientation_groupcoords;
		my $direction;
		
		my $bac_to_insert = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
		my %tpf_lines;
		if ( $self->has_tpf_lines() ) {
			%tpf_lines = %{ $self->get_tpf_lines() };
		}
		my @sorted_tpf_line_numbers = sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive

		#set BAC variables
		$bac_to_insert->set_accession($bac_name);
		if ( $bac[1] < $bac[2] ) {
			$ref_orientation_groupcoords = 'PLUS';  #records the orientation of ref region that aligned to bac
			#$bac_to_insert->set_orientation('PLUS'); #records the orientation of ref region that aligned to bac
			$bac_ref_start = $bac[1];
			$bac_ref_end   = $bac[2];
		}
		elsif ( $bac[1] > $bac[2] ) {#as mummer flips coords for alignments on MINUS strand
			$ref_orientation_groupcoords = 'MINUS';  #records the orientation of ref region that aligned to bac
			#$bac_to_insert->set_orientation('MINUS'); #records the orientation of ref region that aligned to bac
			$bac_ref_start = $bac[2];
			$bac_ref_end   = $bac[1];
		}
		else {
			die	"Error in BAC ref coordinates for BAC $bac_name Start: $bac_ref_start End: $bac_ref_end\n";
		}
		if ( $bac[3] < $bac[4] ) {#query alignment on positive strand
			$bac_query_start        = $bac[3];
			$bac_query_end          = $bac[4];
			$qry_orientation_groupcoords = 'PLUS';
		}
		elsif ( $bac[3] > $bac[4] ) {#query alignment on negative strand
			$bac_query_start = $bac[4];
			$bac_query_end   = $bac[3];
			$qry_orientation_groupcoords = 'MINUS';
		}
		else {
			die	"Error in BAC query coordinates for BAC $bac_name Start: $bac_query_start End: $bac_query_end\n";
		}
		$bac_query_length = $bac[5];
		$direction        = $bac[6]; #direction (+1 if in ref and query align in same orientation, -1 otherwise)
		
		#init vars for placing the BACs in the correct location in TPF file
		my $prev_agp_sequence_start = 0;
		my $prev_agp_sequence_end   = 0;
		my $prev_accession = 'none';
		my $prev_line_key;
		my $bac_is_contained = 0;
		my $agp_sequence_start;
		my $agp_sequence_end;
		my $bac_is_inserted = 0;
		my %gaps_to_resize;  #key will be line number and value will be new size
		my @sorted_gaps_to_resize;
		my %sequences_and_gaps_to_remove; #key is line number and value is delete, can be undef
		my $remove = undef;
		my @rev_sorted_sequences_and_gaps_to_remove;
		my $insert_before_or_after = undef;
		my $insert_line_number     = undef;
		my $insert_line_accession  = undef;
		my $insert_line_strand     = undef;
		my %contained_contigs;    #key will be line number and value will be the contig accession
		my $past_bac = 0;
		my $line_key = 1; # starting from line 1 of TPF

		#add BAC coordinates to AGP info (NOT saved or output to STDOUT or file, maybe write to inserted.scaffolds.agp??)
		my %add_scaffold_agp_coords;
		$add_scaffold_agp_coords{'start'} = $bac_ref_start;
		$add_scaffold_agp_coords{'end'}   = $bac_ref_end;
		if ( $ref_orientation_groupcoords eq 'PLUS' ) {
			$add_scaffold_agp_coords{'orientation'} = '+';
		}
		elsif ( $ref_orientation_groupcoords eq 'MINUS' ) {
			$add_scaffold_agp_coords{'orientation'} = '-';
		}
		else {
			die "No orientation specified for ref for BAC: $bac_name\n";
		}
		$scaffold_agp_coords{$bac_name} = \%add_scaffold_agp_coords;
		
		###############
		#Parse through the TPF line by line and RECORD modifications as required
		###############
		while ( $past_bac == 0 && $line_key <= @sorted_tpf_line_numbers ) {
			#print STDERR "** processing line $line_key\n";
			if ( !exists $tpf_lines{$line_key} ){print STDERR "No TPF line for $line_key\n";}
			
			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
				my $accession = $tpf_lines{$line_key}->get_accession();
				
				#reset remove flag
				$remove = 0;
				
				#skip BACs/accessions already inserted in TPF with no scaffold AGP records
				if ( exists $bac_inserted_accessions{$accession} ){ $line_key++; next; }  
				
				#print STDERR "** processing accession $accession\n";
				my $agp_line_coords_ref = $scaffold_agp_coords{$accession};
				my %line_coords         = %$agp_line_coords_ref;
				$agp_sequence_start = $line_coords{'start'};
				$agp_sequence_end   = $line_coords{'end'};

				#check if past the BAC
				if ( $bac_ref_end < $prev_agp_sequence_start ) {
					$past_bac = 1;
				}

				#check if current contig is contained in the BAC, delete if yes
				if ( $agp_sequence_start >= $bac_ref_start && $agp_sequence_end <= $bac_ref_end ) {
					
					$sequences_and_gaps_to_remove{ $line_key } = 'delete';
					$sequence_accessions_to_remove{ $accession } = 'delete' ; # remember so that accession is not used in contained clause
					#print STDERR "Added $accession to \%sequence_accessions_to_remove\n";
					my $sequence_location = $line_key;
					
					print STDERR "Removing sequence at line $sequence_location for $accession\n";
					$remove = 1 ; #set flag
				}

				#if current BAC is contained in the contig, set CONTAINED and create switch points
				#happens a lot in SL2.50->SL3.0, final sequence will come from BAC but leads to conflicts in GRC alignment
				if ( 	$bac_ref_start >= $agp_sequence_start
					&& $bac_ref_end <= $agp_sequence_end
					&& !$remove ){
					$bac_to_insert->set_contains('CONTAINED');
					$bac_to_insert->set_containing_accession($accession);
					print STDERR "Setting $bac_name contained in $accession\n";
					
					#if BAC is not flush with either end of the WGS contig then
					#set switch point line for transition between BAC and WGS contig
					my ($accession_prefix_orientation, $accession_suffix_orientation, $accession_prefix_last_base, $accession_suffix_first_base);
					if ($tpf_lines{$line_key}->get_orientation() eq 'PLUS'){
						$accession_prefix_orientation = '+';
					}
					elsif($tpf_lines{$line_key}->get_orientation() eq 'MINUS'){
						$accession_prefix_orientation = '-';
					}
					else{
						die "No orientation for $accession. Exiting.. \n";
					}
					if ($ref_orientation_groupcoords eq 'PLUS'){
						$accession_suffix_orientation = '+';
					}
					elsif($ref_orientation_groupcoords eq 'MINUS'){
						$accession_suffix_orientation = '-';
					}
					else{
						die "No orientation for $bac_name. Exiting.. \n";
					}
					
					$accession_prefix_last_base = _get_accession_coordinate ($accession, $bac_ref_start - 1 );
					print STDERR "$bac_name does not align from base 1 on accession $accession. Might be an error\n" if $bac_query_start != 1;
					$accession_suffix_first_base = $bac_query_start;

					print STDERR "$bac_name aligns from base 1 on accession $accession. Not creating switch point\n" if $accession_prefix_last_base == -1;
					if ($accession_prefix_last_base != -1){
						my $sp_prefix_line = Bio::GenomeUpdate::SP::SPLine->new( 
								chromosome => $chromosome,
								accession_prefix => $accession,
								accession_suffix => $bac_name,
								accession_prefix_orientation => $accession_prefix_orientation,
								accession_suffix_orientation => $accession_suffix_orientation,
								accession_prefix_last_base => $accession_prefix_last_base,
								accession_suffix_first_base => $accession_suffix_first_base,
								comment => "BAC $bac_name is contained within WGS contig $accession from previous version. Designates switch point from WGS contig to BAC."
							);
						print STDERR "Adding switch point for transition from $accession to contained $bac_name\n";
						$switch_points->add_line_to_end($sp_prefix_line);					
					}
					
					
					my $temp_orientation = $accession_prefix_orientation;
					$accession_prefix_orientation = $accession_suffix_orientation;
					$accession_suffix_orientation = $temp_orientation;
					
					print STDERR "$bac_name does not align till its end on accession $accession. Might be an error\n" if $bac_query_end != $bac_query_length;
					$accession_prefix_last_base = $bac_query_end;
					$accession_suffix_first_base = _get_accession_coordinate ($accession, $bac_ref_end + 1 );
					
					print STDERR "$bac_name aligns till end on accession $accession. Not creating a switch point\n" if $accession_suffix_first_base > _get_accession_length ($accession);
					if ($accession_suffix_first_base > _get_accession_length ($accession)){
						my $sp_suffix_line = Bio::GenomeUpdate::SP::SPLine->new( 
									chromosome => $chromosome,
									accession_prefix => $bac_name,
									accession_suffix => $accession,
									accession_prefix_orientation => $accession_prefix_orientation,
									accession_suffix_orientation => $accession_suffix_orientation,
									accession_prefix_last_base => $accession_prefix_last_base,
									accession_suffix_first_base => $accession_suffix_first_base,
									comment => "BAC $bac_name is contained within WGS contig $accession from previous version. Designates switch point to from BAC to WGS contig."
								);
						print STDERR "Adding switch point for transition from contained $bac_name to $accession\n";
						$switch_points->add_line_to_end($sp_suffix_line);
					}
				}
				
				#no need to shrink sequence line as GRC aligner will figure out the switch over points
				#set switch point in case BAC has partial overlap with WGS contig on 5' end
				if (	!$remove
					&& $bac_ref_start <= $agp_sequence_start
					&& $bac_ref_end   >  $agp_sequence_start
					&& $bac_ref_end   <  $agp_sequence_end ){ #for all lines incl first line in TPF
				
					#set switch point line for transition between BAC and WGS contig
					my ($accession_prefix_orientation, $accession_suffix_orientation, $accession_prefix_last_base, $accession_suffix_first_base);
					if ($ref_orientation_groupcoords eq 'PLUS'){
						$accession_prefix_orientation = '+';
					}
					elsif($ref_orientation_groupcoords eq 'MINUS'){
						$accession_prefix_orientation = '-';
					}
					else{
						die "No orientation for $bac_name. Exiting.. \n";
					}
					
					if ($tpf_lines{$line_key}->get_orientation() eq 'PLUS'){
						$accession_suffix_orientation = '+';
					}
					elsif($tpf_lines{$line_key}->get_orientation() eq 'MINUS'){
						$accession_suffix_orientation = '-';
					}
					else{
						die "No orientation for $accession. Exiting.. \n";
					}

					print STDERR "$bac_name does not align till its end on accession $accession. Might be an error\n" if $bac_query_end != $bac_query_length;
					$accession_prefix_last_base = $bac_query_end;
					$accession_suffix_first_base = _get_accession_coordinate ($accession, $bac_ref_end + 1 );
					
					my $sp_5prime_line = Bio::GenomeUpdate::SP::SPLine->new( 
								chromosome => $chromosome,
								accession_prefix => $bac_name,
								accession_suffix => $accession,
								accession_prefix_orientation => $accession_prefix_orientation,
								accession_suffix_orientation => $accession_suffix_orientation,
								accession_prefix_last_base => $accession_prefix_last_base,
								accession_suffix_first_base => $accession_suffix_first_base,
								comment => "BAC $bac_name aligns to the 5' end of WGS contig $accession from previous version. Designates switch point from BAC to WGS contig."
							);
					print STDERR "Adding switch point for transition from 5' aligned $bac_name to $accession\n";
					$switch_points->add_line_to_end($sp_5prime_line);
				}
				#set switch points in case BAC has partial overlap with WGS contig on 3' end
				#should handle last sequence line in TPF
				elsif ( $prev_line_key
					&& !$remove
					&& $bac_ref_start >= $agp_sequence_start
					&& $bac_ref_start <  $agp_sequence_end
					&& $bac_ref_end   >  $agp_sequence_end){

					#set switch point line for transition between BAC and WGS contig
					my ($accession_prefix_orientation, $accession_suffix_orientation, $accession_prefix_last_base, $accession_suffix_first_base);
					if ($tpf_lines{$line_key}->get_orientation() eq 'PLUS'){
						$accession_prefix_orientation = '+';
					}
					elsif($tpf_lines{$line_key}->get_orientation() eq 'MINUS'){
						$accession_prefix_orientation = '-';
					}
					else{
						die "No orientation for $accession. Exiting.. \n";
					}
					if ($ref_orientation_groupcoords eq 'PLUS'){
						$accession_suffix_orientation = '+';
					}
					elsif($ref_orientation_groupcoords eq 'MINUS'){
						$accession_suffix_orientation = '-';
					}
					else{
						die "No orientation for $bac_name. Exiting.. \n";
					}

					print STDERR "$bac_name does not align from base 1 on accession $accession. Might be an error\n" if $bac_query_start != 1;
					$accession_prefix_last_base = _get_accession_coordinate ($accession, $bac_ref_start - 1 );
					$accession_suffix_first_base = $bac_query_start;
					
					my $sp_3prime_line = Bio::GenomeUpdate::SP::SPLine->new( 
								chromosome => $chromosome,
								accession_prefix => $accession,
								accession_suffix => $bac_name,
								accession_prefix_orientation => $accession_prefix_orientation,
								accession_suffix_orientation => $accession_suffix_orientation,
								accession_prefix_last_base => $accession_prefix_last_base,
								accession_suffix_first_base => $accession_suffix_first_base,
								comment => "BAC $bac_name aligns to the 3' end of WGS contig $accession from previous version. Designates switch point from WGS contig to BAC."
							);
					print STDERR "Adding switch point for transition from $accession to 3' aligned $bac_name\n";
					$switch_points->add_line_to_end($sp_3prime_line);					
				}
				
				###############
				### processing GAP lines
				###############
				
				#check if gap is spanned by the BAC
				if (   $prev_line_key
					&& $bac_ref_start <= $prev_agp_sequence_end
					&& $bac_ref_end >= $agp_sequence_start
					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
				{
					$sequences_and_gaps_to_remove{ $line_key - 1 } = 'delete'; # subtract one as we want to remove the prev TPF line  
					my $gap_location = $line_key - 1;
					print STDERR "Removing gap at line $gap_location between $accession and $prev_accession w.r.t. original TPF\n";
				}

				#SHRINK gaps when partially spanned by a BAC
				#if BAC start is in the middle of a gap
				# $bac_ref_start < $agp_sequence_start works as there is a offset of 1 even if mummer does not align at N 
				if (   $prev_line_key
					&& $bac_ref_start < $agp_sequence_start
					&& $bac_ref_start > $prev_agp_sequence_end
					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
				{
					if ( ($bac_query_end - $bac_query_start + 1) == $bac_query_length ){ #whole BAC aligned to ref
						$gaps_to_resize{ $line_key - 1 } = $bac_ref_start - $prev_agp_sequence_end;
						my $gap_location = $line_key - 1;
						print STDERR "BAC $bac_name starts within a gap. Resizing gap at line $gap_location between $accession and $prev_accession to ".
							$gaps_to_resize{$gap_location}." bp w.r.t. original TPF\n";
					}
					elsif( ($bac_query_end - $bac_query_start + 1) < $bac_query_length ){ #partial BAC aligned to ref
						print STDERR "Partial alignment for $bac_name query start: $bac_query_start query end: $bac_query_end\n";
						if ( $bac_query_start > 1 ){#fix overhang at begininng
							my $bac_start_overhang =  $bac_query_start - 1;
							
							$gaps_to_resize{ $line_key - 1 } = $bac_ref_start - $bac_start_overhang - $prev_agp_sequence_end;
							my $gap_location = $line_key - 1;
							print STDERR "BAC $bac_name starts within a gap. Resizing gap at line $gap_location between $accession and $prev_accession to ".
								$gaps_to_resize{$gap_location}." bp w.r.t. original TPF\n";
						}
					}
					else{
						print STDERR "Error for aligned region calculation for $bac_name query start: $bac_query_start query end: $bac_query_end\n";
					}
				}
				
				#if BAC end is in the middle of a gap
				#$bac_ref_end < $agp_sequence_end works as there is a offset of 1 even if mummer does not align at N
				if (   $prev_line_key
					&& $bac_ref_end < $agp_sequence_start
					&& $bac_ref_end > $prev_agp_sequence_end
					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
				{
					if ( ($bac_query_end - $bac_query_start + 1) == $bac_query_length ){ #whole BAC aligned to ref
						$gaps_to_resize{ $line_key - 1 } = $agp_sequence_start - $bac_ref_end;
						my $gap_location = $line_key - 1;
						print STDERR "BAC $bac_name ends within a gap. Resizing gap at line $gap_location between $accession and $prev_accession to ".
							$gaps_to_resize{$gap_location}." bp w.r.t. original TPF\n";
					}
					elsif( ($bac_query_end - $bac_query_start + 1) < $bac_query_length ){ #partial BAC aligned to ref
						print STDERR "Partial alignment for $bac_name query start: $bac_query_start query end: $bac_query_end\n";
						my $bac_end_overhang = $bac_query_length - $bac_query_end + 1;
						
						$gaps_to_resize{ $line_key - 1 } = $agp_sequence_start - $bac_ref_end - $bac_end_overhang ;
						my $gap_location = $line_key - 1;
						print STDERR "BAC $bac_name ends within a gap. Resizing gap at line $gap_location between $accession and $prev_accession to ".
							$gaps_to_resize{$gap_location}." bp w.r.t. original TPF\n";
					}
					else{
						print STDERR "Error for aligned region calculation for $bac_name query start: $bac_query_start query end: $bac_query_end\n";
					}
				}
				$prev_line_key  = $line_key;
				$prev_accession = $accession;
				$prev_agp_sequence_start = $agp_sequence_start;
				$prev_agp_sequence_end   = $agp_sequence_end;
			}
			$line_key++;
		}
		
		
		###############
		#EDITING lines already in TPF
		###############
		
		#resizing gap lines first before order is distrupted by deletion
		@sorted_gaps_to_resize     = sort { $a <=> $b } keys %gaps_to_resize;
		foreach my $line_number (@sorted_gaps_to_resize) {
			$tpf_lines{$line_number}->set_gap_size( $gaps_to_resize{$line_number} );
		}
		$self->set_tpf_lines( \%tpf_lines );

		#DELETING sequence lines AND gaps, line number logic needed for gaps as no unique ID for gap lines
		@rev_sorted_sequences_and_gaps_to_remove = sort { $b <=> $a } keys %sequences_and_gaps_to_remove;
		foreach my $line_number (@rev_sorted_sequences_and_gaps_to_remove) {
			$self->delete_line($line_number);
		}
		%tpf_lines = %{ $self->get_tpf_lines() };

		@sorted_tpf_line_numbers = sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive
		
		#set the containing scaffold for the BAC, i.e. the local contig identifier
		#get the TPF line number for insertion before/after
		$line_key = 1;
		while ($bac_is_inserted == 0 && $line_key <= @sorted_tpf_line_numbers + 1 )
		{
			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
				my $accession = $tpf_lines{$line_key}->get_accession();
				
				#skip BACs/accessions already inserted in TPF with no scaffold AGP records
				if ( exists $bac_inserted_accessions{$accession} ){ $line_key++; next; }
				 
				my $agp_line_coords_ref = $scaffold_agp_coords{$accession};
				my %line_coords         = %$agp_line_coords_ref;
				$agp_sequence_start = $line_coords{'start'};
				$agp_sequence_end   = $line_coords{'end'};
				if ( $line_key == 1 ) {          #deal with first one
					if ( $bac_ref_start <= 0 ) {
						$insert_before_or_after = 'before';
						$insert_line_number     = $line_key;
						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier());
						$bac_is_inserted = 1;
					}
				}
				elsif ( $line_key == @sorted_tpf_line_numbers + 1 ) #deal with last one
				{
					if ( $bac_ref_start >= $agp_sequence_start ) {
						$insert_before_or_after = 'after';
						$insert_line_number     = $line_key;
						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier());
						$bac_is_inserted = 1;
					}
				}
				elsif ($bac_ref_start >= $prev_agp_sequence_start
					&& $bac_ref_start < $agp_sequence_start ) #for all other cases
				{
					if ( $bac_ref_start <= $prev_agp_sequence_end ) {
						$insert_before_or_after = 'after';        # insert after prev seq line if BAC starts before prev seq end coordinate
						$insert_line_number     = $prev_line_key; # set line number to prev line
						$bac_to_insert->set_local_contig_identifier($tpf_lines{$prev_line_key}->get_local_contig_identifier());
						$bac_is_inserted = 1;
					}
					elsif ( $bac_ref_start > $prev_agp_sequence_end ) {
						$insert_before_or_after = 'before';  # insert before current seq line if BAC starts before current seq start coordinate
						$insert_line_number     = $line_key; # set line number to current seq
						$bac_to_insert->set_local_contig_identifier($tpf_lines{$line_key}->get_local_contig_identifier());
						$bac_is_inserted = 1;
					}
				}
				$prev_line_key  = $line_key;
				$prev_accession = $accession;
				$prev_agp_sequence_start = $agp_sequence_start;
				$prev_agp_sequence_end   = $agp_sequence_end;
			}
			$line_key++;
		}
		
		
		###############
		#finally INSERTING the BAC TPF line using the hash and array data structures constructed from $self->get_tpf_lines()
		###############
		
		#set the orientation of the BAC/ contig BAC
		#3 levels, BAC aligned to REF generated from TPF
		#get the TPF line accession
		$insert_line_accession = $tpf_lines{$insert_line_number}->get_accession();
		$insert_line_strand    = $scaffold_agp_coords{$insert_line_accession}->{orientation}; # + or -
		
		if ($direction == 1 ){
			$bac_to_insert->set_orientation('PLUS');
			print STDERR "$insert_line_accession in $insert_line_strand orientation in original TPF, Mummer alignment in same direction in ref and qry $bac_name. Setting orientation to PLUS\n";
		}
		elsif ($direction == -1 ){
			$bac_to_insert->set_orientation('MINUS');
			print STDERR "$insert_line_accession in $insert_line_strand orientation in original TPF, Mummer alignment in opposite direction in ref and qry $bac_name. Setting orientation to MINUS\n";
		}
		
		#create line number to accession array and accession to TPF hash for insertions later
		#not using line number + offset logic as it breaks down in complicated cases
		my %accession_tpflines; # accession key = TPF line value in array 
		my @tpfline_accession; # [line number] = accession
		$line_key = 1;
		while ($line_key <= @sorted_tpf_line_numbers )
		{
			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
				my $accession = $tpf_lines{$line_key}->get_accession();
		
				#initialize data structures for insertions
				$tpfline_accession[$line_key] = $accession;
			
				my @temp_arr;
				push @temp_arr, $tpf_lines{$line_key}; #storing value instead of ref
				$accession_tpflines{$accession} = \@temp_arr;
			}
			$line_key++;
		}		

		#record insertion event and modify $insert_line_number if that line was already inserted at, issue #53
		if ( $insert_before_or_after eq 'before'){
			if (! exists $before_insertion_counter{$insert_line_number}){
				$before_insertion_counter{$insert_line_number} = 1 ;
			}
			else{
				my $old_insert_line_number = $insert_line_number;
				$insert_line_number = $insert_line_number - $before_insertion_counter{$insert_line_number};
				$before_insertion_counter{$old_insert_line_number}++ ;
			}
		}
		elsif ( $insert_before_or_after eq 'after'){
			if (! exists $after_insertion_counter{$insert_line_number}){
				$after_insertion_counter{$insert_line_number} = 1 ;
			}
			else{
				my $old_insert_line_number = $insert_line_number;
				$insert_line_number = $insert_line_number + $after_insertion_counter{$insert_line_number};
				$after_insertion_counter{$old_insert_line_number}++ ;
			}
		}

		my $accession = $tpfline_accession[$insert_line_number];
		my $tpf_line_arr_ref = $accession_tpflines{$accession};
		my @temp_arr = @$tpf_line_arr_ref;
		print STDERR "Fetched TPF line $insert_line_number for $accession for insertion\n";

		
		if ($bac_name =~ /^Contig/ ){
			print STDERR "Substituting in BACs for assembled contig $bac_name\n";
			if (!(exists $scaffold_component_contigs{$bac_name}) || !(exists $scaffold_component_contig_directions{$bac_name})){
				print STDERR "$bac_name not found in user supplied ACE file. Exiting ....\n\n"; exit 1;
			}
			my $component_accessions_ref = $scaffold_component_contigs{$bac_name};
			my @component_accessions_arr = @$component_accessions_ref;
			my $component_accession_directions_ref = $scaffold_component_contig_directions{$bac_name};# orientation 
			my @component_accession_directions_arr = @$component_accession_directions_ref;            # +1 for positive strand in contig alignment, -1 if on negative strand
			
			my $component_accessions_count = scalar @component_accessions_arr;
			my $contig_bac_loop_counter;
			
			if ((($bac_to_insert->get_orientation() eq 'PLUS')
				&& ($insert_before_or_after eq 'after'))
				|| 
				(($bac_to_insert->get_orientation() eq 'MINUS')
				&& ($insert_before_or_after eq 'before'))
				){
				#reverse order, flip orientation if MINUS
				$contig_bac_loop_counter = 0;
				while ($contig_bac_loop_counter < $component_accessions_count){#first to last
					my $contig_bac_to_insert = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
					print STDERR "******** inserting ";
					print STDERR $component_accessions_arr[$contig_bac_loop_counter];
					print STDERR "\n";
					#print Dumper $component_accessions_arr[$contig_bac_loop_counter];
					$contig_bac_to_insert->set_accession($component_accessions_arr[$contig_bac_loop_counter]);
					$bac_inserted_accessions{$component_accessions_arr[$contig_bac_loop_counter]} = 'inserted'; #recording accession name
					$contig_bac_to_insert->set_local_contig_identifier($bac_to_insert->get_local_contig_identifier() );
					
					#set contained, do NOT set contained if the accession was already deleted from TPF
					if (($bac_to_insert->has_contains()) &&
						(!exists $sequence_accessions_to_remove{$bac_to_insert->get_containing_accession()})){
						$contig_bac_to_insert->set_contains('CONTAINED');
						$contig_bac_to_insert->set_containing_accession($bac_to_insert->get_containing_accession());
						print STDERR "Added containing clause for BAC ";
						print STDERR $contig_bac_to_insert->get_accession();
						print STDERR " with containing accession ";
						print STDERR $contig_bac_to_insert->get_containing_accession();
						print STDERR " and BAC Contig on PLUS strand\n";
					}
					
					if ($component_accession_directions_arr[$contig_bac_loop_counter] == 1){
						if ($bac_to_insert->get_orientation() eq 'PLUS'){
							$contig_bac_to_insert->set_orientation('PLUS');
						}
						elsif($bac_to_insert->get_orientation() eq 'MINUS'){
							$contig_bac_to_insert->set_orientation('MINUS'); #flip
						}
					}
					elsif ($component_accession_directions_arr[$contig_bac_loop_counter] == -1){
						if ($bac_to_insert->get_orientation() eq 'PLUS'){
							$contig_bac_to_insert->set_orientation('MINUS');
						}
						elsif($bac_to_insert->get_orientation() eq 'MINUS'){
							$contig_bac_to_insert->set_orientation('PLUS'); #flip
						}
					}
					
					#insert 
					if ( $insert_before_or_after eq 'before' ){
						unshift @temp_arr, $contig_bac_to_insert;
						
					}
					elsif ( $insert_before_or_after eq 'after' ){
						push @temp_arr, $contig_bac_to_insert;
					}
					
					print STDERR "Inserted BAC: ";
					print STDERR $component_accessions_arr[$contig_bac_loop_counter];
					print STDERR " for assembled contig $bac_name in reversed order using simple loop ";
					print STDERR "$insert_before_or_after accession $accession\n";
					$contig_bac_loop_counter++;
				}
			}
			elsif ((($bac_to_insert->get_orientation() eq 'PLUS')
				&& ($insert_before_or_after eq 'before'))
				|| 
				(($bac_to_insert->get_orientation() eq 'MINUS')
				&& ($insert_before_or_after eq 'after'))
				){
				#simple order, flip orientation if MINUS
				$contig_bac_loop_counter = $component_accessions_count - 1 ;
				while ($contig_bac_loop_counter >= 0 ){#last to first
					my $contig_bac_to_insert = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
					print STDERR "******** inserting ";
					print STDERR $component_accessions_arr[$contig_bac_loop_counter];
					print STDERR "\n";
					$contig_bac_to_insert->set_accession($component_accessions_arr[$contig_bac_loop_counter]);
					$bac_inserted_accessions{$component_accessions_arr[$contig_bac_loop_counter]} = 'inserted'; #recording accession name
					$contig_bac_to_insert->set_local_contig_identifier($bac_to_insert->get_local_contig_identifier() );
					
					#set contained, do NOT set contained if the accession was already deleted from TPF
					if (($bac_to_insert->has_contains()) &&
						(!exists $sequence_accessions_to_remove{$bac_to_insert->get_containing_accession()})){
						$contig_bac_to_insert->set_contains('CONTAINED');
						$contig_bac_to_insert->set_containing_accession($bac_to_insert->get_containing_accession());
						
						print STDERR "Added containing clause for BAC ";
						print STDERR $contig_bac_to_insert->get_accession();
						if ($bac_to_insert->has_contains()){
							print STDERR " with containing accession ";
							print STDERR $contig_bac_to_insert->get_containing_accession();							
						}
						print STDERR " and BAC Contig on MINUS strand\n";
					}
					
					if ($component_accession_directions_arr[$contig_bac_loop_counter] == 1){
						if ($bac_to_insert->get_orientation() eq 'PLUS'){
							$contig_bac_to_insert->set_orientation('PLUS');
						}
						elsif($bac_to_insert->get_orientation() eq 'MINUS'){
							$contig_bac_to_insert->set_orientation('MINUS'); #flip
						}
					}
					elsif ($component_accession_directions_arr[$contig_bac_loop_counter] == -1){
						if ($bac_to_insert->get_orientation() eq 'PLUS'){
							$contig_bac_to_insert->set_orientation('MINUS');
						}
						elsif($bac_to_insert->get_orientation() eq 'MINUS'){
							$contig_bac_to_insert->set_orientation('PLUS'); #flip
						}
					}
					
					#insert
					if ( $insert_before_or_after eq 'before' ){
						unshift @temp_arr, $contig_bac_to_insert;
						
					}
					elsif ( $insert_before_or_after eq 'after' ){
						push @temp_arr, $contig_bac_to_insert;
					}

					print STDERR "Inserted BAC: ";
					print STDERR $component_accessions_arr[$contig_bac_loop_counter];
					print STDERR " for assembled contig $bac_name in simple order using reverse loop ";
					print STDERR "$insert_before_or_after accession $accession\n";
					$contig_bac_loop_counter--;
				}
			}
			else{
				die "BAC $bac_name not inserted. Exiting ....\n";
			}
		}
		else{#if its a singleton
			$bac_inserted_accessions{$bac_to_insert->get_accession()} = 'inserted'; #recording accession name
			if ( $insert_before_or_after eq 'before' ) {
				unshift @temp_arr, $bac_to_insert;
			}
			elsif ( $insert_before_or_after eq 'after' ) {
				push @temp_arr, $bac_to_insert;
			}
			else {
				die "BAC $bac_name not inserted\n";
			}
			print STDERR "Inserted singleton BAC: $bac_name $insert_before_or_after accession $accession\n";
		}
		
		$accession_tpflines{$accession} = \@temp_arr;
		
		#reconstruct tpf_lines and assign to module data structure
		my %modified_tpf_lines;
		%tpf_lines = %{ $self->get_tpf_lines() };
		@sorted_tpf_line_numbers = sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive
		my $modified_line_key;
		$modified_line_key = $line_key = 1;
		#$self->clear_tpf_lines();#remove old values
		while ($line_key <= @sorted_tpf_line_numbers )
		{
			if (defined $tpfline_accession[$line_key]){#sequence TPF lines
				my $temp_arr_ref = $accession_tpflines{$tpfline_accession[$line_key]};
				my @temp_arr = @$temp_arr_ref;
				
				foreach my $value (@temp_arr){
					$modified_tpf_lines{$modified_line_key++} = $value;
				}
				
			}
			else{#gap TPF lines
				$modified_tpf_lines{$modified_line_key++} = $tpf_lines{$line_key};
			}
			$line_key++;
		}

		$self->set_tpf_lines(\%modified_tpf_lines); #assign back to module data structure
		
		%tpf_lines = %{ $self->get_tpf_lines() };
		#print STDERR $self->get_formatted_tpf();
		
		###############
		#TODO FIX Contig names in switch points to first or last BAC in the contig
		###############
		
		
		###############
	}	
	
	###############
	#remove contained if the accession (typically WGS from prev assembly) was deleted from TPF
	###############
	
	my %lines = %{ $self->get_tpf_lines() };
	my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
	foreach my $line_key (@sorted_line_numbers) {
		if ( $lines{$line_key}->get_line_type() eq "sequence" ) {
			if (($lines{$line_key}->has_contains()) &&
				(exists $sequence_accessions_to_remove{$lines{$line_key}->get_containing_accession()})){
				my $accession_to_delete = $lines{$line_key}->get_containing_accession();
				$lines{$line_key}->clear_contains();
				$lines{$line_key}->clear_containing_accession();
				print STDERR "Removed containing accession $accession_to_delete for BAC ";
				print STDERR $lines{$line_key}->get_accession();
				print STDERR "\n";
			}
		}
	}
	$self->set_tpf_lines(\%lines); #assign back to module data structure
	
	return ($self,$switch_points,$trim_points);
}


sub move_scaffold_before {
	my $self                      = shift;
	my $scaffold_to_insert_before = shift;
	my $scaffold_to_modify        = shift;
	my $previous_accession;
	my %lines;
	my $tpf_copy                     = Bio::GenomeUpdate::TPF->new();
	my $first_line_number            = 0;
	my $last_line_number             = 0;
	my $in_scaffold_range            = 0;
	my $line_number_to_insert_before = 0;
	my $last_line;

	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
	}
	my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
	foreach my $line_key (@sorted_line_numbers) {
		if ( $lines{$line_key}->get_line_type() eq "sequence" ) {
			if ($previous_accession) {
				if ( $lines{$line_key}->get_local_contig_identifier() eq
					$scaffold_to_modify )
				{
					$in_scaffold_range = 1;
					$tpf_copy->add_line_to_end( $lines{$line_key} );
					if ($first_line_number) {
						$last_line_number = $line_key;
					}
					else {
						$first_line_number = $line_key;
						$last_line_number  = $line_key;
					}
				}
				else {
					$in_scaffold_range = 0;
				}
			}
			else {
				$in_scaffold_range = 0;
			}
			$previous_accession = $lines{$line_key}->get_accession();
		}
		elsif ( $lines{$line_key}->get_line_type() eq "gap" ) {
			if ( $in_scaffold_range == 1 ) {
				$tpf_copy->add_line_to_end( $lines{$line_key} );
				$last_line_number = $line_key;
			}
		}
	}

	#remove lines in range
	foreach my $line_key ( reverse( $first_line_number .. $last_line_number ) )
	{
		$self->delete_line($line_key);
	}
	%lines = %{ $self->get_tpf_lines() };
	@sorted_line_numbers = sort { $a <=> $b } keys %lines;
	foreach my $line_key (@sorted_line_numbers) {
		$last_line = $line_key;
		if ( $lines{$line_key}->get_line_type() eq "sequence" ) {
			if ( $line_number_to_insert_before == 0 ) {
				if ( $lines{$line_key}->get_local_contig_identifier() eq
					$scaffold_to_insert_before )
				{
					$line_number_to_insert_before = $line_key;
				}
			}
		}
	}
	my %lines_in_copy = %{ $tpf_copy->get_tpf_lines() };
	my @sorted_line_numbers_in_copy = sort { $a <=> $b } keys %lines_in_copy;
	foreach my $line_key (@sorted_line_numbers_in_copy) {
		$self->insert_line_after( $line_number_to_insert_before - 1,
			$lines_in_copy{$line_key} );
		$line_number_to_insert_before++;
		$last_line++;
	}
}

sub flip_scaffold {
	my $self               = shift;
	my $scaffold_to_modify = shift;
	my $previous_accession;
	my %lines;
	my $tpf_copy                     = Bio::GenomeUpdate::TPF->new();
	my $first_line_number            = 0;
	my $last_line_number             = 0;
	my $in_scaffold_range            = 0;
	my $line_number_to_insert_before = 0;
	my $last_line;

	if ( $self->has_tpf_lines() ) {
		%lines = %{ $self->get_tpf_lines() };
	}
	my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
	foreach my $line_key (@sorted_line_numbers) {
		if ( $lines{$line_key}->get_line_type() eq "sequence" ) {
			if ($previous_accession) {
				if ( $lines{$line_key}->get_local_contig_identifier() eq
					$scaffold_to_modify )
				{
					$in_scaffold_range = 1;
					$tpf_copy->add_line_to_end( $lines{$line_key} );
					if ($first_line_number) {
						$last_line_number = $line_key;
					}
					else {
						$first_line_number = $line_key;
						$last_line_number  = $line_key;
					}
				}
				else {
					$in_scaffold_range = 0;
				}
			}
			else {
				$in_scaffold_range = 0;
			}
			$previous_accession = $lines{$line_key}->get_accession();
		}
		elsif ( $lines{$line_key}->get_line_type() eq "gap" ) {
			if ( $in_scaffold_range == 1 ) {
				$tpf_copy->add_line_to_end( $lines{$line_key} );
				$last_line_number = $line_key;
			}
		}
	}
	$tpf_copy->delete_line($last_line_number);

	#remove lines in range
	foreach
	  my $line_key ( reverse( $first_line_number .. $last_line_number - 1 ) )
	{
		$self->delete_line($line_key);
	}
	$line_number_to_insert_before = $first_line_number - 1;
	my %lines_in_copy               = %{ $tpf_copy->get_tpf_lines() };
	my @sorted_line_numbers_in_copy =
	  reverse( sort { $a <=> $b } keys %lines_in_copy );
	foreach my $line_key (@sorted_line_numbers_in_copy) {
		$self->insert_line_after( $line_number_to_insert_before,
			$lines_in_copy{$line_key} );
		$line_number_to_insert_before++;
		$last_line++;
	}
}

###
1;    #do not remove
###

=pod

=back

=head1 LICENSE

Same as Perl.

=head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>
  Surya Saha        <suryasaha@cornell.edu>

=cut
