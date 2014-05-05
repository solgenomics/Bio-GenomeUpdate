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
							$tpf_sequence_line->set_containing_accesion(
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

=item C<get_summary ()>

Get a string with summary bases covered by sequence and gaps in TPF.

=cut

sub get_summary{
	my $self = shift;
	
	
	return $self;
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

sub get_tpf_with_bacs_inserted {
	my $self           = shift;
	my $bacs_ref       = shift;
	my $agp_coords_ref = shift;
	my @bacs           = @$bacs_ref;
	my %agp_coords     = %$agp_coords_ref;

	#make sure BACs are sorted by position
	foreach my $bac_ref (@bacs) {
		my @bac      = @$bac_ref;
		my $bac_name = $bac[0];
		my $bac_start;
		my $bac_end;
		my $bac_to_insert = Bio::GenomeUpdate::TPF::TPFSequenceLine->new();
		my %tpf_lines;
		if ( $self->has_tpf_lines() ) {
			%tpf_lines = %{ $self->get_tpf_lines() };
		}
		my @sorted_tpf_line_numbers =
		  sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive
		$bac_to_insert->set_accession($bac_name);
		if ( $bac[1] < $bac[2] ) {
			$bac_to_insert->set_orientation('PLUS');
			$bac_start = $bac[1];
			$bac_end   = $bac[2];
		}
		elsif ( $bac[1] > $bac[2] ) {
			$bac_to_insert->set_orientation('MINUS');
			$bac_start = $bac[2];
			$bac_end   = $bac[1];
		}
		else {
			die
"Error in BAC coordinates for BAC $bac_start Start: $bac_start End: $bac_end\n";
		}
		my $prev_agp_start = 0;
		my $prev_agp_end   = 0;
		my $prev_accession = 'none';
		my $prev_line_key;
		my $bac_is_contained = 0;
		my $agp_start;
		my $agp_end;
		my $bac_is_inserted = 0;
		my %gaps_to_resize;  #key will be line number and value will be new size
		my @sorted_gaps_to_resize;
		my %gaps_to_remove;    #key is line number value is undef
		my @rev_sorted_gaps_to_remove;
		my $insert_before_or_after = undef;
		my $insert_line_number     = undef;
		my %contained_contigs
		  ;    #key will be line number and value will be the contig accession
		my $past_bac = 0;
		my $line_key = 1;

		#add BAC coordinates to AGP info (not saved)
		my %add_agp_coords;
		$add_agp_coords{'start'} = $bac_start;
		$add_agp_coords{'end'}   = $bac_end;
		if ( $bac_to_insert->get_orientation() eq 'PLUS' ) {
			$add_agp_coords{'orientation'} = '+';
		}
		elsif ( $bac_to_insert->get_orientation() eq 'MINUS' ) {
			$add_agp_coords{'orientation'} = '-';
		}
		else {
			die "No orientation specified for BAC: $bac_name\n";
		}
		$agp_coords{$bac_name} = \%add_agp_coords;

		while ( $past_bac == 0 && $line_key <= @sorted_tpf_line_numbers + 1 ) {
			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
				my $accession = $tpf_lines{$line_key}->get_accession();
				my $agp_line_coords_ref = $agp_coords{$accession};
				my %line_coords         = %$agp_line_coords_ref;
				$agp_start = $line_coords{'start'};
				$agp_end   = $line_coords{'end'};

				#check if past the BAC
				if ( $bac_end < $prev_agp_start ) {
					$past_bac = 1;
				}

				#check if current contig is contained in the BAC
				if ( $agp_start >= $bac_start && $agp_end <= $bac_end ) {
					$tpf_lines{$line_key}->set_contains('CONTAINED');
					$tpf_lines{$line_key}->set_containing_accession($bac_name);
					$self->set_tpf_lines( \%tpf_lines );

					#$contained_contigs{$line_key}=$bac_name;
				}

				#check if current BAC is contained in the contig
				if ( $bac_start >= $agp_start && $bac_end <= $agp_end ) {
					$bac_to_insert->set_contains('CONTAINED');
					$bac_to_insert->set_containing_accession($accession);
				}

				#check if gap is spanned by the BAC
				if (   $prev_line_key
					&& $bac_start <= $prev_agp_end
					&& $bac_end >= $agp_start
					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
				{
					$gaps_to_remove{ $line_key - 1 } = 'delete';
					my $gap_location = $line_key - 1;
					print STDERR
"Removing gap at line $gap_location between $accession and $prev_accession\n";
				}

				#shrink gaps when partially spanned by a BAC
				if (   $prev_line_key
					&& $bac_start < $agp_start
					&& $bac_start > $prev_agp_end
					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
				{
					$gaps_to_resize{ $line_key - 1 } =
					  $bac_start - $prev_agp_end;
				}
				if (   $prev_line_key
					&& $bac_end < $agp_start
					&& $bac_end > $prev_agp_end
					&& $tpf_lines{ $line_key - 1 }->get_line_type() eq 'gap' )
				{
					$gaps_to_resize{ $line_key - 1 } = $agp_start - $bac_end;
				}
				$prev_line_key  = $line_key;
				$prev_accession = $accession;
				$prev_agp_start = $agp_start;
				$prev_agp_end   = $agp_end;
			}
			$line_key++;
		}

		@sorted_gaps_to_resize     = sort { $a <=> $b } keys %gaps_to_resize;
		@rev_sorted_gaps_to_remove = sort { $b <=> $a } keys %gaps_to_remove;

		foreach my $line_number (@sorted_gaps_to_resize) {
			$tpf_lines{$line_number}
			  ->set_gap_size( $gaps_to_resize{$line_number} );
		}
		$self->set_tpf_lines( \%tpf_lines );
		foreach my $line_number (@rev_sorted_gaps_to_remove) {
			$self->delete_line($line_number);
		}
		%tpf_lines               = %{ $self->get_tpf_lines() };
		@sorted_tpf_line_numbers =
		  sort { $a <=> $b } keys %tpf_lines;    #lines should be consecutive

		$line_key = 1;
		while ($bac_is_inserted == 0
			&& $line_key <= @sorted_tpf_line_numbers + 1 )
		{
			if ( $tpf_lines{$line_key}->get_line_type() eq 'sequence' ) {
				my $accession = $tpf_lines{$line_key}->get_accession();
				my $agp_line_coords_ref = $agp_coords{$accession};
				my %line_coords         = %$agp_line_coords_ref;
				$agp_start = $line_coords{'start'};
				$agp_end   = $line_coords{'end'};
				if ( $line_key == 1 ) {          #deal with first one
					if ( $bac_start <= 0 ) {
						$insert_before_or_after = 'before';
						$insert_line_number     = $line_key;
						$bac_to_insert->set_local_contig_identifier(
							$tpf_lines{$line_key}->get_local_contig_identifier()
						);
						$bac_is_inserted = 1;
					}
				}
				elsif ( $line_key == @sorted_tpf_line_numbers + 1 )
				{                                #deal with last one
					if ( $bac_start >= $agp_start ) {
						$insert_before_or_after = 'after';
						$insert_line_number     = $line_key;
						$bac_to_insert->set_local_contig_identifier(
							$tpf_lines{$line_key}->get_local_contig_identifier()
						);
						$bac_is_inserted = 1;
					}
				}
				elsif ($bac_start >= $prev_agp_start
					&& $bac_start < $agp_start )
				{
					if ( $bac_start <= $prev_agp_end ) {
						$insert_before_or_after = 'after';
						$insert_line_number     = $prev_line_key;
						$bac_to_insert->set_local_contig_identifier(
							$tpf_lines{$prev_line_key}
							  ->get_local_contig_identifier() );
						$bac_is_inserted = 1;
					}
					elsif ( $bac_start > $prev_agp_end ) {
						$insert_before_or_after = 'before';
						$insert_line_number     = $line_key;
						$bac_to_insert->set_local_contig_identifier(
							$tpf_lines{$line_key}->get_local_contig_identifier()
						);
						$bac_is_inserted = 1;
					}
				}
				$prev_line_key  = $line_key;
				$prev_accession = $accession;
				$prev_agp_start = $agp_start;
				$prev_agp_end   = $agp_end;
			}
			$line_key++;
		}

		if ( $insert_before_or_after eq 'before' ) {
			$self->insert_line_before( $insert_line_number, $bac_to_insert );
		}
		elsif ( $insert_before_or_after eq 'after' ) {
			$self->insert_line_after( $insert_line_number, $bac_to_insert );
		}
		else {
			die "BAC $bac_name not inserted\n";
		}
		%tpf_lines = %{ $self->get_tpf_lines() };
	}
	return $self;
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
