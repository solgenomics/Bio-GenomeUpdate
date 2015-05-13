package  Bio::GenomeUpdate::SP;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::SP::SPLine;

=head1 NAME

    SP - Switch point information for NCBI GRC pipeline with instructions used to generate a Accessioned Golden Path (AGP) file

=head1 SYNOPSIS

    my $variable = SwitchPoint->new();

=head1 DESCRIPTION

This class stores information for transitions between TPF components including TPF type, taxonomy and assembly and generates a switch point (SP) file. The switch point file specifies the exact point of transition between two consecutive components in the AGP file. 

=head2 Methods

=over 

=item C<set_taxid ( $taxonomy_id )>

Sets the taxonomy identifier for the SP file, e.g. 4081 for Solanum lycopersicum.

=item C<get_taxid>

Gets the taxonomy identifier for the SP file.

=cut

has 'taxid' => (
	isa     => 'Str',
	is      => 'rw',
	default => '4081',
	required => 1,
	clearer => 'clear_taxid'
);

=item C<set_assembly_group ( $organism_string )>

Sets the assembly_group (required).

=item C<get_assembly_group>

Gets the assembly_group.

=cut

has 'assembly_group' =>
  ( isa => 'Str', is => 'rw', default => 'TGP', required => 1, clearer => 'clear_assembly_group' );

=item C<set_assembly_unit ( $organism_string )>

Sets the assembly_unit (required).

=item C<get_assembly_unit>

Gets the assembly_unit.

=cut

has 'assembly_unit' =>
  ( isa => 'Str', is => 'rw', default => 'Primary', required => 1, clearer => 'clear_assembly_unit' );

=item C<set_tpf_type ( $tpf_type )>

Sets the TPF type. Valid values are chromosome and contig (latter used for unlocalized or unplaced scaffolds, required).

=item C<get_tpf_type >

Gets the TPF type.

=cut

subtype 'SPTPFType', as 'Str', where { $_ eq "chromosome" || $_ eq "contig" },
  message { "The string, $_, was not a valid TPF type. See http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/overlap/ specification link" };

has 'tpf_type' => ( isa => 'SPTPFType', is => 'rw', default => 'chromosome', required => 1, clearer => 'clear_chromosome' );

subtype 'SPLine',
  as 'Bio::GenomeUpdate::SP::SPLine',
  message { "The object was not a SP line" };

has 'sp_lines' => (
	isa       => 'HashRef[SPLine]',
	is        => 'rw',
	predicate => 'has_sp_lines',
	clearer   => 'clear_sp_lines'
);

=item C<add_line_to_end >

Add a SPLine object.

=cut

sub add_line_to_end {
	my $self        = shift;
	my $line_to_add = shift;
	my %lines;
	if ( $self->has_sp_lines() ) {
		%lines = %{ $self->get_sp_lines() };
	}
	my $last_line = $self->get_number_of_lines();
	$lines{ $last_line + 1 } = $line_to_add;#key is just the index or line number
	$self->set_sp_lines( {%lines} );
}

=item C<get_number_of_lines >

Return number of lines in the SPLine object.

=cut

sub get_number_of_lines {
	my $self = shift;
	my %lines;
	if ( $self->has_sp_lines() ) {
		%lines = %{ $self->get_sp_lines() };
		my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
		return $sorted_line_numbers[-1];
	}
	else {
		return 0;
	}
}

=item C<get_formatted_sp >

Return string with all lines in the SPLine object.

=cut

sub get_formatted_sp {
	my $self = shift;
	my %lines;
	my $out_str;

	if ( $self->has_sp_lines() ) {
		%lines = %{ $self->get_sp_lines() };
		my @sorted_line_numbers = sort { $a <=> $b } keys %lines;
		foreach my $line_key (@sorted_line_numbers) {
			$out_str .= $self->get_taxid() . "\t";
			$out_str .= $self->get_assembly_group() . "\t";
			$out_str .= $self->get_assembly_unit() . "\t";
			$out_str .= $lines{$line_key}->get_chromosome() . "\t";
			$out_str .= $self->get_tpf_type() . "\t";
			$out_str .= $lines{$line_key}->get_accession_prefix() . "\t";
			$out_str .= $lines{$line_key}->get_accession_suffix() . "\t";
			$out_str .= $lines{$line_key}->get_accession_prefix_orientation() . "\t";
			$out_str .= $lines{$line_key}->get_accession_suffix_orientation() . "\t";
			$out_str .= $lines{$line_key}->get_accession_prefix_last_base() . "\t";
			$out_str .= $lines{$line_key}->get_accession_suffix_first_base() . "\t";
			$out_str .= $lines{$line_key}->get_comment() . "\n";
		}
	}
	return $out_str;
}

###
1;    #do not remove
###

=pod

=back

=head1 LICENSE

Same as Perl.

=head1 AUTHORS

  Surya Saha        <suryasaha@cornell.edu>

=cut
