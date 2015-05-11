package  Bio::GenomeUpdate::SP;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::SP::SPLine;
use File::Slurp;

=head1 NAME

    SP - Switch point information for NCBI GRC pipeline with instructions used to generate a Accessioned Golden Path (AGP) file

=head1 SYNOPSIS

    my $variable = SwitchPoint->new();

=head1 DESCRIPTION

This class stores information for transitions between TPF components including TPF type, taxonomy and assembly and generates a switch point (SP) file. The switch point file specifies the exact point of transition between two consecutive components in the TPF file. 

=head2 Methods

=over 

=item C<set_taxid ( $taxonomy_id )>

Sets the taxonomy identifier for the SP file, e.g. 4081 for Solanum lycopersicum.

=item C<get_taxid>

Gets the taxonomy identifier for the SP file.

=cut

has 'sp_taxid' => (
	isa     => 'Str',
	is      => 'rw',
	default => '4081',
	required => 1,
	clearer => 'clear_sp_taxid'
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

has 'sp_tpf_type' => ( isa => 'SPTPFType', is => 'rw', default => 'chromosome', required => 1, clearer => 'clear_chromosome' );

subtype 'SPLine',
  as 'Bio::GenomeUpdate::SP::SPLine',
  message { "The object was not a SP line" };

has 'sp_lines' => (
	isa       => 'HashRef[SPLine]',
	is        => 'rw',
	predicate => 'has_sp_lines',
	clearer   => 'clear_sp_lines'
);

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
