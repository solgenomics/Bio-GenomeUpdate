package Bio::GenomeUpdate::AlignmentCoords;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;

=head1 NAME

    AlignmentCoords - Stores alignment coordinates for a reference and query

=head1 SYNOPSIS

    my $variable = AlignmentCoords->new();

=head1 DESCRIPTION

    This class stores alignment coordinates (e.g., from Nucmer/delta-filter/show-coords output) for a query sequence (e.g., from a BAC) relative to reference sequences (e.g., chromosome pseudomolecules).

=head2 Methods

=over 

=cut

has 'reference_start_coord' => (isa => 'Num', is => 'rw', predicate => 'has_reference_start_coord');
has 'reference_end_coord' => (isa => 'Num', is => 'rw', predicate => 'has_reference_end_coord');
has 'query_start_coord' => (isa => 'Num', is => 'rw', predicate => 'has_query_start_coord');
has 'query_end_coord' => (isa => 'Num', is => 'rw', predicate => 'has_query_end_coord');
has 'reference_id' => (isa => 'Str', is => 'rw', predicate => 'has_reference_id');
has 'query_id' => (isa => 'Str', is => 'rw',  predicate => 'has_query_id');

=item C<get_direction>

Returns 1 if the reference and query sequences are aligned in the same direction (strand) and returns -1 if they are reversed.

=cut

sub get_direction {
  my $self = shift;
  if (!($self->has_reference_start_coord() && $self->has_reference_end_coord() && $self->has_query_start_coord() && $self->has_query_end_coord())) {
    return undef;
  }
  if ((($self->get_reference_start_coord() < $self->get_reference_end_coord()) && ($self->get_query_start_coord < $self->get_query_end_coord))
      || (($self->get_reference_start_coord > $self->get_reference_end_coord) && ($self->get_query_start_coord > $self->get_query_end_coord))) {
    return 1;
  } elsif ((($self->get_reference_start_coord < $self->get_reference_end_coord) && ($self->get_query_start_coord > $self->get_query_end_coord))
	   || (($self->get_reference_start_coord > $self->get_reference_end_coord) && ($self->get_query_start_coord < $self->get_query_end_coord))) {
    return -1;
  } else {
    return undef;
  }
}

=item C<get_reference_length>

Returns the length of the aligned reference sequence.

=cut

sub get_reference_length {
  my $self = shift;
  if (!($self->has_reference_start_coord() && $self->has_reference_end_coord())) {
    return undef;
  } else {
    return (($self->has_reference_end_coord() + 1) - $self->has_reference_start_coord());
  }
}

=item C<get_query_length>

Returns the length of the aligned query sequence.

=cut

sub get_query_length {
  my $self = shift;
  if (!($self->has_query_start_coord() && $self->has_query_end_coord())) {
    return undef;
  } else {
    return (($self->has_query_end_coord() + 1) - $self->has_query_start_coord());
  }
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
