package Bio::GenomeUpdate::TPFGapLine;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
#use MooseX::Types;

=head1 NAME

    TPFGapLine - A gap line in a Tiling Path File (TPF)

=head1 SYNOPSIS

    my $variable = TPFGapLine->new();

=head1 DESCRIPTION

    This class stores information for a gap line in a Tiling Path File (TPF).

=head2 Methods

=over 

=item C<set_gap_identifier ( $identifier )>

Sets the string used to identify the line as a gap.  The default identifier is "GAP".

=item C<get_gap_identifier>

Gets the string used to identify the line as a gap.

=cut

has 'gap_identifier' => (isa => 'Str', is => 'rw', default => 'GAP');

=item C<set_gap_type ( $type )>

Sets the type of gap.
TYPE-1: [Deprecated]- was a place-holder for a picked clone.
TYPE-2: clone gap
TYPE-3: contig gap- unable to close using available technology
Biological Gap: If there is a biological gap such as a centromere, etc. 
then use the name rather than type-4. This is a controlled vocabulary:
CENTROMERE
TELOMERE
HETEROCHROMATIN
SHORT-ARM

=item C<get_gap_type>

Gets the gap type

=cut

subtype 'GapType',
  as 'Str',
  where { $_ eq "TYPE-1" || $_ eq "TYPE-2" || $_ eq "TYPE-3" || $_ eq "CENTROMERE"|| $_ eq "TELOMERE" || $_ eq "HETEROCHROMATIN" || $_ eq "SHORT-ARM" },
  message { "The string, $_, was not a valid gap type" };
has 'gap_type' => (isa => 'GapType', is => 'rw', predicate => 'has_gap_type');

=item C<set_gap_size ( $type )>

Sets the estimated size of the gap.  The gap size must be a positive integer.  Used for TYPE-2, TYPE-3 or Biological gap.  Setting a gap size is optional unless the 
gap type is ‘Biological’.

=item C<get_gap_size>

Gets the estimated size of the gap.

=cut

subtype 'PositiveInt',
  as 'Int',
  where { $_ > 0 },
  message { "Int is not larger than 0" };
has 'gap_size' => (isa => 'PositiveInt', is => 'rw', predicate => 'has_gap_size');



=item C<get_gap_methods>

Gets the method(s) used to determine the gap size as an array of strings;

=cut

subtype 'GapMethod',
  as 'Str',
  where {  $_ eq "FISH" ||  $_ eq "OPTICAL MAP" ||  $_ eq "RADIATION HYBRID" ||  $_ eq "PCR" ||  $_ eq "FINGERPRINT" ||  $_ eq "PAIRED ENDS" ||  $_ eq "ALIGN GENUS" ||  $_ eq "ALIGN XGENUS" ||  $_ eq "ALIGN TRNSCPT" },
  message { "The string, $_, was not a valid gap method" };
has 'gap_methods' => (isa => 'ArrayRef[GapMethod]', is => 'rw', predicate => 'has_gap_methods');


=item C<add_gap_method ( $method )>

Adds a gap method.  Valid values are:
FISH
OPTICAL MAP
RADIATION HYBRID
PCR
FINGERPRINT
PAIRED ENDS
ALIGN GENUS
ALIGN XGENUS
ALIGN TRNSCPT

=cut

sub add_gap_method {
  my $self = shift;
  my $method_to_add = shift;
  if ($self->has_gap_methods()) {
    my @methods = @{$self->get_gap_methods()};
    push(@methods, $method_to_add);
    $self->set_gap_methods([@methods]);
  } else {
    $self->set_gap_methods([$method_to_add]);
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
