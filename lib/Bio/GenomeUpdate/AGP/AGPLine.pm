package Bio::GenomeUpdate::AGP::AGPLine;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;

=head1 NAME

    AGPSequenceLine - A line in an Accessioned Golden Path (AGP)

=head1 SYNOPSIS

    my $variable = AGPLine->new();

=head1 DESCRIPTION

    This class stores information for a  line in an Accessioned Golden Path (AGP).

=head2 Methods

=over 

=cut

has 'object_being_assembled' => (isa => 'Str', is => 'rw', predicate => 'has_object_being_assembled');
has 'object_begin' => (isa => 'Int', is => 'rw', predicate => 'has_object_begin');
has 'object_end' => (isa => 'Int', is => 'rw', predicate => 'has_object_end');
has 'line_number' => (isa => 'Int', is => 'rw', predicate => 'has_line_number');
subtype 'ComponentType',
  as 'Str',
  where {$_ eq "A" || $_ eq "D" ||  $_ eq "F" || $_ eq "G" || $_ eq "O" || $_ eq "P" || $_ eq "W"};
has 'component_type' => (isa => 'ComponentType', is => 'rw', predicate => 'has_component_type');
has 'component_id' => (isa => 'Str', is => 'rw', predicate => 'has_component_id');
has 'component_begin' => (isa => 'Int', is => 'rw', predicate => 'has_component_begin');
has 'component_end' => (isa => 'Int', is => 'rw', predicate => 'has_component_end');
subtype 'AGPOrientationType',
  as 'Str',
  where { $_ eq "+" || $_ eq "-" || $_ eq "?" || $_ eq "0" || $_ eq "na"},
  message { "The string, $_, was not a valid orientation type.  Valid types are: + - ? 0 na" };
has 'orientation' => (isa => 'AGPOrientationType', is => 'rw', predicate => 'has_orientation');
sub get_line_type {
  return 'sequence';
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
