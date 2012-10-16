package Bio::GenomeUpdate::TPFSequenceLine;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;

=head1 NAME

    TPFSequenceLine - A sequence line in a Tiling Path File (TPF)

=head1 SYNOPSIS

    my $variable = TPFSequenceLine->new();

=head1 DESCRIPTION

    This class stores information for a sequence line in a Tiling Path File (TPF).

=head2 Methods

=over 

=cut

has 'accession' => (isa => 'Str', is => 'rw', predicate => 'has_accession');
has 'clone_name' => (isa => 'Str', is => 'rw', predicate => 'has_clone_name');
has 'local_contig_identifier' => (isa => 'Str', is => 'rw', predicate => 'has_local_contig_identifier');
subtype 'ContainsType',
  as 'Str',
  where { $_ eq "CONTAINED" || $_ eq "CONTAINED_TURNOUT"  },
  message { "The string, $_, was not CONTAINED or CONTAINED_TURNOUT" };
has 'contains' => (isa => 'ContainsType', is => 'rw', predicate => 'has_contains');
subtype 'OrientationType',
  as 'Str',
  where { $_ eq "PLUS" || $_ eq "MINUS"  },
  message { "The string, $_, was not PLUS or MINUS" };
has 'orientation' => (isa => 'OrientationType', is => 'rw', predicate => 'has_orientation');
has 'containing_accession' => (isa => 'Str', is => 'rw',  predicate => 'has_containing_accession');
has 'containing_clone_name' => (isa => 'Str', is => 'rw',  predicate => 'has_containing_clone_name');
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
