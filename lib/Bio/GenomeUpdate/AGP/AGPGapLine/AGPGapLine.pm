package Bio::GenomeUpdate::AGPGapLine;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;

=head1 NAME

    AGPGapLine - A gap line in an Accessioned Golden Path (AGP)

=head1 SYNOPSIS

    my $variable = AGPGapLine->new();

=head1 DESCRIPTION

    This class stores information for a gap line in an Accessioned Golden Path (AGP).

=head2 Methods

=over 

=cut

has 'object_being_assembled' => (isa => 'Str', is => 'rw', predicate => 'has_object_being_assembled');
has 'object_begin' => (isa => 'Int', is => 'rw', predicate => 'has_object_begin');
has 'object_end' => (isa => 'Int', is => 'rw', predicate => 'has_object_end');
has 'line_number' => (isa => 'Int', is => 'rw', predicate => 'has_line_number');
subtype 'AGPGapComponentType',
  as 'Str',
  where {$_ eq "N" || $_ eq "U"};
has 'component_type' => (isa => 'AGPGapComponentType', is => 'rw', predicate => 'has_component_type');
has 'gap_length' => (isa => 'Int', is => 'rw', predicate => 'has_gap_length');
subtype 'AGPGapType',
  as 'Str',
  where {$_ eq "scaffold" || $_ eq "contig" || $_ eq "centromere" || $_ eq "short_arm" || $_ eq "heterochromatin" || $_ eq "telomere" || $_ eq "repeat" || $_ eq "fragment" || $_ eq "clone"}; #fragment and clone are not valid in v2.0
has 'gap_type' => (isa => 'AGPGapType', is => 'rw', predicate => 'has_gap_type');
subtype 'LinkageType',
  as 'Str',
  where {$_ eq "yes" || $_ eq "no"};
has 'linkage' => (isa => 'LinkageType', is => 'rw', predicate => 'has_linkage');
subtype 'LinkageEvidenceType',
  as 'Str',
  where {$_ eq "na" || $_ eq "paired-ends" || $_ eq "align_genus" || $_ eq "align_xgenus" || $_ eq "align_trnscpt" || $_ eq "within_clone" || $_ eq "clone_contig" || $_ eq "map" || $_ eq "strobe" || $_ eq "unspecified"};
has 'linkage_evidence' => (isa => 'ArrayRef[LinkageEvidenceType]', is => 'rw', predicate => 'has_linkage_evidence');
sub add_linkage_evidence {
  my $self = shift;
  my $evidence_to_add = shift;
  my @evidence;
  if ($self->has_linkage_evidence()) {
    @evidence = @{$self->get_linkage_evidence()};
  }
  push (@evidence, $evidence_to_add);
  $self->set_linkage_evidence([@evidence]);
}
sub get_line_type {
  return 'gap';
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
