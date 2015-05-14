package  Bio::GenomeUpdate::TP::TPLine;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;

=head1 NAME

    TP - Trim point lines for NCBI GRC pipeline with instructions used to generate a Accessioned Golden Path (AGP) file

=head1 SYNOPSIS

    my $tp_line = Bio::GenomeUpdate::TP::TPLine->new( 
					chromosome => 'Un',
					accession => 'accession1',
					accession_prefix_first_or_last_base => 100,
					trim_from_end => 'L',
					comment => 'test line 1 comment some nonsense here and here too');

=head1 DESCRIPTION

This class stores information for transitions between TPF components including chromosome, accessions and coordinates for generating a trim point (TP) file. The trim point file specifies the extent of a component in the AGP file. 

=head2 Methods

=over 

=cut

=item C<set_chromosome ( $chromosome )>

Sets the chromosome. Valid values for Solanum lycopersicum are 1-12 and Un (required).

=item C<get_chromosome >

Gets the chromosome.

=cut

subtype 'TPChromosome',
  as 'Str',
  where { $_ eq '1' || $_ eq '2' || $_ eq '3' || $_ eq '4' || $_ eq '5' || $_ eq '6' || $_ eq '7' || $_ eq '8' || $_ eq '9' || $_ eq '10' || $_ eq '11' || $_ eq '12' || $_ eq "Un" },
  message { "The string, $_, was not a valid chromosome number.  Valid values for Solanum lycopersicum are 1-12 and Un." };
has 'chromosome' => ( isa => 'TPChromosome', is => 'rw', required => 1, clearer => 'clear_chromosome' );

=item C<set_accession ( $accession )>

Sets the accession (required).

=item C<get_accession >

Gets the accession.

=cut

has 'accession' => ( isa => 'Str', is => 'rw', required => 1, clearer => 'clear_accession' );

=item C<set_accession_prefix_first_or_last_base ( $accession_prefix_first_or_last_base )>

Sets the accession_prefix_first_or_last_base. Needs to be a positive integer... duh (required).

=item C<get_accession_prefix_first_or_last_base >

Gets the accession_prefix_first_or_last_base.

=cut

subtype 'TPPositiveInt',
	as 'Int',
	where { $_ > 0 },
	message { "The string, $_, was not a positive coordinate" };
has 'accession_prefix_first_or_last_base' => ( isa => 'TPPositiveInt', is => 'rw', required => 1, clearer => 'clear_accession_prefix_first_or_last_base' );

=item C<set_trim_from_end ( $trim_from_end )>

Sets the trim_from_end. Valid values values are L and H where L: trim bases with values lt accession_prefix_first_or_last_base; H: trim bases with values gt accession_prefix_first_or_last_base. (required).

=item C<get_trim_from_end >

Gets the trim_from_end.

=cut

subtype 'TPTrimDirection',
  as 'Str',
  where { ( $_ eq 'L' ) || ( $_ eq "H" )},
  message { "The string, $_, was not a valid trim direction.  Valid values values are L and H where L: trim bases with values lt accession_prefix_first_or_last_base; H: trim bases with values gt accession_prefix_first_or_last_base." };
has 'trim_from_end' => ( isa => 'TPTrimDirection', is => 'rw', required => 1, clearer => 'clear_trim_from_end' );

=item C<set_comment ( $comment )>

Sets the comment. Needs to be at least 25 chars (required).

=item C<get_comment >

Gets the comment.

=cut

subtype 'TPComment',
  as 'Str',
  where { (length $_) >= 25 },
  message { "The string, $_, was shorter than the minimum length of 25 characters." };
has 'comment' => ( isa => 'TPComment', is => 'rw', required => 1, clearer => 'clear_comment' );

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