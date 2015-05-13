package  Bio::GenomeUpdate::SP::SPLine;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;

=head1 NAME

    SP - Switch point lines for NCBI GRC pipeline with instructions used to generate a Accessioned Golden Path (AGP) file

=head1 SYNOPSIS

    my $sp_line = Bio::GenomeUpdate::SP::SPLine->new( 
						chromosome => 'Un',
						accession_prefix => 'accession1',
						accession_suffix => 'accession2',
						accession_prefix_orientation => '+',
						accession_suffix_orientation => '-',
						accession_prefix_last_base => 100,
						accession_suffix_first_base => 200,
						comment => 'test line 1 comment some nonsense here and here too');

=head1 DESCRIPTION

This class stores information for transitions between TPF components including chromosome, accessions and coordinates for generating a switch point (SP) file. The switch point file specifies the exact point of transition between two consecutive components in the AGP file. 

=head2 Methods

=over 

=cut

=item C<set_chromosome ( $chromosome )>

Sets the chromosome. Valid values for Solanum lycopersicum are 1-12 and Un (required).

=item C<get_chromosome >

Gets the chromosome.

=cut

subtype 'SPChromosome', #only works for genomes with 12 chrs or less.
  as 'Str',
  where { $_ eq '1' || $_ eq '2' || $_ eq '3' || $_ eq '4' || $_ eq '5' || $_ eq '6' || $_ eq '7' || $_ eq '8' || $_ eq '9' || $_ eq '10' || $_ eq '11' || $_ eq '12' || $_ eq "Un" },
  message { "The string, $_, was not a valid chromosome number.  Valid values for Solanum lycopersicum are 1-12 and Un." };
has 'chromosome' => ( isa => 'SPChromosome', is => 'rw', required => 1, clearer => 'clear_chromosome' );

=item C<set_accession_prefix ( $accession_prefix )>

Sets the accession_prefix (required).

=item C<get_accession_prefix >

Gets the accession_prefix.

=cut

has 'accession_prefix' => ( isa => 'Str', is => 'rw', required => 1, clearer => 'clear_accession_prefix' );

=item C<set_accession_suffix ( $accession_suffix )>

Sets the accession_suffix (required).

=item C<get_accession_suffix >

Gets the accession_suffix.

=cut

has 'accession_suffix' => ( isa => 'Str', is => 'rw', required => 1, clearer => 'clear_accession_suffix' );

=item C<set_accession_prefix_orientation ( $accession_prefix_orientation )>

Sets the accession_prefix_orientation. Valid types are: + - ? 0 na (required).

=item C<get_accession_prefix >

Gets the accession_prefix_orientation.

=cut

=item C<set_accession_suffix_orientation ( $accession_suffix_orientation )>

Sets the accession_suffix_orientation. Valid types are: + - ? 0 na (required).

=item C<get_accession_suffix_orientation >

Gets the accession_suffix_orientation.

=cut

subtype 'SPOrientationType',
  as 'Str',
  where { $_ eq "+" || $_ eq "-" || $_ eq "?" || $_ eq "0" || $_ eq "na"},
  message { "The string, $_, was not a valid orientation type.  Valid types are: + - ? 0 na" };
has 'accession_prefix_orientation' => ( isa => 'SPOrientationType', is => 'rw', required => 1, clearer => 'clear_accession_prefix_orientation' );
has 'accession_suffix_orientation' => ( isa => 'SPOrientationType', is => 'rw', required => 1, clearer => 'clear_accession_suffix_orientation' );

=item C<set_accession_prefix_last_base ( $accession_prefix_last_base )>

Sets the accession_prefix_last_base. Needs to be positive integer (required).

=item C<get_accession_suffix_orientation >

Gets the accession_prefix_last_base.

=cut

=item C<set_accession_suffix_first_base ( $accession_suffix_first_base )>

Sets the accession_suffix_first_base. Needs to be positive integer (required).

=item C<get_accession_suffix_first_base >

Gets the accession_suffix_first_base.

=cut

subtype 'PositiveInt',
	as 'Int',
	where { $_ > 0 },
	message { "The string, $_, was not a positive coordinate" };
has 'accession_prefix_last_base' => ( isa => 'PositiveInt', is => 'rw', required => 1, clearer => 'clear_accession_prefix_last_base' );
has 'accession_suffix_first_base' => ( isa => 'PositiveInt', is => 'rw', required => 1, clearer => 'clear_accession_suffix_first_base' );

=item C<set_comment ( $comment )>

Sets the comment. Needs to be at least 25 chars (required).

=item C<get_comment >

Gets the comment.

=cut

subtype 'SPComment',
  as 'Str',
  where { (length $_) >= 25 },
  message { "The string, $_, was shorter than the minimum length of 25 characters." };
has 'comment' => ( isa => 'SPComment', is => 'rw', required => 1, clearer => 'clear_comment' );

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
