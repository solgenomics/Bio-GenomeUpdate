package  Bio::GenomeUpdate::SP::SPLine;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;

use Data::Dumper;#for debugging

=head1 NAME

    SP - Switch point lines for NCBI GRC pipeline with instructions used to generate a Accessioned Golden Path (AGP) file

=head1 SYNOPSIS

    my $variable = SPLine->new();

=head1 DESCRIPTION

This class stores information for transitions between TPF components including chromosome, accessions and coordinates for generating a switch point (SP) file. The switch point file specifies the exact point of transition between two consecutive components in the TPF file. 

=head2 Methods

=over 

=cut

=item C<set_chromosome ( $chromosome )>

Sets the chromosome. Valid values for Solanum lycopersicum are 1-12 and Un (required).

=item C<get_chromosome >

Gets the chromosome.

=cut

has 'chromosome' => ( isa => 'Str', is => 'rw', required => 1, clearer => 'clear_chromosome' );

has 'accession_prefix' => ( isa => 'Str', is => 'rw', required => 1, clearer => 'clear_accession_prefix' );
has 'accession_suffix' => ( isa => 'Str', is => 'rw', required => 1, clearer => 'clear_accession_suffix' );

subtype 'SPOrientationType',
  as 'Str',
  where { $_ eq "+" || $_ eq "-" || $_ eq "?" || $_ eq "0" || $_ eq "na"},
  message { "The string, $_, was not a valid orientation type.  Valid types are: + - ? 0 na" };

has 'accession_prefix_orientation' => ( isa => 'SPOrientationType', is => 'rw', required => 1, clearer => 'clear_accession_prefix_orientation' );
has 'accession_suffix_orientation' => ( isa => 'SPOrientationType', is => 'rw', required => 1, clearer => 'clear_accession_suffix_orientation' );

has 'accession_prefix_last_base' => ( isa => 'Int', is => 'rw', required => 1, clearer => 'clear_accession_prefix_last_base' );
has 'accession_suffix_first_base' => ( isa => 'Int', is => 'rw', required => 1, clearer => 'clear_accession_suffix_first_base' );

subtype 'SPComment',
  as 'Str',
  where { scalar $_ >= 25 },
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
