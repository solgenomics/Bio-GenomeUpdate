package  Bio::GenomeUpdate::TP::TPLine;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;

use Data::Dumper;#for debugging

=head1 NAME

    TP - Trim point lines for NCBI GRC pipeline with instructions used to generate a Accessioned Golden Path (AGP) file

=head1 SYNOPSIS

    my $variable = TPLine->new();

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
  where { ( $_ >= 1 && $_ <=12 ) || ( $_ eq "Un" )},#does NOT work. need to do -> if int check 1-12, if str check Un
  message { "The string, $_, was not a valid chromosome number.  Valid values for Solanum lycopersicum are 1-12 and Un." };
has 'chromosome' => ( isa => 'TPChromosome', is => 'rw', required => 1, clearer => 'clear_chromosome' );

has 'accession' => ( isa => 'Str', is => 'rw', required => 1, clearer => 'clear_accession' );

subtype 'PositiveInt',
	as 'Int',
	where { $_ > 0 },
	message { "The string, $_, was not a positive coordinate" };
has 'accession_prefix_first_or_last_base' => ( isa => 'PositiveInt', is => 'rw', required => 1, clearer => 'clear_accession_prefix_first_or_last_base' );

subtype 'TPComment',
  as 'Str',
  where { (scalar $_) >= 25 },#does not work!!
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