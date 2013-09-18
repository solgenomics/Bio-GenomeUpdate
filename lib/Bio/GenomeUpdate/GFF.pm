package Bio::GenomeUpdate::GFF;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GFF3::LowLevel::Parser;

=head1 NAME

    GFF - Generic Feature Format (GFF) contains the annotation for a genome.

=head1 SYNOPSIS

    my $variable = Bio::GenomeUpdate::GFF->new();

=head1 DESCRIPTION

    This class stores Generic Feature Format (GFF) information including coordinates, strand and source. It reads in old and new Accessioned Golden Path (AGP) files and prints a GFF with updated coordinates. 

=head2 Methods

=over

=item C<parse_gff ( $gff_file_content )>







###
1;				#do not remove
###

=back

=head1 LICENSE

    Same as Perl.

=head1 AUTHORS

Surya Saha <suryasaha@cornell.edu , @SahaSurya>   

=cut
