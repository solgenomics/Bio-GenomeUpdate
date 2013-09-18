package Bio::GenomeUpdate::GFF::GFFRearrange;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::AGP;

=head1 NAME

    GFFRearrange - Generic Feature Format (GFF) Rearrange modifies the coordinates of a GFF file.

=head1 SYNOPSIS

    my $variable = Bio::GenomeUpdate::GFF::GFFRearrange->new();

=head1 DESCRIPTION

    This class modifies Generic Feature Format (GFF) coordinates using old and new Accessioned Golden Path (AGP) files. 

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
