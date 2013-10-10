package Bio::GenomeUpdate::GFF;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GFF3::LowLevel qw (gff3_parse_feature  gff3_format_feature);
use Bio::GenomeUpdate::GFF::GFFRearrange;

=head1 NAME

    GFF - Generic Feature Format (GFF) contains the annotation for a genome.

=head1 SYNOPSIS

    my $variable = Bio::GenomeUpdate::GFF->new();

=head1 DESCRIPTION

    This class stores Generic Feature Format (GFF) information including coordinates, strand and source. It reads in old and new Accessioned Golden Path (AGP) files and prints a GFF with updated coordinates. GFF features than span scaffolds are not handled and written out to errors.gff3 file. 


=head2 Methods

=over

=cut

has 'comment_lines' => (
	isa       => 'ArrayRef[Str]',
	is        => 'rw',
	predicate => 'has_comment_lines',
	clearer   => 'clear_comment_lines'
);

=item C<add_comment_line ( @comment_string )>

Adds a comment line (optional).

=cut
sub add_comment_line {
	my $self           = shift;
	my $comment_to_add = shift;
	my @lines;
	if ( $self->has_comment_lines() ) {
		@lines = @{ $self->get_comment_lines() };
	}
	push( @lines, $comment_to_add );
	$self->set_comment_lines( [@lines] );
}

has 'gff_lines' => (
	isa       => 'ArrayRef[HashRef]',
	is        => 'rw',
	predicate => 'has_gff_lines',
	clearer   => 'clear_gff_lines'
);

=item C<add_gff_line ( $HashRef )>

Adds hash of GFF3 line.

=cut

sub add_gff_line {
	my $self = shift;
	my $line_to_add = shift;
	my @lines;
	if ( $self->has_gff_lines() ) {
		@lines = @{ $self->get_gff_lines() };
	}
	push( @lines, $line_to_add );
	$self->set_gff_lines( [@lines] );
}

=item C<parse_gff ( $gff_file_content )>

Reads in the input GFF3 using Bio::GFF3::LowLevel::Parser.

=cut

sub parse_gff{
	my $self = shift;
	my $gff_file_content = shift;
	$self->clear_comment_lines();
	$self->clear_gff_lines();
	
	my @lines = split ( /\n/, $gff_file_content);
	my $line_count         = scalar(@lines);
	my $line_counter       = 0;
	
	foreach my $line (@lines) {
		$line_counter++;
		chomp($line);

		if ( $line =~ m/^#/ ) {
			$self->add_comment_line($line);
			next;
		}
		else{
			my $gff_line_hash = gff3_parse_feature($line);
			$self->add_gff_line($gff_line_hash);
			
		}
		print STDERR "\rParsing GFF3 line ". $line_counter . " of "
				  . $line_count;
	}
	print STDERR "\nParsed $line_count GFF3 lines\n";
	return $self;
}

=item C<get_formatted_gff ()>

Returns formatted GFF3.  

=cut

sub get_formatted_gff{
	my $self = shift;
	my $out_str;
	
	if ( $self->has_comment_lines() ) {
		foreach my $comment ( @{ $self->get_comment_lines() } ) {
			$out_str .= $comment . "\n";
		}
	}
	
	if ( $self-> has_gff_lines()){
		foreach my $gff_line_hash ( @{ $self->get_gff_lines()} ){
			$out_str .= gff3_format_feature($gff_line_hash);
		}
	}
	return $out_str;
}


=item C<get_reordered_coordinates ( $agp_old, $agp_new )>

Returns a hash of updated coordinates. Only accepts AGP objects now. Calls GFFRearrange object for reordering.  

=cut
sub get_reordered_coordinates{
	my $self = shift;
	my $agp_old = shift;
	my $agp_new = shift;
	
	#reset current line number if already processed once
	$agp_old->set_current_agp_line_number(1);
	$agp_new->set_current_agp_line_number(1);
	
	my $gff_rearrange_obj = Bio::GenomeUpdate::GFF::GFFRearrange->new();
	my %coordinates = $gff_rearrange_obj->reordered_coordinates_AGP( $agp_old, $agp_new); 
	
	return %coordinates;
}

=item C<get_flipped_coordinates ( $agp_old, $agp_new )>

Returns a hash of flipped coordinates (0 or 1). Only accepts AGP objects now. Calls GFFRearrange object.  

=cut
sub get_flipped_coordinates{
	my $self = shift;
	my $agp_old = shift;
	my $agp_new = shift;
	
	#reset current line number if already processed once
	$agp_old->set_current_agp_line_number(1);
	$agp_new->set_current_agp_line_number(1);

	my $gff_rearrange_obj = Bio::GenomeUpdate::GFF::GFFRearrange->new();
	my %flipped = $gff_rearrange_obj->flipped_coordinates_AGP( $agp_old, $agp_new); 
	
	return %flipped;
}

=item C<remap_coordinates ( $agp_old, $agp_new )>

Updates GFF coordinates according to mapping function in GFFRearrange::updated_coordinates_strand_AGP routine. GFF features than span scaffolds are not handled and written out to errors.gff3 file.  

=cut
sub remap_coordinates{
	my $self = shift;
	my $agp_old = shift;
	my $agp_new = shift;
	
	if ( $self-> has_gff_lines()){
		my @lines = @{ $self->get_gff_lines()};
		$self->clear_gff_lines();
		my $errors = '';
		
		foreach my $gff_line_hash ( @lines ){
			my $start = $gff_line_hash->{'start'};
			my $end = $gff_line_hash->{'end'};
			my $strand = $gff_line_hash->{'strand'};
			
			my $gff_rearrange_obj = Bio::GenomeUpdate::GFF::GFFRearrange->new();
			my ($nstart, $nend, $nstrand) = $gff_rearrange_obj->updated_coordinates_strand_AGP( $start, $end, $strand, $agp_old, $agp_new);
			
			if (($nstart == 0) && ($nend == 0) && ($nstrand == 0)){
				$errors .= gff3_format_feature($gff_line_hash);
			}
			else{
				$gff_line_hash->{'start'} = $nstart;
				$gff_line_hash->{'end'} = $nend;
				$gff_line_hash->{'strand'} = $nstrand;
				
				#add back to $self
				$self->add_gff_line($gff_line_hash);
			}
		}
		
		#print GFF lines where feature spanned scaffolds	
		if ($errors ne ''){
			open(EGFF,">errors.gff3");
			print EGFF $errors;
			close(EGFF);
		}
	}	
	return $self;
}
=item C<remap_coordinates_hash ( %coordinates, %flipped )>

Updates GFF coordinates according to mapping in %coordinates hash. Deprecated as GFF features than span scaffolds are not handled separately and may be erroneous. Please double check. Use remap_coordinates()

=cut
sub remap_coordinates_hash{

	my $self = shift;
	my $coordinates = shift;
	my $flipped = shift;
	
	if ( $self-> has_gff_lines()){
		my @lines = @{ $self->get_gff_lines()};
		$self->clear_gff_lines();
		
		foreach my $gff_line_hash ( @lines ){
			my $start = $gff_line_hash->{'start'};
			my $end = $gff_line_hash->{'end'};
			my $strand = $gff_line_hash->{'strand'};
			
			#modify strand
			if ( $flipped->{$start} ){
				if ( $strand eq '+'){
					$gff_line_hash->{'strand'} = '-';
				}
				elsif( $strand eq '-'){
					$gff_line_hash->{'strand'} = '+';
				}
			}
			#modify coords
			if ( $flipped->{$start} ){
				$gff_line_hash->{'start'} = $coordinates->{$end};
				$gff_line_hash->{'end'} = $coordinates->{$start};
			}
			else{
				$gff_line_hash->{'start'} = $coordinates->{$start};
				$gff_line_hash->{'end'} = $coordinates->{$end};
			}
			
			#add back to $self
			$self->add_gff_line($gff_line_hash);
		}
	}	
	return $self;
}

###
1;				#do not remove
###


=back

=head1 LICENSE

Same as Perl.

=head1 AUTHORS

Surya Saha <suryasaha@cornell.edu , @SahaSurya>   

=cut
