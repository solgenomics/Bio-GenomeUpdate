package Bio::GenomeUpdate::GFF::GFFRearrange;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::AGP;
use Bio::GenomeUpdate::GFF;
use Scalar::Util 'looks_like_number';

=head1 NAME

    GFFRearrange - Generic Feature Format (GFF) Rearrange modifies the coordinates of a GFF file. 

=head1 SYNOPSIS

    my $variable = Bio::GenomeUpdate::GFF::GFFRearrange->new();

=head1 DESCRIPTION

    This class modifies Generic Feature Format (GFF) coordinates using old and new Accessioned Golden Path (AGP) files. It does not currently handle Tiling Path Files (TPF). It does NOT handle cases where component sizes have changed in new AGP. It only handles changes in gap sizes and component flips. GFF features than span scaffolds or map to gaps or map outside the chr or GFF strand = 0 are not handled and written out to errors.gff3 and log to error.messages files. 

=head2 Methods

=over

=item C<reordered_coordinates_AGP ( $agp_old, $agp_new)>

Returns a hash of updated coordinates. Uses the component names to match positions. Will need mapping function if component names are different. Names should be same in case of accessions for contigs.

=cut

sub reordered_coordinates_AGP {
	my $self    = shift;
	my $agp_old = shift;
	my $agp_new = shift;
	my ( %coordinates, %obj_old_start, %obj_old_end, %comp_old_or );

	#get coords from old AGP
	while ( my $agp_line = $agp_old->get_next_agp_line() ) {

		next if ( ref($agp_line) eq 'Bio::GenomeUpdate::AGP::AGPGapLine' );

		my $agp_line_comp          = $agp_line->get_component_id();
		my $agp_line_obj_old_start = $agp_line->get_object_begin();
		my $agp_line_obj_old_end   = $agp_line->get_object_end();
		my $agp_line_comp_old_or   = $agp_line->get_orientation();

		$obj_old_start{$agp_line_comp} = $agp_line_obj_old_start;
		$obj_old_end{$agp_line_comp}   = $agp_line_obj_old_end;

		#set 0 as + for comparison with new AGP
		if ( looks_like_number($agp_line_comp_old_or) ) {
			if ( $agp_line_comp_old_or == 0 ) {
				$comp_old_or{$agp_line_comp} = '+';
			}
		} elsif (    ( $agp_line_comp_old_or eq '+' )
				  || ( $agp_line_comp_old_or eq '-' ) )
		{
			$comp_old_or{$agp_line_comp} = $agp_line_comp_old_or;
		} elsif (    ( $agp_line_comp_old_or eq '?' )
				  || ( $agp_line_comp_old_or eq 'na' ) )
		{
			$comp_old_or{$agp_line_comp} = '+';
		}

		for my $base ( $agp_line_obj_old_start .. $agp_line_obj_old_end ) {
			die "$agp_line_comp already has $base recorded, aborting.. "
			  if exists $coordinates{$base};
			$coordinates{$base} = 'X'
			  ; #use X to post-check for unordered coords, will cause datatype error
		}
	}

	#get coords from new AGP
	while ( my $agp_line = $agp_new->get_next_agp_line() ) {

		next if ( ref($agp_line) eq 'Bio::GenomeUpdate::AGP::AGPGapLine' );

		my $agp_line_comp          = $agp_line->get_component_id();
		my $agp_line_obj_new_start = $agp_line->get_object_begin();
		my $agp_line_obj_new_end   = $agp_line->get_object_end();
		my $agp_line_comp_new_or;

		#set 0,na,? as + for comparison with old AGP
		if ( looks_like_number( $agp_line->get_orientation() ) ) {
			if ( $agp_line->get_orientation() == 0 ) {
				$agp_line_comp_new_or = '+';
			}
		} elsif (    ( $agp_line->get_orientation() eq '+' )
				  || ( $agp_line->get_orientation() eq '-' ) )
		{
			$agp_line_comp_new_or = $agp_line->get_orientation();
		} elsif (    ( $agp_line->get_orientation() eq '?' )
				  || ( $agp_line->get_orientation() eq 'na' ) )
		{
			$agp_line_comp_new_or = '+';
		}

		#error check
		die "$agp_line_comp not found in old AGP, aborting.."
		  if ( !exists $obj_old_start{$agp_line_comp} );

		#same position and strand
		if (    ( $obj_old_start{$agp_line_comp} == $agp_line_obj_new_start )
			 && ( $obj_old_end{$agp_line_comp} == $agp_line_obj_new_end )
			 && ( $comp_old_or{$agp_line_comp} eq $agp_line_comp_new_or ) )
		{

			for my $base (
				$obj_old_start{$agp_line_comp} .. $obj_old_end{$agp_line_comp} )
			{
				$coordinates{$base} = $base;
			}
		}

		#same strand, moved downstream
		elsif (    ( $obj_old_start{$agp_line_comp} < $agp_line_obj_new_start )
				&& ( $obj_old_end{$agp_line_comp} < $agp_line_obj_new_end )
				&& ( $comp_old_or{$agp_line_comp} eq $agp_line_comp_new_or ) )
		{

			for my $base (
				$obj_old_start{$agp_line_comp} .. $obj_old_end{$agp_line_comp} )
			{
				$coordinates{$base} = $base +
				  ( $agp_line_obj_new_start - $obj_old_start{$agp_line_comp} );
			}
		}

		#same strand, moved upstream
		elsif (    ( $obj_old_start{$agp_line_comp} > $agp_line_obj_new_start )
				&& ( $obj_old_end{$agp_line_comp} > $agp_line_obj_new_end )
				&& ( $comp_old_or{$agp_line_comp} eq $agp_line_comp_new_or ) )
		{

			for my $base (
				$obj_old_start{$agp_line_comp} .. $obj_old_end{$agp_line_comp} )
			{
				$coordinates{$base} = $base -
				  ( $obj_old_start{$agp_line_comp} - $agp_line_obj_new_start );
			}
		}

		#diff strand, flipped, old start = new end
		elsif (    ( $obj_old_start{$agp_line_comp} == $agp_line_obj_new_start )
				&& ( $obj_old_end{$agp_line_comp} == $agp_line_obj_new_end )
				&& ( $comp_old_or{$agp_line_comp} ne $agp_line_comp_new_or ) )
		{

			my $counter = 0;
			for my $base (
				$obj_old_start{$agp_line_comp} .. $obj_old_end{$agp_line_comp} )
			{
				$coordinates{$base} = $agp_line_obj_new_end - $counter;
				$counter++;
			}

			#err check
			die "Problem in assigning coords for flipped $agp_line_comp"
			  if ( ( $counter - 1 ) !=
				   $agp_line_obj_new_end - $agp_line_obj_new_start );
		}

		#diff strand, start, end
		elsif (    ( $obj_old_start{$agp_line_comp} != $agp_line_obj_new_start )
				&& ( $obj_old_end{$agp_line_comp} != $agp_line_obj_new_end )
				&& ( $comp_old_or{$agp_line_comp} ne $agp_line_comp_new_or ) )
		{

			my $counter = 0;
			for my $base (
				$obj_old_start{$agp_line_comp} .. $obj_old_end{$agp_line_comp} )
			{
				$coordinates{$base} = $agp_line_obj_new_end - $counter;
				$counter++;
			}

			#err check
			die "Problem in assigning coords for flipped/moved $agp_line_comp"
			  if ( ( $counter - 1 ) !=
				   $agp_line_obj_new_end - $agp_line_obj_new_start );
		} else {
			die "This should not happen!";
		}
	}
	return %coordinates;
}

=item C<flipped_coordinates_AGP ( $agp_old, $agp_new)>

Returns a hash of coordinates that are flipped(0 or 1). Uses the component names to match positions. Will need mapping function if component names are different. Names should be same in case of accessions for contigs.

=cut

sub flipped_coordinates_AGP {
	my $self    = shift;
	my $agp_old = shift;
	my $agp_new = shift;
	my ( %flipped, %obj_old_start, %obj_old_end, %comp_old_or );

	#get coords from old AGP
	while ( my $agp_line = $agp_old->get_next_agp_line() ) {

		next if ( ref($agp_line) eq 'Bio::GenomeUpdate::AGP::AGPGapLine' );

		my $agp_line_comp          = $agp_line->get_component_id();
		my $agp_line_obj_old_start = $agp_line->get_object_begin();
		my $agp_line_obj_old_end   = $agp_line->get_object_end();
		my $agp_line_comp_old_or   = $agp_line->get_orientation();

		$obj_old_start{$agp_line_comp} = $agp_line_obj_old_start;
		$obj_old_end{$agp_line_comp}   = $agp_line_obj_old_end;
		$comp_old_or{$agp_line_comp}   = $agp_line_comp_old_or;

		for my $base ( $agp_line_obj_old_start .. $agp_line_obj_old_end ) {
			die "$agp_line_comp already has $base recorded, aborting.. "
			  if exists $flipped{$base};
			$flipped{$base} = 0;    #use 0 for not flipped
		}
	}

	#get coords from new AGP
	while ( my $agp_line = $agp_new->get_next_agp_line() ) {

		next if ( ref($agp_line) eq 'Bio::GenomeUpdate::AGP::AGPGapLine' );

		my $agp_line_comp          = $agp_line->get_component_id();
		my $agp_line_obj_new_start = $agp_line->get_object_begin();
		my $agp_line_obj_new_end   = $agp_line->get_object_end();
		my $agp_line_comp_new_or   = $agp_line->get_orientation();

		#error check
		die "$agp_line_comp not found in old AGP, aborting.."
		  if ( !exists $obj_old_start{$agp_line_comp} );

		#diff strand, flipped, old start = new end
		if (    ( $obj_old_start{$agp_line_comp} == $agp_line_obj_new_start )
			 && ( $obj_old_end{$agp_line_comp} == $agp_line_obj_new_end )
			 && ( $comp_old_or{$agp_line_comp} ne $agp_line_comp_new_or ) )
		{

			for my $base (
				$obj_old_start{$agp_line_comp} .. $obj_old_end{$agp_line_comp} )
			{
				$flipped{$base} = 1;
			}
		}

		#diff strand, start, end
		elsif (    ( $obj_old_start{$agp_line_comp} != $agp_line_obj_new_start )
				&& ( $obj_old_end{$agp_line_comp} != $agp_line_obj_new_end )
				&& ( $comp_old_or{$agp_line_comp} ne $agp_line_comp_new_or ) )
		{

			for my $base (
				$obj_old_start{$agp_line_comp} .. $obj_old_end{$agp_line_comp} )
			{
				$flipped{$base} = 1;
			}
		} else {

			#all cases where component did not flip
		}
	}
	return %flipped;
}

=item C<updated_coordinates_strand_AGP ( $start, $end, $strand, $agp_old, $agp_new, $gff_file_name)>

Returns new coordinates(int) and strand(char) wrt new AGP. Sets 0 to + as strand in old AGP for comparison purposes. Uses the component names to match positions. Will need mapping function if component names are different. Names should be same in case of accessions for contigs. It will return all 0's if GFF feature spans scaffolds and all 1's if feature maps to a gap or outside the chr and 100 for all other errors. 

=cut

sub updated_coordinates_strand_AGP {
	my $self          = shift;
	my $start         = shift;
	my $end           = shift;
	my $strand        = shift;
	my $agp_old       = shift;
	my $agp_new       = shift;
	my $gff_file_name = shift;
	my (
		 %obj_old_start, %obj_old_end, %comp_old_or,
		 %obj_new_start, %obj_new_end, %comp_new_or
	);
	my ( $nstart, $nend, $nstrand );
	my $errors = '';

	#reset current line number if already processed once
	$agp_old->set_current_agp_line_number(1);
	$agp_new->set_current_agp_line_number(1);

	#get coords from old AGP
	while ( my $agp_line = $agp_old->get_next_agp_line() ) {

		next if ( ref($agp_line) eq 'Bio::GenomeUpdate::AGP::AGPGapLine' );

		my $agp_line_comp          = $agp_line->get_component_id();
		my $agp_line_obj_old_start = $agp_line->get_object_begin();
		my $agp_line_obj_old_end   = $agp_line->get_object_end();
		my $agp_line_comp_old_or   = $agp_line->get_orientation();

		$obj_old_start{$agp_line_comp} = $agp_line_obj_old_start;
		$obj_old_end{$agp_line_comp}   = $agp_line_obj_old_end;

		#set 0,?,na as + for comparison with new AGP
		if ( looks_like_number($agp_line_comp_old_or) ) {
			if ( $agp_line_comp_old_or == 0 ) {
				$comp_old_or{$agp_line_comp} = '+';
			}
		} elsif (    ( $agp_line_comp_old_or eq '+' )
				  || ( $agp_line_comp_old_or eq '-' ) )
		{
			$comp_old_or{$agp_line_comp} = $agp_line_comp_old_or;
		} elsif (    ( $agp_line_comp_old_or eq '?' )
				  || ( $agp_line_comp_old_or eq 'na' ) )
		{
			$comp_old_or{$agp_line_comp} = '+';
		}
	}

	#get coords from new AGP
	while ( my $agp_line = $agp_new->get_next_agp_line() ) {

		next if ( ref($agp_line) eq 'Bio::GenomeUpdate::AGP::AGPGapLine' );

		my $agp_line_comp          = $agp_line->get_component_id();
		my $agp_line_obj_new_start = $agp_line->get_object_begin();
		my $agp_line_obj_new_end   = $agp_line->get_object_end();
		my $agp_line_comp_new_or   = $agp_line->get_orientation();

		$obj_new_start{$agp_line_comp} = $agp_line_obj_new_start;
		$obj_new_end{$agp_line_comp}   = $agp_line_obj_new_end;

		#set 0,na,? as + for comparison with old AGP
		if ( looks_like_number($agp_line_comp_new_or) ) {
			if ( $agp_line_comp_new_or == 0 ) {
				$comp_new_or{$agp_line_comp} = '+';
			}
		} elsif (    ( $agp_line_comp_new_or eq '+' )
				  || ( $agp_line_comp_new_or eq '-' ) )
		{
			$comp_new_or{$agp_line_comp} = $agp_line_comp_new_or;
		} elsif (    ( $agp_line_comp_new_or eq '?' )
				  || ( $agp_line_comp_new_or eq 'na' ) )
		{
			$comp_new_or{$agp_line_comp} = '+';
		}
	}

	#get component from old AGP
	my $start_component = $self->get_component_AGP( $start, $agp_old );
	my $end_component   = $self->get_component_AGP( $end,   $agp_old );

	#in case start/end base was in gap or outside chr
	if ( ( $start_component eq 'NA' ) || ( $end_component eq 'NA' ) ) {
		print STDERR "Component not found.\nStart: ", $start, ' Component: ',
		  $start_component, "\nEnd: ", $end, ' Component: ',
		  $end_component, "\n";
		$errors .=
		    "Component not found.\nStart: " 
		  . $start
		  . ' Component: '
		  . $start_component
		  . "\nEnd: "
		  . $end
		  . ' Component: '
		  . $end_component . "\n";
		$nstart  = 1;
		$nend    = 1;
		$nstrand = 1;
	}

	#Require start component == end component
	elsif ( $start_component ne $end_component ) {
		print STDERR "Diff component for start and stop.\nStart: ", $start,
		  ' Component: ',
		  $start_component, "\nEnd: ", $end, ' Component: ',
		  $end_component, "\n";
		$errors .=
		    "Diff component for start and stop.\nStart: " 
		  . $start
		  . ' Component: '
		  . $start_component
		  . "\nEnd: "
		  . $end
		  . ' Component: '
		  . $end_component . "\n";
		$nstart  = 0;
		$nend    = 0;
		$nstrand = 0;
	}

	#same position and strand
	elsif (
		 (
		   $obj_old_start{$start_component} == $obj_new_start{$start_component}
		 )
		 && ( $obj_old_end{$start_component} == $obj_new_end{$start_component} )
		 && ( $comp_old_or{$start_component} eq $comp_new_or{$start_component} )
	  )
	{
		$nstart  = $start;
		$nend    = $end;
		$nstrand = $strand;
	}

	#same strand, moved downstream
	elsif (
		 ( $obj_old_start{$start_component} < $obj_new_start{$start_component} )
		 && ( $obj_old_end{$start_component} < $obj_new_end{$start_component} )
		 && ( $comp_old_or{$start_component} eq $comp_new_or{$start_component} )
	  )
	{
		$nstart =
		  $start +
		  ( $obj_new_start{$start_component} -
			$obj_old_start{$start_component} );
		$nend =
		  $end +
		  ( $obj_new_start{$start_component} -
			$obj_old_start{$start_component} );
		$nstrand = $strand;
	}

	#same strand, moved upstream
	elsif (
		 ( $obj_old_start{$start_component} > $obj_new_start{$start_component} )
		 && ( $obj_old_end{$start_component} > $obj_new_end{$start_component} )
		 && ( $comp_old_or{$start_component} eq $comp_new_or{$start_component} )
	  )
	{
		$nstart =
		  $start -
		  ( $obj_old_start{$start_component} -
			$obj_new_start{$start_component} );
		$nend =
		  $end -
		  ( $obj_old_start{$start_component} -
			$obj_new_start{$start_component} );
		$nstrand = $strand;
	}

	#diff strand i.e. flipped
	elsif (
		  ( $comp_old_or{$start_component} ne $comp_new_or{$start_component} ) )
	{
		if ( $strand eq '+' ) {
			$nstart =
			  $obj_new_end{$start_component} -
			  ( ( $end - $start ) +
				( $start - $obj_old_start{$start_component} ) );
			$nend =
			  $obj_new_end{$start_component} -
			  ( $start - $obj_old_start{$start_component} );
			$nstrand = '-';
		} elsif ( $strand eq '-' ) {
			$nstart =
			  $obj_new_start{$start_component} +
			  ( ( $obj_old_end{$start_component} - $end ) );
			$nend =
			  $obj_new_start{$start_component} +
			  ( ( $obj_old_end{$start_component} - $end ) + ( $end - $start ) );
			$nstrand = '+';
		}
	} else {
		die "This should not happen!";
	}

	#For all other errors like GFF feature strand = 0, can have 9 error codes
	if (    ( !defined($nstart) )
		 || ( !defined($nend) )
		 || ( !defined($nstrand) ) )
	{
		print STDERR "Other error for GFF feature.\nStart: ", $start,
		  '  End: ', $end, '  Strand: ', $strand, "\n";
		$errors .=
		    "Other error for GFF feature.\nStart: " 
		  . $start
		  . '  End: '
		  . $end
		  . '  Strand: '
		  . $strand . "\n";
		$nstart  = 1;
		$nend    = 0;
		$nstrand = 0;
	}

	#print ERR messages
	if ( $errors ne '' ) {
		open( EGFF, ">>${gff_file_name}.error.messages" );
		print EGFF $errors;
		close(EGFF);
	}
	return ( $nstart, $nend, $nstrand );
}

=item C<get_component_AGP ( $base, $agp)>

Returns component name given a base and associated AGP file. It returns 'NA' when the base maps to a gap or outside the chr.

=cut

sub get_component_AGP {
	my $self = shift;
	my $base = shift;
	my $agp  = shift;
	my $component;
	my $found = 0;

	#ERR
	#print "get_component_AGP called for $base\n";

	my $current_agp_line_number = $agp->get_current_agp_line_number();
	$agp->set_current_agp_line_number(1);

	while ( my $agp_line = $agp->get_next_agp_line() ) {

		next if ( ref($agp_line) eq 'Bio::GenomeUpdate::AGP::AGPGapLine' );

		my $agp_line_comp      = $agp_line->get_component_id();
		my $agp_line_obj_start = $agp_line->get_object_begin();
		my $agp_line_obj_end   = $agp_line->get_object_end();

		if (    ( $base >= $agp_line_obj_start )
			 && ( $base <= $agp_line_obj_end ) )
		{
			$component = $agp_line_comp;
			$found     = 1;
			last;
		}
	}
	$agp->set_current_agp_line_number($current_agp_line_number);

	#return NA for base in gap or outside chr
	if ( $found == 0 ) {
		$component = 'NA';
	}
	return ($component);

	#die "No component found containing $base. Exiting..." if($found == 0);
}

###
1;    #do not remove
###

=back

=head1 LICENSE

    Same as Perl.

=head1 AUTHORS

    Surya Saha <suryasaha@cornell.edu , @SahaSurya>   

=cut
