package Bio::GenomeUpdate::AlignmentCoordsGroup;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Bio::GenomeUpdate::AlignmentCoords;

=head1 NAME

    AlignmentCoordsGroup - Stores an array of AlignmentCoords objects and reports information

=head1 SYNOPSIS

    my $variable = AlignmentCoordsGroup->new();

=head1 DESCRIPTION

    This package performs calculations on an array of objects containing alignment coordinates (e.g., from Nucmer/delta-filter/show-coords output) for a single query sequence (e.g., from a BAC) relative to multiple reference sequences (e.g., chromosome pseudomolecules).

=cut

has 'array_of_alignment_coords' => (
						  isa => 'ArrayRef[Bio::GenomeUpdate::AlignmentCoords]',
						  is  => 'rw',
						  predicate => 'has_array_of_alignment_coords'
);

=head2 Methods

=over

=item C<get_reference_ids>

Returns a hash of reference IDs and the corresponding number of alignment clusters (hits) for each reference ID.

=cut

sub get_reference_ids_and_hit_counts {
	my $self = shift;
	my %seen_ids_and_counts;
	foreach my $align_coords ( @{ $self->get_array_of_alignment_coords() } ) {
		if ( $seen_ids_and_counts{ $align_coords->get_reference_id() } ) {
			$seen_ids_and_counts{ $align_coords->get_reference_id() }++;
		}
		else {
			$seen_ids_and_counts{ $align_coords->get_reference_id() } = 1;
		}
	}
	return %seen_ids_and_counts;
}

=item C<get_count_of_reference_sequence_ids>

Returns the number of reference sequence IDs that have alignment clusters (hits) with the query.

=cut

sub get_count_of_reference_sequence_ids {
	my $self    = shift;
	my %ref_ids = $self->get_reference_ids_and_hit_counts();
	return scalar( keys %ref_ids );
}

=item C<includes_reference_id ( $reference_id )>

Returns 1 if specified reference ID is included in the group of alignment clusters, and returns 0 if not.

=cut

sub includes_reference_id {
	my $self            = shift;
	my $ref_id_to_check = shift;
	my %ref_ids         = $self->get_reference_ids_and_hit_counts();
	if ( $ref_ids{$ref_id_to_check} ) {
		return 1;
	}
	else {
		return 0;
	}
}

=item C<order_alignment_clusters_on_each_reference_sequence>

Calculates and returns a hash of reference IDs and corresponding arrays of alignment clusters sorted by reference sequence coordinates.

=cut

sub order_alignment_clusters_on_each_reference_sequence {
	my $self = shift;
	my %clusters;
	foreach my $align_coords ( @{ $self->get_array_of_alignment_coords() } ) {
		if ( !$clusters{ $align_coords->get_reference_id() } ) {

			#if a key in clusters doesn't exist, create one as an empty array
			$clusters{ $align_coords->get_reference_id() } = ();
		}
		push(
			  @{ $clusters{ $align_coords->get_reference_id() } },
			  $align_coords
		);

	}
	foreach my $ref_seq_key ( keys %clusters ) {
		my @sorted_clusters = sort {
			$a->get_reference_start_coord() <=> $b->get_reference_start_coord()
		} @{ $clusters{$ref_seq_key} };

		#print "$ref_seq_key:";
		#@{$clusters{$ref_seq_key}} = @sorted_clusters;
		$clusters{$ref_seq_key} = \@sorted_clusters;
		foreach my $c (@sorted_clusters) {

			#print $c->get_reference_start_coord().",";
		}

		#print "\n";
	}
	return %clusters;
}

=item C<group_alignment_clusters ( $gap_allowed )>

Calculates and returns a hash of reference IDs and corresponding arrays of arrays of alignment clusters grouped by proximity on the reference sequence.  The proximity grouping is done using the specified length of an allowed gap between aligned clusters in the reference sequence ($gap allowed). Not implemented: an adjustment based on a percent ($gap_percent) of the length of the gap between the aligned clusters on the query sequence. 

=cut

sub group_alignment_clusters {
	my $self        = shift;
	my $gap_allowed = shift;

	#hash of
	my %clusters =
	  $self->order_alignment_clusters_on_each_reference_sequence()
	  ; #Calculates and returns a hash of reference IDs and corresponding arrays of alignment clusters sorted by reference sequence coordinates.
	my %cluster_groups;
	foreach my $ref_seq_key ( keys %clusters ) {
		$cluster_groups{$ref_seq_key} = ();
		my @group;
		my $end_of_previous_ref;
		my $end_of_previous_query;
		foreach my $align_coords ( @{ $clusters{$ref_seq_key} } ) {
			if ($end_of_previous_ref) {

				#my $start_of_curr_query;
				#$start_of_curr_query = $align_coords->get_query_start_coord();
				if (
					 (
					   $align_coords->get_reference_start_coord() -
					   ( $end_of_previous_ref + 1 )
					 ) < $gap_allowed
				  )
				{

#	  print STDERR $ref_seq_key.",".$align_coords->get_reference_start_coord().",",$end_of_previous_ref,",",$gap_allowed,"\n";
					push( @group, $align_coords );
				}
				else {
					push( @{ $cluster_groups{$ref_seq_key} }, [@group] );
					@group = ();
					push( @group, $align_coords );
				}
			}
			else {
				push( @group, $align_coords );
			}
			$end_of_previous_ref = $align_coords->get_reference_end_coord();

			#$end_of_previous_query = $align_coords->get_query_end_coord();
		}
		push( @{ $cluster_groups{$ref_seq_key} }, [@group] )
		  ;    # store last cluster group
	}
	return %cluster_groups;
}

=item C<get_cluster_groups_sorted_by_length ( $gap_allowed )>

Calculates and returns an array of arrays containing the proximity-grouped alignment clusters sorted by longest to shortest length of non-overlapping sequence covered by alignment clusters.  The proximity grouping is done using the specified length of an allowed gap between aligned clusters in the reference sequence ($gap allowed). Not implemented: plus an adjustment based on a percent ($gap_percent) of the length of the gap between the aligned clusters on the query sequence. 

=cut

sub get_cluster_groups_sorted_by_length {
	my $self           = shift;
	my $gap_allowed    = shift;
	my %cluster_groups = $self->group_alignment_clusters($gap_allowed);
	my @groups;

	foreach my $ref_seq_key ( keys %cluster_groups ) {
		foreach my $align_coords_groups ( @{ $cluster_groups{$ref_seq_key} } ) {
			my @align_group = @$align_coords_groups;
			push( @groups, [@align_group] );
		}
	}
	my @sorted_groups = sort cluster_group_length_sort @groups;
	return @sorted_groups;

	sub cluster_group_length_sort {
		my $a_length = calculate_length( [ @{$a} ] );
		my $b_length = calculate_length( [ @{$b} ] );
		if ( $a_length > $b_length ) {
			return -1;

		}
		elsif ( $a_length == $b_length ) {
			return 0;

		}
		elsif ( $a_length < $b_length ) {
			return 1;
		}
		else {
			return 0;
		}

		sub calculate_length {
			my $cluster_ref   = shift;
			my @cluster_group = @{$cluster_ref};
			my $length        = 0;
			my $prev_end_coord;
			foreach my $align_coords (@cluster_group) {
				$length +=
				  ( ( $align_coords->get_reference_end_coord() + 1 ) -
					$align_coords->get_reference_start_coord() );
				if ($prev_end_coord) {
					if ( $align_coords->get_reference_start_coord() <=
						 $prev_end_coord )
					{
						$length -=
						  ( ( $prev_end_coord + 1 ) -
							$align_coords->get_reference_start_coord() )
						  ; #subtract overlap between alignment clusters from total alignment length
					}
				}
				$prev_end_coord = $align_coords->get_reference_end_coord();
			}
			return $length;
		}
	}
}

=item C<get_id_coords_and_direction_of_longest_alignment_cluster_group ( $gap_allowed )>

Returns IDs, start and end coordinates, total aligned sequence length, co-linear alignment flag and direction for the longest proximity-grouped alignment clusters sorted by longest to shortest length of non-overlapping sequence covered by alignment clusters.  The proximity grouping is done using the specified length of an allowed gap between aligned clusters in the reference sequence ($gap allowed). 

=cut

sub get_id_coords_and_direction_of_longest_alignment_cluster_group {
	my $self        = shift;
	my $gap_allowed = shift;

#returns an array of arrays containing the
#proximity-grouped alignment clusters sorted by longest to shortest
#length of non-overlapping sequence covered by alignment clusters.
#aligned_coords are ordered acc to ref chr. see order_alignment_clusters_on_each_reference_sequence()
	my @returnclusters = $self->get_cluster_groups_sorted_by_length($gap_allowed);
	my $first_group    = shift(@returnclusters);    #first is largest

	#my @second_largest_cluster_group;
	#if (defined(@{$returnclusters[1]})){
	#@second_largest_cluster_group = @{$returnclusters[1]};
	#}
	my @largest_cluster_group = @{$first_group};
	my $ref_start = $largest_cluster_group[0]->get_reference_start_coord();
	my $ref_end   = $largest_cluster_group[-1]->get_reference_end_coord();
	my $ref_id    = $largest_cluster_group[0]->get_reference_id();
	my $query_id  = $largest_cluster_group[0]->get_query_id();
	my $query_start;
	my $query_end;
	my $sequence_aligned_in_clusters = 0;
	my $prev_end_coord;
	my $overlap_allowed = 100;    #between tiles in cluster
	my $is_overlapping  = 0;

	foreach my $align_coords (@largest_cluster_group)
	{                             #parse through alignments in cluster
		$sequence_aligned_in_clusters +=
		  ( ( $align_coords->get_reference_end_coord() + 1 ) -
			$align_coords->get_reference_start_coord() );
		if ($prev_end_coord) {    #if this is not the first hit
			if ( $align_coords->get_reference_start_coord() <= $prev_end_coord )
			{
				$sequence_aligned_in_clusters -=
				  ( ( $prev_end_coord + 1 ) -
					$align_coords->get_reference_start_coord() )
				  ; #subtract overlap between alignment clusters from total alignment length
			}
			if (
				(
				   $align_coords->get_reference_start_coord() + $overlap_allowed
				) <= $prev_end_coord
			  )
			{
				$is_overlapping = 1;
			}
		}
		$prev_end_coord = $align_coords->get_reference_end_coord();
	}

	if ( $largest_cluster_group[0]->get_query_start_coord() <
		 $largest_cluster_group[-1]->get_query_end_coord() )
	{
		$query_start = $largest_cluster_group[0]->get_query_start_coord();
		$query_end   = $largest_cluster_group[-1]->get_query_end_coord();
	}
	else {
		$query_end   = $largest_cluster_group[0]->get_query_start_coord();
		$query_start = $largest_cluster_group[-1]->get_query_end_coord();
	}

	# if ($largest_cluster_group[0]->get_direction() == 1) {
	#   $query_start = $largest_cluster_group[0]->get_query_start_coord();
	# } else {
	#   $query_end = $largest_cluster_group[0]->get_query_start_coord();
	# }
	# if ($largest_cluster_group[-1]->get_direction() == 1) {
	#   $query_end = $largest_cluster_group[-1]->get_query_end_coord();
	# } else {
	#   $query_start = $largest_cluster_group[-1]->get_query_end_coord();
	# }

	my $alternates;

	# if exists @returnclusters?
	foreach my $cluster_group (@returnclusters) {
		my @current_cluster_group               = @{$cluster_group};
		my $sequence_aligned_in_current_cluster = 0;
		my $curr_ref_start                      =
		  $current_cluster_group[0]->get_reference_start_coord();
		my $curr_ref_end =
		  $current_cluster_group[-1]->get_reference_end_coord();
		my $curr_ref_id = $current_cluster_group[0]->get_reference_id();
		my $prev_end_coord_of_current_cluster;
		foreach my $align_coords (@current_cluster_group) {
			$sequence_aligned_in_current_cluster +=
			  ( ( $align_coords->get_reference_end_coord() + 1 ) -
				$align_coords->get_reference_start_coord() );
			if ($prev_end_coord_of_current_cluster) {
				if ( $align_coords->get_reference_start_coord() <=
					 $prev_end_coord_of_current_cluster )
				{
					$sequence_aligned_in_current_cluster -=
					  ( ( $prev_end_coord_of_current_cluster + 1 ) -
						$align_coords->get_reference_start_coord() )
					  ; #subtract overlap between alignment clusters from total alignment length
				}
			}
			$prev_end_coord_of_current_cluster =
			  $align_coords->get_reference_end_coord();
		}
		if ( defined($alternates) ) {
			$alternates =
			    $alternates . ", "
			  . $curr_ref_id . ":"
			  . $sequence_aligned_in_current_cluster . ":"
			  . $curr_ref_start . ":"
			  . $curr_ref_end;
		}
		else {
			$alternates =
			    $curr_ref_id . ":"
			  . $sequence_aligned_in_current_cluster . ":"
			  . $curr_ref_start . ":"
			  . $curr_ref_end;
		}
	}

	if ( !defined($alternates) ) {
		$alternates = "";
	}

	my $size_of_next_largest_match = 0;

	my $second_group = shift(@returnclusters);
	if ($second_group) {
		my @second_largest_cluster_group = @{$second_group};
		$size_of_next_largest_match =
		  ( $second_largest_cluster_group[-1]->get_reference_end_coord() -
			$second_largest_cluster_group[0]->get_reference_start_coord() );
	}

	my $direction_check =
	  check_direction_of_cluster_group( \@largest_cluster_group );
	my $colinear_order_check =
	  check_colinear_order_of_cluster_group( \@largest_cluster_group );
	return (
			 $ref_id,                       $query_id,
			 $ref_start,                    $ref_end,
			 $query_start,                  $query_end,
			 $sequence_aligned_in_clusters, $direction_check,
			 $colinear_order_check,         $is_overlapping,
			 $size_of_next_largest_match,   $alternates
	);

=item C<check_direction_of_cluster_group ( @largest_cluster_group )>

Returns 0 for if all align_coord's in cluster are not in same orientation.    

=cut

	sub check_direction_of_cluster_group {
		my $cluster_group = shift;
		my $direction;
		foreach my $align_coords (@$cluster_group) {
			if ($direction) {
				if ( $direction != $align_coords->get_direction() ) {
					print STDERR "Mixed orientation for "
					  . $align_coords->get_query_id() . "\n";
					return 0;
				}
			}
			else {
				$direction = $align_coords->get_direction();
			}
		}
		return $direction;
	}

=item C<check_colinear_order_of_cluster_group ( @largest_cluster_group )>

Returns 1 for error condition, i.e. if query start coordinate of a align_coord is out of order with prev one in the cluster. Note that alignments are ordered by ref start coordinate.   

=cut

	sub check_colinear_order_of_cluster_group {

		#aligned_coords are ordered acc to ref chr
		#see order_alignment_clusters_on_each_reference_sequence()
		my $cluster_group = shift;
		my ( $query_start, $colinear_order );
		$colinear_order = 0;
		foreach my $align_coords (@$cluster_group) {
			if ($query_start){
				#cluster is ordered acc to ref coords so checking if query starts are in order
				#order is affected by strand as nucmer records alignments on -ive strand in opposite order
				#presuming all query sequences are aligned in same orientation
				if ((($query_start >= $align_coords->get_query_start_coord())
					&& ($align_coords->get_query_strand() == 1 ))
					||($query_start <= $align_coords->get_query_start_coord())
					&& ($align_coords->get_query_strand() == -1 )) {
					#print STDERR "******** query_start: $query_start\talign_coords->get_query_start_coord(): ",$align_coords->get_query_start_coord(),"\n";
					print STDERR "Query start error for "
					  . $align_coords->get_query_id() . "\n";
					return 1;
				}
			}
			else {
				#print STDERR "\n******** align_coords->get_query_start_coord(): ",$align_coords->get_query_start_coord(),"\n";
				$query_start = $align_coords->get_query_start_coord();
			}
		}
		return $colinear_order;
	}

}

###
1;    #do not remove
###

=back

=head1 LICENSE

    Same as Perl.

=head1 AUTHORS

    Jeremy D. Edwards <jde22@cornell.edu>
    Surya Saha <suryasaha at cornell.edu, @SahaSurya>   

=cut

