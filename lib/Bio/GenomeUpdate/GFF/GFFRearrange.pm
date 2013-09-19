package Bio::GenomeUpdate::GFF::GFFRearrange;
use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Bio::GenomeUpdate::AGP;
use Bio::GenomeUpdate::GFF;

=head1 NAME

    GFFRearrange - Generic Feature Format (GFF) Rearrange modifies the coordinates of a GFF file.

=head1 SYNOPSIS

    my $variable = Bio::GenomeUpdate::GFF::GFFRearrange->new();

=head1 DESCRIPTION

    This class modifies Generic Feature Format (GFF) coordinates using old and new Accessioned Golden Path (AGP) files. It does not currently handle Tiling Path Files (TPF). 

=head2 Methods

=over

=item C<reorder_coordinates_AGP ( $agp_old, $agp_new)>

Returns a hash of updated coordinates. Uses the component names to match positions. Will need mapping function if component names are different. Names should be same in case of accessions for contigs.

=cut

sub reorder_coordinates_AGP{
	my $self = shift;
	my $agp_old = shift;
	my $agp_new = shift;
	my (%coordinates, %obj_old_start, %obj_old_end, %comp_old_or);
	
	#get coords from old AGP
	while ( my $agp_line = $agp_old->get_next_agp_line()){
		my $agp_line_comp = $agp_line->get_component_id();
		my $agp_line_obj_old_start = $agp_line->get_object_begin();
		my $agp_line_obj_old_end = $agp_line->get_object_end();
		my $agp_line_comp_old_or = $agp_line->get_orientation();
		
		$obj_old_start{$agp_line_comp} = $agp_line_obj_old_start;
		$obj_old_end{$agp_line_comp} = $agp_line_obj_old_end;
		$comp_old_or{$agp_line_comp} = $agp_line_comp_old_or;
		
		for my $base ($agp_line_obj_old_start..$agp_line_obj_old_end){
			die "$agp_line_comp already has $base recorded, aborting.. " if exists $coordinates{$base};
			$coordinates{$base}='X';#use X to post-check for unordered coords, will cause datatype error 
		}
	}
	
	#get coords from new AGP
	while ( my $agp_line = $agp_new->get_next_agp_line()){
		my $agp_line_comp = $agp_line->get_component_id();
		my $agp_line_obj_new_start = $agp_line->get_object_begin();
		my $agp_line_obj_new_end = $agp_line->get_object_end();
		my $agp_line_comp_new_or = $agp_line->get_orientation();
		
		#error check
		die "$agp_line_comp not found in old AGP, aborting.." if (!exists $obj_old_start{$agp_line_comp});
		
		#same position and strand
		if ( ($obj_old_start{$agp_line_comp} == $agp_line_obj_new_start) &&
			 ($obj_old_end{$agp_line_comp} == $agp_line_obj_new_end) &&
			 ($comp_old_or{$agp_line_comp} == $agp_line_comp_new_or)){
			
			for my $base ($obj_old_start{$agp_line_comp}..$obj_old_end{$agp_line_comp}){
				$coordinates{$base} = $base;
			}
		}
		#same strand, moved downstream
		elsif ( ($obj_old_start{$agp_line_comp} < $agp_line_obj_new_start) &&
			 ($obj_old_end{$agp_line_comp} < $agp_line_obj_new_end) &&
			 ($comp_old_or{$agp_line_comp} == $agp_line_comp_new_or)){

			for my $base ($obj_old_start{$agp_line_comp}..$obj_old_end{$agp_line_comp}){
				$coordinates{$base} = $base + ($agp_line_obj_new_start - $obj_old_start{$agp_line_comp});
			}
		}
		#same strand, moved upstream
		elsif ( ($obj_old_start{$agp_line_comp} > $agp_line_obj_new_start) &&
			 ($obj_old_end{$agp_line_comp} > $agp_line_obj_new_end) &&
			 ($comp_old_or{$agp_line_comp} == $agp_line_comp_new_or)){

			for my $base ($obj_old_start{$agp_line_comp}..$obj_old_end{$agp_line_comp}){
				$coordinates{$base} = $base - ($obj_old_start{$agp_line_comp} - $agp_line_obj_new_start);
			}
		}
		#diff strand, flipped, old start = new end
		elsif ( ($obj_old_start{$agp_line_comp} == $agp_line_obj_new_start) &&
			 ($obj_old_end{$agp_line_comp} == $agp_line_obj_new_end) &&
			 ($comp_old_or{$agp_line_comp} != $agp_line_comp_new_or)){
			 
			my $counter = 0;
			for my $base ($obj_old_start{$agp_line_comp}..$obj_old_end{$agp_line_comp}){
				$coordinates{$base} = $agp_line_obj_new_end - $counter;
				$counter++;
			}
			#err check
			die "Problem in assigning coords for flipped $agp_line_comp" if ( $counter != $agp_line_obj_new_end - $agp_line_obj_new_start);
		}
		#diff strand, start, end
		elsif ( ($obj_old_start{$agp_line_comp} != $agp_line_obj_new_start) &&
			 ($obj_old_end{$agp_line_comp} != $agp_line_obj_new_end) &&
			 ($comp_old_or{$agp_line_comp} != $agp_line_comp_new_or)){
			 
			my $counter = 0;
			for my $base ($obj_old_start{$agp_line_comp}..$obj_old_end{$agp_line_comp}){
				$coordinates{$base} = $agp_line_obj_new_end - $counter;
				$counter++;
			}
			#err check
			die "Problem in assigning coords for flipped/moved $agp_line_comp" if ( $counter != $agp_line_obj_new_end - $agp_line_obj_new_start);
		}
		else{
			die "This should not happen!";
		}
	}
}

=item C<flipped_coordinates_AGP ( $agp_old, $agp_new)>

Returns a hash of coordinates that are flipped(0 or 1). Uses the component names to match positions. Will need mapping function if component names are different. Names should be same in case of accessions for contigs.

=cut
sub flipped_coordinates_AGP{
	my $self = shift;
	my $agp_old = shift;
	my $agp_new = shift;
	my (%flipped, %obj_old_start, %obj_old_end, %comp_old_or);
	
	#get coords from old AGP
	while ( my $agp_line = $agp_old->get_next_agp_line()){
		my $agp_line_comp = $agp_line->get_component_id();
		my $agp_line_obj_old_start = $agp_line->get_object_begin();
		my $agp_line_obj_old_end = $agp_line->get_object_end();
		my $agp_line_comp_old_or = $agp_line->get_orientation();
		
		$obj_old_start{$agp_line_comp} = $agp_line_obj_old_start;
		$obj_old_end{$agp_line_comp} = $agp_line_obj_old_end;
		$comp_old_or{$agp_line_comp} = $agp_line_comp_old_or;
		
		for my $base ($agp_line_obj_old_start..$agp_line_obj_old_end){
			die "$agp_line_comp already has $base recorded, aborting.. " if exists $flipped{$base};
			$flipped{$base}=0;#use 0 for not flipped 
		}
	}
	
	#get coords from new AGP
	while ( my $agp_line = $agp_new->get_next_agp_line()){
		my $agp_line_comp = $agp_line->get_component_id();
		my $agp_line_obj_new_start = $agp_line->get_object_begin();
		my $agp_line_obj_new_end = $agp_line->get_object_end();
		my $agp_line_comp_new_or = $agp_line->get_orientation();
		
		#error check
		die "$agp_line_comp not found in old AGP, aborting.." if (!exists $obj_old_start{$agp_line_comp});
		
		#diff strand, flipped, old start = new end
		if ( ($obj_old_start{$agp_line_comp} == $agp_line_obj_new_start) &&
			 ($obj_old_end{$agp_line_comp} == $agp_line_obj_new_end) &&
			 ($comp_old_or{$agp_line_comp} != $agp_line_comp_new_or)){
			 
			for my $base ($obj_old_start{$agp_line_comp}..$obj_old_end{$agp_line_comp}){
				$flipped{$base} = 1;
			}
		}
		#diff strand, start, end
		elsif ( ($obj_old_start{$agp_line_comp} != $agp_line_obj_new_start) &&
			 ($obj_old_end{$agp_line_comp} != $agp_line_obj_new_end) &&
			 ($comp_old_or{$agp_line_comp} != $agp_line_comp_new_or)){
			 
			for my $base ($obj_old_start{$agp_line_comp}..$obj_old_end{$agp_line_comp}){
				$flipped{$base} = 1;
			}
		}
		else{
			die "This should not happen!";
		}
	}
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
