package Utilities;

use strict;
use warnings;

use Moose;
use MooseX::FollowPBP;
use Moose::Util::TypeConstraints;
use Proc::ProcessTable;

=head1 NAME

    Utilities - Module with general utilities. 

=head1 SYNOPSIS

    my $variable = Utilities->new();

=head1 DESCRIPTION

    This class contains general utilities for performance metrics. 

=head2 Methods

=over

=item C<run_time ()>

Prints run time details.

=cut

sub run_time {
	my ( $user_t, $system_t, $cuser_t, $csystem_t );
	( $user_t, $system_t, $cuser_t, $csystem_t ) = times;
	print STDERR "System time for process: $system_t\n";
	print STDERR "User time for process: $user_t\n\n";
}

=item C<mem_used ()>

Prints runtime details.

=cut

sub mem_used {
	my ( $i, $t );
	$t = new Proc::ProcessTable;
	foreach my $got ( @{ $t->table } ) {
		next if not $got->pid eq $$;
		$i = $got->size;
	}
	print STDERR "Process id=", $$, "\n";
	print "Memory used(MB)=", $i / 1024 / 1024, "\n";
}

###
1;   #do not remove
###

=pod

=back

=head1 LICENSE

    Same as Perl.

=head1 AUTHORS

    Surya Saha <suryasaha@cornell.edu , @SahaSurya>   

=cut

