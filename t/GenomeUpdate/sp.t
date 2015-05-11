#!/usr/bin/perl

=head1 NAME

tpf.t
A test for Bio::GenomeUpdate::SP class

=cut

=head1 SYNOPSIS

perl sp.t


=head1 DESCRIPTION

Creating new SPLines and adding to SP object. Checking if returned object is correctly formatted.

=head2 Author

Surya Saha        <suryasaha@cornell.edu>

=cut

use strict;
use warnings;
use autodie;

use Test::More tests => 1;
BEGIN {use_ok( 'Bio::GenomeUpdate::SP' ); }
require_ok( 'Bio::GenomeUpdate::SP::SPLine' );

#create a sp line
ok(my $sp_line1 = Bio::GenomeUpdate::SP::SPLine->new( 
	chromosome => 'Un',
	accession_prefix => 'accession1',
	accession_suffix => 'accession2',
	accession_prefix_orientation => '+',
	accession_suffix_orientation => '-',
	accession_prefix_last_base => 100,
	accession_suffix_first_base => 200,
	comment => 'test line 1 comment'
));
#create another sp line
ok(my $sp_line2 = Bio::GenomeUpdate::SP::SPLine->new( 
	chromosome => 1,
	accession_prefix => 'accession3',
	accession_suffix => 'accession4',
	accession_prefix_orientation => '-',
	accession_suffix_orientation => '-',
	accession_prefix_last_base => 300,
	accession_suffix_first_base => 100,
	comment => 'test line 1 comment'
));

#get various values and test


#add to SP object


#print out SP object and compare to string


#modify values for #sp_line1 and add to SP object


#print out SP object and compare to string
