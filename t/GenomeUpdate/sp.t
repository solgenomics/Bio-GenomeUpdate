#!/usr/bin/perl

=head1 NAME

sp.t - A test for Bio::GenomeUpdate::SP class

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

use Test::More;
use Test::Exception;

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
	comment => 'test line 1 comment some nonsense here and here too'
), 'created $sp_line1');
#create another sp line
ok(my $sp_line2 = Bio::GenomeUpdate::SP::SPLine->new( 
	chromosome => 1,
	accession_prefix => 'accession3',
	accession_suffix => 'accession4',
	accession_prefix_orientation => '-',
	accession_suffix_orientation => '-',
	accession_prefix_last_base => 300,
	accession_suffix_first_base => 100,
	comment => 'test line 2 comment some nonsense here and here too'
),'created $sp_line2');

#get various values 
is($sp_line1->get_chromosome(), 'Un', 'Got correct chr from $sp_line1');
is($sp_line1->get_accession_prefix(), 'accession1', 'Got correct accession_prefix from $sp_line1');
is($sp_line1->get_accession_suffix(), 'accession2', 'Got correct accession_suffix from $sp_line1');
is($sp_line1->get_accession_prefix_orientation(), '+', 'Got correct accession_prefix_orientation from $sp_line1');
is($sp_line1->get_accession_suffix_orientation(), '-', 'Got correct accession_suffix_orientation from $sp_line1');
is($sp_line1->get_accession_prefix_last_base(), 100, 'Got correct accession_prefix_last_base from $sp_line1');
is($sp_line1->get_accession_suffix_first_base(), 200, 'Got correct accession_suffix_first_base from $sp_line1');
is($sp_line1->get_comment(), 'test line 1 comment some nonsense here and here too', 'Got correct comment from $sp_line1');

#test for type checks for SP line
throws_ok{$sp_line1->set_chromosome(13)} qr/not a valid chromosome number/, "invalid chr exception caught"; 
throws_ok{$sp_line1->set_accession_prefix_orientation('3')} qr/not a valid orientation type/, "invalid orientation exception caught";
throws_ok{$sp_line1->set_accession_prefix_last_base(-1)} qr/not a positive coordinate/, "negative coordinate exception caught";
throws_ok{$sp_line1->set_accession_suffix_first_base('A')} qr/not a positive coordinate/, "character coordinate exception caught";
throws_ok{$sp_line1->set_comment('short comment')} qr/shorter than the minimum length/, "short comment exception caught";

#add to SP object
ok(my $sp = Bio::GenomeUpdate::SP->new(
	taxid => '001',
	assembly_group => 'TGP',
	assembly_unit => 'Primary',
	tpf_type => 'chromosome',
), 'Created sp object');
ok ($sp->add_line_to_end($sp_line1), 'added $sp_line1');
ok ($sp->add_line_to_end($sp_line2), 'added $sp_line2');

#test for type checks for SP object
throws_ok{$sp->set_tpf_type('blah blah')} qr/not a valid TPF type/, "invalid TPF type exception caught";

#print out SP object and compare to string
ok(my $sp_file_string = $sp->get_formatted_sp(), 'got switch point file content' );
my $expected_sp_file_string = q(001	TGP	Primary	Un	chromosome	accession1	accession2	+	-	100	200	test line 1 comment some nonsense here and here too
001	TGP	Primary	1	chromosome	accession3	accession4	-	-	300	100	test line 2 comment some nonsense here and here too
);
ok( $sp_file_string eq $expected_sp_file_string, 'switch point file content looks good ' );

#modify values for #sp_line1 and add to SP object
ok($sp_line1->set_accession_prefix_last_base(110), 'Set accession_prefix_last_base in $sp_line1');
ok($sp_line1->set_accession_suffix_first_base(210), 'Set accession_suffix_first_base in $sp_line1');
ok ($sp->add_line_to_end($sp_line1), 'added modified $sp_line1');

#print out SP object and compare to string
ok($sp_file_string = $sp->get_formatted_sp(), 'got switch point file content again' );
$expected_sp_file_string = q(001	TGP	Primary	Un	chromosome	accession1	accession2	+	-	110	210	test line 1 comment some nonsense here and here too
001	TGP	Primary	1	chromosome	accession3	accession4	-	-	300	100	test line 2 comment some nonsense here and here too
001	TGP	Primary	Un	chromosome	accession1	accession2	+	-	110	210	test line 1 comment some nonsense here and here too
); #$sp_line1 also changes as its a hashref
ok( $sp_file_string eq $expected_sp_file_string, 'modified switch point file content looks good' );

#no more tests
done_testing();