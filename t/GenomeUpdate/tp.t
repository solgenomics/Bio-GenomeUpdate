#!/usr/bin/perl

=head1 NAME

tp.t - A test for Bio::GenomeUpdate::TP class

=cut

=head1 SYNOPSIS

perl tp.t


=head1 DESCRIPTION

Creating new TPLines and adding to TP object. Checking if returned object is correctly formatted.

=head2 Author

Surya Saha        <suryasaha@cornell.edu>

=cut

use strict;
use warnings;
use autodie;

use Test::More;
use Test::Exception;

BEGIN {use_ok( 'Bio::GenomeUpdate::TP' ); }
require_ok( 'Bio::GenomeUpdate::TP::TPLine' );

#create a tp line
ok(my $tp_line1 = Bio::GenomeUpdate::TP::TPLine->new( 
	chromosome => 'Un',
	accession => 'accession1',
	accession_prefix_first_or_last_base => 100,
	trim_from_end => 'L',
	comment => 'test line 1 comment some nonsense here and here too'
), 'created $tp_line1');
#create another sp line
ok(my $tp_line2 = Bio::GenomeUpdate::TP::TPLine->new( 
	chromosome => 1,
	accession => 'accession2',
	accession_prefix_first_or_last_base => 200,
	trim_from_end => 'H',
	comment => 'test line 2 comment some nonsense here and here too'
), 'created $tp_line2');

#get various values 
is($tp_line1->get_chromosome(), 'Un', 'Got correct chr from $tp_line1');
is($tp_line1->get_accession(), 'accession1', 'Got correct accession from $tp_line1');
is($tp_line1->get_accession_prefix_first_or_last_base(), 100, 'Got correct accession_prefix_first_or_last_base from $tp_line1');
is($tp_line1->get_trim_from_end(), 'L', 'Got correct trim_from_end from $tp_line1');
is($tp_line1->get_comment(), 'test line 1 comment some nonsense here and here too', 'Got correct comment from $tp_line1');

#test for type checks
throws_ok{$tp_line1->set_chromosome(13)} qr/not a valid chromosome number/, "invalid int chr exception caught";
throws_ok{$tp_line1->set_chromosome('MyChr')} qr/not a valid chromosome number/, "invalid char chr exception caught";
throws_ok{$tp_line1->set_accession_prefix_first_or_last_base(-1)} qr/not a positive coordinate/, "invalid coordinate exception caught";
throws_ok{$tp_line1->set_trim_from_end(1)} qr/not a valid trim direction/, "invalid int trim_from_end exception caught";
throws_ok{$tp_line1->set_trim_from_end('A')} qr/not a valid trim direction/, "invalid char trim_from_end exception caught";
throws_ok{$tp_line1->set_comment('short comment')} qr/shorter than the minimum length/, "short comment exception caught";

#add to TP object
ok(my $tp = Bio::GenomeUpdate::TP->new(
	taxid => '001',
	assembly_group => 'TGP',
	assembly_unit => 'Primary',
	tpf_type => 'chromosome',
), 'Created tp object');
ok ($tp->add_line_to_end($tp_line1), 'added $tp_line1');
ok ($tp->add_line_to_end($tp_line2), 'added $tp_line2');

#test for type checks for TP object
throws_ok{$tp->set_tpf_type('blah blah')} qr/not a valid TPF type/, "invalid TPF type exception caught";

#print out TP object and compare to string
ok(my $tp_file_string = $tp->get_formatted_tp(), 'got trim point file content' );
my $expected_tp_file_string = q(001	TGP	Primary	Un	chromosome	accession1	100	L	test line 1 comment some nonsense here and here too
001	TGP	Primary	1	chromosome	accession2	200	H	test line 2 comment some nonsense here and here too
);
ok( $tp_file_string eq $expected_tp_file_string, 'trim point file content looks good ' );

#modify values for $tp_line1 and add to SP object
ok($tp_line1->set_accession_prefix_first_or_last_base(110), 'Re-set accession_prefix_first_or_last_base in $tp_line1');
ok($tp_line1->set_accession('accession3'), 'Re-set accession in $tp_line1');
ok ($tp->add_line_to_end($tp_line1), 'added modified $tp_line1');

#print out TP object and compare to string
ok($tp_file_string = $tp->get_formatted_tp(), 'got trim point file content again' );
$expected_tp_file_string = q(001	TGP	Primary	Un	chromosome	accession3	110	L	test line 1 comment some nonsense here and here too
001	TGP	Primary	1	chromosome	accession2	200	H	test line 2 comment some nonsense here and here too
001	TGP	Primary	Un	chromosome	accession3	110	L	test line 1 comment some nonsense here and here too
); #$tp_line1 also changes as its a hashref
ok( $tp_file_string eq $expected_tp_file_string, 'modified trim point file content looks good' );

#no more tests
done_testing();
