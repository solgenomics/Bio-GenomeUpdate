#!/usr/local/bin/perl
#Author: Mirella Flores

=head1 NAME
compare_seq.pl
=head1 SYNOPSIS
compare_seq.pl
=cut

#use strict;
use Data::Dumper qw(Dumper);
use List::Util 'max';

#use warnings;
open OUT, ">list.txt" or die;

open ($sfh, "<ITAG3.0_proteins_cDNA.fasta");

my %hNew;
my @list;

while (<$sfh>) {
        chomp;
        if($_ =~ /^>(.+)/){
            $id = trim($1);
        }else{
            $hNew{$id} .= trim($_);
        }
    }

open ($sfh2, "<ITAG2.4_transcript.fasta");

my %hOld;

while (<$sfh2>) {
    chomp;
        if($_ =~ /^>(.+)/){
            $id = trim($1);
        }else{
            $hOld{$id} .= trim($_);
        }
    }

foreach my $gene (keys %hNew){

    if ($hNew{$gene} eq $hOld{$gene}) {
        print OUT $gene ."\n";
    }
}

sub trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };

