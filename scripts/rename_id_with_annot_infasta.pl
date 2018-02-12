
#!/usr/local/bin/perl
#Author: Mirella Flores

=head1 NAME
rename_id_with_annot_infasta.pl

=head1 DESCRIPTION
Cutomized for the following inputs files.
    File with description:
    Solyc00g005000.2    Solyc00g005000.3    Note=Eukaryotic aspartyl protease family protein (AHRD V3.3 *** 
    Fasta File:
    mRNA:Solyc00g005000.3.1 (joined) (translated)
=cut

#use strict;
use Data::Dumper qw(Dumper);
use List::Util 'max';

#use warnings;
 
open OUT, ">ITAG3.0_proteins_full_desc.fasta" or die;
open sOUT, ">ITAG3.0_proteins_.fasta" or die;


### To get Function
my %hashdesc;

open FILE3, "ITAG3.0_annotationlist.v2.txt" or die;

while (my $line3=<FILE3>) {   
    # chomp;
    (my $prot,my $protid, my $desc) = split /\t/, $line3;  
    $desc =~ s/Note=//g;
    (my $sdesc) = split /\;/, $desc;  
      if (length(trim($protid)) > 0){
	   	$hashdesc{trim($protid.".1")}{desc}   = $desc;
        $hashdesc{trim($protid.".1")}{sdesc}   = trim($sdesc);
	}
}


open ($sfh, "<ITAG3.0_proteins.fasta");

while (<$sfh>) {
    ## if we find a header line ...
    if (/^\>(.*)/) {
        chomp;
        ## write the previous sequence before continuing with this one
        unless ($first) {
			(my $id) = split /\ /, $header; 
            $id =~ s/mRNA://g;
      	    print OUT ">" . $id . " ". $hashdesc{$id}{desc} ;
			print OUT  $seq;
            print sOUT ">" . $id . " ". $hashdesc{$id}{sdesc} ;
            print sOUT  "\n" . $seq;           
            ## reset the sequence
            $seq = '';
        }

        $first = 0;
        $header = $1;

    ## else we've found a sequence line
    } else {
        ## skip it if it is just whitespace
        next if (/^\s*$/);

        ## record this portion of the sequence
        $seq .= $_;
    }
}
          (my $id) = split /\ /, $header; 
          $id =~ s/mRNA://g;
      	    print OUT ">" . $id . " ". $hashdesc{$id}{desc} ;
			print OUT  $seq;
            print sOUT ">" . $id . " ". $hashdesc{$id}{sdesc} ;
            print sOUT  "\n" . $seq;           

sub trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_; } ;




