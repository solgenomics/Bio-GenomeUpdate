
#!/usr/local/bin/perl
#Author: Mirella Flores

=head1 NAME
get_functional_annot.pl
=head1 SYNOPSIS
get_functional_annot.pl
=cut

#use strict;
use Data::Dumper qw(Dumper);
use List::Util 'max';

#use warnings;
 
open OUT, ">merged_list.txt" or die;


### To get Function
my %hashahrd;
my %hashncbi;
my %hashlocus;
my %hashinterPfam;
my %hashinterGO;
my %hashversion;


open FILE1, "solyc_curated_loci.txt" or die;

while (my $line1=<FILE1>) {   
    # chomp;
    (my $protid, my $desc) = split /\t/, $line1;  
	   if (length(trim($protid)) > 0){
	   	$hashlocus{trim($protid)}{desc}   = trim($desc);  #Solyc06g064520          (1-4)-beta-mannan endohydrolase precursor
	}
}

open FILE2, "ncbi_annot.txt" or die;

while (my $line2=<FILE2>) {   
    # chomp;
    (my $protidc,my $ort,my $protid, my $rate, my $desc) = split /\t/, $line2;  
	   if (length(trim($protid)) > 0){
	   	$hashncbi{trim($protid)}{desc}   = trim($desc);   #Solyc03g111050        CONSTANS interacting protein 5
	}
}


open FILE3, "ahrd_ITAG_filtered.txt" or die;

while (my $line3=<FILE3>) {   
    # chomp;
    (my $protid, my $desc, my $desc1) = split /\t/, $line3;  
    $protid =~ s/\.[0-9]//g;
    (my $db,my $ortid,my $geneid) = split /\|/, $desc1;  
    if (length(trim($geneid))<1) { $geneid= $db;}

	   if (length(trim($protid)) > 0){
	   	$hashahrd{trim($protid)}{desc}   = trim($desc)." ". trim($geneid).")";   #Solyc05g015840.2.1      Squamosa promoter binding protein NtabSPL13 (AHRD V3.3 *** A0A125SZN9_TOBAC)
	}
}

open FILE4, "ITAG3.0_proteins.go_interpro.txt" or die;

while (my $line4=<FILE4>) {   
    # chomp;
    (my $protid, my $interpro, my $gos) = split /\t/, $line4;  
    $gos =~ s/\|/,/g; 
    $protid =~ s/mRNA://g;   
    $protid =~ s/\.[0-9]//g;

	   if (length(trim($protid)) > 0){
	   	$hashinterGO{trim($protid)}{interpro}   = trim($interpro);   # mRNA:Solyc00g017910.1.1 IPR002504       GO:0003951|GO:0006741|GO:0008152
	   	$hashinterGO{trim($protid)}{go}   = trim($gos);  
	   }
}
#print Dumper \%hashinterGO;

open FILE5, "ITAG3.0_proteins.pfam_interpro.txt" or die;

while (my $line5=<FILE5>) {   
    # chomp;
    (my $protid,my $interpro,my $q, my $pfam) = split /\t/, $line5;  
    $protid =~ s/mRNA://g; 
        $protid =~ s/\.[0-9]//g;

	   if (length(trim($protid)) > 0){
	   	$hashinterPfam{trim($protid)}{interpro}   = trim($interpro);  
	   	$hashinterPfam{trim($protid)}{pfam}   = trim($pfam);  # mRNA:Solyc00g041190.1.1 IPR009081       Pfam    PF00550
	   }
}
 

open FILE7, "list4versioning.txt" or die;

while (my $line7=<FILE7>) {   
    # chomp;
    (my $protid, my $protidold, my $protidnew) = split /\t/, $line7;  
	   if (length(trim($protid)) > 0){
			$hashversion{trim($protid)}{old}   = trim($protidold); 
		   	$hashversion{trim($protid)}{new}   = trim($protidnew);  #Solyc06g064520          (1-4)-beta-mannan endohydrolase precursor
	}
}

open FILE6, "listSolyc.txt" or die;
my $output; my $i,$j,$k,$l;
while (my $line6=<FILE6>) {   
    # chomp;
    (my $id1,my $aed,my $note) = split /\t/, $line6;  
    $id1 = trim($id1);
    $id =  $id1 ;
     $id1 =~ s/\.1\.1/\.1/g;
      $id1 =~ s/\.2\.1/\.2/g;
        $note = trim($note);
    $id =~ s/\.[0-9]//g;

    if (trim($hashversion{$id}{old}) ne '') {
    	    $hashversion{$id}{old} =~ s/\.1\.1/\.1/g;
      		$hashversion{$id}{old} =~ s/\.2\.1/\.2/g;
		$output =  $hashversion{$id}{old}."\t".$hashversion{$id}{new}."\t". $note ;
	}
	else {
		$output = $id1 ."\t".$id1 ."\t" . $note ;
	}
    if (trim($hashlocus{$id}{desc}) ne '') {
		$output = $output .  $hashlocus{$id}{desc} ; $j++;
    }
    elsif (trim($hashncbi{$id}{desc}) ne '') {
		$output = $output .  $hashncbi{$id}{desc} ;  $i++; #print $id."\t";
    }
    elsif (trim($hashahrd{$id}{desc}) ne '') {
		$output = $output .  $hashahrd{$id}{desc} ; $k++;
	}
	else{
$l++;
	}
		if (trim($hashinterGO{$id}{interpro}) ne ''){
			#$output = $output . " contains Interpro domain(s) " . $hashinterGO{$id}{interpro} .
			$output = $output . ";Dbxref=InterPro:". $hashinterGO{$id}{interpro} ;
		}
		elsif (trim($hashinterPfam{$id}{interpro}) ne ''){
			$output = $output . ";Dbxref=InterPro:". $hashinterPfam{$id}{interpro} ;
		}

		if (trim($hashinterPfam{$id}{pfam}) ne ''){
			$output = $output . ",Pfam:" . $hashinterPfam{$id}{pfam} ;
		}
		if (trim($hashinterGO{$id}{go}) ne ''){
			$output = $output . ";Ontology_term=" . trim($hashinterGO{$id}{go}) ;
		}

print OUT $output . "\n"; 
}
print $i."\t".$j."\t".$k."\t".$l."\n";

sub trim { my $s = shift; $s =~ s/^\s+|\s+$//g; return $s };







