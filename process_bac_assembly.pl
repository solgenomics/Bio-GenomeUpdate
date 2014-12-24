#!/usr/bin/perl

=head1 NAME

process_bac_assembly.pl

=head1 SYNOPSIS

process_bac_assembly.pl -f [ACE file] -o [output directory] 

=head1 COMMAND-LINE OPTIONS

 -f  ACE file from Phrap assembly (required)
 -o  Output directory
 -h  Help

=cut

use strict;
use warnings;








#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script analyzes a ACE file from a Phrap assembly of BACs. Contigs with 1 read (BAC) are written to singleton_BACs.fas. Contigs with multiple reads (BACs) are written to contigs_BACs.fas and the corresponding meta-data to contigs_BACs.txt. 

    Usage:
      process_bac_assembly.pl -f [ACE file] -o [output directory]
      
    Flags:

         -f  ACE file from Phrap assembly (required)
         -o  Output directory 
         -h  Help

EOF
	exit(1);
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

__END__