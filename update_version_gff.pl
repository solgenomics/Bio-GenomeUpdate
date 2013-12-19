#!/usr/bin/perl

=head1 NAME

update_version_gff.pl

=head1 SYNOPSIS

update_coordinates_gff.pl -i [old GFF file] -o [output GFF file] <other params>

=head1 COMMAND-LINE OPTIONS

 -i  old GFF file (required)
 -o  output GFF file with new versions (required)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;
use Bio::GenomeUpdate::AGP;
use Bio::GenomeUpdate::GFF;
use Utilities;
use File::Slurp;
use Getopt::Std;

our ( $opt_i, $opt_o, $opt_d, $opt_h );
getopts('o:n:g:m:c:d:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_i || !$opt_o ) {
	print
"\nOld GFF3 and New GFF3 are required. 
See help below\n\n\n";
	help();
}

#get input files
my $gff_input_file = $opt_i;
my $input_gff      = read_file($gff_input_file)
  or die "Could not open GFF input file: $gff_input_file\n";
my $gff_output_file = $opt_o;

my $util = Utilities->new();
if ($opt_d){ 
	print STDERR "Params parsed..\n";
	$util->mem_used();
	$util->run_time();
}

#object for gff
my $gff = Bio::GenomeUpdate::GFF->new();
$gff->parse_gff($input_gff, $gff_input_file); 
if ($opt_d){ 
	print STDERR "Files read..\n";
	$util->mem_used();
	$util->run_time();
}

#update version(s)
#TODO



if ($opt_d){ 
	print STDERR "Version updated..\n";
	$util->mem_used();
	$util->run_time();
}

#print the GFF ($gff_output_file)
my $new_gff = $gff->get_formatted_gff();
unless(open(OGFF,">$gff_output_file")){print STDERR "Cannot open $gff_output_file\n"; exit 1;}
print OGFF $new_gff;
close(OGFF);
if ($opt_d){ 
	print STDERR "GFF written..\n";
	$util->mem_used();
	$util->run_time();
}

#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script updates the version of features in a Generic Feature Format version 3 (GFF)
     
    NOTE:


    Usage:
      update_coordinates_gff.pl -i [old GFF file] -o [output GFF file] <other params>
      
    Flags:

 -i  old GFF file (required)
 -o  output GFF file with new versions (required)
 -d  debugging messages (1 or 0)
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

# GFF
SL2.40ch00	SL2.40_assembly	supercontig	21801721	21803721	.	+	.	ID=SL2.40sc03931;Parent=SL2.40ch00;Name=SL2.40sc03931;Target=SL2.40sc03931 1 2001 +;reliably_oriented=0
SL2.40ch00	SL2.40_assembly	contig	21801721	21803721	.	+	.	ID=SL2.40ct09637;Name=SL2.40ct09637;Parent=SL2.40sc03931;Target=SL2.40ct09637 1 2001 +;scaffold_reliably_oriented=0;reliably_oriented=1
SL2.40ch00	SL2.40_assembly	remark	21803722	21803821	.	+	.	Name=contig_gap;Note=type: unknown_gap%2C description: a gap between clone contigs (also called a "layout gap")
SL2.40ch00	SL2.40_assembly	supercontig	21803822	21805821	.	+	.	ID=SL2.40sc04627;Parent=SL2.40ch00;Name=SL2.40sc04627;Target=SL2.40sc04627 1 2000 +;reliably_oriented=0
SL2.40ch00	SL2.40_assembly	contig	21803822	21805821	.	+	.	ID=SL2.40ct16352;Name=SL2.40ct16352;Parent=SL2.40sc04627;Target=SL2.40ct16352 1 2000 +;scaffold_reliably_oriented=0;reliably_oriented=1
SL2.40ch01	SL2.40_assembly	ultracontig	1	90304244	.	.	.	ID=SL2.40ch01;Name=SL2.40ch01
SL2.40ch01	SL2.40_assembly	supercontig	1	32987597	.	-	.	ID=SL2.40sc04133;Parent=SL2.40ch01;Name=SL2.40sc04133;Target=SL2.40sc04133 1 32987597 +;reliably_oriented=1
SL2.40ch01	SL2.40_assembly	contig	1	3610	.	-	.	ID=SL2.40ct13037;Name=SL2.40ct13037;Parent=SL2.40sc04133;Target=SL2.40ct13037 1 3610 +;scaffold_reliably_oriented=1;reliably_oriented=1
SL2.40ch01	SL2.40_assembly	contig	3864	12848	.	-	.	ID=SL2.40ct13036;Name=SL2.40ct13036;Parent=SL2.40sc04133;Target=SL2.40ct13036 1 8985 +;scaffold_reliably_oriented=1;reliably_oriented=1
SL2.40ch01	SL2.40_assembly	contig	12869	15910	.	-	.	ID=SL2.40ct13035;Name=SL2.40ct13035;Parent=SL2.40sc04133;Target=SL2.40ct13035 1 3042 +;scaffold_reliably_oriented=1;reliably_oriented=1
SL2.40ch01	SL2.40_assembly	contig	16487	76523	.	-	.	ID=SL2.40ct13034;Name=SL2.40ct13034;Parent=SL2.40sc04133;Target=SL2.40ct13034 1 60037 +;scaffold_reliably_oriented=1;reliably_oriented=1
SL2.40ch01	SL2.40_assembly	contig	77166	85600	.	-	.	ID=SL2.40ct13033;Name=SL2.40ct13033;Parent=SL2.40sc04133;Target=SL2.40ct13033 1 8435 +;scaffold_reliably_oriented=1;reliably_oriented=1
