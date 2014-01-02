#!/usr/bin/perl

=head1 NAME

validate_gff.pl

=head1 SYNOPSIS

validate_gff.pl -g [old GFF file] -f [old Fasta file] -u [updated GFF file] -n [new Fasta file] <other params>

=head1 COMMAND-LINE OPTIONS

 -g  old GFF file (required)
 -f  old Fasta file (required)
 -u  Updated GFF file output by update_coordinates_gff.pl (required)
 -n  New Fasta file based on new scaffold AGP (required)
 -d  debugging messages (1 or 0)
 -h  Help

=cut


#TODO: Write GFF to fasta, read them in and then compare the files using Bio::Seq 

use strict;
use warnings;

use Bio::GenomeUpdate::GFF;
use Bio::SeqIO;
use Utilities;
use File::Slurp;
use Getopt::Std;

our ( $opt_g, $opt_f, $opt_u, $opt_n, $opt_d, $opt_h );
getopts('g:f:u:n:d:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_g || !$opt_f || !$opt_u || !$opt_n ) {
	print
"\nOld GFF3/Fasta and updated GFF3/Fasta are required. 
See help below\n\n\n";
	help();
}

#prep input data
my $gff_old_file = $opt_g;
my $old_gff      = read_file($gff_old_file)
  or die "Could not open old GFF file: $gff_old_file\n";
my $fasta_old_file;
if ( -e $opt_f ){
	$fasta_old_file = $opt_f;
}
else{
	die "Could not open old Fasta file: $opt_f\n";
}
my $gff_updated_file = $opt_u;
my $updated_gff      = read_file($gff_updated_file)
  or die "Could not open updated GFF file: $gff_updated_file\n";
my $fasta_updated_file;
if ( -e $opt_n ){
	$fasta_updated_file = $opt_n;
}
else{
	die "Could not open old Fasta file: $opt_n\n";
}


my $util = Utilities->new();
if ($opt_d){ 
	print STDERR "Params parsed..\n";
	$util->mem_used();
	$util->run_time();
}

#Read in gff and fasta files
my $old_gff_obj = Bio::GenomeUpdate::GFF->new();
$old_gff_obj->parse_gff($old_gff, $gff_old_file);
$old_gff_obj->write_fasta('old_gff.fas', $fasta_old_file);

my $updated_gff_obj = Bio::GenomeUpdate::GFF->new();
$updated_gff_obj->parse_gff($updated_gff, $gff_updated_file);
$update_gff_obj->write_fasta('updated_gff.fas', $fasta_updated_file);
 
if ($opt_d){ 
	print STDERR "Files read..\n";
	$util->mem_used();
	$util->run_time();
}

#Read in both Fasta files and compare
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

     This script checks if the sequences produced using the updated GFF file match the sequences produced by the old GFF file. GFF files should be in Generic Feature Format version 3.
     
    NOTE:


    Usage:
      validate_gff.pl -g [old GFF file] -f [old Fasta file] -u [updated GFF file] -n [new Fasta file] <other params>
      
    Flags:

 -g  old GFF file (required)
 -f  old Fasta file (required)
 -u  Updated GFF file output by update_coordinates_gff.pl (required)
 -n  New Fasta file based on new scaffold AGP (required)
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
