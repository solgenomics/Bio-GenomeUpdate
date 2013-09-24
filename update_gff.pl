#!/usr/bin/perl

=head1 NAME

update_gff.pl

=head1 SYNOPSIS

update_gff.pl -o [old AGP file] -n [new AGP file] -g [GFF file] -m [output GFF file]

=head1 COMMAND-LINE OPTIONS

 -o  old AGP file (required)
 -n  new scaffold AGP file with updated coordinates (required)
 -g  GFF3 file based on old AGP to update to new AGP file (required)
 -m  output mapped GFF file (optional)
 -d  debugging messages (1 or 0)
 -h  Help

=cut

use strict;
use warnings;
use Bio::GenomeUpdate::AGP;
use Bio::GenomeUpdate::GFF;
use File::Slurp;
use Getopt::Std;
use Proc::ProcessTable;

our ( $opt_o, $opt_n, $opt_g, $opt_m, $opt_d, $opt_h );
getopts('o:n:g:m:d:h');
if ($opt_h) {
	help();
	exit;
}
if ( !$opt_o || !$opt_n || !$opt_g ) {
	print
"\nOld AGP, New AGP and GFF3 based on old AGP are required. 
See help below\n\n\n";
	help();
}

#get input files
my $old_agp_input_file = $opt_o;
my $input_old_agp      = read_file($old_agp_input_file)
  or die "Could not open old AGP input file: $old_agp_input_file\n";
my $new_agp_input_file = $opt_n;
my $input_new_agp      = read_file($new_agp_input_file)
  or die "Could not open new AGP input file: $new_agp_input_file\n";
my $gff_input_file = $opt_g;
my $input_gff      = read_file($gff_input_file)
  or die "Could not open GFF input file: $gff_input_file\n";
  
my $gff_output_file;
#set output file
if ($opt_m) {
	$gff_output_file = $opt_m;
}
if ($opt_d){ 
	print STDERR "Params parsed..\n";
	mem_used();
}

#object for old AGP
my $old_agp = Bio::GenomeUpdate::AGP->new();
$old_agp->parse_agp($input_old_agp);

#object for new AGP
my $new_agp = Bio::GenomeUpdate::AGP->new();
$new_agp->parse_agp($input_new_agp);

#object for gff
my $gff = Bio::GenomeUpdate::GFF->new();
$gff->parse_gff($input_gff); 
if ($opt_d){ 
	print STDERR "Files read..\n";
	mem_used();
}

#get coordinates mapped from old AGP to new AGP space
my %coords = $gff->get_reordered_coordinates($old_agp,$new_agp);
my %flips = $gff->get_flipped_coordinates($old_agp,$new_agp);
if ($opt_d){ 
	print STDERR "Hashes populated..\n";
	mem_used();
}

#remapping the GFF
$gff->remap_coordinates(\%coords,\%flips);
if ($opt_d){ 
	print STDERR "Coords remapped..\n";
	mem_used();
}

#print the GFF ($gff_output_file)
my $new_gff = $gff->get_formatted_gff();
unless(open(OGFF,">$gff_output_file")){print STDERR "Cannot open $gff_output_file\n"; exit 1;}
print OGFF $new_gff;
close(OGFF);
if ($opt_d){ 
	print STDERR "GFF written..\n";
	mem_used();
}

#----------------------------------------------------------------------------

sub help {
	print STDERR <<EOF;
  $0:

    Description:

     This script updates the coordinates of features in a Generic Feature Format version 3 (GFF)
     using a new Accessioned Golden Path (AGP) file.

    Usage:
      update_gff.pl -o [old AGP file] -n [new AGP file] -g [GFF file] -m [output GFF file]
      
    Flags:

 -o  old AGP file (required)
 -n  new scaffold  AGP file with updated coordinates (required)
 -g  GFF3 file based on old AGP to update to new AGP file (required)
 -m  output mapped GFF file (optional)
 -h  help

EOF
	exit(1);
}

sub mem_used{
	my ($i,$t); 
	$t = new Proc::ProcessTable;
	foreach my $got ( @{$t->table} ) {
		next if not $got->pid eq $$; $i=$got->size;
	}
	print STDERR "Process id=",$$,"\n"; print "Memory used(MB)=", $i/1024/1024, "\n";
}

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Surya Saha <suryasaha@cornell.edu , @SahaSurya>

=cut

__END__

# AGP old
SL2.40ch00	21803722	21803821	6262	U	100	contig	no
SL2.40ch00	21803822	21805821	6263	W	SL2.40sc04627	1	2000	0
SL2.40ch01	1	32987597	1	W	SL2.40sc04133	1	32987597	-
SL2.40ch01	32987598	32987697	2	U	100	contig	no
SL2.40ch01	32987698	35621724	3	W	SL2.40sc03666	1	2634027	0
SL2.40ch01	35621725	35621824	4	U	100	contig	no

# AGP new
SL2.40ch00	21801721	21803721	6261	W	SL2.40sc03931	1	2001	0
SL2.40ch00	21803722	21803821	6262	U	100	contig	no
SL2.40ch00	21803822	21805821	6263	W	SL2.40sc04627	1	2000	0
SL2.40ch01	1	32987597	1	F	SL2.40SC04133	1	32987597	+
SL2.40ch01	32987598	35267597	2	N	2280000	contig	no	
SL2.40ch01	35267598	36990194	3	F	SL2.40SC04191	1	1722597	+
SL2.40ch01	36990195	39120194	4	N	2130000	contig	no	
SL2.40ch01	39120195	41754221	5	F	SL2.40SC03666	1	2634027	+
SL2.40ch01	41754222	42324221	6	N	570000	contig	no	

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
