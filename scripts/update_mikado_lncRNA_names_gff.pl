#!/usr/bin/perl

=head1 NAME

update_mikado_lncRNA_names_gff.pl

=head1 SYNOPSIS

update_mikado_lncRNA_names_gff.pl -i [Mikado GFF file]

=head1 COMMAND-LINE OPTIONS

 -i  Mikado GFF [required]
 -h  Help

=cut

use strict;
use warnings;
use File::Slurp;
use Getopt::Std;

our ( $opt_i, $opt_h );
getopts('i:h');
if ($opt_h) {
    help();
    exit;
}
if ( !$opt_i ) {
    print "\nMikado GFF3 is required.
See help below\n\n\n";
    help();
}

#get input files
my $old_gff_input_file = $opt_i;
my $input_old_gff      = read_file($old_gff_input_file)
  or die "Could not open Mikado GFF input file: $old_gff_input_file\n";

my $new_id_output_file;
$new_id_output_file = $old_gff_input_file . '_new-ids.names';

my @lines        = split( /\n/, $input_old_gff );
my $line_count   = scalar(@lines);
my $line_counter = 0;
my $gene_flag    = 0;
my @gene_gff_line_arr;
my $current_mRNA_Solycid;
my $prev_mRNA_Solycid;
my $new_id_output           = '';
my $outofrange_gene_counter = 0;
my %mRNA_Solycid_hash;

#new
#my $species      = 'Solyc'; tomato ITAG
my $species      = 'Dcitr'; # Dcitr OGSv3
my $prefix       = 'r';
my $gene_counter = -1;
my $lnc_RNA_counter;
my $lnc_RNA_exon_counter;

foreach my $line (@lines) {
    $line_counter++;
    chomp($line);

    if ( $line =~ m/^#/ ) {
				print STDOUT $line . "\n";
				next;
    }

    # print STDERR "\rParsing GFF3 line "
    #   . $line_counter . " of "
    #   . $line_count . "\n";
    if ( $line =~ /\tgene\t/ ) {
        $lnc_RNA_counter      = 0;
        if ( $gene_counter == -1 ){
          $gene_counter++; # increment for first gene
        }
        else{
          $gene_counter = $gene_counter + 10; # increment for all other genes

        }

        my @line_arr = split( "\t", $line );
        $line_arr[1] =
          # 'mikado_ITAG';    #using source to reflect ITAG/eugene fed into maker
        'mikado_OGSv3';    #using source to reflect mikado/OGSv3 fed into maker
        my $chr = $line_arr[0];
        # $chr =~ s/^SL4\.0ch//; # tomato
        $chr =~ s/^DC3\.0sc//; # Dcitr
        my $gene_count   = sprintf( "%06d", $gene_counter );
        my $gene_version = 1;
        my $gene_id      = $species . $chr . $prefix . $gene_count;

        my $new_attr =
            'ID=gene:'
          . $gene_id . '.'
          . $gene_version
          . ';Alias='
          . $gene_id
          . ';Name='
          . $gene_id . '.'
          . $gene_version . "\n";

        for ( 0 .. 7 ) {
            print STDOUT $line_arr[$_] . "\t";
        }
        print STDOUT $new_attr;
    }
    elsif ( $line =~ /\tlncRNA\t/ ) {
				$lnc_RNA_exon_counter = 1;
				$lnc_RNA_counter++; # incrementing for next isoform, if any

        my @line_arr = split( "\t", $line );
        $line_arr[1] =
          # 'mikado_ITAG';    #using source to reflect ITAG/eugene fed into maker
        'mikado_OGSv3';    #using source to reflect mikado/OGSv3 fed into maker
        $line_arr[2] =
          'lnc_RNA'
          ; #using SO compatible long non-coding term http://www.sequenceontology.org/browser/current_svn/term/SO:0001877
        my $chr = $line_arr[0];
        # $chr =~ s/^SL4\.0ch//; #tomato
        $chr =~ s/^DC3\.0sc//; # Dcitr
        my $gene_count   = sprintf( "%06d", $gene_counter );
        my $gene_version = 1;
        my $gene_id      = $species . $chr . $prefix . $gene_count;

        my $new_attr =
            'ID=lnc_RNA:'
          . $gene_id . '.'
          . $gene_version . '.'
          . $lnc_RNA_counter
          . ';Name='
          . $gene_id . '.'
          . $gene_version . '.'
          . $lnc_RNA_counter
          . ';Parent=gene:'
          . $gene_id . '.'
          . $gene_version . "\n";

        for ( 0 .. 7 ) {
            print STDOUT $line_arr[$_] . "\t";
        }
        print STDOUT $new_attr;
    }
    elsif ( $line =~ /\texon\t/ ) {
        my @line_arr = split( "\t", $line );
        # $line_arr[1] =
        #   'mikado_ITAG';    #using source to reflect ITAG/eugene fed into maker for tomato
        $line_arr[1] =
          'mikado_OGSv3';    #using source to reflect ITAG/eugene fed into maker for OGSv3
        my $chr = $line_arr[0];
        # $chr =~ s/^SL4\.0ch//; #tomato
        $chr =~ s/^DC3\.0sc//; # Dcitr
        my $gene_count   = sprintf( "%06d", $gene_counter );
        my $gene_version = 1;
        my $gene_id      = $species . $chr . $prefix . $gene_count;

        my $new_attr =
            'ID=exon:'
          . $gene_id . '.'
          . $gene_version . '.'
					. $lnc_RNA_counter . '.'
          . $lnc_RNA_exon_counter
          . ';Name='
          . $gene_id . '.'
          . $gene_version . '.'
					. $lnc_RNA_counter . '.'
          . $lnc_RNA_exon_counter
          . ';Parent=lnc_RNA:'
          . $gene_id . '.'
          . $gene_version . '.'
          . $lnc_RNA_counter . "\n";

        $lnc_RNA_exon_counter++;    # incrementing for next isoform, if any

        for ( 0 .. 7 ) {
            print STDOUT $line_arr[$_] . "\t";
        }
        print STDOUT $new_attr;
    }
}

#----------------------------------------------------------------------------

sub help {
    print STDERR <<EOF;
  $0:

    Description:

    Add gene, lnc_RNA and exon ids in the ITAG convention to Mikado lnc_RNA predictions while discarding old ids. The mapping of old ids to new ids is written to a new-ids.names file. Alternative transcripts for a gene are allowed. Non-coding genes follow the Solyc00r000000 convention. Using a gap of 10 between successive ids.


    Usage:
      update_mikado_lncRNA_names_gff.pl

    Flags:

     -i  Mikado GFF file (required)
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


Input
SL4.0ch00	.	gene	381860	398979	.	+	.	ID=gene1796
SL4.0ch00	.	mRNA	381860	398979	.	+	.	ID=mRNA2051;Parent=gene1796
SL4.0ch00	Mikado_loci	exon	381860	381923	.	+	.	Parent=mRNA2051
SL4.0ch00	Mikado_loci	exon	392183	394836	.	+	.	Parent=mRNA2051
SL4.0ch00	Mikado_loci	exon	395898	398979	.	+	.	Parent=mRNA2051
###
SL4.0ch00	.	gene	546387	548638	.	-	.	ID=gene1759
SL4.0ch00	.	mRNA	546387	548638	.	-	.	ID=mRNA2009;Parent=gene1759
SL4.0ch00	Mikado_loci	exon	546387	546860	.	-	.	Parent=mRNA2009
SL4.0ch00	Mikado_loci	exon	547762	548199	.	-	.	Parent=mRNA2009
SL4.0ch00	Mikado_loci	exon	548575	548638	.	-	.	Parent=mRNA2009

Output
SL4.0ch00	mikado_ITAG	gene	381860	398979	.	+	.	ID=gene:Solyc00r000000.1;Alias=Solyc00r000000;Name=Solyc00r000000.1
SL4.0ch00	mikado_ITAG	lnc_RNA	381860	398979	.	+	.	ID=lnc_RNA:Solyc00r000000.1.1;Parent=gene:Solyc00r000000.1;Name=Solyc00r000000.1.1
SL4.0ch00	mikado_ITAG	exon	381860	381923	.	+	.	ID=exon:Solyc00r000000.1.1.1;Parent=lnc_RNA:Solyc00r000000.1.1
SL4.0ch00	mikado_ITAG	exon	392183	394836	.	+	.	ID=exon:Solyc00r000000.1.1.2;Parent=lnc_RNA:Solyc00r000000.1.1
SL4.0ch00	mikado_ITAG	exon	395898	398979	.	+	.	ID=exon:Solyc00r000000.1.1.3;Parent=lnc_RNA:Solyc00r000000.1.1
###
SL4.0ch00	mikado_ITAG	gene	546387	548638	.	-	.	ID=gene:Solyc00r000010.1;Alias=Solyc00r000010;Name=Solyc00r000010.1
SL4.0ch00	mikado_ITAG	lnc_RNA	546387	548638	.	-	.	ID=lnc_RNA:Solyc00r000010.1.1;Parent=gene:Solyc00r000010.1;Name=Solyc00r000010.1.1
SL4.0ch00	mikado_ITAG	exon	546387	546860	.	-	.	ID=exon:Solyc00r000010.1.1.1;Parent=lnc_RNA:Solyc00r000010.1.1
SL4.0ch00	mikado_ITAG	exon	547762	548199	.	-	.	ID=exon:Solyc00r000010.1.1.2;Parent=lnc_RNA:Solyc00r000010.1.2
SL4.0ch00	mikado_ITAG	exon	548575	548638	.	-	.	ID=exon:Solyc00r000010.1.1.3;Parent=lnc_RNA:Solyc00r000010.1.3
