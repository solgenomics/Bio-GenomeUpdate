#!/usr/bin/perl

=head1 NAME
copy_updated_coordinates_to_vcf.pl
=head1 SYNOPSIS
copy_updated_coordinates_to_vcf.pl [updated pseudo GFF file] [old VCF file] [FASTA file] [chrom_label]
=cut

use warnings;
use strict;

my $GFF = $ARGV[0];
my $VCF = $ARGV[1];
my $FASTA = $ARGV[2];
my $chrom_label = $ARGV[3];

my $build_name = $chrom_label;
$build_name =~ s/(.*)ch$/$1/;

$VCF =~ s/.pseudo.gff(.*)/.vcf$1/;

my $new_vcf = $VCF;
$new_vcf =~ s/(.*)/updated_vcfs\/$1/;

open(MAPPEDGFF, "<", $GFF) || die "Can't open gff file $GFF";

open(OLDVCF, "<", $VCF) || die "Can't open vcf file $VCF";

open(my $NEWVCF, ">", $new_vcf) || die "Can't create new VCF file $new_vcf";

my $data_line_counter = 0;

for (<OLDVCF>) {
    chomp;
    if (m/^##reference=/) {
	    print $NEWVCF "##reference=$FASTA\n";
    }
    elsif (m/^#/) {
	    print $NEWVCF $_ . "\n";
    }
    else {
      my ($chrom, $position, $id, $ref, $alt, $qual, $filter, $info, @extra) = split /\t/;
      chomp (my $gff_string = <MAPPEDGFF>);
      my @values = split "\t", $gff_string;
      my $new_position = $values[3];
      my $new_orientation = $values[6];
      $chrom =~ s/^(.*)([0-9][0-9])$/$chrom_label$2/;

	    if ($new_orientation eq '-') { #if scaffold has been flipped
        my ($new_ref, $new_alt, @new_alts);

        if ($info =~ m/^INDEL/) { #if indel, fix position, get new leading ref base, and calculate new ref and alt seqs (reverse complements of old seqs)

          my $indel_length = length $ref;
          $new_position = $new_position - $indel_length;
          my $new_ref_base = &get_ref_base_by_position($FASTA, $chrom, $new_position);

          $new_ref = &find_reverse_compliment($ref, $new_ref_base);

          my @old_seqs = split ",", $alt;
          my @new_seqs;
          foreach my $seq (@old_seqs) {
            push @new_seqs, &find_reverse_compliment($seq, $new_ref_base);
          }
          $new_alt = join ",", @new_seqs;

          print STDERR "ref $ref with alt $alt is an INDEL\n";
          print STDERR "new leading base at position $new_position not affected by indel event is $new_ref_base\n";
          print STDERR "corrected to ref $new_ref and alt $new_alt\n";

	      } else { #if a SNP, find simple complement
          $new_ref = &replace_with_complementary_base($ref);

	  my @old_alts = split ",", $alt;
	  foreach my $base(@old_alts) {
	      push @new_alts, &replace_with_complementary_base($base);
	  }
	  $new_alt = join ",", @new_alts;
        }

        print $NEWVCF join "\t", ($chrom, $new_position, $id, $new_ref, $new_alt, $qual, $filter, $info, @extra);
        print $NEWVCF "\n";

      } else { # if not flipped

	       print $NEWVCF join "\t", ($chrom, $new_position, $id, $ref, $alt, $qual, $filter, $info, @extra);
         print $NEWVCF "\n";
      }
    }
  }

# sort new vcf file to ensure newly mapped coordinates are in order

system ("grep '^#' $new_vcf > $new_vcf.sorted");
system ("grep -v '^#' $new_vcf | sort -V >> $new_vcf.sorted");
system ("mv $new_vcf.sorted $new_vcf");

print "New vcf file $new_vcf with updated coordinates successfully created.\n";

sub replace_with_complementary_base {
    if ($_[0] eq 'A') {
	return 'T';
    }
    elsif ($_[0] eq 'T') {
	return 'A';
    }
    elsif ($_[0] eq 'C') {
	return 'G';
    }
    elsif ($_[0] eq 'G') {
	return 'C';
    }
    else {
    	return $_[0];  # to handle Ns
    }
}

sub find_reverse_compliment {
  my $old_seq = $_[0];
  my $new_ref_base = $_[1];
  my (@new_bases, $new_seq);

  foreach my $base (split '', $old_seq) {
    push @new_bases, &replace_with_complementary_base($base);
  }
  @new_bases = reverse @new_bases;
  my $old_leading_ref_base = pop @new_bases;
  unshift @new_bases, $new_ref_base;
  $new_seq = join "", @new_bases;
  return $new_seq;
}

sub get_ref_base_by_position {
  my $fasta = $_[0];
  my $chrom = $_[1];
  my $pos = $_[2];
  my @results = `samtools faidx $fasta $chrom:$pos-$pos`;
  my $ref_base = $results[1];
  chomp $ref_base;
  return $ref_base;
}
