#!/usr/bin/env perl
# CallSam: Call bases from a sequence alignment/map
# Author: Lee Katz <lkatz@cdc.gov>
# Run with no options for usage help

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw/sum min max/;
use Bio::Perl;
use File::Basename;
use threads;
use Thread::Queue;
use FindBin;
use lib "$FindBin::FindBin/../lib";
use CallSam qw/logmsg/;

$0=fileparse($0);
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help basename));
  die usage() if($$settings{help} || !@ARGV);

  my @VCF=@ARGV;
  die usage() if(!@VCF);

  vcfToFasta(\@VCF,$settings);
}

sub vcfToFasta{
  my($VCF,$settings)=@_;

  for my $vcf(@$VCF){
    my $fasta=vcfToFastaEntry($vcf,$settings);
    print $fasta;
    logmsg $vcf;
  }
  return 1;
}

sub vcfToFastaEntry{
  my($vcf,$settings)=@_;
  my @vcfHeader=qw(chrom pos id ref alt qual filter info);
  my @alt;
  open(VCF,$vcf) or die "ERROR: could not open $vcf: $!";
  while(<VCF>){
    next if(/^#/);
    chomp;
    my @F=split /\t/;
    push(@alt,$F[4]);
  }
  close VCF;
  my $defline=">$vcf\n";
    $defline=">".basename($vcf)."\n" if($$settings{basename});
  my $entry=$defline.join("",@alt)."\n";
  return $entry;
}

# TODO erase this sub from here and put it into the lib just in case
sub readVcf{
  my($vcf,$settings)=@_;
  my %vcfContent;
  my @vcfHeader=qw(chrom pos id ref alt qual filter info);
  open(VCF,$vcf) or die "ERROR: could not open $vcf: $!";
  while(<VCF>){
    next if(/^#/);
    chomp;
    my @F=split /\t/;
    my %line;
    @line{@vcfHeader}=@F;
    $vcfContent{$F[0]}{$F[1]}=\%line;
  }
  close VCF;
  return \%vcfContent;
}

sub usage{
  "Converts several VCFs to a multiple sequence alignment
  Usage: $0 *.vcf > out.aln.fasta
    -b to include only the basename of the filename in the fasta entries
  VCF files must be all-sites files.
  "
}
