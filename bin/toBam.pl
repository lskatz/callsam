#!/usr/bin/env perl
# creates a sorted and indexed bam file from
# any Sequence Alignment/Map file.
# Author: Lee Katz <lkatz@cdc.gov>

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/../lib";
use CallSam qw/logmsg mktempdir is_sam is_bam is_cram/;

exit main();
sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s keep reference=s));
  die usage() if($$settings{help} || @ARGV<1);
  $$settings{tempdir}||=CallSam::mktempdir();
  
  my $file=\@ARGV;

  logmsg "Transforming files to bam";
  my $bam=toBam($file,$settings);
  logmsg "Sorting bam files";
  my $sorted=sortBam($bam,$settings);
  logmsg "Indexing bam files";
  indexBam($sorted,$settings);
  logmsg "Merging bam files";
  my $merged=mergeBam($sorted,$settings);
  system("cat '$merged'");
  die if $?;
  return 0;
}

sub toBam{
  my($file,$settings)=@_;
  my @bam;
  for my $f(@$file){
    my $bam;
    if(is_sam($f,$settings)){
      $bam=samToBam($f,$settings);
    } elsif(is_bam($f,$settings)){
      $bam=$f;
    } elsif(is_cram($f,$settings)){
      $bam=cramToBam($f,$settings);
    } else {
      die "ERROR: I do not understand the format of $_";
    }
    push(@bam,$bam);
  }
  return \@bam;
}

# converts one sam file to bam
sub samToBam{
  my($file,$settings)=@_;
  system("samtools view -bSh '$file' > '$file.bam'");
  die if $?;
  return "$file.bam";
}

# sorts an array of bams
sub sortBam{
  my($file,$settings)=@_;
  my @sorted;
  for(@$file){
    system("samtools sort '$_' '$_.sorted'");
    die if $?;
    push(@sorted,"$_.sorted.bam");
  }
  return \@sorted;
}

# indexes an array of sorted bams
sub indexBam{
  my($file,$settings)=@_;
  my @indexed;
  for(@$file){
    system("samtools index '$_'");
    die if $?;
    push(@indexed,$_);
  }
  return \@indexed;
}

# merges an array of sorted bams
sub mergeBam{
  my($file,$settings)=@_;
  my $merged="$$settings{tempdir}/merged.bam";
  return $$file[0] if(@$file==1);
  my $bamString="'".join("' '",@$file)."'";
  my $command="samtools merge '$merged' $bamString";
  #logmsg $command;
  system($command);
  die if $?;
  return $merged;
}

sub usage{
  "Transforms sam, bam, ... to a sorted and indexed bam file
  If multiple files are given, then they are merged with samtools merge
  Usage: $0 file.sam > file.bam
  -r reference.assembly.fasta
  -t tempdir/
  --keep to keep all intermediate files
  "
}
