#!/usr/bin/env perl
# Indexes a VCF using callsam's sortVcf; bgzip; and tabix
# Author: Lee Katz

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use File::Basename;
use File::Copy qw/copy move/;

$0=fileparse $0;
sub logmsg{print STDERR "$0: @_\n";}

exit(main());
sub main{
  my $settings={};
  GetOptions($settings,qw(help force));
  die usage() if($$settings{help});
  die "ERROR: need input file\n". usage() if(!@ARGV);

  my $file=shift(@ARGV);

  # File derivatives
  my($sorted,$compressed,$indexFile)=(
    "$file.sorted.tmp",
    "$file.gz",
    "$file.gz.tbi",
  );
  # File checking
  if(!-e $file){
    die "ERROR: could not find $file, but I did find $compressed (already compressed?)\n" if(-f $compressed);
    die "ERROR: could not find $file\n".usage();
  }

  # Check if output files already exist
  for ($sorted,$compressed,$indexFile){
    last if($$settings{force});
    die "ERROR: $_ already exists and so this program will not continue. Use --force to force." if(-e $_);
  }

  # sort the VCF
  system("sortVcf.sh '$file' > '$sorted'");
  die "ERROR: could not sort $file into $sorted with sortVcf.sh" if $?;
  move($sorted,$file) or die $!;

  # compress the VCF
  system("bgzip -f '$file'");
  die "ERROR with bgzip on $file" if $?;

  # index the VCF
  system("tabix '$compressed'");
  die "ERROR with tabix" if $?;

  return 0;
}

sub usage{
  "Sort and index a VCF
  Usage: $0 in.vcf # produces in.vcf.gz*
  --force Overwrite any existing file
  "
}
