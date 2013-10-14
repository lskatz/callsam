#!/usr/bin/env perl

package CallSam;
require 5.005;
my ($VERSION) = ('$Id$' =~ /,v\s+(\d+\S+)/o);

use strict;
use List::Util qw(min max sum reduce shuffle);
use File::Basename;
use File::Spec;
use File::Copy;
use File::Temp ('tempdir');
use Data::Dumper;

use Exporter;
our @ISA = "Exporter";
our @methods = qw(logmsg fullPathToExec mktempdir is_sam is_bam);
our %EXPORT_TAGS = (all => [@methods]);
Exporter::export_ok_tags('all');

local $SIG{'__DIE__'} = sub { my $e = $_[0]; $e =~ s/(at [^\s]+? line \d+\.$)/\nStopped $1/; die("$0: ".(caller(1))[3].": ".$e); };

sub logmsg {my $FH = *STDERR; print $FH "$0: ".(caller(1))[3].": @_\n";}

sub mktempdir(;$) {
	my ($settings) = @_;
	my $tempdir_path = File::Spec->join(File::Spec->tmpdir(), (split("::",(caller(1))[3]))[1].".$$.XXXXX");
	my $tempdir = tempdir($tempdir_path, CLEANUP => !($$settings{keep}));
	return $tempdir;
}

# If argument is an executable in the current path, returns the full path to it, otherwise returns undef
# arguments: warn_on_error
sub fullPathToExec($;$) {
	my ($executable,$settings) = @_;
	my $fullpath;
	for ("", split(/:/, $ENV{PATH})) {
    my $path=$_."/".$executable;
		if (-x $path && -f $path) { $fullpath = File::Spec->rel2abs($path); last; }
	}
  if(! -x $fullpath){
	  my $errStr="Error finding full path to executable ($executable)";
    warn $errStr if($$settings{warn_on_error});
    die $errStr if(!$$settings{warn_on_error});
  }
	return $fullpath;
}

sub getNumCPUs() {
	my $num_cpus;
	open(IN, '<', '/proc/cpuinfo'); while (<IN>) { /processor\s*\:\s*\d+/ or next; $num_cpus++; } close IN;
	return $num_cpus || 1;
}

sub is_sam{
  my($file,$settings)=@_;
  # A sam file ends with .sam and has at least 11 columns, after any headers
  my($name,$path,$suffix)=fileparse($file,qw(.sam));
  return 0 if($suffix !~/sam$/i);

  my $is_sam=1;
  open(TESTFILE,$file) or die "Could not open file $file for reading:$!";
  while(<TESTFILE>){
    s/^\s+|\s+$//g; # trim
    if(/^$/){       # empty line
      next;
    } elsif (/^@/){ # a header
      next;
    } else {
      my @F=split /\t/;
      $is_sam=0 if(@F<11);
      last;
    }
  }
  close TESTFILE;
  return $is_sam;
}

sub is_bam{
  my($file,$settings)=@_;
  $$settings{tempdir} || die "ERROR: tempdir was not set in settings";
  my $tmpSam="$$settings{tempdir}/tmp.sam";

  # A bam file ends with .bam.
  # It works with samtools view and has at least 11 columns, after any headers
  my($name,$path,$suffix)=fileparse $file;
  return 0 if($suffix !~/\.bam$/i);
  
  system("samtools view '$file' | head -9999 > '$tmpSam'");
  die if $?;
  
  return is_sam($tmpSam,$settings);
}

1;

