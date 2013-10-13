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

$0=fileparse($0);
sub logmsg{$|++; print STDERR "$0: @_\n"; $|--;}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help min-coverage=i min-frequency=s reference=s numcpus=i unsorted variants-only mpileupxopts=s debug));
  die usage() if($$settings{help});
  $$settings{'min-coverage'}||=10;
  $$settings{'min-frequency'}||=0.75;
  $$settings{'reference'} || logmsg("Warning: reference not given");
  $$settings{numcpus}||=1;
  $$settings{mpileupxopts}||="-q 1";
  my ($file)=@ARGV;
  die "ERROR: need input file\n".usage() if(!$file);
  
  # If there is only 1 cpu, then the output can be unsorted streaming.
  if($$settings{numcpus} < 2){
    $$settings{numcpus}=1;
    if(!$$settings{unsorted}){
      logmsg "Num cpus are 1, so there are no race conditions. So, I am setting the option 'unsorted' so that the output is streaming.";
      $$settings{unsorted}=1;
    }
  }

  # get all the reference bases into a hash
  my $refBase=readReference($settings);
  printHeaders($settings);
  my $numpositions=bamToVcf($file,$refBase,$settings);

  logmsg "Done. $numpositions positions were analyzed.";

  return 0;
}

sub readReference{
  my($settings)=@_;
  my %seq;

  return \%seq if(!$$settings{reference});
  die "Could not locate the reference $$settings{reference}\n".usage() if(!-f $$settings{reference});

  my $in=Bio::SeqIO->new(-file=>$$settings{reference});
  while(my $seq=$in->next_seq){
    # undef the zero position to make this 1-based
    $seq{$seq->id}=[undef,split(//,$seq->seq)];
  }
  return \%seq;
}

sub printHeaders{
  my($settings)=@_;
  my($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
  my $date=sprintf("%04d%02d%02d",$year+1900,$mon+1,$mday);

  # simply put all the headers into an array to process later
  my @header=split/\s*\n+\s*/,qq(
      fileformat=VCFv4.2
      fileDate=$date
      source=$0
      INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
      INFO=<ID=AC,Number=1,Type=Integer,Description="allele count in genotypes, for each ALT allele, in the same order as listed">
  );
  push(@header,"reference=".$$settings{reference}) if($$settings{reference});
  
  # clean up and process the headers
  @header=grep(!/^\s*$/,@header);
  for(@header){
    $_='##'.$_;
  }
  # add the final header with one hash tag
  push(@header,join("\t",'#CHROM',qw(POS ID REF ALT QUAL FILTER INFO)));
  # finally print it out
  print $_."\n" for(@header);
  # return the number of headers
  return scalar(@header);
}
sub bamToVcf{
  my($file,$refBase,$settings)=@_;

  my $Q=Thread::Queue->new;
  my $printQueue=Thread::Queue->new;
  my @thr;
  $thr[$_]=threads->new(\&pileupWorker,$Q,$printQueue,$refBase,$settings) for(0..$$settings{numcpus}-1);
  my $printer=threads->new(\&vcfPrinter,$printQueue,$settings);

  # Run mpileup to show the pileup at each position.
  # Feed each line of mpileup to the threads that analyze them.
  my $fp;
  my $numPositions=0;
  my $command="samtools mpileup $$settings{mpileupxopts} -O -s '$file'";
  logmsg "\n  $command";
  open($fp,"$command | ") or die "Could not open $file with samtools mpileup:$!";
  # Because it is cpu-heavy to enqueue once each line,
  # Save a bunch of lines in an array before enqueuing them.
  my @buffer;
  while(<$fp>){
    $numPositions++;
    push(@buffer,$_);
    if($numPositions % 10000 == 0){
      $Q->enqueue(@buffer);
      @buffer=();
    }
    last if($$settings{debug} && $numPositions>9999);
  }
  $Q->enqueue(@buffer); # enqueue any remaining lines
  close $fp;

  # wrap up the threads
  $Q->enqueue(undef) for(@thr);
  for(@thr){
    $_->join;
  }
  $printQueue->enqueue(undef);
  $printer->join;

  return $numPositions;
}

sub pileupWorker{
  my($Q,$printQ,$refBase,$settings)=@_;
  my @bamField=qw(contig pos basecall depth dna qual mappingQual readPos);
  while(defined(my $line=$Q->dequeue)){
    chomp $line;
    # %F and @F have all the mpileup fields.
    # These values will be parsed to make VCF output fields.
    my @F=split /\t/, $line;
    my %F;
    @F{@bamField}=@F;
    $F{info}={DP=>$F{depth}}; # put the depth into the info field so that it is displayed correctly.
    # Turn the DNA cigar line to an array.
    $F{dnaArr}=parseDnaCigar(\%F,$refBase,$settings);
    # Use the DNA array and other %F fields to make a consensus base.
    my ($basecall,$passFail,$qual)=findConsensus(\%F,$settings);
    # Find the reference base in the complex hash. A dot if not found.
    my $ref=$$refBase{$F{contig}}[$F{pos}] || '.'; 
    # A samtools-style identifier for the appropriate VCF field.
    my $ID=$F{contig}.':'.$F{pos};

    # Use the info hash to generate the VCF info field
    my $info="";
    while(my($key,$value)=each(%{$F{info}})){
      $info.="$key=$value;";
    } 
    # chop off that semicolon
    $info=~s/;$//;
    # I wonder if substr($info,0,-1) would be faster to remove the semicolon
    #$info=substr($info,0,-1);

    # if the user only wants variants and it is not a variant site, then skip it
    next if($$settings{'variants-only'} && ($ref eq $basecall || $basecall eq 'N' || $ref eq 'N'));
    $printQ->enqueue( join("\t",$F{contig},$F{pos},$ID,$ref,$basecall,$qual,$passFail,$info)."\n" );
  }
}

# The threads feed into the vcfPrinter.
# This subroutine should be the only thing that prints to stdout
sub vcfPrinter{
  my($Q,$settings)=@_;

  if($$settings{unsorted}){
    $|++;
    while(defined(my $line=$Q->dequeue)){
      print $line;
    }
    return;
  }

  my @unsorted;
  while(defined(my $line=$Q->dequeue)){
    push(@unsorted,$line);
  }
  return if($$settings{unsorted});
  
  my @sorted=sort {
    my($contigA,$posA)=split /\t/,$a;
    my($contigB,$posB)=split /\t/,$b;
    return $contigA cmp $contigB if($contigA ne $contigB);
    $posA <=> $posB;
  } @unsorted;
  print $_ for(@sorted);
}  

# finds the consensus for a position
sub findConsensus{
  my($F,$settings)=@_;
  my $passFail=""; # for the filter field in the VCF
  # min depth requirement
  $passFail.="d$$F{depth};" if($$F{depth} < $$settings{'min-coverage'});

  my $dna=uc($$F{dna});

  # find counts
  my %nt;
  $nt{$_}=0 for(("A".."Z"),'*'); # start the counts at zero
  for (@{$$F{dnaArr}}){
    $nt{uc($_)}++;
  }

  # Sort the counts and find the majority
  my @majorityNt=sort {
    return $nt{$b}<=>$nt{$a}; 
  } keys(%nt);
  my $winner=$majorityNt[0];

  # alter the hash to show the Allele Count
  $$F{info}{AC}=$nt{$winner};

  # Majority consensus requirement
  my $frequency=sprintf("%0.2f",$nt{$winner}/$$F{depth});
  $passFail.="freq$frequency;" if($frequency < $$settings{'min-frequency'});

  # set the pass/fail field correctly
  if($passFail){
    $winner="N";
    $passFail=~s/;+$//;
  } else {
    $passFail="PASS";
  }

  # Make some kind of score for the SNP
  my $score=0;
  my @qual=map(ord($_)-33,split(//,$$F{qual}));
  my @baq =map(ord($_)-33,split(//,$$F{mappingQual}));
  for(my $i=0;$i<$$F{depth};$i++){
    my $weight=$qual[$i]*$baq[$i];
    if($$F{dnaArr}[$i] eq $winner){
      $score+=$weight;
    } else {
      $score-=$weight;
    }
  }
  # transform the score
  # TODO maybe put some rationale and thinking into this formula
  $score=sprintf("%0.2f",$score/sum(@qual,@baq));

  # If the score is too small, then the base call is ambiguous and doesn't pass the filter
  # TODO: test what the score theshold should be 
  if($score < 0){
    $passFail="score:$score;mostlikely:$winner";
    $winner='N';
  }

  return ($winner,$passFail,$score);
}

# Turn a cigar string into a meaningful array of bases
sub parseDnaCigar{
  my($bamField,$refBase,$settings)=@_;
  my $cigar=$$bamField{dna};
  $cigar=uc($cigar); 
  
  my $length=length($cigar);
  # Parse the cigar.
  # ^Xn is the start of a read with base quality with a nucleotide.
  # n$ is the end of the read with its nucleotide.
  # . is a match.
  # , is the reverse-strand match.
  my @base=();
  for(my $i=0;$i<$length;$i++){
    my $x=substr($cigar,$i,1);
    my $nt="";

    # If this is the beginning of a read,
    # then get the mapping quality and advance to the next nt.
    if($x eq '^'){
      $i++;
      my $mappingQuality=substr($cigar,$i,1);
      $i++;
      $x=substr($cigar,$i,1);
    } 
    
    # If this is the end of a read, 
    # then mark it and move on.
    if ($x eq '$'){
      my $startEnd=1;
      next;
    }

    # If this is an insertion, 
    # then the format is +Nn... where N is the number inserted
    # and n... is the nucleotide(s).
    if ($x eq '+'){
      $i++;
      die "Insertion shown in mpileup, but the length was not given" if(substr($cigar,$i,1)!~/(\d+)/);
      my $lengthOfInsertion=$1;
      my $digitsLen=length($lengthOfInsertion);

      $i+=$digitsLen; # advance to the actual insertion
      $x=substr($cigar,$i,$lengthOfInsertion); # the insertion
      $i+=$lengthOfInsertion; # advance past the insertion
    }
    
    # If this is not the beginning or end of a read,
    # then grab the nt. This would also grab deletions.
    if ($x eq '.' || $x eq ','){
      my $pos=$$bamField{pos};
      my $contig=$$bamField{contig};
      $x=$$refBase{$contig}{$pos};
      #die Dumper [$x,\@base,$bamField];
    } else {
      $nt=$x;
    }
    push(@base,$nt);
  }
  return \@base;
}

sub usage{
  "Creates a vcf from a sorted bam file.
  Usage: $0 file.sorted.bam > out.vcf
  --min-coverage 10 Min depth at a position
  --min-frequency 0.75 Min needed for majority
  -ref reference.fasta (optional)
  --numcpus 1
  --unsorted Produces streaming output, but unsorted due to thread race conditions
  --variants-only Do not print invariant sites
  -mpileup '-q 1' Send options to mpileup 'samtools mpileup' for additional help
  "
}
