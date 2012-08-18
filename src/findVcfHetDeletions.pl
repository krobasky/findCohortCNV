#!/usr/local/bin/perl
#
# Licensed under a Creative Commons Attribution 3.0 Unported License.
#

#
# Purpose:
#
# Given a series of vcf files, one per sample, use the relative
# coverage in the population to find heterozygous deletions in
# specific samples. This application uses a state-machine to identify
# gaps in coverage (homozygous deletions) along with the het-deletions
# in a single pass. Heuristics are implemented to accommodate for
# up-stream data analysis software anomalies.
#
# Notes:
#
# This program uses configuration by convention; project directory
# structure and configuration file names are implied.
#
# Inputs:
#
#  <root> : root directory
#  <proj-name> : name of project for this instance of runs
#  <root>/proj/<proj-name>/config/samples.cfg : space-delimited list of integer sample IDs
#  <root>/proj/<proj-name>/config/primers.bed : indicates amplicon start and end loci
#  <root>/proj/<proj-name>/config/codons.txt : tab-delimited codon annotations, format:
#     <chr> <start-loci-list> <end-loci-list> <gene-name> <refseq-name>
# e.g.: chr1	40538381,40542513,40544231,40546068,40555081,40562786,	40539855,40542585,40544330,40546159,40555184,40563142,	PPT1	NM_001142604
#
#  <root>/src/findVcfHetDeletions-config.pl : application defaults
#  <root>/config/config.pl: system-level configuration
#  <root>/proj/<proj-name>/config/config.pl : project-specific configuration (over-rides system config)
#  <root>/proj/<proj-name>/data/<sample-i>/<proj-name>-realigned-targetCalls.vcf: 
#                 VCF-formatted file containing all calls and coverage on all target loci, both reference and variant 
#                  (e.g., GATK can be used with parameters like the following to create this file:
#                    java -jar $gatkdir/GenomeAnalysisTK.jar -I realigned.bam -R ref.fa -T UnifiedGenotyper  
#                               -D dbsnp_132.hg19.vcf   --excludeIntervals $projdir/config/exclude.intervals 
#                               -o <root>/proj/<proj-name>/data/<sample-i>/<proj-name>-realigned-targetCalls.vcf
#                               -glm BOTH -gt_mode DISCOVERY 
#                               -out_mode EMIT_ALL_SITES
#                               -deletions -1 -metrics metrics.txt -mbq 0 -minIndelCnt 4 -dcov 500000
# 

# $debug = 1;

$ARGC = $#ARGV + 1;
if($ARGC != 2 ) {
    print "usage: findVcfHetDeletions.pl <root> <proj-name>\n";
    print "       e.g.: ./findVcfHetDeletions.pl /home/username/findCohortCNV ncl\n";
    die;
}

$rootdir = $ARGV[0];
$proj = $ARGV[1];

$projdir = "$rootdir/proj/$proj";
$samplesfile = "$projdir/config/samples.cfg"; 
$targetfile = "$projdir/config/exclude.intervals";

require "$rootdir/src/findVcfHetDeletions-config.pl";
require "$rootdir/src/findVcfHetDeletions-util.pl";

require "$rootdir/config/config.pl";
require "$projdir/config/config.pl";

open( SMP, "<", "$samplesfile" )  or die $!;
# get sample filenames
while( <SMP> ) {
    chomp($_); 
    @samples = (@samples, split( '\s', $_ ));
}
close( SMP ) or die $!;
for $i (0..$#samples) {
    # open files
    open($vcfFiles[$i], "<", "$projdir/data/$samples[$i]/$proj-realigned-targetCalls.vcf" ) or die($!." $i:$projdir/data/$samples[$i]/$proj-realigned-targetCalls.vcf\n");
}

$primersfile = "$rootdir/proj/$proj/config/primers.bed";
($pchr_ref, $pstart_ref, $pend_ref, $ptag_ref ) = getPrimers($primersfile);

use Switch 'Perl5', 'Perl6';

# A CNV is a run of coverage that is multiple of the median.
# However, runs can span chromosome and interval borders,
# so more than just a simple 'CNV/noCNV' state machine is necessary
# for splitting CNVs out of runs that cross these boundaries

print "0-based coordinates for CNVs, interval ends marked with '>>' and '<<'\n";
for $sample (0..$#vcfFiles) {
    $state[$sample] = $stateNotCNV;
    $cnvstart[$sample] = 0;$cnvend[$sample] = 0;$cnvchr[$sample] = "";
    $cnvtype[$sample] = "";
    $hiccup[$sample] = 0;
    $chrLast[$sample] = "";
    $locusLast[$sample] = 0;
}

# From config files:
print "deletion threshold = $hetdelMedianThreshold\n";
print "normalized control reads threshold = $hetdelMinMedianThreshold\n"; 
print "normalized sample reads threshold = $hetdelMinNormCovThreshold\n"; 

while (!eof($vcfFiles[0])) { # check for eof on any file; all should reach eof at the same time

    # read in all the samples for the current locus for the purpose of normalization
  for $sample (0..$#vcfFiles) {
      $fp = $vcfFiles[$sample];
      $line = <$fp>;
      chomp($line);
      if(substr($line,0,1) eq "#") { next; } # skip comments

      $chrLast[$sample] = $chr[$sample];
      $locusLast[$sample] = $locus[$sample];
      @inFlds = split( '\t', $line);
      $chr[$sample] = $inFlds[0];
      $locus[$sample] = $inFlds[1];
      $info[$sample] = $inFlds[7];

      # gatk sometimes prints double lines for same loci when using EMIT_ALL_LOCI with indel detection;
      # just take the first line in these cases.
      # Skip lines in this sample until you get the next location (but use the coverage from the first, non-indel call)
      while($chrLast[$sample] == $chr[$sample] && $locusLast[$sample] == $locus[$sample]) { 
	  if($printDuplicateLines) {
	      print "DUP[$sample]($chr[$sample]:$locus[$sample]) $line\n\n";
	  }
	  $line = <$fp>; 
	  @inFlds = split( '\t', $line);
	  $chr[$sample] = $inFlds[0];
	  $locus[$sample] = $inFlds[1];
	  $info[$sample] = $inFlds[7];
      }


      if($printStateChanges[$sample]) {
	  print "$sample*($stateStrings[$state[$sample]]) read: $line";
      }

      # don't let 'hiccups' cause transition from CNV
      $hiccup[$sample] = 0;
      if($chrLast[$sample] == $chr[$sample] && 
	 ($locus[$sample] - $locusLast[$sample]) > 1 &&
	 ($locus[$sample] - $locusLast[$sample]) < $hiccupThreshold) {
	  $hiccup[$sample] = ($locus[$sample] - $locusLast[$sample]);
      }

      $cov[$sample] = parseCoverage(\@inFlds);
      $normCov[$sample] = $cov[$sample] * $fracMax[$sample]; # fracMax is in the project's config file
      
  }


  if(substr($line,0,1) eq "#") { next; } # skip comments

  # sort coverage values to get median
  #  @sortedNormCov = sort{$a <=> $b} @normCov;
  #  $median = $sortedNormCov[int($#normCov/2)];
  $median = $normCov[$controlIdx]; # xxx just use control as 'valid' for now


  # perform state transitions
  for $sample (0..$#vcfFiles) {

      # check that the current line is in synch with the control:
      if($sample != $controlIdx && 
	 ($chr[$sample] != $chr[$controlIdx] || $locus[$sample] != $locus[$controlIdx]) ) { 
	  print "WARNING: $sample($chr[$sample]:$locus[$sample]) out of synch with $controlIdx($chr[$controlIdx]:$locus[$controlIdx]) \n"; 
      }

      if($printStateChanges[$sample]) {
	  $ret = inCNV($info[$sample], $median, $normCov[$sample], $sample, $chr[$sample], $locus[$sample], $type[$sample]);
	  print "$sample ($chr[$sample]:$locus[$sample]) --> [$ret] ($info[$sample]) \n";
      }

      switch ($state[$sample]) {
	  case($stateNotCNV) {
	      $ret = inCNV($info[$sample], $median, $normCov[$sample], $sample, $chr[$sample], $locus[$sample], $type[$sample]);
	      if($ret ne "no") {
		  ($cnvchr[$sample], $cnvstart[$sample], $cnvend[$sample], $state[$sample], $type[$sample]) = 
		      goStateOpenCNV($chr[$sample], $locus[$sample], $ret);
		  if($debug == 1) {
		      print ">> to openCNV ( [$sample] 6: $inFlds[6], 7: $inFlds[7], 8: $inFlds[8]) \n";
		  }
	      }

	  }
	  
	  case($stateOpenCNV) {
	      $ret = inCNV($info[$sample], $median, $normCov[$sample], $sample, $chr[$sample], $locus[$sample], $type[$sample]);
	      if($ret eq "no") {
		  $state[$sample] = 
		      goStateCloseCNV($cnvchr[$sample], $cnvstart[$sample], $cnvend[$sample], $pchr_ref, $pstart_ref, $pend_ref, $ptag_ref, $sample,$type[$sample]);
	      } else {
		  if(isContiguous($chr[$sample],$locus[$sample],$cnvchr[$sample],$cnvend[$sample],$hiccup[$sample])) {
		      ($cnvend[$sample], $state[$sample]) = 
			  goStateCNV($cnvend[$sample], $incr);
		      if($debug == 1) {
			  print ">>current-> [$sample] $chr[$sample]:$locus[$sample]\n";
			  print ">> to CNV [$cnvchr[$sample]:$cnvstart[$sample]-$cnvend[$sample]]\n";
		      }
		  } else {
		      ($cnvchr[$sample], $cnvstart[$sample], $cnvend[$sample], $state[$sample]) = 
			  goStateCloseOpenCNV($cnvchr[$sample], $cnvstart[$sample], $cnvend[$sample], $pchr_ref, $pstart_ref, $pend_ref, $ptag_ref,$sample,$type[$sample]);
		      if($debug == 1) {
			  print "[$sample] $cnvchr[$sample] $cnvstart[$sample]-$cnvend[$sample] [".($cnvend[$sample] - $cnvstart[$sample])." bp]\n";
		      }
		  }
	      }
	  }
	  
	  case($stateCloseCNV) {
	      $ret = inCNV($info[$sample], $median, $normCov[$sample], $sample, $chr[$sample], $locus[$sample], $type[$sample]);
	      if($ret eq "no") {
		  $state[$sample] = 
		      goStateNotCNV();
	      } else {
		  ($cnvchr[$sample], $cnvstart[$sample], $cnvend[$sample], $state[$sample]) = 
		      goStateOpenCNV($chr[$sample], $locus[$sample]);
	      }
	  }
	  
	  case($stateCNV) {
	      $ret = inCNV($info[$sample], $median, $normCov[$sample], $sample, $chr[$sample], $locus[$sample], $type[$sample]);

	      if($ret eq "no") {
		  $state[$sample] = 
		      goStateCloseCNV($cnvchr[$sample],$cnvstart[$sample],$cnvend[$sample],$pchr_ref,$pstart_ref,$pend_ref,$ptag_ref,$sample,$type[$sample]);
		  if($debug == 1) {
		      print "[$sample] $cnvchr[$sample] $cnvstart[$sample]-$cnvend[$sample] [".($cnvend[$sample] - $cnvstart[$sample])." bp]\n";
		  }
	      } else {
		  if($debug == 1) {
		      print ">> [$sample] IN gap [$cnvchr[$sample]:$cnvstart[$sample]-$cnvend[$sample]], $chr[$sample]:$locus[$sample] == $cnvend[$sample]+1\n";
		  }

		  if(isContiguous($chr[$sample],$locus[$sample],$cnvchr[$sample],$cnvend[$sample],$hiccup[$sample])) {
		      ($cnvend[$sample], $state[$sample]) = 
			  goStateCNV($cnvend[$sample], $incr);
		  } else {
		      ($cnvchr[$sample], $cnvstart[$sample], $cnvend[$sample], $state[$sample]) = 
			  goStateCloseOpenCNV($cnvchr[$sample],$cnvstart[$sample],$cnvend[$sample],$pchr_ref,$pstart_ref,$pend_ref,$ptag_ref,$sample,$type[$sample]);
		      if($debug == 1) {
			  print " [$sample] $chr[$sample] $cnvstart[$sample]-$cnvend[$sample] [".($cnvend[$sample] - $cnvstart[$sample])." bp]\n";
		      }
		  }
	      }
	  }
	  
	  case($stateCloseOpenCNV) {
	      $ret = inCNV($info[$sample], $median, $normCov[$sample], $sample, $chr[$sample], $locus[$sample], $type[$sample]);
	      if($ret eq "no") {
		  $state[$sample] = 
		      goStateCloseCNV($cnvchr[$sample], $cnvstart[$sample], $cnvend[$sample], $pchr_ref, $pstart_ref, $pend_ref, $ptag_ref,$sample,$type[$sample]);
		  if($debug == 1) {
		      print "[$sample]  $cnvchr[$sample] $cnvstart[$sample]-$cnvend[$sample] [".($cnvend[$sample] - $cnvstart[$sample])." bp]\n";
		  }
	      } else {
		if(isContiguous($chr[$sample],$locus[$sample],$cnvchr[$sample],$cnvend[$sample],$hiccup[$sample])) {
		    ($cnvend[$sample], $state[$sample]) = 
			goStateCNV($cnvend[$sample], $incr);
		} else {
		    ($cnvchr[$sample], $cnvstart[$sample], $cnvend[$sample], $state[$sample]) = 
			goStateCloseOpenCNV($cnvchr[$sample],$cnvstart[$sample],$cnvend[$sample],$pchr_ref,$pstart_ref,$pend_ref,$ptag_ref,$sample,$type[$sample]);
		    if($debug == 1) {
			print "[$sample] $cnvchr[$sample] $cnvstart[$sample]-$cnvend[$sample] [".($cnvend[$sample] - $cnvstart[$sample])." bp]\n";
		    }

		}
	      }
	  }
	  
	  else { die "($state[$sample])Error\n"; }
      } # end state machine

  } # for each sample
} # for each line (in any file)
