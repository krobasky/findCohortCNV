# Licensed under a Creative Commons Attribution 3.0 Unported License.


sub inCNV() {

    $info = $_[0]; $median = $_[1]; $cov = $_[2]; 
    $sample = $_[3]; $chr = $_[4]; $pos = $_[5];
    $cnvtype = $_[6];
    @inFlds = @$inFlds_ref;

    if($info eq '.') {
	return "homdel";
    } elsif($median > $hetdelMinMedianThreshold  && $cov > $hetdelMinNormCovThreshold) {
	$fracMedian = int(100.0 * ($median - $cov)/$median);
	if($fracMedian > $hetdelMedianThreshold) {
	    if($printAllHetDelLoci == 1) {
		print $sample."\t$chr:$pos "; 
		print "--$cov\t$fracMedian\t$median\n";
	    }
	    return "hetdel";
	}
    }

    return "no";
}

sub isContiguous() {
    $chr = $_[0]; $locus = $_[1], 
    $cnvchr = $_[2]; $cnvend = $_[3]; $hiccup = $_[4]; 

    # handle hiccups
    $incr = $hiccup > 1 ? $hiccup : 1;
    return ($chr eq $cnvchr && $locus == $cnvend + $incr) ;
}


# -----State transitions:----

sub goStateOpenCNV() {
    $chr = $_[0]; $locus = $_[1]; $cnvType = $_[2];

    return ($chr, $locus-1, $locus, $stateOpenCNV, $cnvType);
}

sub goStateCloseCNV() {
    $chr = $_[0]; $start = $_[1]; $end = $_[2]; 
    $pchr_ref = $_[3]; $pstart_ref = $_[4]; $pend_ref = $_[5]; $ptag_ref = $_[6]; $sample = $_[7]; $cnvtype = $_[8];
    
    printCNV($chr, $start, $end, $pchr_ref, $pstart_ref, $pend_ref, $ptag_ref, $sample, $cnvtype);

    return $stateCloseCNV;
}

sub goStateCNV() {
    $end = $_[0]; $incr = $_[1];

    $end += $incr;
    
    return ($end, $stateCNV);
}

sub goStateCloseOpenCNV() {
    $chr = $_[0]; $start = $_[1]; $end = $_[2]; 
    $pchr_ref = $_[3]; $pstart_ref = $_[4]; $pend_ref = $_[5]; $ptag_ref = $_[6]; $sample = $_[7]; $cnvtype = $_[8];

    printCNV($chr, $start, $end, $pchr_ref, $pstart_ref, $pend_ref, $ptag_ref, $sample, $cnvtype ); 

    return ($chr, $locus-1, $locus, $stateCloseOpenCNV);
}

sub goStateNotCNV() {
    return $stateNotCNV;
}

#-----end state transitions -----


sub getPrimers {
    $primersfile = $_[0];

    my @tchr, @tstart, @tend, @ttag;

    @tchr=(); @tstart=(); @tend=(); @ttag=();

    open( PRIMERS, "$primersfile" ) or die ($!." $primersfile");
# read in primer interval coordinates:
    $i = 0;
    while (<PRIMERS>) {
	($tchr[$i], $tstart[$i], $tend[$i], $ttag[$i]) = split('\t', $_);
	chomp($ttag[$i]);
	$i++;
    }

    return (\@tchr, \@tstart, \@tend, \@ttag);
}

# return an array of indices
# representing all the primers that contain some portion of the interval.
sub getTiles {
    $chr = $_[0];    $start = $_[1];     $end = $_[2];
    $tchr_ref = $_[3]; $tstart_ref = $_[4]; $tend_ref = $_[5]; $ttag_ref = $_[6];
    @tchr = @$tchr_ref; @tstart = @$tstart_ref; @tend = @$tend_ref; @ttag_ref = @$ttag_ref;

    $j = 0;
    for $i (0..$#tstart) {
	if($tchr[$i] == $chr &&
	   ( $tstart[$i] <= $start && $start <= $tend[$i] ||
	     $tstart[$i] <= $end   && $end   <= $tend[$i] ||
	     $start <= $tstart[$i] && $tend[$i] <= $end  )
	   ) {
	    $p[$j] = $i;
	    $j++;
	}
    }

    return @p;
}

# Reformat how gaps are printed here:
sub printCNV {
    $cnvchr = $_[0]; $cnvstart = $_[1]; $cnvend = $_[2];
    $tchr_ref = $_[3]; $tstart_ref = $_[4]; $tend_ref = $_[5]; $ttag_ref = $_[6];
    $sample = $_[7]; $cnvtype = $_[8];


    # start is reported here as 0-based, but its not that way in vcf file
    $len = ($cnvend - $cnvstart);

    $label = $len."bp";

    ($genes_ref,$exons_ref) = whichCodons($cnvchr, $cnvstart+1, $cnvend);
    @genes=@$genes_ref; @exons=@$exons_ref;
    if($#exons > 0 ) {
	$codonStr = "$genes[0]($exons[0])";
	for $i (1..$#exons) {
	    $codonStr = $codonStr . ",$genes[$i]($exons[$i])";
	}
    } else {
	$codonStr = "NOT CODING";
    }
    $label = $label . " {$codonStr}";

    # Primers are tiled sometimes, so this locus may belong to more than one
    @targets = getTiles($cnvchr, $cnvstart+1, $cnvend,
			$tchr_ref, $tstart_ref, $tend_ref, $ttag_ref);


    for $i (0..$#targets) {
	if($cnvstart+1 == $tstart[$targets[$i]]) {
	    $label = $label . ", [3p--".@$ttag_ref[$targets[$i]]."]";
	} 
	if ($cnvend == $tend[$targets[$i]]) {
	    $label = $label . ", [5p--".@$ttag_ref[$targets[$i]]."]";
	}
	$label = $label .",".@$ttag_ref[$targets[$i]];
    }


    print "$samples[$sample]\t$cnvchr\t$cnvstart\t$cnvend\t$cnvtype\t$label\n";
}


sub getIntervalEnds {
    $targetfile = $_[0];  $printTargets = $_[1];
    
#
# Get the interval ends from excluded intervals ROD:
#
    
# INCLUDED intervals are in between last excluded interval and next excluded interval:
    $i = 0;
    $lastchr = "";
    open( TARGETS, "<", "$targetfile" )  or die $!;
    while( <TARGETS> ) {
	chomp($_);
	($chr, $range) = split(':', $_);
	($start, $end) = split('-',$range);
	if($chr ne $lastchr) { # current excluded interval starts a new chromosome, so there's no interval to include
	    next;
	}
	$intervalEnds[$i][0] = $chr;
	$intervalEnds[$i][1] = $start; # the end of the interval to include = the start of the interval to exclude
	$i++;
    } continue {
	$lastchr = $chr;
	$lastend = $end;
    }
    
    if($printTargets == 1) {
	for $i (0..$#intervalEnds) {
	    print $intervalEnds[$i][0].":".$intervalEnds[$i][1]."\n";
	}
    }
    close( TARGETS ) or die $!;
    return @intervalEnds;
}


sub parseCoverage() {
    $inFlds_ref = $_[0];
    @flds = @$inFlds_ref;

    $data = $flds[9];
    chomp($data);
    if($data eq "./."){
	$cov = 0;
    } else {
	$labels = $flds[8];
	@flds = split(':',$labels);
	if($flds[1] eq "DP") {
	    @flds = split(':',$data);
	    $cov = $flds[1];
	} elsif ($flds[2] eq "DP") {
	    @flds = split(':',$data);
	    $cov = $flds[2];
	} else {
	    print "ERROR - $sample\t$chr[$sample]:$pos[$sample]\tCan't find DP\n";
	}
    }
    return $cov;
}


sub getCodons {

    @chr=(); @starts=(); @ends=(); @name=(); @geneName=();

# read in codon interval coordinates:
    $i = 0;
$codonsfile = "$rootdir/proj/$proj/config/codons.txt";
open(CODONS, "<", $codonsfile ) or die ($!." $codonsfile");


    while (<CODONS>) {
	($chr[$i], $starts[$i], $ends[$i], $name[$i], $geneName[$i]) = split('\t', $_);
	chomp($geneName[$i]);
	$i++;
    }

    return (\@chr, \@starts, \@ends, \@name, \@geneName);
}

sub whichCodons {
    $chr = $_[0]; $start = $_[1]; $end = $_[2]; $codonsfp = $_[3];

    ($cchr_ref, $cstarts_ref, $cends_ref, $cname_ref, $cgeneName_ref ) = getCodons();
    @cchr = @$cchr_ref; @cstarts = @$cstarts_ref; @cends = @$cends_ref; @name = @$cname_ref; @cgeneName = @$cgeneName_ref; 

    @genes=(); @exons=();
    $exon = 0; $str = ""; $i = 0;
    for $gene (0..$#cstarts) {
	if($cchr[$gene] eq $chr) {
	    @starts = split(',', $cstarts[$gene]);
	    @ends = split(",",$cends[$gene]);
	    for $exon (0..$#starts) {
		if($starts[$exon] <= $start && $start <= $ends[$exon] ||
		   $starts[$exon] <= $end && $end <= $ends[$exon] ||
		   $start <= $starts[$exon] && $ends[$exon] <= $end   ) {

		    $genes[$i] = $geneName[$gene];
		    $exons[$i] = $exon+1;
		    $i++
		}
	    }
	}
    }
    return (\@genes, \@exons);
}

return 1;
