#
# Licensed under a Creative Commons Attribution 3.0 Unported License.
#

# application-level defaults can be overriden by system and project config files at the system and project levels, respectively.

# Every loci in the VCF should be called either reference or variant.
$hiccupThreshold = 5; # Occassionally, GATK will inexplicably skip a locus; this is deemed a "hiccup"
$hetdelMedianThreshold = 35.0; # 40% catches more (less stringent), 60% catches fewer hetdels
$hetdelMinMedianThreshold = 20.0; # median has at least this many reads to be counted
$hetdelMinNormCovThreshold = 5.0; # normalized coverage has at least this many reads to be counted
$hetdelMinDist = 2; # Tolerance for how many non-hetdels in a row wil be counted before closing the run

# $printAllHetDelLoci = 1;
# $printDuplicateLines = 1;

# define the state constants
$stateNotCNV = 0;
$stateCNV = 1;
$stateCloseCNV = 2;
$stateOpenCNV = 3;
$stateCloseOpenCNV = 4;
@stateStrings = ("NotCNV", "CNV", "Close", "Open", "CloseOpen");


return 1;
