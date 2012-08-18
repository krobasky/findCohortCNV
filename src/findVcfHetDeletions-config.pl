#
# Licensed under a Creative Commons Attribution 3.0 Unported License.
#

# application-level defaults can be overriden by system and project config files at the system and project levels, respectively.
$hiccupThreshold = 4; # can be overridden in both system and project configfiles
$hetdelMedianThreshold = 35.0;
$hetdelMinMedianThreshold = 20.0; 
$hetdelMinNormCovThreshold = 5.0;
$hetdelMinDist = 2;

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
