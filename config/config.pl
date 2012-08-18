#
# Licensed under a Creative Commons Attribution 3.0 Unported License.
#

# These are system-level variable values that can be over-ridden in project level config files on a per-project basis

$hetdelMedianThreshold = 45.0; # 40% catches more (less stringent), 60% catches fewer hetdels
$hetdelMinMedianThreshold = 500.0; # median has at least this many reads to be counted
$hetdelMinNormCovThreshold = 5.0; # normalized coverage has at least this many reads to be counted
$hetdelMinDist = 2; # Tolerance for how many non-hetdels in a row wil be counted before closing the run

