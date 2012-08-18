# hetdel calling

$hetdelMedianThreshold = 40.0; # 40% catches more (less stringent), 60% catches fewer hetdels
$hetdelMinMedianThreshold = 0.0; # median has at least this many reads to be counted
$hetdelMinNormCovThreshold = 0.0; # normalized coverage has at least this many reads to be counted
$hetdelMinDist = 1; # Tolerance for how many non-hetdels in a row wil be counted before closing the run

@fracMax = (1.0,1.0,1.0);
# @fracMax = (0.76,0.76,1.0);



return 1;
