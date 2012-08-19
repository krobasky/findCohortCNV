# hetdel calling

$hetdelMedianThreshold = 50.0; # 40% catches more (less stringent), 60% catches fewer hetdels
$hetdelMinMedianThreshold = 300.0; # median has at least this many reads to be counted
$hetdelMinNormCovThreshold = 5.0; # normalized coverage has at least this many reads to be counted
$hetdelMinDist = 100; # Tolerance for how many non-hetdels in a row wil be counted before closing the run

# @fracMax = (1.0,1.0,1.0);
@fracMax = (1.315,1.315,1.0);



return 1;
