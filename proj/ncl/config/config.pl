# fraction of total mapped reads per sample
# in same order as samples in samples.cfg
# given numReads = number of mapped reads (taken from "#   Unique Alignment: numReads" in $proj.stat.txt)
# for each sample x, fracMax[x] = max(numReads[0..totSamples]) / numReads[x]

# nclUnfiltered:
# max = 44,892,778 [sample 3]
#
#sample  all mapped
#1   17,523,480
#2   15,844,728
#3   44,892,778
#4   27,030,573
#5   24,085,406
#6   23,235,728
#7   24,080,095
#12  18,610,069

@fracMax = (2.5618643100571,2.8332943298238,1.0,1.6608148854262,1.8638995747051,1.9320581649088,1.8643106682096,2.4122843391929);
$controlIdx = 7; # control is 12, or 7

$hetdelMedianThreshold = 50.0; # 40% catches more (less stringent), 60% catches fewer hetdels
$hetdelMinMedianThreshold = 500.0; # median has at least this many reads to be counted
$hetdelMinNormCovThreshold = 5.0; # normalized coverage has at least this many reads to be counted
$hetdelMinDist = 2; # Tolerance for how many non-hetdels in a row wil be counted before closing the run

return 1;
