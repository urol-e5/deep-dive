# The following array must contain a list of proTRAC output folders.
# Listed proTRAC runs must base on the same genome. Otherwise,
# merging cluster coordinates makes no sense!

@proTRAC_folders=
	(
	"proTRAC_sRNA-ACR-140_preproc.map.weighted-10000-1000-b-0_2024y6m20d14h55m16s",
    "proTRAC_sRNA-ACR-145_preproc.map.weighted-10000-1000-b-0_2024y6m20d15h1m25s",
    "proTRAC_sRNA-ACR-150_preproc.map.weighted-10000-1000-b-0_2024y6m20d15h8m8s",
    "proTRAC_sRNA-ACR-173_preproc.map.weighted-10000-1000-b-0_2024y6m20d15h13m58s",
    "proTRAC_sRNA-ACR-178_preproc.map.weighted-10000-1000-b-0_2024y6m20d15h19m33s",
	);


$min_dist=10000; # This is the minimum distance between two independent piRNA clusters. Clusters closer to each other will be merged.
$output_table="merged_cluster_coordinates.txt";
$|=1;

%clusters=();
foreach$folder(@proTRAC_folders)
	{
	open(TABLE,"$folder/results.table")||print"\nCannot open $folder/results.table.\n$!\n\n";
	while(<TABLE>)
		{
		if($_=~/^Cluster/)
			{
			$_=~s/Location: [^\t]+//;
			$loc=$&;
			$loc=~s/Location: //;
			$_=~s/Coordinates: [^\t]+//;
			$coord=$&;
			$coord=~s/Coordinates: //;
			@coord=split('-',$coord);
			foreach$p($coord[0]..$coord[1])
				{
				$clusters{$loc}{$p}++;
				}
			print"\n$loc -> $coord[0] .. $coord[1]";
			}
		}
	close TABLE;
	}

open(MERGE,">$output_table");
$prev_loc="";
$cluster_id=0;
foreach$loc(sort{$a cmp $b}keys%clusters)
	{
	$prev_p=-$min_dist;
	$cluster_id++;
	$start=0;
	foreach$p(sort{$a<=>$b}keys%{$clusters{$loc}})
		{
		if($start==0)
			{
			print MERGE"$cluster_id\t$loc\t$p";
			$start=1;
			}
		if($p>($prev_p+$min_dist)&&$start==1&&$prev_p!=-$min_dist)
			{
			print MERGE"\t$prev_p\n";
			$cluster_id++;
			$start=0;
			}
		$prev_p=$p;
		}
	print MERGE"\t$prev_p\n";
	}
close MERGE;
exit;