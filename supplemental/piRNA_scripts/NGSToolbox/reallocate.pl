#!/usr/bin/perl
use Getopt::Long;
$|=1;
$script=$0;
$script=~s/^.+[\\\/]//;

$information="
 ================================ reallocate ==================================
 VERSION 1.1                                  LAST MODIFIED: 24. September 2015

 SCOPE:
 This script will process map files in order to reallocate read counts of
 multiple mapping sequences according to the transcription rate of genomic
 loci based on uniquely mapping reads. Map files must be in ELAND format
 and can be created using sRNAmapper which is provided along with the proTRAC
 software. reallocate will output a modified map file that contains two
 additional columns that refer to i) number of genomic hits of a sequence
 and ii) read counts that are assigned to this locus. proTRAC 2.0.5 and later
 versions accept this format and utilize this information for cluster
 prediction. Generally, using reallocate may result in a higher amount of
 sequence reads that can be assigned to predicted piRNA clusters and may
 also result in a higher number of predicted piRNA clusters.
 
 ------------------------------------------------------------------------------
 
 USAGE:
 reallocate requires five arguments. Start the script with the following
 command: 
  perl $script mapfile[file] perimeter[bp] resolution[bp] shape[b,e,l,n]
  cutoff[value]
 
 <mapfile> refers to your input file which should be the output of the
 sRNAmapper script that comes along with proTRAC.
 <perimeter> refers to the up/downstream region that will be considered for
 calculation of transcription rate (e.g. 10000).
 <resolution> refers to the space between the centers of the windows for
 that transcription rates will be calculated (e.g. 1000).
 <shape> refers to the shape of the function that determines how reads
 are weighted according to the distance from the center of the window.
 Valid values are <b> (bell-shaped), <e> (exponential), <l> (linear) or
 <n> (no weighting).
 <cutoff> refers to the threshold value for the minimum number of assigned
 sequence reads for a locus to appear in the output file. Use -1 if you
 want to output all loci. Use 0 if you want to reject loci with no allocated
 sequence reads. Use higher thresholds as desired.
 
 The processed map will be saved as e.g
 mapfile.map.weighted-10000-1000-b-0
 (mapfile.weighted-perimeter-resolution-shape-cutoff)

 ------------------------------------------------------------------------------

 OPTIONS:
 -h
 -help
             Will print this information.

 ------------------------------------------------------------------------------

 (c) David Rosenkranz
 Institute of Anthropology, small RNA group
 Johannes Gutenberg University Mainz, Germany
 Contact: rosenkranz\@uni-mainz.de
 http://www.smallRNAgroup-mainz.de/
 
 ##############################################################################



";

GetOptions
	(
	"h"=>\$print_info,
	"help"=>\$print_info,
	);

if($print_info)
	{
	print$information;
	exit;
	}

###   parse arguments   ###
$unweighted_map=$ARGV[0];
$orbit=$ARGV[1];
$resolution=$ARGV[2];
$sw_shape=$ARGV[3];
$cutoff=$ARGV[4];

unless(@ARGV==5)
	{
	print"\nERROR! Invalid number of arguments.\n$information";exit;
	}
unless(-e $unweighted_map)
	{
	print"\nERROR! $unweighted_map does not exist.\n";exit;
	}
if($orbit!~/^\d+$/)
	{
	print"\nERROR! Value of orbit is '$orbit' but must be an integer.\n";exit;
	}
if($resolution!~/^\d+$/)
	{
	print"\nERROR! Value of resolution is '$resolution' but must be an integer.\n";exit;
	}
if($resolution>=$orbit)
	{
	print"\nERROR! Value for perimeter ($orbit) must exceed value for resolution ($resolution).\n";exit;
	}
if($sw_shape!~/^[beln]$/)
	{
	print"\nERROR! Value for weighting of reads within the orbit must be [b,e,l,n].\nb:bell-shaped, e:exponential, l:linear, n:no weighting\n";exit;
	}
if($cutoff!~/^-?(\d+\.)?\d$/)
	{
	print"\nERROR! Cutoff must be a numeric value. E.g. -1, 0, 0.5 etc.\n";exit;
	}

$|=1;
###   calculate weighting factors according to distance from hit   ###
@factors=();
foreach$up_down(0..$orbit)
	{
	if($sw_shape eq 'l')
		{
		# linear regression
		push(@factors,($orbit-$up_down)/$orbit);	
		}
	elsif($sw_shape eq 'e')
		{
		# exponential regression
		push(@factors,(($orbit-$up_down)/$orbit)**2);
		}
	elsif($sw_shape eq 'b')
		{
		# bell-shaped regression
		push(@factors,1-((cos(((($orbit-$up_down)/$orbit)*3.1415926535897932384626433832795))+1)/2));	
		}
	elsif($sw_shape eq 'n')
		{
		# no regression
		push(@factors,1);
		}
	}
@reverse_factors=reverse@factors;
shift@factors;
@factors=(@reverse_factors,@factors);


###   read map file and save genomic hits and read counts for each sequence   ###
%hits=();
%reads=();
%loc2id=();
%id2loc=();
$id=0;
$hits=0;
$prev_loc="";
print"\nSave genomic hits and reads for each sequence...";
open(UNWEIGHTED,$unweighted_map)||die print"\nCould not open input file $unweighted_map.\n$!\n\n";
while(<UNWEIGHTED>)
	{
	if($_!~/^trans_id/&&$_!~/^#/)
		{
		$hits++;
		@d=split("\t",$_);
		if($prev_loc ne $d[0])
			{
			$id++;
			$loc2id{$d[0]}=$id;
			$id2loc{$id}=$d[0];
			print" $d[0]";
			}
		$hits{$d[4]}++;
		$reads{$d[4]}=$d[3];
		$prev_loc=$d[0];
		}
	}
close UNWEIGHTED;


###   load location and coordinates of uniquely mapping sequences to memory   ###
open(UNWEIGHTED,$unweighted_map)||die print"\nCould not open input file $unweighted_map.\n$!\n\n";
$c=0;
$t0=time;
$p="0.00 %";
print"\nSave coordinates of uniquely mapping sequences... $p";
while(<UNWEIGHTED>)
	{
	$t=time;
	if($t-$t0>2)
		{
		foreach(1..length$p)
			{
			print"\b";
			}
		$perc=(int((($c/$hits)*10000)+0.5))/100;
		$p="$perc %  ";
		print$p;
		$t0=time;
		}
	
	$_=~s/\s*$//;
	if($_!~/^trans_id/&&$_!~/^#/)
		{
		$c++;
		@d=split("\t",$_);
		if($hits{$d[4]}==1)
			{
			$umap{$loc2id{$d[0]}}{$d[6]}{$d[1]}=$reads{$d[4]};
			}
		}
	}
close UNWEIGHTED;
foreach(1..length$p)
	{
	print"\b";
	}
print"100.0 %";


###   calculate and save scores for local transcription of sites corresponding to multiple-mapping sequences   ###
$c=0;
$t0=time;
$p="0.00 %";
%scores=();
%total_score=();
print"\nCalculate total transcription scores... $p";
open(UNWEIGHTED,$unweighted_map);
$last_loc="";
$last_coord=0;
while(<UNWEIGHTED>)
	{
	$t=time;
	if($t-$t0>2)
		{
		foreach(1..length$p)
			{
			print"\b";
			}
		$perc=(int((($c/$hits)*10000)+0.5))/100;
		$p="$perc %  ";
		print$p;
		$t0=time;
		}
	
	$_=~s/\s*$//;
	if($_!~/^trans_id/&&$_!~/^#/)
		{
		$c++;
		@d=split("\t",$_);
		if($hits{$d[4]}>1)
			{
			# calculate sw transcription only if distance to last sw exceeds defined resolution
			if($d[0]ne$last_loc||($last_coord+$resolution)<=$d[1])
				{
				$i=-1;
				$score=0;
				foreach$up_down(-($orbit-($resolution/2))..$orbit+($resolution/2))
					{
					$i++;
					if($umap{$loc2id{$d[0]}}{$d[6]}{$d[1]+$up_down})
						{
						$score+=$umap{$loc2id{$d[0]}}{$d[6]}{$d[1]+$up_down}*$factors[$i];
						}
					}
				$last_loc=$d[0];
				$last_coord=$d[1];
				}
			$scores{$d[4]}.="$score ";
			$total_score{$d[4]}+=$score;
			}
		}
	}
close UNWEIGHTED;
foreach(1..length$p)
	{
	print"\b";
	}
print"100.0 %";


###   create weighted map -> add genomic hits and locus reads  ###
$c=0;
$t0=time;
$p="0.00 %";
print"\nCreate weighted map... $p";
open(UNWEIGHTED,$unweighted_map)||die print"\nCould not open input file $unweighted_map.\n$!\n\n";
open(WEIGHTED,">$unweighted_map.weighted-$orbit-$resolution-$sw_shape-$cutoff")||die print"\nCould not create output file $unweighted_map.weighted-$orbit-$resolution-$sw_shape-$cutoff.\n$!\n\n";
while(<UNWEIGHTED>)
	{
	$t=time;
	if($t-$t0>2)
		{
		foreach(1..length$p)
			{
			print"\b";
			}
		$perc=(int((($c/$hits)*10000)+0.5))/100;
		$p="$perc %  ";
		print$p;
		$t0=time;
		}
	
	$_=~s/\s*$//;
	if($_!~/^trans_id/&&$_!~/^#/)
		{
		$c++;
		@d=split("\t",$_);
		if($hits{$d[4]}>1)
			{
			$scores{$d[4]}=~s/[^ ]+ //;
			$score=$&;
			$score=~s/ //;
			if($total_score{$d[4]}>0)
				{
				$locus_reads=$d[3]*$score/$total_score{$d[4]};
				}
			else
				{
				$locus_reads=$d[3]/$hits{$d[4]};
				}
			if($locus_reads>$cutoff)
				{
				print WEIGHTED"$_\t$hits{$d[4]}\t$locus_reads\n";
				}
			}
		else
			{
			print WEIGHTED"$_\t$hits{$d[4]}\t$d[3]\n";
			}
		}
	else
		{
		print WEIGHTED$_;
		}
	}
close WEIGHTED;
close UNWEIGHTED;
foreach(1..length$p)
	{
	print"\b";
	}
print"100.0 %\nProcessing finished.";
exit;