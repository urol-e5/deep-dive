#!/usr/bin/perl
use Getopt::Long;
use Parallel::ForkManager;
use Scalar::Util qw( looks_like_number );
use Math::Random::Discrete;

# info
$script=$0;
$script=~s/^.+[\\\/]//;
$info="
=================================== ppMeter ===================================
VERSION: 0.4
LAST MODIFIED: 30. April 2018

Usage: perl $script -i input.file [-option value]

When testing for significantly different ping-pong signatures between two
conditions we strongly recommend to use at least 2 biological replicates per
condition to account for variance introduced during RNA extraction and/or
NGS library construction.

Command line options:
-------------------------------------------------------------------------------
[s] = string, e.g. a file name.
[i] = integer

-h OR -help            Will print this information.
-i OR -input [s]       sRNA map file in SAM or ELAND3 (sRNAmapper) format. Use
                       this option multiple times for multiple input files:
                       -i input_1.map -i input_2.map [...]. Input files must
                       be sorted by chromosomal position in ascending order.
-f OR -format [s]      Specify the input format. Allowed values are 'SAM' and
                       'ELAND3'. This is only required if on of the input files
                       contains less than 1000 hits.
-o OR -output [s]      Name of the output file (default='PPmeter_results.txt').
-t OR -threads [i]     Number of maximum allowed parallel threads.
-b OR -bootstraps [i]  Number of bootstrap pseudo-replicate datasets that will
                       be created for each input dataset (default=100).
-d OR -depth [i]       Number of sequence reads for each pseudo-replicate
                       (default=1000000).
-c OR -compare         Compare ping-pong signatures of pseudo-replicates across
                       different input files. This will perform a statistical
                       test (Mann-Whitney U test) to check whether the ping-
                       pong signatures are different between two input files.
                       By default all pairwise comparisons will be conducted.
                       Use the option -g1 AND -g2 to build two groups of input
                       files (e.g. each group comprising biological
                       replicates). Groups do not have to comprise the same
                       number of input files.
-g1 OR -group1 [s]     Specify files that will constitite group 1 for
                       statistical comparison of ping-pong signatures. This
                       option can only be used together with -c AND -g2.
-g2 OR -group2 [s]     Specify files that will constitite group 2 for
                       statistical comparison of ping-pong signatures. This
                       option can only be used together with -c AND -g1.
-less_memory           Will use a different algorithm for bootstrapping that
                       uses less memory but might be considerably slower.
-more_output           Temporary map files for bootstrap replicates and ping-
                       pong results for each temporary map file will not be
                       removed. Sequences producing 10 nt 5' overlaps will
                       be saved for each bootstrap replicate in a separate
                       FASTA file together with information on sequence read
                       length distribution and positional nucleotide
                       composition in an additional text file.
-silent                Less output to STDOUT.


Contact:
David Rosenkranz
Institute Organismic and Molecular Evolution, small RNA group
Johannes Gutenberg University Mainz, Germany
Web: www.smallRNAgroup.uni-mainz.de
Email: rosenkranz\@uni-mainz.de

===============================================================================
";

# defaults and options
$|=1;
$t0=time;
ZSCORE();
GetOptions
	(
	"h|help"=>\$print_info,
	"i|input=s"=>\@files,
	"f|format=s"=>\$file_format,
	"o|output=s"=>\$output_file,
	"t|threads=s"=>\$threads,
	"b|bootstraps=i"=>\$bootstraps,
	"d|depth=i"=>\$resampling_depth,
	"c|compare"=>\$compare_with_MWU,
	"g1|group1=s"=>\@group1,
	"g2|group2=s"=>\@group2,
	"less_memory"=>\$less_memory,
	"more_output"=>\$more_output,
	"silent"=>\$silent,
	);
if(!$threads){$threads=1;}
if(!$bootstraps){$bootstraps=100;}
if(!$resampling_depth){$resampling_depth=1000000;}
if(!$output_file){$output_file="PPmeter_results.txt";}
if(@files==0)
	{
	$print_info=1;
	}
if($print_info){print$info;exit;}

# check options and files
%input_files=();
foreach$file(@files)
	{
	if(!-e$file)
		{
		print"\nERROR. $file does not exist.\n\n";
		exit;
		}
	$input_files{$file}=1;
	}
if($compare_with_MWU)
	{
	$all_vs_all==0;
	$group1_vs_group2=0;
	if(@group1>0&&@group2>0)
		{
		print"\nFiles in group 1:";
		foreach$f1(@group1)
			{
			print"\n -$f1";
			if(!$input_files{$f1})
				{
				print"\nERROR. $f1 (group 1) is not among the input files.\n\n";
				exit;
				}
			}
		print"\nFiles in group 2:";
		foreach$f2(@group2)
			{
			print"\n -$f2";
			if(!$input_files{$f2})
				{
				print"\nERROR. $f2 (group 2) is not among the input files.\n\n";
				exit;
				}
			}
		$group1_vs_group2=1;
		print"\n\n";
		}
	elsif(@group1==0&&@group2==0)
		{
		$all_vs_all=1;
		}
	else
		{
		if(@group1>0)
			{
			print"\nERROR. Group 2 comprises no files.\n\n";
			exit;
			}
		else
			{
			print"\nERROR. Group 1 comprises no files.\n\n";
			exit;
			}
		}
	}

# do not use identical file names
%filenames_used=();
foreach$file(@files)
	{
	$file_parsed=$file;
	$file_parsed=~s/^.+[\\\/]//;
	if($filenames_used{$file_parsed})
		{
		print"\nDo not use multiple files with identical names (in different folders)!\n$filenames_used{$file_parsed}\n$file\n\n";
		exit;
		}
	$filenames_used{$file_parsed}=$file;
	}

# load datasets to memory
$max_length=0;
foreach$file(@files)
	{
	$file_parsed=$file;
	$file_parsed=~s/^.+[\\\/]//;
	
	if(!$silent){print"\nRead $file...";}
	@dataset=();
	$total_reads=0;
	%processed=();
	%hits_per_seq=();
	
	# check file format (SAM or ELAND3, sorted, read counts available)
	if(!$file_format)
		{
		$sam_lines=0;
		$eland3_lines=0;
		$tested_lines=0;
		open(IN,$file)||die print"\nUnable to read map file $file.\n$!\n\n";
		$prev_ref="";
		$prev_coord=-1;
		%listed_locations_or_scaffolds=();
		while(<IN>)
			{
			next if($_=~/^trans_id/||$_=~/^[#\@]/);
			$_=~s/\s*$//;
			@d=split("\t",$_);
			if(looks_like_number($d[0])&&$d[3]=~/^\d+$/&&$d[9]=~/^[ATGCUN]+$/)
				{
				$sam_lines++;
				$hit_without_1st_colum=$_;
				$hit_without_1st_colum=~s/^[^\t]+\t//;
				if($prev_ref eq $d[2])
					{
					if($prev_coord>$d[3])
						{
						print"\nERROR. Input file $map is not sorted according to chromosome coordinates in ascending order.\n\n";
						exit;
						}
					}
				elsif($listed_locations_or_scaffolds{$d[2]}==1)
					{
					print"\nERROR. Input file $map is not sorted according to chromosome coordinates in ascending order.\n\n";
					exit;
					}
				if($hit_without_1st_colum eq $prev_hit_without_1st_colum)
					{
					print"
ERROR. Input file $map is not made from collapsed FASTA/FASTQ file.

Uncollapsed:
>read1
TTTCGAATCGCGATCGAGGAGATCGATC
>read2
TTTCGAATCGCGATCGAGGAGATCGATC
>read3
ACCGATGGCTAGACTCGATACGC

Collapsed: (headers will refer to sequence read counts)
>2
TTTCGAATCGCGATCGAGGAGATCGATC
>1
ACCGATGGCTAGACTCGATACGC

You can use the collapse tool from the NGS TOOLBOX to convert your sequence
file into a collapsed sequence file. You can download NGS TOOLBOX here:
http://www.smallrnagroup.uni-mainz.de/software/TBr2.zip
Once you have created a collapsed sequence file, repeat mapping with
the mapper of your choice and run PPmeter with the new map file.\n\n";
					exit;
					}
				$prev_ref=$d[2];
				$prev_coord=$d[3];
				$listed_locations_or_scaffolds{$prev_ref}=1;
				$prev_hit_without_1st_colum=$hit_without_1st_colum;
				}
			elsif($d[1]=~/^\d+$/&&$d[2]=~/^[ATGCUN]+$/&&looks_like_number($d[3])&&$d[4]=~/^[ATGCUN]+$/&&$d[5]=~/^\d+$/&&$d[6]=~/^[\+-]$/)
				{
				$eland3_lines++;
				if($prev_ref eq $d[0])
					{
					if($prev_coord>$d[1])
						{
						print"\nERROR. Input file $map is not sorted according to chromosome coordinates in ascending order.\n\n";
						exit;
						}
					}
				elsif($listed_locations_or_scaffolds{$d[0]}==1)
					{
					print"\nERROR. Input file $map is not sorted according to chromosome coordinates in ascending order.\n\n";
					exit;
					}
				$prev_ref=$d[0];
				$prev_coord=$d[1];
				$listed_locations_or_scaffolds{$prev_ref}=1;
				}
			else
				{
				print"
ERROR. Input file $map is not valid SAM or ELAND3 format. Note that FASTA
headers in the file that was used to create SAM or ELAND3 files must refer
to read counts like this:
>1
TACCCGTAGCTACGTCTCTTAGC
>17
TTTAGCCCCGATTCAGACTCGACT
[...]
Doing so, read counts will appear in SAM files column 1 and ELAND3 files
column 4. You can create collapsed FASTA files using the collapse tool from
the NGS toolbox which you can download here:
http://www.smallrnagroup.uni-mainz.de/software/TBr2.zip.
";
				exit;
				}
			$tested_lines++;
			if($tested_lines==1000)
				{
				if($sam_lines==1000)
					{
					$file_format='SAM';
					last;
					}
				elsif($eland3_lines==1000)
					{
					$file_format='ELAND3';
					last;
					}
				else
					{
					print"\nERROR. Unknown input file format.\n\n";
					exit;
					}
				}
			}
		}
	if($tested_lines<1000&&!$file_format)
		{
		print"\nERROR. Cannot check file format for $file. Less than 100 hits). Use the option -format to specify input file format (SAM or ELAND3).\n\n";
		exit;
		}
	close IN;
	if(!$silent){print" (=$file_format format)"};
	
	# save flags that refer to '-' or '+' strand in SAM files
	if($file_format eq'SAM')
		{foreach$bit1(0,1)
		{foreach$bit2(0,2)
		{foreach$bit3(0,4)
		{foreach$bit4(0,8)
		{foreach$bit5(0,16)
		{foreach$bit6(0,32)
		{foreach$bit7(0,64)
		{foreach$bit8(0,128)
		{foreach$bit9(0,256)
		{foreach$bit10(0,512)
		{foreach$bit11(0,1024)
		{foreach$bit12(0,2048)
			{
			$flag=$bit1+$bit2+$bit3+$bit4+$bit5+$bit6+$bit7+$bit8+$bit9+$bit10+$bit11+$bit12;
			if($bit3==4)
				{
				$sam_mapped{$flag}=0;
				}
			else
				{
				$sam_mapped{$flag}=1;
				}
			if($bit5==16)
				{
				$sam_strand{$flag}='-';
				}
			else
				{
				$sam_strand{$flag}='+';
				}}}}}}}}}}}}}}
	
	open(IN,$file)||die print"\nUnable to read map file $file.\n$!\n\n";
	while(<IN>)
		{
		next if($_=~/^trans_id/||$_=~/^[#\@]/); # skip SeqMap, sRNAmapper and SAM header lines
		$_=~s/\s*$//;
		@d=split("\t",$_);
		if($file_format eq'SAM')
			{
			next if($sam_mapped{$d[1]}!=1);
			$seq=$d[9];
			$reads=$d[0];
			}
		elsif($file_format eq'ELAND3')
			{
			$seq=$d[4];
			$reads=$d[3];
			}
		if(length$seq>$max_length) # save max. sequence length (-> maximum overlap)
			{
			$max_length=length$seq;
			}
		if(!$processed{$seq})
			{
			# load reads into array that represents their probabilities to be drawn during bootstrapping
			if(!$less_memory)
				{
				foreach$i(1..$reads)
					{
					push(@dataset,$seq);
					}
				}
			$total_reads+=$reads;
			$processed{$seq}=1;
			}
		$hits_per_seq{$seq}++;
		$reads_per_seq{$seq}=$reads;
		}
	close IN;
	undef%processed;
	$reads_per_file{$file_parsed}=$total_reads;
	if(!$silent){print" done.\nSaved $total_reads reads from input file.";}
	
	# create arrays for memory saving bootstrapping
	if($less_memory)
		{
		if(!$silent){print"\nGenerate arrays for memory efficient bootstrapping...";}
		@seq_to_draw=();
		@prob_to_draw=();
		foreach$seq(keys%reads_per_seq)
			{
			push(@seq_to_draw,$seq);
			push(@prob_to_draw,$reads_per_seq{$seq});
			}
		$math_random_discrete=Math::Random::Discrete->new
			(
			\@prob_to_draw,
			\@seq_to_draw,
			);
		if(!$silent){print" done.";}
		}
	undef%reads_per_seq;
	
	# start bootstrapping
	if(!$silent){print"\nBuild $bootstraps bootstrap pseudo-replicates with $resampling_depth reads each.";}
	$pm=Parallel::ForkManager->new($threads);
	foreach$bs(1..$bootstraps)
		{
		$pm->start and next;
		$t_bs=time;
		if(!$silent){print"\n  -Start bootstrap procedure $bs.";}

		# create pseudo-replicate
		$t0=time;
		%bs_dataset=();
		srand($bs); # make sure to use different series of random numbers for each bootstrap procedure
		if(!$less_memory)
			{
			foreach(1..$resampling_depth)
				{
				$bs_dataset{$dataset[int(rand(scalar(@dataset)))]}++;
				}
			}
		else
			{
			foreach(1..$resampling_depth)
				{
				$draw_seq=$math_random_discrete->rand;
				$bs_dataset{$draw_seq}++;
				}
			}
		$tused=time-$t0;
		if(!$silent){print"\n  -Bootstraped $resampling_depth reads (bootstrap:$bs) in $tused second(s).";}
		
		# print pseudo-replicate dataset to temporary map file
		$t0=time;
		open(TEMPIN,">$file_parsed.$bs.temp")||die print"\nUnable to create temporary file $file.$bs.temp.\n$!\n\n";
		open(ORIG,$file)||die print"\nUnable to open file $file.\n$!\n\n";
		while(<ORIG>)
			{
			$_=~s/\s*$//;
			@d=split("\t",$_);
			if($file_format eq'SAM'&&$bs_dataset{$d[9]}>0)
				{
				print TEMPIN"$d[2]\t$d[3]\t-\t$bs_dataset{$d[9]}\t$d[9]\t-\t$sam_strand{$d[1]}\n";
				}
			elsif($file_format eq'ELAND3'&&$bs_dataset{$d[4]}>0)
				{
				# columns 3  and 6 from eland3 are not necessary for pp calculation ($d[2] -> -, $d[6] -> -)
				print TEMPIN"$d[0]\t$d[1]\t-\t$bs_dataset{$d[4]}\t$d[4]\t-\t$d[6]\n";
				}
			}
		close ORIG;
		close TEMPIN;
		$tused=time-$t0;
		if(!$silent){print"\n  -Created temporary map file (bootstrap:$bs) in $tused second(s).";}
		
		# test for ping-pong signature in pseudo replicate
		$t0=time;
		@sw=();
		%overlap=();
		open(TEMPIN,"$file_parsed.$bs.temp")||die print"\nUnable to read temporary file $file.$bs.temp.\n$!\n\n";
		%pingpong_sequences=();
		while(<TEMPIN>)
			{
			@d=split("\t",$_);
			if($d[0]ne$prev_loc) # empty sliding window when location changes
				{
				@sw=();
				}
			if($d[6]=~/-/) # add hit to sliding window if current hit is on the minus strand
				{
				push(@sw,[$d[1],$d[4]]);
				while(1) # check if hit(s) must be removed from sliding window
					{
					if(@sw>1&&$d[1]>($sw[0][0])+$max_length)
						{
						shift@sw;
						}
					else
						{
						last;
						}
					}
				}
			else # count overlaps in the sliding window if current hit is on the plus strand
				{
				foreach$i(0..@sw-1)
					{
					if($sw[$i][0]+(length$sw[$i][1])>$d[1]&&$sw[$i][0]+(length$sw[$i][1])<$d[1]+length$d[4]) # check whether hits overlap
						{
						$overlap{$sw[$i][0]+(length$sw[$i][1])-$d[1]}+=(($bs_dataset{$d[4]}/$hits_per_seq{$d[4]})*($bs_dataset{$sw[$i][1]}/$hits_per_seq{$sw[$i][1]}));
						if(($sw[$i][0]+(length$sw[$i][1])-$d[1])==10)
							{
							$pingpong_sequences{$d[4]}=$bs_dataset{$d[4]};
							$pingpong_sequences{$sw[$i][1]}=$bs_dataset{$sw[$i][1]};
							}
						}
					}
				}
			$prev_loc=$d[0];
			}
		close TEMPIN;
		if(!$more_output)
			{
			unlink"$file_parsed.$bs.temp";
			open(PP_PLAYERS,">pp_players.$bs.fas");
			foreach$seq(sort{$pingpong_sequences{$b}<=>$pingpong_sequences{$a}}keys%pingpong_sequences)
				{
				print PP_PLAYERS">$pingpong_sequences{$seq}\n$seq\n";
				}
			close PP_PLAYERS;
			$ba_in="pp_players.$bs.fas";
			$ba_out="pp_players.$bs.info";
			basic_analysis();
			}
		$tused=time-$t0;
		if(!$silent){print"\n  -Checked pp signature for temporary map file (bootstrap:$bs) in $tused second(s).";}
		
		# save results to temporary file (is necessary because information gets lost when this thread is finished), calculate ping-pong Z-score according to Zhang et al. 2012 Mol Cell 44(4):572-584
		open(PP,">$file_parsed.$bs.pp")||die print"\nUnable to create temporary output file $file_parsed.$bs.pp.\n$!\n\n";
		foreach$overlap(1..$max_length)
			{
			if(!$overlap{$overlap})
				{
				$overlap{$overlap}=0;
				}
			print PP"$overlap\t$overlap{$overlap}\n";
			
			# calculate average background
			if($overlap>=1&&$overlap<=9)
				{
				$average_background+=$overlap{$overlap};
				}
			elsif($overlap>=11&&$overlap<=20)
				{
				$average_background+=$overlap{$overlap};
				}
			}
		$average_background=$average_background/19;
		
		# calculate variance
		$variance=0;
		foreach$overlap(1..$max_length)
			{
			if($overlap>=1&&$overlap<=9)
				{
				$variance+=($overlap{$overlap}-$average_background)**2;
				}
			elsif($overlap>=11&&$overlap<=20)
				{
				$variance+=($overlap{$overlap}-$average_background)**2;
				}
			}
		$variance=$variance/19;
		$stddeviation=sqrt($variance);
		
		#calculate ppZ-score
		if($stddeviation>0)
			{
			$Zscore=($overlap{10}-$average_background)/$stddeviation;
			}
		else
			{
			if($overlap{10}>0)
				{
				$Zscore="infinite";
				}
			else
				{
				$Zscore=0;
				}
			}
		$pingpong_reads=0;
		$pingpong_sequences=keys%pingpong_sequences;
		foreach$seq(keys%pingpong_sequences)
			{
			$pingpong_reads+=$pingpong_sequences{$seq};
			}
		print PP"ppZScore:\t$Zscore\nppSequences:\t$pingpong_sequences\nppReads:\t$pingpong_reads";
		close PP;
		undef%overlap;
		undef%pingpong_sequences;
		$tused=time-$t_bs;
		if(!$silent){print"\n  -Finished bootstrap procedure $bs in $tused second(s).";}
		$pm->finish;
		}
	undef%bs_dataset;
	undef%hits_per_seq;
	$pm->wait_all_children();
	}

# print results to output file
open(OUT,">$output_file")||die print"\nUnable to create output file $output_file.\n$!\n\n";
foreach$file(@files)
	{
	%avg=();
	%dev=();
	$file_parsed=$file;
	$file_parsed=~s/^.+[\\\/]//;
	
	print OUT"Ping-pong signatures for $bootstraps pseudo-replicate datasets ($resampling_depth reads each) based on $file_parsed ($reads_per_file{$file_parsed} reads)\n";
	foreach$i(1..$max_length-1,'ppZScore','ppSequences','ppReads')
		{
		print OUT"\t$i";
		}
	
	foreach$bs(1..$bootstraps)
		{
		print OUT"\n$bs";
		open(PP,"$file_parsed.$bs.pp")||die print"\nUnable to read temporary results file $file_parsed.$bs.pp.\n$!\n\n";
		while(<PP>)
			{
			$_=~s/\s*$//;
			@d=split("\t",$_);
			next if($d[0]==$max_length-1); # this is not counted as 5' overlap
			$d[0]=~s/://; # remove colom from ppZScore/ppSequences/ppReads
			if($d[1]ne'infinite') # ignore ppZScores = 'infinite'
				{
				$d[1]=0+$d[1];
				$avg{$d[0]}+=$d[1];
				}
			push(@{$all_val{$file_parsed}{$d[0]}},$d[1]);
			print OUT"\t$d[1]";
			}
		close PP;
		if(!$more_output)
			{
			unlink"$file_parsed.$bs.pp"||print"\nUnable to remove temporary results file $file_parsed.$bs.pp.\n$!\n\n";
			}
		}
	
	# calculate average and std-deviation
	print OUT"\n\nAverage (AVG) ping-pong signature and standard deviation (DEV) per [nt] 5' overlap (last 3 columns refer to avg. ppZScore, avg. ppSequeces and avg. ppReads incl. standard deviation).\n";
	foreach$i(1..$max_length-1,'ppZScore','ppSequences','ppReads')
		{
		print OUT"\t$i";
		$avg{$i}=$avg{$i}/$bootstraps;
		foreach$val(@{$all_val{$file_parsed}{$i}})
			{
			$dev{$i}+=($avg{$i}-$val)**2;	
			}
		$dev{$i}=sqrt($dev{$i}/$bootstraps);
		}
	
	print OUT"\nAVG";
	foreach$overlap(1..$max_length-1,'ppZScore','ppSequences','ppReads')
		{
		print OUT"\t$avg{$overlap}";
		}
	print OUT"\nDEV";
	foreach$overlap(1..$max_length-1,'ppZScore','ppSequences','ppReads')
		{
		print OUT"\t$dev{$overlap}";
		}
	print OUT"\n\n\n# # # # # # # # # # # # # # # # # # # # # # # # # # # # # #\n\n";
	undef%avg;
	undef%dev;
	}

# compare number of ping-pong pairs (score for 10bp 5' overlap) and ppZ-scores (10bp 5' overlap over background): Mann-Whitney U test.
if($compare_with_MWU)
	{
	print OUT"Testing for significantly different ping-pong signatures (Mann-Whitney U-test)\n";
	%comparisons=();
	if($all_vs_all==1)
		{
		@files1=@files;
		@files2=@files;
		$file_names_group1="";
		$file_names_group2="";
		}
	elsif($group1_vs_group2==1)
		{
		# collect all values from group 1
		foreach$file_group1(@group1)
			{
			$file_group1_parsed=$file_group1;
			$file_group1_parsed=~s/^.+[\\\/]//;
			foreach$val(@{$all_val{$file_group1_parsed}{10}})
				{
				push(@{$all_val{"group1"}{10}},$val);
				}
			}
		# collect all values from group 2
		foreach$file_group2(@group2)
			{
			$file_group2_parsed=$file_group2;
			$file_group2_parsed=~s/^.+[\\\/]//;
			foreach$val(@{$all_val{$file_group2_parsed}{10}})
				{
				push(@{$all_val{"group2"}{10}},$val);
				}
			}
		$file_names_group1=join(', ',@group1);
		$file_names_group2=join(', ',@group2);
		@files1=("group1 ($file_names_group1)");
		@files2=("group2 ($file_names_group2)");
		$file_names_group1=" (".$file_names_group1.")";
		$file_names_group2=" (".$file_names_group2.")";
		}
	foreach$file1(@files1)
		{
		$file1_parsed=$file1;
		$file1_parsed=~s/^.+[\\\/]//;
		foreach$file2(@files2)
			{
			if($group1_vs_group2==1)
				{
				$file1_parsed="group1";
				$file2_parsed="group2";
				}
			else
				{
				$file2_parsed=$file2;
				$file2_parsed=~s/^.+[\\\/]//;
				next if($file1_parsed eq $file2_parsed);
				next if($comparisons{"$file2_parsed $file1_parsed"});
				$comparisons{"$file1_parsed $file2_parsed"}=1;
				}
			
			# check occurence of a specific value and calculate ranks
			%val_occ=();
			%val_rank=();
			foreach$val(@{$all_val{$file1_parsed}{10}},@{$all_val{$file2_parsed}{10}})
				{
				$val_occ{$val}++;
				}
			$rank=0;
			foreach$val(sort{$a<=>$b}keys%val_occ)
				{
				$rank+=1+(0.5*($val_occ{$val}-1));
				$val_rank{$val}=$rank;
				}
			foreach$val(@{$all_val{$file1_parsed}{10}})
				{
				$sum_ranks_file1+=$val_rank{$val};
				}
			foreach$val(@{$all_val{$file2_parsed}{10}})
				{
				$sum_ranks_file2+=$val_rank{$val};
				}
			$n1=@{$all_val{$file1_parsed}{10}};
			$n2=@{$all_val{$file2_parsed}{10}};
			if($sum_ranks_file1<$sum_ranks_file2)
				{
				($n1,$n2)=($n2,$n1);
				($sum_ranks_file1,$sum_ranks_file2)=($sum_ranks_file2,$sum_ranks_file1);
				}
			
			# calculate Z
			$U=$n1*$n2+(($n1*($n1+1))/2)-$sum_ranks_file1;
			$EU=($n1*$n2)/2;
			$corr=0;
			foreach$val(keys%val_occ)
				{
				$corr+=(($val_occ{$val}**3)-$val_occ{$val})/(($n1+$n2)*($n1+$n2-1))
				}
			$sig_corr=sqrt((($n1*$n2)/12)*($n1+$n2+1-$corr));
			if($sig_corr!=0)
				{
				$Z=abs((int(((($U-$EU)/$sig_corr)*100)+0.5))/100);
				}
			else
				{
				$Z="INF (division by zero)";
				}
			
			if($Z<=3.72)
				{
				$p="p=$Z2p{$Z}";
				}
			else
				{
				$p="p<0.0001";
				}
			print OUT"\n\n\n+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + +\n\n${file1_parsed}$file_names_group1\n   vs.\n${file2_parsed}$file_names_group2\n-> Scores for 10nt 5' overlaps are different with significance level: $p (Z=$Z).";
			print OUT"\n\nSorted scores and ranks for $file1_parsed";
			foreach$val(sort{$b<=>$a}@{$all_val{$file1_parsed}{10}})
				{
				print OUT"\n$val\t$val_rank{$val}";
				}
			print OUT"\n\nSorted scores and ranks for $file2_parsed";
			foreach$val(sort{$b<=>$a}@{$all_val{$file2_parsed}{10}})
				{
				print OUT"\n$val\t$val_rank{$val}";
				}
			}
		}
	}
close OUT;

if(!$silent)
	{
	$computation_time=time-$t0;
	$days=int($computation_time/86400);
	$hours=int(($computation_time-$days*86400)/3600);
	$minutes=int((($computation_time-$days*86400)-($hours*3600))/60);
	$seconds=$computation_time-($days*86400)-($hours*3600)-($minutes*60);
	print"\n\nComputation finished.\nElapsed time: $days days, $hours hours, $minutes minutes, $seconds seconds.\n\n";
	}
exit;

sub ZSCORE
	{
	%Ztable=
		(
		"0.0"=>"0.5000 0.5040 0.5080 0.5120 0.5160 0.5199 0.5239 0.5279 0.5319 0.5359",
		"0.1"=>"0.5398 0.5438 0.5478 0.5517 0.5557 0.5596 0.5636 0.5675 0.5714 0.5753",
		"0.2"=>"0.5793 0.5832 0.5871 0.5910 0.5948 0.5987 0.6026 0.6064 0.6103 0.6141",
		"0.3"=>"0.6179 0.6217 0.6255 0.6293 0.6331 0.6368 0.6406 0.6443 0.6480 0.6517",
		"0.4"=>"0.6554 0.6591 0.6628 0.6664 0.6700 0.6736 0.6772 0.6808 0.6844 0.6879",
		"0.5"=>"0.6915 0.6950 0.6985 0.7019 0.7054 0.7088 0.7123 0.7157 0.7190 0.7224",
		"0.6"=>"0.7257 0.7291 0.7324 0.7357 0.7389 0.7422 0.7454 0.7486 0.7517 0.7549",
		"0.7"=>"0.7580 0.7611 0.7642 0.7673 0.7704 0.7734 0.7764 0.7794 0.7823 0.7852",
		"0.8"=>"0.7881 0.7910 0.7939 0.7967 0.7995 0.8023 0.8051 0.8078 0.8106 0.8133",
		"0.9"=>"0.8159 0.8186 0.8212 0.8238 0.8264 0.8289 0.8315 0.8340 0.8365 0.8389",
		"1.0"=>"0.8413 0.8438 0.8461 0.8485 0.8508 0.8531 0.8554 0.8577 0.8599 0.8621",
		"1.1"=>"0.8643 0.8665 0.8686 0.8708 0.8729 0.8749 0.8770 0.8790 0.8810 0.8830",
		"1.2"=>"0.8849 0.8869 0.8888 0.8907 0.8925 0.8944 0.8962 0.8980 0.8997 0.9015",
		"1.3"=>"0.9032 0.9049 0.9066 0.9082 0.9099 0.9115 0.9131 0.9147 0.9162 0.9177",
		"1.4"=>"0.9192 0.9207 0.9222 0.9236 0.9251 0.9265 0.9279 0.9292 0.9306 0.9319",
		"1.5"=>"0.9332 0.9345 0.9357 0.9370 0.9382 0.9394 0.9406 0.9418 0.9429 0.9441",
		"1.6"=>"0.9452 0.9463 0.9474 0.9484 0.9495 0.9505 0.9515 0.9525 0.9535 0.9545",
		"1.7"=>"0.9554 0.9564 0.9573 0.9582 0.9591 0.9599 0.9608 0.9616 0.9625 0.9633",
		"1.8"=>"0.9641 0.9649 0.9656 0.9664 0.9671 0.9678 0.9686 0.9693 0.9699 0.9706",
		"1.9"=>"0.9713 0.9719 0.9726 0.9732 0.9738 0.9744 0.9750 0.9756 0.9761 0.9767",
		"2.0"=>"0.9772 0.9778 0.9783 0.9788 0.9793 0.9798 0.9803 0.9808 0.9812 0.9817",
		"2.1"=>"0.9821 0.9826 0.9830 0.9834 0.9838 0.9842 0.9846 0.9850 0.9854 0.9857",
		"2.2"=>"0.9861 0.9864 0.9868 0.9871 0.9875 0.9878 0.9881 0.9884 0.9887 0.9890",
		"2.3"=>"0.9893 0.9896 0.9898 0.9901 0.9904 0.9906 0.9909 0.9911 0.9913 0.9916",
		"2.4"=>"0.9918 0.9920 0.9922 0.9925 0.9927 0.9929 0.9931 0.9932 0.9934 0.9936",
		"2.5"=>"0.9938 0.9940 0.9941 0.9943 0.9945 0.9946 0.9948 0.9949 0.9951 0.9952",
		"2.6"=>"0.9953 0.9955 0.9956 0.9957 0.9959 0.9960 0.9961 0.9962 0.9963 0.9964",
		"2.7"=>"0.9965 0.9966 0.9967 0.9968 0.9969 0.9970 0.9971 0.9972 0.9973 0.9974",
		"2.8"=>"0.9974 0.9975 0.9976 0.9977 0.9977 0.9978 0.9979 0.9979 0.9980 0.9981",
		"2.9"=>"0.9981 0.9982 0.9982 0.9983 0.9984 0.9984 0.9985 0.9985 0.9986 0.9986",
		"3.0"=>"0.9987 0.9987 0.9987 0.9988 0.9988 0.9989 0.9989 0.9989 0.9990 0.9990",
		"3.1"=>"0.9990 0.9991 0.9991 0.9991 0.9992 0.9992 0.9992 0.9992 0.9993 0.9993",
		"3.2"=>"0.9993 0.9993 0.9994 0.9994 0.9994 0.9994 0.9994 0.9995 0.9995 0.9995",
		"3.3"=>"0.9995 0.9995 0.9995 0.9996 0.9996 0.9996 0.9996 0.9996 0.9996 0.9997",
		"3.4"=>"0.9997 0.9997 0.9997 0.9997 0.9997 0.9997 0.9997 0.9997 0.9997 0.9998",
		"3.5"=>"0.9998 0.9998 0.9998 0.9998 0.9998 0.9998 0.9998 0.9998 0.9998 0.9998",
		"3.6"=>"0.9998 0.9998 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999 0.9999",
		"3.7"=>"0.9999 0.9999 0.9999",
		);
	foreach$Z(keys%Ztable)
		{
		$digit=-1;
		@p=split(' ',$Ztable{$Z});
		foreach$p(@p)
			{
			$digit++;
			$Z2p{$Z.$digit}=(int(((1-$p)*10000)+0.5))/10000;
			}
		}
	}

sub basic_analysis
	{
	open(BA_IN,$ba_in)||print"\nERROR! Unable to analyze output file $ba_in.\n$!\n\n";
	open(BA_OUT,">$ba_out")||print"\nERROR! Unable to create output file $ba_out.\n$!\n\n";
	%length=();
	%comp=();
	%nt=();
	while(<BA_IN>)
		{
		$_=~s/\s*$//;
		if($_=~s/^>//)
			{
			$reads=$_;
			}
		else
			{
			$length{length$_}+=$reads;
			foreach$p(1..length$_)
				{
				$comp{$p}{substr($_,$p-1,1)}+=$reads;
				$nt{substr($_,$p-1,1)}=1;
				}
			}
		}
	close IN;
	
	@size=keys%length;
	@size=sort{$a<=>$b}@size;
	@nt=keys%nt;
	@nt=sort@nt;
	print BA_OUT"LENGTH DISTRIBUTION\nlength\treads";
	foreach$l($size[0]..$size[-1])
		{
		print BA_OUT"\n$l\t$length{$l}";
		}

	print BA_OUT"\n\nPOSITIONAL NUCLEOTIDE COMPOSITION\npos.";
	foreach$nt(@nt)
		{
		print BA_OUT"\t$nt";
		}
	foreach$p(1..$size[-1])
		{
		print BA_OUT"\n$p";
		foreach$nt(@nt)
			{
			if(!$comp{$p}{$nt})
				{
				$comp{$p}{$nt}=0;
				}
			print BA_OUT"\t$comp{$p}{$nt}";
			}
		}
	undef%length;
	undef%comp;
	undef%nt;
	close BA_IN;
	close BA_OUT;
	}