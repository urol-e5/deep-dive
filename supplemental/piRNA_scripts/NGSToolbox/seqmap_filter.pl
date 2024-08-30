$input=$ARGV[0]; # input map files must be generated using SeqMap (https://doi.org/10.1093/bioinformatics/btn429)
$max_int=$ARGV[1]; # maximum number of allowed internal mismatches.
$max_tail=$ARGV[2]; # maximum number of allowed non-template 3' nucleotides.
$|=1;

if(@ARGV!=3)
	{
	print
"
USAGE:
perl seqmap_filter.pl mapfile max_int_mm max_tail

where mapfile is a map file in ELAND3 format produced
with SeqMap (with option /output_all_matches), max_int_mm
is an integer giving the maximum number of allowed internal
mismatches and max_tail is an integer giving the maximum
number of allowed non-template 3' nucleotides. The output
file will be named according to your input file name with
extension .filtered (eg. mapfile.filtered).

A standard procedure could be composed of the following two
steps:

seqmap 3 sRNAdata.fas genome.fas mapfile /output_all_matches

perl seqmap_filter.pl mapfile 1 2

";exit;
	};

open(SEQMAP,$input)||die print$!;
open(OUT,">$input.filtered")||die print$!;
print"\nSave best alignments...";
$head=<SEQMAP>;
while(<SEQMAP>)
	{
	@d=split("\t",$_);
	if(!$min_mm{$d[4]}||$min_mm{$d[4]}>$d[5])
		{
		$min_mm{$d[4]}=$d[5];
		}
	}
close SEQMAP;
print" ok.";

print"\nCheck best alignments and output the valid ones...";
open(SEQMAP,$input)||die print$!;
$head=<SEQMAP>;
while(<SEQMAP>)
	{
	@d=split("\t",$_);
	if($d[5]==$min_mm{$d[4]})
		{
		# total mm is low, doesn't matter were mm occurs
		if($d[5]<=$max_int)
			{
			print OUT$_;
			}
		# check if max. internal mismatch is not exceeded
		elsif($d[5]<=$max_int+$max_tail)
			{
			# calculate mm for 3' tail
			$tail_mm=0;
			foreach$i(-$max_tail..-1)
				{
				if(substr($d[4],$i,1) ne substr($d[2],$i,1))
					{
					$tail_mm=-$i;
					last;
					}
				}
			# calculate internal mm
			$int_mm=0;
			foreach$i(1..(length$d[4])-$tail_mm)
				{
				if(substr($d[4],$i-1,1) ne substr($d[2],$i-1,1))
					{
					$int_mm++;
					}
				}
			if($int_mm<=$max_int&&($int_mm+$tail_mm)<=($max_int+$max_tail)&&($int_mm+$tail_mm)==$min_mm{$d[4]})
				{
				print OUT$_;
				}
			}
		}
	}
close SEQMAP;
close OUT;
print" ok.";
unlink$input;
exit;