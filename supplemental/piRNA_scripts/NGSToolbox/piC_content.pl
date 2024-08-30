#!/usr/bin/perl
use strict;
use warnings;

# Constants
$|=1; #Autoflush
# Variables

# Program name
print("\n--- $0 ---\n");
# Collect command line arguments
my $USAGE = "perl $0 <genome_locs.txt> <repeatmasker.out> <gene_set.gtf>\n";
unless ($ARGV[0] && $ARGV[1] && $ARGV[2]) {
	die("\nUsage: $USAGE\n");
}
# Input files
my $genome_loc_file = $ARGV[0];
my $repeatmask_file = $ARGV[1];
my $genome_gtf_file = $ARGV[2];

# Get pic loci
my $pic_data = get_pic_data($genome_loc_file);
# Get repeat loci
my $rep_data = get_repeatmask_data($repeatmask_file);
# Get gff gene and cds locations
my($gen_locs,$cds_locs) = get_gff_gene_locs($genome_gtf_file);

# Get repeat share for each pic
my $rep_share_per_pic = get_pic_repeat_shares($pic_data,$rep_data);
# Get repeat bps for each pic
my $rep_bps_per_pic = get_pic_repeat_bps($pic_data,$rep_data);
# Get gene share for each pic
my $gen_share_per_pic = get_pic_gene_shares($pic_data,$gen_locs,$cds_locs);
# Get gene bps for each pic
my $gen_bps_per_pic = get_pic_gene_bps($pic_data,$gen_locs,$cds_locs);
# Get repeat divergence per pic
my $rep_div_per_pic = get_pic_repeat_divergence($pic_data,$rep_data);

# Output results
my $out = open_outfile('pic_content.txt');
print($out "row_id\tpic_chr\tpic_beg\tpic_end\tpic_len\tTEshare\tCDSshare\tTEbp_per_pic\tCDSbp_per_pic\tTEdiv_per_pic\n");
foreach my $row_id (sort keys %{$pic_data}) {
	if ($row_id eq '#') { print($out "\n") and next }
	my $pic_chr = $pic_data->{$row_id}->[1];
	my $pic_beg = $pic_data->{$row_id}->[2];
	my $pic_end = $pic_data->{$row_id}->[3];
	my $pic_len = $pic_end-$pic_beg+1;
	print($out "row_$row_id\t$pic_chr\t$pic_beg\t$pic_end\t$pic_len\t");
	print($out "$rep_share_per_pic->{$row_id}\t$gen_share_per_pic->{$row_id}\t");
	print($out "$rep_bps_per_pic->{$row_id}\t$gen_bps_per_pic->{$row_id}\t");
	print($out "$rep_div_per_pic->{$row_id}\n");
}

exit;

################################# subroutines #################################

sub get_pic_data {
	# Take pic expression file name
	my($expression_file) = @_;
	# Variables
	my %pic_data = ();
	my $row_id = 0;
	# Get file data
	my @expression_data = get_file_data_array($expression_file);
	# Parse pic expression file
	foreach my $line (@expression_data) {
		# Line starts with pic id
		if ($line =~ /^\d+/) {
			$row_id++;
			# Get line data
			my @d = split(/\t/,$line);
			@{$pic_data{$row_id}} = @d;
		}
	}
	return \%pic_data;
}

sub get_repeatmask_data {
	# Take repeatmasker file name
	my($repeatmask_file) = @_;
	# Storage variable
	my %rep_data = ();
	# Get file data
	my @repeatmask_data = get_file_data_array($repeatmask_file);
	# Parse repeatmasker file
	foreach my $line (@repeatmask_data) {
		# Line starts with sw score
		if ($line =~ /^\s*\d+/) {
			# Remove leading whitespace
			$line =~ s/^\s*//;
			# Get line data
			my @d = split(/\s+/,$line);
			my $chr = $d[4];
			push(@{$rep_data{$chr}},\@d);
		}
	}
	return \%rep_data;
}

sub get_gff_gene_locs {
	# Take infile name and keyword
	my($gff_file) = @_;
	# Get file data
	my $gff = open_infile($gff_file);
	# Initialize variables
	my %gff_cds_locs = ();
	my %gff_gen_locs = ();
	my $gene = '';
	my $tran = '';
	# Parse gff file
	while (my $line = <$gff>) {
		# Results line
		if ($line =~ /^[^#]\w+/) {
			# Get info
			my @d = split(/\t/, $line);
			my $chr  = $d[0];
			my $type = $d[2];
			my $beg  = $d[3];
			my $end  = $d[4];
			my $info = $d[8];
			# Gene line
			if ($type =~ /gene/i) {
				# Get gene id
				if ($info =~ /gene_id\s\"([^\"]+)\";/) {
					($gene) = ($info =~ /gene_id\s\"([^\"]+)\";/);
				} elsif ($info =~ /^ID=maker/) {
					($gene) = ($info =~ /ID=([^;:]+)/);
				}
				# Save entries
				$gff_gen_locs{$chr}{$gene} = \@d;
			}
			# CDS line
			elsif ($type =~ /cds/i) {
				# Get transcript id
				if ($info =~ /transcript_id\s\"([^\"]+)\";/) {
					($tran) = ($info =~ /transcript_id\s\"([^\"]+)\";/);
				} elsif ($info =~ /^ID=maker/) {
					($tran) = ($info =~ /ID=([^;:]+)/);
				}
				# Save entries
				push(@{$gff_cds_locs{$gene}{$tran}},\@d);
			}
		}
	}
	return \(%gff_gen_locs,%gff_cds_locs);
}

sub get_pic_repeat_shares {
	# Take pic and repeat data
	my($pic_data,$rep_data) = @_;
	# Storage variable
	my %rep_share_per_pic = ();
	# Go through each pic
	foreach my $row_id (sort keys %{$pic_data}) {
		# Get pic coordinates
		my $pic_chr = $pic_data->{$row_id}->[1];
		my $pic_beg = $pic_data->{$row_id}->[2];
		my $pic_end = $pic_data->{$row_id}->[3];
		my $pic_len = $pic_end-$pic_beg+1;
		# Variable for non-redundant repeat positions
		my %rep_positions = ();
		# Get pic repeats
		foreach my $rep (@{$rep_data->{$pic_chr}}) {
			# Get repeat coordinates
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			my $rep_name = $rep->[9];
			my $rep_class = $rep->[10];
			# Check if repeat lies in pic
			if ($rep_end >= $pic_beg && $rep_beg <= $pic_end) {
				# Count repeat bps
				foreach my $pos ($rep_beg..$rep_end) {
					if ($pos >= $pic_beg && $pos <= $pic_end) {
						$rep_positions{$pos} = 1;
					}
				}
			}
		}
		# Get sum rep shares per pic
		$rep_share_per_pic{$row_id} = 0;
		foreach my $pos (keys %rep_positions) {
			$rep_share_per_pic{$row_id} += 1/$pic_len*100;
		}
	}
	return \%rep_share_per_pic;
}

sub get_pic_repeat_bps {
	# Take pic and repeat data
	my($pic_data,$rep_data) = @_;
	# Storage variable
	my %rep_bps_per_pic = ();
	# Go through each pic
	foreach my $row_id (sort keys %{$pic_data}) {
		# Get pic coordinates
		my $pic_chr = $pic_data->{$row_id}->[1];
		my $pic_beg = $pic_data->{$row_id}->[2];
		my $pic_end = $pic_data->{$row_id}->[3];
		my $pic_len = $pic_end-$pic_beg+1;
		# Variable for non-redundant repeat positions
		my %rep_positions = ();
		# Get pic repeats
		foreach my $rep (@{$rep_data->{$pic_chr}}) {
			# Get repeat coordinates
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			my $rep_name = $rep->[9];
			my $rep_class = $rep->[10];
			# Check if repeat lies in pic
			if ($rep_end >= $pic_beg && $rep_beg <= $pic_end) {
				# Count repeat bps
				foreach my $pos ($rep_beg..$rep_end) {
					if ($pos >= $pic_beg && $pos <= $pic_end) {
						$rep_positions{$pos} = 1;
					}
				}
			}
		}
		# Get sum rep shares per pic
		$rep_bps_per_pic{$row_id} = 0;
		foreach my $pos (keys %rep_positions) {
			$rep_bps_per_pic{$row_id} += 1;
		}
	}
	return \%rep_bps_per_pic;
}

sub get_pic_repeat_divergence {
	# Take pic and repeat data
	my($pic_data,$rep_data) = @_;
	# Storage variable
	my %rep_div_per_pic = ();
	# Go through each pic
	foreach my $row_id (sort keys %{$pic_data}) {
		# Get pic coordinates
		my $pic_chr = $pic_data->{$row_id}->[1];
		my $pic_beg = $pic_data->{$row_id}->[2];
		my $pic_end = $pic_data->{$row_id}->[3];
		my $pic_len = $pic_end-$pic_beg+1;
		# Variable for non-redundant repeat positions
		my @div_sum = ();
		# Get pic repeats
		foreach my $rep (@{$rep_data->{$pic_chr}}) {
			# Get repeat coordinates
			my $rep_div = $rep->[1];
			my $rep_beg = $rep->[5];
			my $rep_end = $rep->[6];
			my $rep_name = $rep->[9];
			my $rep_class = $rep->[10];
			# Check if repeat lies in pic
			if ($rep_end >= $pic_beg && $rep_beg <= $pic_end) {
				# Count repeat bps
				push(@div_sum,$rep_div);
			}
		}
		# Get sum rep shares per pic
		$rep_div_per_pic{$row_id} = 0;
		foreach my $div (@div_sum) {
			$rep_div_per_pic{$row_id} += $div/scalar(@div_sum);
		}
	}
	return \%rep_div_per_pic;
}

sub get_pic_gene_shares {
	# Take pic and repeat data
	my($pic_data,$gen_locs,$cds_locs) = @_;
	# Storage variable
	my %gen_share_per_pic = ();
	# Go through each pic
	foreach my $row_id (sort keys %{$pic_data}) {
		# Get pic coordinates
		my $pic_chr = $pic_data->{$row_id}->[1];
		my $pic_beg = $pic_data->{$row_id}->[2];
		my $pic_end = $pic_data->{$row_id}->[3];
		my $pic_len = $pic_end-$pic_beg+1;
		# Variable for non-redundant cds positions
		my %cds_positions = ();
		# Get pic genes
		foreach my $gene (sort {$gen_locs->{$pic_chr}->{$a}->[3] <=> $gen_locs->{$pic_chr}->{$b}->[3]} keys %{$gen_locs->{$pic_chr}}) {
			# Skip if gene lies outside of pic
			my $gen_beg = $gen_locs->{$pic_chr}->{$gene}->[3];
			my $gen_end = $gen_locs->{$pic_chr}->{$gene}->[4];
			unless ($gen_end >= $pic_beg && $gen_beg <= $pic_end) {next}
			# Go through each transcript
			foreach my $tran (sort keys %{$cds_locs->{$gene}}) {
				# Go through each exon
				foreach my $exon (@{$cds_locs->{$gene}->{$tran}}) {
					# Get exon coordinates
					my $exon_beg = $exon->[3];
					my $exon_end = $exon->[4];
					# Count repeat bps
					foreach my $pos ($exon_beg..$exon_end) {
						if ($pos >= $pic_beg && $pos <= $pic_end) {
							$cds_positions{$pos} = 1;
						}
					}
				}
			}
		}
		# Get sum rep shares per pic
		my $cds_count = 0;
		$gen_share_per_pic{$row_id} = 0;
		foreach my $pos (keys %cds_positions) {
			$gen_share_per_pic{$row_id} += 1/$pic_len*100;
			$cds_count++;
		}
	}
	return \%gen_share_per_pic;
}

sub get_pic_gene_bps {
	# Take pic and repeat data
	my($pic_data,$gen_locs,$cds_locs) = @_;
	# Storage variable
	my %gen_bps_per_pic = ();
	# Go through each pic
	foreach my $row_id (sort keys %{$pic_data}) {
		# Get pic coordinates
		my $pic_chr = $pic_data->{$row_id}->[1];
		my $pic_beg = $pic_data->{$row_id}->[2];
		my $pic_end = $pic_data->{$row_id}->[3];
		my $pic_len = $pic_end-$pic_beg+1;
		# Variable for non-redundant cds positions
		my %cds_positions = ();
		# Get pic genes
		foreach my $gene (sort {$gen_locs->{$pic_chr}->{$a}->[3] <=> $gen_locs->{$pic_chr}->{$b}->[3]} keys %{$gen_locs->{$pic_chr}}) {
			# Skip if gene lies outside of pic
			my $gen_beg = $gen_locs->{$pic_chr}->{$gene}->[3];
			my $gen_end = $gen_locs->{$pic_chr}->{$gene}->[4];
			unless ($gen_end >= $pic_beg && $gen_beg <= $pic_end) {next}
			# Go through each transcript
			foreach my $tran (sort keys %{$cds_locs->{$gene}}) {
				# Go through each exon
				foreach my $exon (@{$cds_locs->{$gene}->{$tran}}) {
					# Get exon coordinates
					my $exon_beg = $exon->[3];
					my $exon_end = $exon->[4];
					# Count repeat bps
					foreach my $pos ($exon_beg..$exon_end) {
						if ($pos >= $pic_beg && $pos <= $pic_end) {
							$cds_positions{$pos} = 1;
						}
					}
				}
			}
		}
		# Get sum rep shares per pic
		$gen_bps_per_pic{$row_id} = 0;
		foreach my $pos (keys %cds_positions) {
			$gen_bps_per_pic{$row_id} += 1;
		}
	}
	return \%gen_bps_per_pic;
}

sub open_infile {
	# Take input file name
    my($file) = @_;
    # Open input file
    my $fh;
    if ($file =~ /.gz$/) {
		open($fh, "gunzip -c $file |") or die("Cannot open file '$file': $!\n");
	} else {
    	open($fh, '<', $file) or die("Cannot open file '$file': $!\n");
    }
    # Return filehandle
    return $fh;
}

sub open_outfile {
	# Take output file name
    my($file) = @_;
    # Open output file
    open(my $fh, '>', $file) or die("Cannot open file '$file': $!\n");
    # Return filehandle
    return $fh;
}

sub get_file_data_array {
	# Take input file name
    my($file,$ref_opt) = @_;
    my @filedata = ();
    $ref_opt = 0 unless $ref_opt;
	# Open input file
    my $fh = open_infile($file);
	# Extract lines and save in array
    while (my $line = <$fh>) {
    	$line =~ s/\s+$//; #better chomp
    	push(@filedata, $line);
    }
	# Close file
    close($fh) or die("Unable to close: $!\n");
	# Return array containing file data
    if ($ref_opt) {
    	return \@filedata;
    } else {
    	return @filedata;
    }
}