#!/usr/bin/perl
use strict;
use warnings;

if (@ARGV != 4)
{
	die "Usage: perl compress_allSamples.pl <input> input_directory out_directory script_directory file_manifest.txt \n";
}

unless( $ARGV[3] =~ "file_manifest.txt") { die "Not a file manifest file: $ARGV[3]\n"; }

open(FM, $ARGV[3]) or die "Can't open file $ARGV[3] $!\n";

my $in_dir=$ARGV[0];
my $out_dir=$ARGV[1];
my $script_dir=$ARGV[2];

mkdir($out_dir);

my ($line, $barcode, $sample, $input, $output, $resp);

while($line=<FM>)
{
	#example input line
	#Platform type	center	platform	level	sample	barcode	file name
	#miRNASeq	BCGSC	50	3	TCGA-18-5592-01	TCGA-18-5592-01A-01T-1634-13	TCGA-18-5592-01A-01T-1634-13.isoform.quantification.txt
	unless( $line =~ m"miRNASeq	BCGSC	(.*?)	(.*?)	(.*?)	(.*?)	(.*?\.isoform\.quantification\.txt)") { next; }	
	$barcode=$3;
	$sample=$5;
	$input="${in_dir}/${sample}";
	$output="${out_dir}/${barcode}_miRNAcounts.txt";
	$resp = `perl "${script_dir}/compress_1sample_miRNAcounts.pl"  "$input" "$output"`;
}
close(FM)
