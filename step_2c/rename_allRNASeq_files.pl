#!/usr/bin/perl
use warnings;
use strict;
use File::Copy qw(copy);
use File::Find;

if (@ARGV != 3)
{
	die "Usage: perl rename_allRNASeq_files.pl <input> input_directory out_directory file_manifest.txt\n";
}

unless( $ARGV[2] =~ "file_manifest.txt" ) { die "Not a file manifest file? $ARGV[2]\n"; }

open(FM, $ARGV[2]) or die "Can't open file $ARGV[2] $!\n";

my $in_dir = $ARGV[0];
my $out_dir = $ARGV[1];
#my $out_dir = "TCGA_isoforms_normalized_results";
mkdir($out_dir);

my ($line, $barcode, $sample, %sample2barcode);

while($line=<FM>)
{
	#Example input line that this script cares for:
	#RNASeqV2	UNC	58	3	TCGA-L4-A4E5-01	TCGA-L4-A4E5-01A-11R-A24X-07	unc.edu.00481f66-1aad-4885-9162-bf802b1ed7f5.2063859.rsem.isoforms.noramlized_results
	unless( $line =~ m"RNASeqV2	UNC	(.*?)	(.*?)	(.*?)	(.*?)	(.*?\.rsem.isoforms.results)") { next; }
	$barcode=$3;
	$sample=$5;

	$sample2barcode{ $sample } = $barcode;
}
close(FM);

find(\&rename_file, $in_dir);

sub rename_file
{
	my $file = $_;
	my $new_file_name;
	if (exists $sample2barcode{$file})
	{
		$barcode = $sample2barcode{$file};
		$new_file_name = "${out_dir}/${barcode}.rsem.isoforms.results";
		copy $file, $new_file_name;
	}
}
