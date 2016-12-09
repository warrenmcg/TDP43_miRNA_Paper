#!/usr/bin/perl
use warnings;
use strict;

if (@ARGV != 2)
{
	die "Usage: perl compress_1sample_miRNAcounts.pl <input> TCGA_quantification_file <output> output_file_name\n";
}

my $QNT=$ARGV[0];
unless( $QNT =~ m"(.*?)\.isoform\.quantification\.txt" ) { die "Not a TCGA quantification file? $QNT\n"; }
my $f_name=$1;

open(IN, "$QNT") or die "Can't open file $QNT $!\n";
open(OUT, ">$ARGV[1]") or die "Can't open file $ARGV[1] $!\n";


my $counter = 0; #counter for sample/miRNA array

my ($line, $miRNA, $genom, $rawReads, $rpm, $mirBase, $index, %miRNAcounts, @reads);

while($line=<IN>)
{
	#Example input line:
	#hsa-let-7a-1	hg19:9:96938242-96938265:+	1	0.994335	N	mature,MIMAT0000062
	unless( $line =~ m"(.*?)	(.*?)	(.*?)	(.*?)	(.*?)	(.*?),(.*)" ) { next; }
	$miRNA=$1; #miRNA name
	#$genom=$2; #genomic coordinates
	$rawReads=$3; #raw reads
	$rpm=$4; #reads per million
	$mirBase=$7; #miRBase accession number
	
	unless( exists $miRNAcounts{$mirBase} ) 
	{ 
		$miRNAcounts{ $mirBase } = $counter; 
		$counter++;
	}
	
	$index=$miRNAcounts{$mirBase};
    
	if(defined $reads[$index])
	{
		$reads[$index][1] += $rawReads;
		$reads[$index][2] += $rpm;
	}
	else { $reads[$index] = [$miRNA, $rawReads, $rpm]; }
}
close(IN);

print OUT "mirBase Accession	miRNA name	raw read count	RPM\n";	
foreach my $i (sort keys %miRNAcounts)
{
	$index=$miRNAcounts{$i};
	print OUT "$i	$reads[$index][0]	$reads[$index][1]	$reads[$index][2]\n";
}
close(OUT);


