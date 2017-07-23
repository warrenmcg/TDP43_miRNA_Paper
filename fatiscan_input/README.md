# Fatiscan Input Files

These files were produced by running steps 3 and 6. Here is a description of each file:

## (LUAD|LUSC)_mrna_old_transcriptSignedLogQ_Rank.txt

This file is produced by step 3. It ranks each transcript using the "Adjusted Rank" metric described in the paper.
Specifically:
+ if the transcript had a base mean of 30 or less, its rank = log10(base mean) + |log2 fold change|
+ otherwise, rank = log10(base mean) + |log2 fold change| + log10(adjusted p-value)

The rank was then given the same sign as the direction of fold change in lung cancer versus control.

This metric was chosen to prevent transcripts with low abundanced expression from dominating the ranking by fold change,
and to otherwise emphasize the transcripts with the most expression, greatest fold change, and greatest statistical significance.

## (LUAD|LUSC)_inverseGMTmirbaseAccExtended.txt

This file lists each miRNA-mRNA interaction predicted by ProMISe (i.e. the output ProMISe probability was greater than 0) in all
samples in the dataset.
