# FatiGO Input Files

The files in this directory are the input files for FatiGO, produced by **step 7** and used for **step 8**.

It is simply a list of the unique gene symbols that are targets of at least one miRNA that was significant after Fatiscan,
and were differentially expressed as identified by DESeq2.

The lists are split into whether the miRNA is up-regulated or down-regulated. Only genes that are changing in the opposite direction
are considered here. Future work would look at genes that are changing in the same direction as the miRNA (because the miRNA
positively regulates the gene expression, for example through transcriptional activation, or the miRNA is regulating the gene in a
"noise-buffering" capacity, rather than in a "shut-off" capacity).

FatiGO needs to be run separately four times for each list.
