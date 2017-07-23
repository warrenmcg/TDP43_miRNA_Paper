# FatiGO Output Files

These are all of the output files produced by FatiGO, directly downloaded from Babelomics v4 website.

There is a directory for each dataset (`LUAD`, `LUSC`), and only considering down-regulated or up-regulated genes
(`downGene` and `upGene` respectively).

The relevant files are the `*txt` files for each of the databases listed in `Annotations/pathway_databases.txt`:
GO Biological Process (`go_biological_process_3_9.txt`), Reactome (`reactome.txt`), KEGG (`kegg.txt`), and Biocarta (`biocarta.txt`).

Here is a description of the columns in each of these files (it's the same for all of them):
Here are the meanings of the columns in that file:
+ **term**: the ID of the pathway
+ **term_size**: the number of mRNA members of the pathway as defined by `*.annot` file
+ **term_size_in_genome**: the number of mRNA members actually found in the rank list (this may be higher because some genes were not expressed)
+ **list1_positives**: the number of mRNA members found above the cut-off as defined by Fatiscan's heuristic
+ **list1_negatives**: the number of mRNAs above the cut-off that are not members of the pathway
+ **list1_percentage**: the percentage of mRNAs above the cut-off that are members of the pathway
+ **list2_positives**: the number of mRNA members found below the cut-off as defined by Fatiscan's heuristic
+ **list2_negatives**: the number of mRNAs below the cut-off that are not members of the pathway
+ **list2_percentage**: the percentage of mRNAs below the cut-off that are members of the pathway
+ **list1_positive_ids**: a list of the UCSC IDs for the mRNA members of the pathway above the cut-off
+ **list2_positive_ids**: a list of the UCSC IDs for the mRNA members of the pathway below the cut-off
+ **odds_ratio_log**: the natural log of the odds ratio: (list1 %) / (list2 %)
+ **pvalue**: the p-value as determined by the hypergeometric test
+ **adj_pvalue**: the adjusted p-value after Benjamini-Hochberg correction

## Reference

FatiGO paper: [10.1093/bioinformatics/btg455](https://dx.doi.org/10.1093/bioinformatics/btg455)
FatiGO+ paper: [10.1093/nar/gkm260](https://dx.doi.org/10.1093/nar/gkm260)
Babelomics paper: [10.1093/nar/gkq388](https://dx.doi.org/10.1093/nar/gkq388)
