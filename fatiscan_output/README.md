# Fatiscan Output Files

The files in this directory were downloaded directly from the Babelomics website.

One job was run for each dataset (`LUAD`, `LUSC`), found in their respective directories.

In order for **step 8** to work correctly, these need to be moved to the `LUAD_data` and `LUSC_data`
directories.

The main file used for downstream steps is `your_annotation.txt` which contains the results from
`Fatiscan` for all of the miRNAs.

Here are the meanings of the columns in that file:
+ **term**: the mirBase (v21) accession number of the miRNA
+ **term_size**: the number of mRNA targets for the miRNA identified by the ProMISe part of the pipeline
+ **term_size_in_genome**: the number of mRNA targets actually found in the rank list (it's the same in this case)
+ **list1_positives**: the number of mRNA targets found above the cut-off as defined by Fatiscan's heuristic
+ **list1_negatives**: the number of mRNAs above the cut-off that are not targeted by the miRNA
+ **list1_percentage**: the percentage of mRNAs above the cut-off that are targets of the miRNA
+ **list2_positives**: the number of mRNA targets found below the cut-off as defined by Fatiscan's heuristic
+ **list2_negatives**: the number of mRNAs below the cut-off that are not targeted by the miRNA
+ **list2_percentage**: the percentage of mRNAs below the cut-off that are targets of the miRNA
+ **list1_positive_ids**: a list of the UCSC IDs for the mRNA targets of the miRNA above the cut-off
+ **list2_positive_ids**: a list of the UCSC IDs for the mRNA targets of the miRNA below the cut-off
+ **odds_ratio_log**: the natural log of the odds ratio: (list1 %) / (list2 %)
+ **pvalue**: the p-value as determined by the hypergeometric test
+ **adj_pvalue**: the adjusted p-value after Benjamini-Hochberg correction

## References

Fatiscan paper: [10.1186/1471-2105-8-114](https://dx.doi.org/10.1186/1471-2105-8-114)
Babelomics paper: [10.1093/nar/gkq388](https://dx.doi.org/10.1093/nar/gkq388)
