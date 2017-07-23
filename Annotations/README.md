# Annotation Files

Here is a description of the files in this directory.

## TDP-43-regulated_miRNAlist.txt

This files contains a list of mature human miRNAs that are putatively regulated by TDP-43.
There are 83+2 miRNAs listed here. The miRNAs on the list met at least one of two criteria
for inclusion. These two criteria are the following:
1. Experimental evidence of TDP-43 interacting with the miRNA (see Figures 1 and 4 from the paper).
2. At least one of the following:
+ Expression was affected in at least two cell lines (see Supplemental Table 2 from the paper).
+ Arm ratio or isomiR expression patterns were affected in at least two cell line (see Supplemental Tables 3-5 from the paper).

**Note**: this list was manually curated in the Spring of 2014 by the first author of the paper, Xiaowei Chen.
There may have been miRNAs that should have been included in this list but were missed. These missing miRNAs
do not change the conclusions in the paper: there are multiple miRNAs putatively regulated by TDP-43 that are
either already associated with cancer, or are predicted to be associated with lung cancer based on the pipeline.

## aliases_list.txt

This lists aliases for gene targets of TDP-43-regulated miRNAs. The goal of this file is to convert
older or non-standard names for these genes to the standard HGNC, modern gene symbols. This is required for step 9
to produce the tables for the Hive Plots.

## pathway_databases.txt

This lists the different databases used in FatiGO. This is used to color the nodes and edges of the pathway 
aspect of the Hive Plot. This is used in step 9.


## pathways_id2names.txt

This contains all of the pathways for the different pathway databases (Gene Ontology, KEGG, Biocarta, and Reactome).
Here are the columns:
+ **category**: which database is this pathway are part of?
+ **id**: what is the ID of this pathway?
+ **name**: what is the human readable name for the pathway?

This is used in step 9 to generate the edge and node tables for the Hive Plots.
