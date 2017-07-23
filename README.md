# Steps to reproduce this analysis

This is a step-by-step guide on how to reproduce exactly the analysis performed in the paper.

__Note #1__: To run all steps, you can type `bash run_all_steps.sh` at the command line while inside
This directory. This script presumes that you've already run Fatiscan and FatiGO. If you do not have the `fatiscan_output`
and `fatigo_output` directories, the script currently breaks at step 7 or step 9, respectively.

__Note #2__: Each bash script for each step assumes that they are run from their directory. The
scripts will break otherwise.

__Note #3__: The annotations and the TCGA data are not included with this repository because of their
size (several GB total). This pipeline assumes that you do step 1a and 1b.

## Step 0
### install necessary R packages (DESeq2, ProMISe) and python packages (BCBio and HTSeq)
* run `bash step_0.sh` to install the necessary packages (use `sudo bash step_0.sh` to install system-wide)
	- this presumes that R is already installed: to run this particular analysis, you need to install R-3.1.3
	- If you need to install it, follow these instructions:
		* download R-3.1.3 from here: https://cran.r-project.org/src/base/R-3/R-3.1.3.tar.gz
		* run `tar -xzf R-3.1.3.tar.gz` in the directory with the downloaded tarball.
		* navigate into the `R-3.1.3` directory, and then run either of the following set of commands:
			- if you have no other versions of R, run: `./configure` and then `make`, then `make install`
			- if you have other versions of R, run `./configure --prefix=/some/other/path` before running `make` and `make install`
		* be sure to use the right path to Rscript in the following steps

**Note**: you may need to change the top line of the `install_packages.R` script to reflect the proper path to Rscript. If there is an error, you can check the path by typing `which Rscript` at the command line.

**Note**: you *will* need to change the directory paths at the beginning of each bash script currently to run this on your own machine.

## Step 1a
### Download all samples with clinical data, RNA-SeqV2, and miRNA-Seq data available.
1. LUAD: 441 (19 control samples, 422 tumor samples)
	* use `luad_file_manifest.txt` with the GDC data transfer tool to download the same files yourself
2. LUSC: 367 (37 control samples, 330 tumor sample)
	* use `lusc_file_manifest.txt` with the GDC data transfer tool to download the same files yourself

The easiest way to do this is to follow the steps in the `tcga_file_manifests` directory using the
GDC data transfer tool. Click [here](https://github.com/warrenmcg/TDP43_miRNA_Paper/tree/master/tcga_file_manifests#how-to-download-tcga-data) to follow those steps.

If you download the data manually, please keep in mind the following note:

__NOTE__: For the downstream steps to work without modification, you must rename the top level directories for this data
as `LUAD_data` and `LUSC_data`, respectively, and place these directories in this repository.
If you wish to use different names, or place directories elsewhere on your computer, you need to modify the
`LUAD_DIR` and `LUSC_DIR` variables at the top of the following scripts: step_2c.sh, step_3.sh, step_4.sh, step_5.sh, step_7.sh, and step_9.sh

## Step 1b
### Download annotations used for TCGA analysis 
* RNA: UCSC, hg19, June 2011; miRNA: miRBase v21 and miRanda 08/2010 release
* run `step_1b.sh` to do two things:
	- download UCSC tables, the miRNA annotations and initial miRanda miRNA-target predictions
	- generate the table that matches pathway IDs to intelligible pathway names
**NOTE**: this download takes up approximately 880MB of space

## Step 2
### Prep the UCSC annotations to facilitate downstream analyses
* run `step_2a.sh` to prepare the UCSC annotations to facilitate downstream analyses
**Output**: a combined kgXref table to facilitate mapping UCSC IDs to gene names
		
## Step 2b
### Prep the miRNA annotations to facilitate downstream analyses
* run `step_2b.sh`
**Outputs**: 
	- a list of all miRNAs (and accessions) for both the target matrix and miRBase v21
	- a list of all targets
	- a matrix of interactions (1 if yes, 0 if no)

## Step 2c
### Prep the TCGA input files for downstream analyses
* run `step_2c.sh`
* this automatically runs on both LUAD and LUSCS directories
**Outputs**:
	- miRNA count matrix ready for step 3 is in `(LUAD|LUSC)_data` directory
	- miRNA counts ready for 4 are in `(LUAD|LUSC)_data/miRNASeq/BCGSC__50/Level_3/compressed_miRNAcounts` directory
	- RNA normalized isoform counts ready for steps 3 and 4 are in `(LUAD|LUSC)_data/RNASeqV2/UNC__58/Level_3/TCGA_isoform_normalized_results` directory
	
## Step 3
### Run DESeq on the miRNA expression values, as well as the RNA isoform values
* run `step_3.sh`
**Output**: DESeq2 results are in `(LUAD|LUSC)_data/deseq_results` directories.
The *Rank.txt file is the one used in step 6 for the Fatiscan analysis.

## Step 4
### Run ProMISe
* run `step_4.sh`
**Output**: the promise RData will be stored in a directory called `output`

## Step 5
### Identify updated miRNA-mRNA predictions for each dataset
* run `step_5.sh` to process the ProMISe results
**Output**: this produces a GMT gene sets file with transcripts, and an extended gene set file.
The latter is used for step 6, along with the rank file generated by step 3.

## Step 6
### Run Babelomics Fatiscan on DESeq2 results of using promise predictions
*This requires manual manipulation of files*
1) Make an account on Babelomics v4 (http://v4.babelomics.org)
2) Upload the extended miRNA-target gene set file with data type `Annotation > Extended Annotation`
3) Upload the DESeq2 rank list of transcripts with data type `ID List > Ranked`
4) Choose `Functional Analysis`, then `Gene Set Analysis` under `Set Enrichment Analysis`
5) Choose the rank list for input data; choose `fatiscan` for test, and choose `two-tailed` 
6) choose `your own annotations` and select the extended list

Precomputed results from Fatiscan are found in the `fatiscan_output` directory.

## Step 7
### Run post-Babelomics processing on the results
1) first download the results from the Babelomics website. On the bottom of the results page, click `Download Job`
2) move the resulting job folder to the `(LUAD|LUSC)_data` directory and rename it `fatiscan_output`
3) run `step_7.sh`
**Output**: a `post_babelomics_results` directory with processed results for up-regulated and down-regulated genes separately
	- in these directories (specifically, the `downregulatedGenes/up-regulated_miRNAs/` and `upregulatedGenes/down-regulated_miRNAs`), you'll find `*FatigoInput.txt` lists of the genes we'll use as input for FatiGO in step 8.

## Step 8
### Run Babelomics FatiGO on gene-level DESeq2 results
*Like step 6, this also requires manual input of data*
1) Upload the up-regulated and down-regulated list of mRNA targets as `ID List > Gene`
2) Choose `Functional Analysis` then `FatiGO` under `Single Enrichment Analysis` 
3) Choose `Id list vs Rest of Genome`
4) Choose the up- or down-regulated list (run separate analysis for each)
5) `Remove duplicates?` should be set to `Remove from list2 those appearing in list1 (complementary list)`
6) For databases, choose `human` for organism, and check boxes for `GO biological process`, `GO cellular component`, `GO molecular function`, `KEGG pathways`, `Reactome`, and `Biocarta`.
**Note**: For this step, if you do both LUSC and LUAD, you'll run four separate jobs.

Precomputed results from this step are found in the `fatigo_output` directory
## Step 9
### Generate the Hive Plots
1) first download the results from the Babelomics website. Like before, you click on `Download job`
2) in the `(LUAD|LUSC)_data` directory, create a directory called `fatigo_out`
3) move the resulting job folder to the `fatigo_out` and rename it either `upGene` or `downGene`, depending on whether the FatiGO job focused on up-regulated or down-regulated genes.
4) run `step_9.sh`
**Output**: PDF files of the hive plots. 

**Note**: You'll have to manually add in annotations afterward using the tables and Adobe, Preview, or Powerpoint.

**Note 2**: If there are any issues with reproducing this analysis, please submit an issue on Github.
