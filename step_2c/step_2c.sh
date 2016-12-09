#!/bin/bash

CURRENT_DIR=$(pwd)
SCRIPT_DIR='${CURRENT_DIR}/..'
ANNO_DIR='${CURRENT_DIR}/../Annotations'
## This script assumes that the LUAD and LUSC data directories
## are in the same directory as this repository, and called
## "LUAD_data" and "LUSC_data".
## If the data files are in a different folder, you need to change the following 
## two lines to the correct directories
LUAD_DIR='${CURRENT_DIR}/../LUAD_data'
LUSC_DIR='${CURRENT_DIR}/../LUSC_data'
RNA_DIR=RNASeqV2/UNC__58/Level_3
MIR_DIR=miRNASeq/BCGSC__50/Level_3

#### miRNA pre-processing ###
## Step 1: compress all of the isoform quantification files to count up all reads for each mature miRNA ###
echo "Pre-processing miRNA data ..."
SCRIPT=${SCRIPT_DIR}/step_2c/compress_allSamples.pl
perl "$SCRIPT" "${LUAD_DIR}/${MIR_DIR}" "${LUAD_DIR}/${MIR_DIR}/compressed_miRNAcounts" "${SCRIPT_DIR}/step_2c" "${LUAD_DIR}/file_manifest.txt"
perl "$SCRIPT" "${LUSC_DIR}/${MIR_DIR}" "${LUSC_DIR}/${MIR_DIR}/compressed_miRNAcounts" "${SCRIPT_DIR}/step_2c" "${LUSC_DIR}/file_manifest.txt"

## Step 2: prepare the miRNA counts for DESeq2 ##
echo "Preparing the miRNA counts for DESeq2"
SCRIPT=${SCRIPT_DIR}/step_2c/prepMiRNACounts4DESeq2.R
MIR_FILE=${ANNO_DIR}/mirna_databases/hg19_predictions_S_miRNAlist.txt
Rscript "$SCRIPT" -m "$MIR_FILE" -c "${LUAD_DIR}/${MIR_DIR}/compressed_miRNAcounts" -t "${LUAD_DIR}/${MIR_DIR}/compressed_miRNAcounts" \
	-o "${LUAD_DIR}/LUAD_miRNA"
Rscript "$SCRIPT" -m "$MIR_FILE" -c "${LUSC_DIR}/${MIR_DIR}/compressed_miRNAcounts" -t "${LUSC_DIR}/${MIR_DIR}/compressed_miRNAcounts" \
	-o "${LUSC_DIR}/LUSC_miRNA"

### RNA-Seq pre-processing ###
## Step 1: convert all RNA-Seq files to have TCGA sample numbers#
echo "Pre-processing the isoform counts for DESeq2"
SCRIPT=${SCRIPT_DIR}/step_2c/rename_allRNASeq_files.pl
perl "$SCRIPT" "${LUAD_DIR}/${RNA_DIR}" "${LUAD_DIR}/${RNA_DIR}/TCGA_isoform_results" "${LUAD_DIR}/file_manifest.txt"
perl "$SCRIPT" "${LUSC_DIR}/${RNA_DIR}" "${LUSC_DIR}/${RNA_DIR}/TCGA_isoform_results" "${LUSC_DIR}/file_manifest.txt"

SCRIPT=${SCRIPT_DIR}/step_2c/rename_allRNASeq_norm_files.pl
perl "$SCRIPT" "${LUAD_DIR}/${RNA_DIR}" "${LUAD_DIR}/${RNA_DIR}/TCGA_isoform_normalized_results" "${LUAD_DIR}/file_manifest.txt"
perl "$SCRIPT" "${LUSC_DIR}/${RNA_DIR}" "${LUSC_DIR}/${RNA_DIR}/TCGA_isoform_normalized_results" "${LUSC_DIR}/file_manifest.txt"
