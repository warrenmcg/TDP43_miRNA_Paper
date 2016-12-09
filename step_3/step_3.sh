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

RNA_DIR=RNASeqV2/UNC__58/Level_3/TCGA_isoform_normalized_results
#RNA_DIR=RNASeqV2/UNC__58/Level_3/TCGA_isoform_results
LUAD_FILE=${LUAD_DIR}/LUAD_miRNA_allmirna_mirnaCounts.txt
LUSC_FILE=${LUSC_DIR}/LUSC_miRNA_allmirna_mirnaCounts.txt

if [[ ! -d "${LUAD_DIR}/deseq_results" ]]; then
	echo "creating a DESeq2 results folder for LUAD data"
	mkdir "${LUAD_DIR}/deseq_results"
fi
if [[ ! -d "${LUSC_DIR}/deseq_results" ]]; then
	echo "creating a DESeq2 results folder for LUSC data"
	mkdir "${LUSC_DIR}/deseq_results"
fi

### STEP 1: run DESeq2 on the miRNA counts ###
#echo "Running DESeq2 on the miRNA counts"
#SCRIPT=${SCRIPT_DIR}/step_3/DESeq_miRNACounts_TCGA_version.R
#Rscript "$SCRIPT" -m "$LUAD_FILE" -o "${LUAD_DIR}/deseq_results/LUAD_miRNA"
#Rscript "$SCRIPT" -m "$LUSC_FILE" -o "${LUSC_DIR}/deseq_results/LUSC_miRNA"

### STEP 2: run DESeq2 on the RNA isoform counts (which includes preprocessing) ###
echo "Running DESeq2 on the mRNA isoform counts"
SCRIPT=${SCRIPT_DIR}/step_3/prep_isoformCounts4DESeq2.R
~/warren/Software/R-3.1.3/bin/Rscript "$SCRIPT" -c "${LUAD_DIR}/${RNA_DIR}" -t "${LUAD_DIR}/${RNA_DIR}" -o "${LUAD_DIR}/deseq_results/LUAD_mrna_old"
~/warren/Software/R-3.1.3/bin/Rscript "$SCRIPT" -c "${LUSC_DIR}/${RNA_DIR}" -t "${LUSC_DIR}/${RNA_DIR}" -o "${LUSC_DIR}/deseq_results/LUSC_mrna_old"
echo "Step 3 complete ..."
