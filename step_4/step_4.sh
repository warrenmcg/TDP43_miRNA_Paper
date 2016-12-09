#!/bin/bash

CURRENT_DIR=$(pwd)
SCRIPT_DIR='${CURRENT_DIR}/..'
ANNO_DIR='${CURRENT_DIR}/../Annotations'
PREDICTIONS='${ANNO_DIR}/mirna_databases/hg19_predictions_S_targetMatrix.txt'
## This script assumes that the LUAD and LUSC data directories
## are in the same directory as this repository, and called
## "LUAD_data" and "LUSC_data".
## If the data files are in a different folder, you need to change the following
## two lines to the correct directories
LUAD_DIR='${CURRENT_DIR}/../LUAD_data'
LUSC_DIR='${CURRENT_DIR}/../LUSC_data'

LUAD_FILE=${LUAD_DIR}/LUAD_miRNA_allmirna_mirnaCounts.txt
LUSC_FILE=${LUSC_DIR}/LUSC_miRNA_allmirna_mirnaCounts.txt
RNA_DIR=RNASeqV2/UNC__58/Level_3/TCGA_isoform_normalized_results
MIR_DIR=miRNASeq/BCGSC__50/Level_3/compressed_miRNAcounts
SCRIPT=${SCRIPT_DIR}/step_4/runPromise_AllSamples.R

CWD=`pwd`
cd "$LUAD_DIR"
Rscript "$SCRIPT" "file_manifest.txt" "$PREDICTIONS" "${MIR_DIR}" "${RNA_DIR}"

cd "$LUSC_DIR"
Rscript "$SCRIPT" "file_manifest.txt" "$PREDICTIONS" "${MIR_DIR}" "${RNA_DIR}"

cd "$CWD"
