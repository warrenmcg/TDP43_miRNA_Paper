#!/bin/bash

CURRENT_DIR=$(pwd)
SCRIPT_DIR='${CURRENT_DIR}/..'
ANNO_DIR='${CURRENT_DIR}/../Annotations'
UCSC_FILE=${ANNO_DIR}/ucsc_tables/ucsc_kgXref_combinedTable.txt
TDP43_FILE=${ANNO_DIR}/TDP-43-regulated_miRNAList_updated.txt
MIRNA=${ANNO_DIR}/mirna_databases/hs_mirna_mature.txt
## This script assumes that the LUAD and LUSC data directories
## are in the same directory as this repository, and called
## "LUAD_data" and "LUSC_data".
## If the data files are in a different folder, you need to change the following
## two lines to the correct directories
LUAD_DIR='${CURRENT_DIR}/../LUAD_data'
LUSC_DIR='${CURRENT_DIR}/../LUSC_data'

LUAD_INPUT=LUAD_old_fatiscan/your_annotation.txt
LUSC_INPUT=LUSC_old_fatiscan/your_annotation.txt

R_PIPELINE=${SCRIPT_DIR}/step_7/modifiedSPIA_pipeline_a.R
BASH_PIPELINE=${SCRIPT_DIR}/step_7/modifiedSPIA_pipeline_b.sh

BABEL_DIR=post_babelomics_results

# Process the LUAD results
OLD_DIR=`pwd`
cd "$LUAD_DIR"
if [[ ! -d "$BABEL_DIR" ]]; then
	mkdir $BABEL_DIR
fi
OUT_DIR="${LUAD_DIR}/${BABEL_DIR}/processed_results"
PREFIX=LUAD
#Rscript "$R_PIPELINE" -b "${LUAD_DIR}/fatiscan_output/${LUAD_INPUT}" \
Rscript "$R_PIPELINE" -b "${SCRIPT_DIR}/fatiscan_output/${LUAD_INPUT}" \
	-d "${LUAD_DIR}/deseq_results/${PREFIX}_miRNA_miRNA_allResults.txt" \
	-m "$MIRNA" -t "$TDP43_FILE" -x ${BABEL_DIR}/$PREFIX

bash "$BASH_PIPELINE" -i $BABEL_DIR -s "${SCRIPT_DIR}/step_7" -o "$OUT_DIR" \
	-p $PREFIX -d "$UCSC_FILE" -r "${LUAD_DIR}/deseq_results/${PREFIX}_mrna_old_transcriptSignedLogQ_Rank.txt" \
	-t "$TDP43_FILE"
#	-p $PREFIX -d "$UCSC_FILE" -r "${LUAD_DIR}/deseq_results/${PREFIX}_mrna_norm_transcriptThresholdRank.txt" \

cd "$LUSC_DIR"
if [[ ! -d "$BABEL_DIR" ]]; then
	mkdir $BABEL_DIR
fi
OUT_DIR="${LUSC_DIR}/${BABEL_DIR}/processed_results"
PREFIX=LUSC
#Rscript "$R_PIPELINE" -b "${LUSC_DIR}/fatiscan_output/${LUSC_INPUT}" \
Rscript "$R_PIPELINE" -b "${SCRIPT_DIR}/fatiscan_output/${LUSC_INPUT}" \
	-d "${LUSC_DIR}/deseq_results/${PREFIX}_miRNA_miRNA_allResults.txt" \
	-m "$MIRNA" -t "$TDP43_FILE" -x ${BABEL_DIR}/$PREFIX

bash "$BASH_PIPELINE" -i $BABEL_DIR -s "${SCRIPT_DIR}/step_7" -o "$OUT_DIR" \
	-p $PREFIX -d "$UCSC_FILE" -r "${LUSC_DIR}/deseq_results/${PREFIX}_mrna_old_transcriptSignedLogQ_Rank.txt" \
	-t "$TDP43_FILE"
#	-p $PREFIX -d "$UCSC_FILE" -r "${LUSC_DIR}/deseq_results/${PREFIX}_mrna_norm_transcriptThresholdRank.txt" \

cd "$OLD_DIR"
