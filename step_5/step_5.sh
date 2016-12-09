#!/bin/bash

CURRENT_DIR=$(pwd)
SCRIPT_DIR='${CURRENT_DIR}/..'
SCRIPT1=${SCRIPT_DIR}/step_5/extract_promise_miRNA_UnionGeneSets_4.0.R
SCRIPT2=${SCRIPT_DIR}/step_5/extract_extendedTranscriptListFromGeneSets.py
ANNO_DIR='${CURRENT_DIR}/../Annotations'
ANNOS='${ANNO_DIR}/mirna_databases/hg19_predictions_S_miRNAlist.txt'
## This script assumes that the LUAD and LUSC data directories
## are in the same directory as this repository, and called
## "LUAD_data" and "LUSC_data".
## If the data files are in a different folder, you need to change the following
## two lines to the correct directories
LUAD_DIR='${CURRENT_DIR}/../LUAD_data'
LUSC_DIR='${CURRENT_DIR}/../LUSC_data'
OUT_DIR=output
PREFIX=LUAD
SUFFIX=all_0_UnionGeneSetsWithCounts.gmt

CWD=`pwd`
cd "${LUAD_DIR}"
Rscript "$SCRIPT1" -d "${LUAD_DIR}/${OUT_DIR}" -x "${LUAD_DIR}/${PREFIX}" -a "$ANNOS"
python "$SCRIPT2" "${LUAD_DIR}/${PREFIX}_${SUFFIX}" "${LUAD_DIR}/${PREFIX}"

cd "${LUSC_DIR}"
PREFIX=LUSC
Rscript "$SCRIPT1" -d "${LUSC_DIR}/${OUT_DIR}" -x "${LUSC_DIR}/${PREFIX}" -a "$ANNOS"
python "$SCRIPT2" "${LUSC_DIR}/${PREFIX}_${SUFFIX}" "${LUSC_DIR}/${PREFIX}"

cd "$CWD"
