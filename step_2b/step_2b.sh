#!/bin/bash

CURRENT_DIR=$(pwd)
SCRIPT_DIR='${CURRENT_DIR}/..'
MIR_DIR='${CURRENT_DIR}/../Annotations/mirna_databases'

python "${SCRIPT_DIR}/step_2b/extract_miRANDA_targets.py" "${MIR_DIR}/human_predictions_S_C_aug2010.txt" "${MIR_DIR}/human_predictions_S_0_aug2010.txt" "${MIR_DIR}/hg19_predictions_S"
