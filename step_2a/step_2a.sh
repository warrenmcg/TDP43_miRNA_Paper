#!/bin/bash

# This assumes that you have completed step 0 and step 1

### You should change the Variables below to your own
CURRENT_DIR=$(pwd)
SCRIPT_DIR='${CURRENT_DIR}/..'
UCSC_DIR='${CURRENT_DIR}/Annotations/ucsc_tables'
###

python "${SCRIPT_DIR}/step_2a/combine_ucscTables.py" "${UCSC_DIR}/kgXref.txt" "${UCSC_DIR}/kgXrefOld6.txt" "${UCSC_DIR}/kgXrefOld5.txt" "${UCSC_DIR}/ucsc_kgXref_combinedTable.txt"
