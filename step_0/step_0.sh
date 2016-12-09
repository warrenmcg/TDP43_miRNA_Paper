#!/bin/bash

CURRENT_DIR=$(pwd)
SCRIPT_DIR='${CURRENT_DIR}/..'

Rscript ${SCRIPT_DIR}/step_0/install_packages.R

pip install bcbio-gff
pip install HTSeq
