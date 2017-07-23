#!/bin/bash

## THIS SCRIPT ASSUMES YOU HAVE GDC Data Transfer Tool Downloaded
## and it is on your $PATH

CURRENT_DIR=$(pwd)
MANIFEST_DIR='${CURRENT_DIR}/../tcga_file_manifests'
LUAD_manifest='${MANIFEST_DIR}/LUAD_new_file_manifest.txt'
LUSC_manifest='${MANIFEST_DIR}/LUSC_new_file_manifest.txt'

GDC_CLIENT_PATH=gdc-client

if [[ ! -d "${CURRENT_DIR}/../LUAD_data" ]]; then
        mkdir '${CURRENT_DIR}/../LUAD_data'
fi

if [[ ! -d "${CURRENT_DIR}/../LUSC_data" ]]; then
        mkdir '${CURRENT_DIR}/../LUSC_data'
fi

'$GDC_CLIENT_PATH' download -m $LUAD_MANIFEST -d '${CURRENT_DIR}/../LUAD_data'
'$GDC_CLIENT_PATH' download -m $LUSC_MANIFEST -d '${CURRENT_DIR}/../LUSC_data'
