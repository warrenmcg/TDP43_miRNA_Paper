#!/bin/bash

DIR=`pwd`

ANNO_DIR='${DIR}/../Annotations'
SCRIPT_DIR='${DIR}/..'

Rscript "${SCRIPT_DIR}/step_1b/generate_pathway_names.R"
mv "${SCRIPT_DIR}/step_1b/pathway_ids2names.txt" "${ANNO_DIR}/"

cd "$ANNO_DIR"
if [[ ! -d "ucsc_tables" ]]; then
	mkdir ucsc_tables
fi
cd ucsc_tables
wget ftp://hgdownload.cse.ucsc.edu/apache/htdocs/goldenPath/hg19/database/kgXref.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/apache/htdocs/goldenPath/hg19/database/kgXrefOld5.txt.gz
wget ftp://hgdownload.cse.ucsc.edu/apache/htdocs/goldenPath/hg19/database/kgXrefOld6.txt.gz
gunzip kgXref*
cd ..

if [[ ! -d "mirna_databases" ]]; then
	mkdir mirna_databases
fi
cd mirna_databases
wget http://cbio.mskcc.org/microrna_data/human_predictions_S_C_aug2010.txt.gz
wget http://cbio.mskcc.org/microrna_data/human_predictions_S_0_aug2010.txt.gz
wget ftp://mirbase.org/pub/mirbase/21/database_files/mirna_mature.txt.gz
gunzip *gz

# make a new file with just human miRNA names, aliases, and miRBase accessions
echo -e "name\talias\taccession" > hs_mirna_mature.txt
cat mirna_mature.txt | cut -f 2,3,4 | grep "^hsa" >> hs_mirna_mature.txt
cd "$DIR"
