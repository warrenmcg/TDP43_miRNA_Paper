#!/bin/bash

CURRENT_DIR=$(pwd)
SCRIPT_DIR='${CURRENT_DIR}/..'
## This script assumes that the LUAD and LUSC data directories
## are in the same directory as this repository, and called
## "LUAD_data" and "LUSC_data".
## If the data files are in a different folder, you need to change the following
## two lines to the correct directories
LUAD_DIR='${CURRENT_DIR}/../LUAD_data'
LUSC_DIR='${CURRENT_DIR}/../LUSC_data'

ANNO_DIR='${CURRENT_DIR}/../Annotations'
ALIAS_FILE=${ANNO_DIR}/aliases_list.txt
PATH_FILE=${ANNO_DIR}/pathway_ids2names.txt
DB_FILE=${ANNO_DIR}/pathway_databases.txt
NUM_DB=`wc -l "$DB_FILE" | cut -f 1 -d ' '`

OUT_DIR=hive_plot_tables

UP_SUFFIX=upGenes_downMirnaTDP43SigGeneList.txt
DOWN_SUFFIX=downGenes_upMirnaTDP43SigGeneList.txt
MIRNA_SUFFIX=miRNA_miRNA_allResults.txt
MRNA_SUFFIX=mrna_norm_transcripts_allResults.txt

DESEQ_DIR=deseq_results/
UP_DIR=post_babelomics_results/processed_results/upregulatedGenes/down-regulated_miRNAs/
DOWN_DIR=post_babelomics_results/processed_results/downregulatedGenes/up-regulated_miRNAs/
FATIGO_UP=fatigo_out/upGene
FATIGO_DOWN=fatigo_out/downGene

COMBINE=${SCRIPT_DIR}/step_9/combineInteractionsAndFatiGOprocessHits.py
HIVE_TABLE=${SCRIPT_DIR}/step_9/createHivePlotTable.py
HIVE_DOT=${SCRIPT_DIR}/step_9/convertHivePlotTable2DOT_extension.py
HIVE_R=${SCRIPT_DIR}/step_9/generateHivePlots.R

OLD_DIR=`pwd`

####LUAD####
cd "$LUAD_DIR"
if [[ ! -d "$OUT_DIR" ]]; then
	mkdir $OUT_DIR
	mkdir $OUT_DIR/initial_list
	mkdir $OUT_DIR/intermediate_table
	mkdir $OUT_DIR/hive_plot
	mkdir $OUT_DIR/hive_plot/up
	mkdir $OUT_DIR/hive_plot/down
fi

PREFIX=LUAD
MIRNA_FILE=${LUAD_DIR}/${DESEQ_DIR}/${PREFIX}_${MIRNA_SUFFIX}
MRNA_FILE=${LUAD_DIR}/${DESEQ_DIR}/${PREFIX}_${MRNA_SUFFIX}

####PROCESS upregulated tables####
echo "Processing up-regulated genes in LUAD"
GENE_FILE=${LUAD_DIR}/${UP_DIR}/${PREFIX}_${UP_SUFFIX}
UP_FILES=""
DONE=false
until $DONE; do
	read || DONE=true
	[[ ! $REPLY ]] && continue
	database=$REPLY
	FATIGO_FILE=${LUAD_DIR}/${FATIGO_UP}/${database}.txt
	python "$COMBINE" "$GENE_FILE" "$PATH_FILE" "$ALIAS_FILE" "$FATIGO_FILE" \
		${OUT_DIR}/initial_list/${PREFIX}_up_${database}
	UP_FILES="$UP_FILES ${OUT_DIR}/initial_list/${PREFIX}_up_${database}_List.txt"
done < "$DB_FILE"

echo "Converting FatiGO processed results into Hive Table"
python "$HIVE_TABLE" $NUM_DB "${MIRNA_FILE}" "${MRNA_FILE}" "$ALIAS_FILE" \
	$UP_FILES ${OUT_DIR}/intermediate_table/${PREFIX}_up

echo "converting Hive Table to Dot File"
python "$HIVE_DOT" ${OUT_DIR}/intermediate_table/${PREFIX}_up_nodeList.txt \
	${OUT_DIR}/intermediate_table/${PREFIX}_up_edgeList.txt \
	${OUT_DIR}/hive_plot/up/${PREFIX}_up

####PROCESS downregulated tables####
echo "Processing down-regulated genes in LUAD"
GENE_FILE=${LUAD_DIR}/${DOWN_DIR}/${PREFIX}_${DOWN_SUFFIX}
DOWN_FILES=""
DONE=false
until $DONE; do
	read || DONE=true
	[[ ! $REPLY ]] && continue
	database=$REPLY
	FATIGO_FILE=${LUAD_DIR}/${FATIGO_DOWN}/${database}.txt
	python "$COMBINE" "$GENE_FILE" "$PATH_FILE" "$ALIAS_FILE" "$FATIGO_FILE" \
		${OUT_DIR}/initial_list/${PREFIX}_down_${database}
	DOWN_FILES="$DOWN_FILES ${OUT_DIR}/initial_list/${PREFIX}_down_${database}_List.txt"
done < "$DB_FILE"

echo "Converting FatiGO processed results into Hive Table"
python "$HIVE_TABLE" $NUM_DB "${MIRNA_FILE}" "${MRNA_FILE}" "$ALIAS_FILE" \
	$DOWN_FILES ${OUT_DIR}/intermediate_table/${PREFIX}_down

echo "converting Hive Table to Dot File"
python "$HIVE_DOT" ${OUT_DIR}/intermediate_table/${PREFIX}_down_nodeList.txt \
	${OUT_DIR}/intermediate_table/${PREFIX}_down_edgeList.txt \
	${OUT_DIR}/hive_plot/down/${PREFIX}_down

#### Generate HivePlots ####
echo "Creating the HivePlots"
Rscript "$HIVE_R" -d ${OUT_DIR}/hive_plot/down \
		  -u ${OUT_DIR}/hive_plot/up \
		  -p ${OUT_DIR}/hive_plot/${PREFIX}

##### LUSC #####
echo "Processing up-regulated genes in LUSC"
cd "$LUSC_DIR"
if [[ ! -d "$OUT_DIR" ]]; then
	mkdir $OUT_DIR
	mkdir $OUT_DIR/initial_list
	mkdir $OUT_DIR/intermediate_table
	mkdir $OUT_DIR/hive_plot
	mkdir $OUT_DIR/hive_plot/up
	mkdir $OUT_DIR/hive_plot/down
fi

PREFIX=LUSC
MIRNA_FILE=${LUSC_DIR}/${DESEQ_DIR}/${PREFIX}_${MIRNA_SUFFIX}
MRNA_FILE=${LUSC_DIR}/${DESEQ_DIR}/${PREFIX}_${MRNA_SUFFIX}

####PROCESS upregulated tables####
GENE_FILE=${LUSC_DIR}/${UP_DIR}/${PREFIX}_${UP_SUFFIX}
UP_FILES=""
DONE=false
until $DONE; do
	read || DONE=true
	[[ ! $REPLY ]] && continue
	database=$REPLY
	FATIGO_FILE=${LUSC_DIR}/${FATIGO_UP}/${database}.txt
	python "$COMBINE" "$GENE_FILE" "$PATH_FILE" "$ALIAS_FILE" "$FATIGO_FILE" \
		${OUT_DIR}/initial_list/${PREFIX}_up_${database}
	UP_FILES="$UP_FILES ${OUT_DIR}/initial_list/${PREFIX}_up_${database}_List.txt"
done < "$DB_FILE"

echo "Converting FatiGO processed results into Hive Table"
python "$HIVE_TABLE" $NUM_DB "${MIRNA_FILE}" "${MRNA_FILE}" "$ALIAS_FILE" \
	$UP_FILES ${OUT_DIR}/intermediate_table/${PREFIX}_up

echo "converting Hive Table to Dot File"
python "$HIVE_DOT" ${OUT_DIR}/intermediate_table/${PREFIX}_up_nodeList.txt \
	${OUT_DIR}/intermediate_table/${PREFIX}_up_edgeList.txt \
	${OUT_DIR}/hive_plot/up/${PREFIX}_up

####PROCESS downregulated tables####
echo "Processing down-regulated genes in LUSC"
GENE_FILE=${LUSC_DIR}/${DOWN_DIR}/${PREFIX}_${DOWN_SUFFIX}
DOWN_FILES=""
DONE=false
until $DONE; do
	read || DONE=true
	[[ ! $REPLY ]] && continue
	database=$REPLY
	FATIGO_FILE=${LUSC_DIR}/${FATIGO_DOWN}/${database}.txt
	python "$COMBINE" "$GENE_FILE" "$PATH_FILE" "$ALIAS_FILE" "$FATIGO_FILE" \
		${OUT_DIR}/initial_list/${PREFIX}_down_${database}
	DOWN_FILES="$DOWN_FILES ${OUT_DIR}/initial_list/${PREFIX}_down_${database}_List.txt"
done < "$DB_FILE"

echo "Converting FatiGO processed results into Hive Table"
python "$HIVE_TABLE" $NUM_DB "${MIRNA_FILE}" "${MRNA_FILE}" "$ALIAS_FILE" \
	$DOWN_FILES ${OUT_DIR}/intermediate_table/${PREFIX}_down

echo "converting Hive Table to Dot File"
python "$HIVE_DOT" ${OUT_DIR}/intermediate_table/${PREFIX}_down_nodeList.txt \
	${OUT_DIR}/intermediate_table/${PREFIX}_down_edgeList.txt \
	${OUT_DIR}/hive_plot/down/${PREFIX}_down

#### Generate HivePlots ####	
#echo "Creating the HivePlots"
Rscript "$HIVE_R" -d ${OUT_DIR}/hive_plot/down \
		   -u ${OUT_DIR}/hive_plot/up \
		   -p ${OUT_DIR}/hive_plot/${PREFIX}

echo "Step 9 Complete"
cd "$OLD_DIR"
