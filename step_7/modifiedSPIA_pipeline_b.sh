#!/bin/bash

set -e

usage () {
echo "

Usage: star_pipeline.sh [-OPTION VALUE]

Options:

    -h  Display this help
    -i  full path to directory with the intermediate post-babelomic files
    -s  full path to scripts
    -o  output directory for the post-babelomic files
    -p  prefix for the output files
    -d  full path to UCSC database file to annotate genes
    -r  full path to the DESeq2 Rank file
    -t  full path to the list of TDP-43 regulated miRNAs (expect mature names, e.g. hsa-miR-1-3p)
"
   exit 0
}

# if no arguments are provided, display help
[[ ! $1 ]] && { usage; }

# initialize variables
IN_DIR=""
SCRIPT_DIR=""
DATABASE=""
OUT_DIR=""
PREFIX=""
RANK_FILE=""
TDP43_FILE=""

while getopts ":i:d:o:p:r:s:t:h" opt; do
   case $opt in
   h ) usage ; exit;;
   i ) IN_DIR="$OPTARG" ;;
   s ) SCRIPT_DIR="$OPTARG" ;;
   d ) DATABASE="$OPTARG" ;;
   o ) OUT_DIR="$OPTARG" ;;
   p ) PREFIX="$OPTARG" ;;
   r ) RANK_FILE="$OPTARG" ;;
   t ) TDP43_FILE="$OPTARG" ;;
   \?) echo "Unknown parameter: $OPTARG"
       usage ;;
   esac
done

REMOVE="${SCRIPT_DIR}/removeBabelomicHitList.py"
EXTRACT="${SCRIPT_DIR}/extract_BabelomicSigGeneLists_TCGAversion.py"
TDP43="${SCRIPT_DIR}/extract_TDP43miRNAsFromBabelomicsGeneList.py"

if [[ ! -d "${OUT_DIR}" ]]; then
    mkdir "${OUT_DIR}"
    mkdir "${OUT_DIR}/downregulatedGenes"
    mkdir "${OUT_DIR}/downregulatedGenes/all_miRNAs"
    mkdir "${OUT_DIR}/downregulatedGenes/down-regulated_miRNAs"
    mkdir "${OUT_DIR}/downregulatedGenes/up-regulated_miRNAs"

    mkdir "${OUT_DIR}/upregulatedGenes"
    mkdir "${OUT_DIR}/upregulatedGenes/all_miRNAs"
    mkdir "${OUT_DIR}/upregulatedGenes/down-regulated_miRNAs"
    mkdir "${OUT_DIR}/upregulatedGenes/up-regulated_miRNAs"

fi

cd "${IN_DIR}"
mv *downGenes_all* "${OUT_DIR}/downregulatedGenes/all_miRNAs/"
mv *downGenes_down* "${OUT_DIR}/downregulatedGenes/down-regulated_miRNAs/"
mv *downGenes_up* "${OUT_DIR}/downregulatedGenes/up-regulated_miRNAs/"
mv *upGenes_all* "${OUT_DIR}/upregulatedGenes/all_miRNAs/"
mv *upGenes_down* "${OUT_DIR}/upregulatedGenes/down-regulated_miRNAs/"
mv *upGenes_up* "${OUT_DIR}/upregulatedGenes/up-regulated_miRNAs/"

cd "${OUT_DIR}/downregulatedGenes/all_miRNAs"
python "$REMOVE" \
	*_downGenes_allResults.txt \
	${PREFIX}_downGenes_allSimplified.txt

cd "${OUT_DIR}/downregulatedGenes/up-regulated_miRNAs"
python "$REMOVE" \
	*_downGenes_upMirnaSigResultsFDR0.05.txt \
	${PREFIX}_downGenes_upMirnaSigSimplified.txt
python "$EXTRACT" \
	*_downGenes_upMirnaSigResultsFDR0.05.txt \
	"$DATABASE" "$RANK_FILE" "down" \
	${PREFIX}_downGenes_upMirnaSigGeneList.txt
python "$TDP43" \
	${PREFIX}_downGenes_upMirnaSigGeneList.txt \
	"$TDP43_FILE" \
	${PREFIX}_downGenes_upMirnaTDP43SigGeneList.txt
cut -f 9 ${PREFIX}_downGenes_upMirnaTDP43SigGeneList.txt | tail -n +2 | \
	sort | uniq > ${PREFIX}_downGenes_upMirnaTDP43FatigoInput.txt

cd "${OUT_DIR}/downregulatedGenes/down-regulated_miRNAs"
python "$REMOVE" \
	*_downGenes_downMirnaSigResultsFDR0.05.txt \
	${PREFIX}_downGenes_downMirnaSigSimplified.txt
python "$EXTRACT" \
	*_downGenes_downMirnaSigResultsFDR0.05.txt \
	"$DATABASE" "$RANK_FILE" "down" \
	${PREFIX}_downGenes_downMirnaSigGeneList.txt
python "$TDP43" \
	${PREFIX}_downGenes_downMirnaSigGeneList.txt \
	"$TDP43_FILE" \
	${PREFIX}_downGenes_downMirnaTDP43SigGeneList.txt

cd "${OUT_DIR}/upregulatedGenes/up-regulated_miRNAs"
python "$REMOVE" \
	*_upGenes_upMirnaSigResultsFDR0.05.txt \
	${PREFIX}_upGenes_upMirnaSigSimplified.txt
python "$EXTRACT" \
	*_upGenes_upMirnaSigResultsFDR0.05.txt \
	"$DATABASE" "$RANK_FILE" "up" \
	${PREFIX}_upGenes_upMirnaSigGeneList.txt
python "$TDP43" \
	${PREFIX}_upGenes_upMirnaSigGeneList.txt \
	"$TDP43_FILE" \
	${PREFIX}_upGenes_upMirnaTDP43SigGeneList.txt


cd "${OUT_DIR}/upregulatedGenes/down-regulated_miRNAs"
python "$REMOVE" \
	*_upGenes_downMirnaSigResultsFDR0.05.txt \
	${PREFIX}_upGenes_downMirnaSigSimplified.txt
python "$EXTRACT" \
	*_upGenes_downMirnaSigResultsFDR0.05.txt \
	"$DATABASE" "$RANK_FILE" "up" \
	${PREFIX}_upGenes_downMirnaSigGeneList.txt
python "$TDP43" \
	${PREFIX}_upGenes_downMirnaSigGeneList.txt \
	"$TDP43_FILE" \
	${PREFIX}_upGenes_downMirnaTDP43SigGeneList.txt
cut -f 9 ${PREFIX}_upGenes_downMirnaTDP43SigGeneList.txt | tail -n +2 | \
	sort | uniq > ${PREFIX}_upGenes_downMirnaTDP43FatigoInput.txt

