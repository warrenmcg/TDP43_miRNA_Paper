#!/usr/bin/python

import re
import sys

if '-h' in sys.argv:
	print "This script is to extract targets of miRNAs from a miRANDA database."
	print "Syntax for using this script is extract_miRANDA_targets.py [options] prediction_file1 prediction_file2 output_file_prefix"
	print "\nOptions:"
	print "\t-h\tprints this help message and exits"
	print "\nCreated by Warren McGee, 2014"
	sys.exit()
elif len(sys.argv) < 4:
	print "Error: syntax for using this script is extract_miRANDA_targets.py [options] prediction_file1 prediction_file2 output_file_prefix"
	sys.exit()

miranda_file1 = sys.argv[1]
miranda_file2 = sys.argv[2]
output_prefix = sys.argv[3]

pairs_dict = {}
miRNA = {}
genes = {}
print "Processing the first miRanda database file"
with open(miranda_file1, 'r') as mir1:
	for line in mir1.readlines():
		lineList = line.split("\t")
		if re.match("#", lineList[0]):
			continue
 
		mirbase_acc = lineList[0]
		miRNA_name = lineList[1]
		gene_name = lineList[3]
		ucsc_ID = lineList[4]
		ensembl_ID = lineList[5]
		key = (miRNA_name, gene_name)
		if key in pairs_dict:
			pairs_dict[key] += 1
		else:
			pairs_dict[key] = 1
		
		key = mirbase_acc
		if key not in miRNA:
			miRNA[key] = miRNA_name
		
		key = ucsc_ID
		if key not in genes:
			genes[key] = [ensembl_ID, gene_name]

print "Processing the second miRanda database file"
with open(miranda_file2, 'r') as mir2:
	for line in mir2.readlines():
		lineList = line.split("\t")
		if re.match("#", lineList[0]):
			continue
 
		mirbase_acc = lineList[0]
		miRNA_name = lineList[1]
		gene_name = lineList[3]
		ucsc_ID = lineList[4]
		ensembl_ID = lineList[5]
		key = (mirbase_acc, ucsc_ID)
		if key in pairs_dict:
			pairs_dict[key] += 1
		else:
			pairs_dict[key] = 1
		
		key = mirbase_acc
		if key not in miRNA:
			miRNA[key] = miRNA_name
		
		key = ucsc_ID
		if key not in genes:
			genes[key] = [ensembl_ID, gene_name]

print "Creating the target pair matrix"
matrix_file = "_".join([output_prefix, "targetMatrix.txt"])
with open(matrix_file, 'w') as matrix:
	sep = "\t"
	line = ""
	matrix.write(sep)
	for key in sorted(miRNA):
		line = sep.join([line,key])
	matrix.write(line)
	matrix.write("\n")
	line = ""
	for gene_key in sorted(genes):
		line = str(gene_key)
		for miRNA_key in sorted(miRNA):
			key = (miRNA_key, gene_key)
			if key in pairs_dict:
				value = pairs_dict[key]
			else:
				value = 0		
			line = sep.join([line,str(value)])
		matrix.write(line)
		matrix.write("\n")

print "Creating the list of miRNAs"
miRNA_file = "_".join([output_prefix, "miRNAlist.txt"])
with open(miRNA_file, 'w') as miRNAlist:
	sep = "\t"
	line = sep.join(["miRBaseAccession", "miRNA_name"]) + "\n"
	miRNAlist.write(line)
	for key in sorted(miRNA):
		line = sep.join([str(key), miRNA[key]]) + "\n"
		miRNAlist.write(line)

print "Creating the list of mRNA targets"
mRNA_file = "_".join([output_prefix, "transcriptList.txt"])
with open(mRNA_file, 'w') as geneList:
	sep = "\t"
	line = sep.join(["UCSC ID", "Ensembl ID", "gene Name"]) + "\n"
	geneList.write(line)
	for key in sorted(genes):
		line = sep.join([str(key), genes[key][0], genes[key][1]]) + "\n"
		geneList.write(line)
