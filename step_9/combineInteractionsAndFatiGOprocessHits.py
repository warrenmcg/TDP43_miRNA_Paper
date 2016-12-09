#!/opt/loca/bin/python
import sys, re, itertools, optparse
from math import *

optParser = optparse.OptionParser( 
   
   usage = "python %prog [options] <TDP-43-regulated extended Gene List> " + 
   			"<TermID2Name file> <Gene Alias file> <FatiGO output file> <out_prefix>",
   			
   description=
      "This script takes output from the extract_TDP43miRNAsFromBabelomicsGeneList.py script " +
      "and outputs two files: one file contains an extended list of all miRNA-mRNA interactions " +
      "with each term annotated for that gene (gene without annotations are ignored), and a second " +
      "file with an edge list for input into Cytoscape.",
      
   epilog = 
      "Written by Warren McGee (warren-mcgee@fsm.northwestern.edu), " +
      "Jane Wu Lab, Northwestern University, Chicago, USA. (c) 2014. " + 
      "Released under the terms of the GNU General Public License v3." )

(opts, args) = optParser.parse_args()


if len( args ) != 5:
   sys.stderr.write( sys.argv[0] + ": Error: Please provide five arguments.\n" )
   sys.stderr.write( "  Call with '-h' to get usage information.\n" )
   sys.exit( 1 )

in_file = args[0]
id2name_file = args[1]
alias_file = args[2]
fatigo_file = args[3]
out_table = args[4] + "_List.txt"
out_edge = args[4] + "_edgeList.txt"

mirna_dict = {}
ucsc_dict = {}
alias_dict = {}
annot_dict = {}
term_dict = {}
edge_dict = {}

with open(in_file, 'r') as input:
	header = input.readline()
	for line in input.readlines():
		lineList = line.strip("\n").split("\t")
		mir_acc = lineList[0]
		name = lineList[1]
		pDE = lineList[3]
		pFatiscan = lineList[4]
		pGFDR = lineList[5]
		ucsc_id = lineList[6]
		logp = lineList[7]
		gene_name = lineList[8]
		gene_desc = lineList[9]
		
		if mir_acc in mirna_dict:
			mirna_dict[mir_acc][4].append(ucsc_id)
			ucsc_dict[ucsc_id] = [gene_name, gene_desc, logp]
		else:
			mirna_dict[mir_acc] = [name, pDE, pFatiscan, pGFDR, [ucsc_id]]
			ucsc_dict[ucsc_id] = [gene_name, gene_desc, logp]
		
with open(id2name_file, 'r') as idname:
	for line in idname.readlines():
		lineList = line.strip("\n").split("\t")
		category = lineList[0]
		id = lineList[1]
		name = lineList[2]
		term_dict[id] = [name, category]

with open(alias_file, 'r') as alias:
	for line in alias.readlines():
		lineList = line.strip("\n").split("\t")
		name = lineList[2]
		alias = lineList[1]
		if alias != "":
			alias_dict[name] = alias

with open(fatigo_file, 'r') as fatigo:
	header = fatigo.readline()
	for line in fatigo.readlines():
		lineList = line.strip("\n").split("\t")
		id = lineList[0]
		size = lineList[2]
		percent = lineList[5]
		geneList = lineList[-5].split(",")
		OR = lineList[-3]
		adjPval = lineList[-1]
		if float(adjPval) > 0.05:
			continue
		
		term_dict[id] = term_dict[id] + [size, percent, OR, adjPval]
		for gene in geneList:
			if gene in annot_dict:
				annot_dict[gene].append(id)
			else:
				annot_dict[gene] = [id]

with open(out_table, 'w') as output:
	headerList = ["miRBase Accession", "mature miRNA name", "pDESeq2", "pFatiscan", 
		"pCombinedFDR","UCSC ID","DESeq2 Signed LogP", "Gene Name", "Gene Description", 
		"Term Category","Term ID","Term Name","Term Size","Percent","Log Odds Ratio","pFatiGO"]
	header = "\t".join(headerList) + "\n"
	output.write(header)
	for mir_acc in mirna_dict:
		mir_name, pDE, pFatiscan, pGFDR, ucsc_list = mirna_dict[mir_acc]
		for ucsc_id in ucsc_list:
			gene_name, desc, logp = ucsc_dict[ucsc_id]
			if gene_name in alias_dict:
				gene_name = alias_dict[name]
			
			if gene_name in annot_dict:
				term_list = annot_dict[gene_name]
			else:
				continue
			
			if gene_name in edge_dict:
				edge_dict[gene_name].add(mir_name)
			else:
				edge_dict[gene_name] = set([mir_name])
				
			for term_id in term_list:
				term_name, category, size, percent, OR, adjPval = term_dict[term_id]
				edge_dict[gene_name].add(term_name)
				lineList = [mir_acc, mir_name, pDE, pFatiscan, pGFDR, ucsc_id, logp, 
					gene_name, desc, category, term_id, term_name, size, percent, OR, adjPval]
				line = "\t".join(lineList) + "\n"
				output.write(line)

with open(out_edge, 'w') as output:
	headerList = ["outgoing Node", "incoming Node", "interaction"]
	header = "\t".join(headerList) + "\n"
	output.write(header)
	for gene_name in edge_dict:
		set = edge_dict[gene_name]
		for item in set:
			if re.match("hsa-",item):
				category = "miRNA"
				lineList = [item, gene_name, category]
			else:
				category = "term"
				lineList = [gene_name, item, category]
			newLine = "\t".join(lineList) + "\n"
			output.write(newLine)