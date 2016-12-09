#!/opt/loca/bin/python
import sys, re, itertools, optparse
from math import *

optParser = optparse.OptionParser( 
   
   usage = "python %prog [options] <number of tables to combine> <miRNA DESeq Table> " +
   	"<transcript DESeq table> <alias table> <That number of tables> <out_prefix>",
   
   description=
      "This script takes output tables from the combineInteractionsAndFatiGOprocessHits.py script " +
      "and outputs a table that is suitable for the online Hive Plot tool or for further processing " +
      "into a DOT table for HiveR.",
     
   epilog = 
      "Written by Warren McGee (warren-mcgee@fsm.northwestern.edu), " +
      "Jane Wu Lab, Northwestern University, Chicago, USA. (c) 2014. " + 
      "Released under the terms of the GNU General Public License v3." )

(opts, args) = optParser.parse_args()

num_tables = int(args[0])
tables = args[4:-1]
if(len(tables) < num_tables):
    tables = tables[0].split(" ")[1:]

if len( args ) != num_tables + 5 and len(tables) != num_tables:
   sys.stderr.write( sys.argv[0] + ": Error: Please provide the number of tables as you specified.\n" )
   sys.stderr.write( "  Call with '-h' to get usage information.\n" )
   sys.exit( 1 )

node_out = args[-1] + "_nodeList.txt"
edge_out = args[-1] + "_edgeList.txt"
mirna_file = args[1]
transcript_file = args[2]
alias_file = args[3]

alias_dict = {}
deseq_dict = {}
mirna_dict = {}
gene_dict = {}
term_dict = {}
edge_dict = {}
max_mirRank = max_mirEdgeRank = max_mirLogp = min_mirLogp = 0
max_geneRank = max_geneLogp = min_geneLogp = 0
max_termRank = 0


def processDeseq(fileName):
	with open(fileName, 'r') as input:
		header = input.readline()
		for line in input.readlines():
			lineList = line.strip("\n").split("\t")
			id = lineList[0]
			baseMean = float(lineList[1]) + 1
			if lineList[2] == "NA":
				log2FC = 0
			else:
				log2FC = float(lineList[2])
			if lineList[-1] == "NA":
				padj = 1
			else:
				padj = float(lineList[-1])
				
			
			if baseMean < 30:
				if log2FC < 0:
					overall_rank = log2FC - log10(baseMean)
				else:
					overall_rank = log2FC + log10(baseMean)
			else:
				if log2FC < 0:
					overall_rank = log2FC + log10(padj) - log10(baseMean)
				else:
					overall_rank = log2FC + abs(log10(padj)) + log10(baseMean)
			deseq_dict[id] = overall_rank

def processTable(fileName):
	global max_mirLogp, min_mirLogp, max_geneLogp, min_geneLogp
	global max_mirRank, max_mirEdgeRank, max_geneRank, max_termRank
	
	with open(fileName, 'r') as input:
		header = input.readline()
		for line in input.readlines():
			
			lineList = line.strip("\n").split("\t")
			mir_acc, mir_name, pDE, pFatiscan, pCombined = lineList[:5]
			ucsc_id, ucsc_logp, gene_name, gene_desc = lineList[5:9]
			category, term_id, term_name, term_size, percent, log_odds, pFatigo = lineList[9:]
			
			mirna_rank = deseq_dict[mir_acc]
			mirna_edgeRank = -log10( float( pFatiscan ) )
			ucsc_rank = deseq_dict[ucsc_id]
			ucsc_fdr = str( 10**( -abs( float( ucsc_logp ) ) ) )
			term_rank = str( -log10( float( pFatigo ) ) )
			
			# calculate the logP of the miRNA's differential expression
			if mirna_rank < 0:
				mir_logp = log10(float(pDE))
			else:
				mir_logp = -log10(float(pDE))

			# add the miRNA, gene, and/or pathway/GO term if it's not already been stored
			if mir_acc not in mirna_dict:
				mirna_dict[mir_acc] = [mir_name, 0, 0, str(mirna_edgeRank), pCombined, mir_logp, mirna_rank]
			if term_id not in term_dict:
				term_dict[term_id] = [term_name, 0, 0, category, pFatigo, term_rank]
			if ucsc_id not in gene_dict:
				gene_dict[ucsc_id] = [gene_name, 0, 0, ucsc_fdr, ucsc_logp, ucsc_rank]
			
			# add the edge if it's not already been stored; increment the edge count
			if mir_acc in edge_dict:
				if ucsc_id not in edge_dict[mir_acc]:
					edge_dict[mir_acc].append(ucsc_id)
					mirna_dict[mir_acc][1] += 1
					gene_dict[ucsc_id][2] += 1
			else:
				edge_dict[mir_acc] = [ucsc_id]
				mirna_dict[mir_acc][1] += 1
				gene_dict[ucsc_id][2] += 1
				
			if term_id in edge_dict:
				if ucsc_id not in edge_dict[term_id]:
					edge_dict[term_id].append(ucsc_id)
					gene_dict[ucsc_id][1] += 1
					term_dict[term_id][2] += 1
			else:
				edge_dict[term_id] = [ucsc_id]
				gene_dict[ucsc_id][1] += 1
				term_dict[term_id][2] += 1
			
			if gene_name in alias_dict:
				gene_name = alias_dict[gene_name]
			
			# update the max values if necessary (use to normalize the colors and node height)
			if abs(mirna_rank) > max_mirRank:
				max_mirRank = abs(mirna_rank)
						
			if float(mir_logp) > max_mirLogp:
				max_mirLogp = float(mir_logp)
			elif float(mir_logp) < min_mirLogp:
				min_mirLogp = float(mir_logp)
			
			if mirna_edgeRank > max_mirEdgeRank:
				max_mirEdgeRank = mirna_edgeRank
			
			if abs(ucsc_rank) > max_geneRank:
				max_geneRank = abs(ucsc_rank)
			
			if float(ucsc_logp) > max_geneLogp:
				max_geneLogp = float(ucsc_logp)
			elif float(ucsc_logp) < min_geneLogp:
				min_geneLogp = float(ucsc_logp)
			
			if float(term_rank) > max_termRank:
				max_termRank = float(term_rank)
			
def returnMirNodeColor(mir_acc):
	mir_logp = float( mirna_dict[mir_acc][-2] )
	mir_rank = float( mirna_dict[mir_acc][-1] )
	if mir_logp < 0:
		red = int(ceil(255*(1- mir_logp / min_mirLogp)))
		green = int(ceil(255*(1- mir_logp / min_mirLogp)))
		blue = 255
		alpha = int(ceil(150*(abs(mir_rank) / max_mirRank))) + 50
	else:
		red = 255
		green = int(ceil(255*(1 - mir_logp / max_mirLogp)))
		blue = int(ceil(255*(1 - mir_logp / max_mirLogp)))
		alpha = int(ceil(150*(mir_rank / max_mirRank))) + 50

	color = '#%02x%02x%02x%02x' % tuple([red, green, blue, alpha])
	if re.search("-", color):
		print >> sys.stderr, mir_acc + " with mir_rank/max_rank " + str([mir_rank, max_mirRank]) + \
			" and mir_logp/max/min_logp " + str([mir_logp, max_mirLogp, min_mirLogp]) + " produced the invalid hex color " + color
		sys.exit(1)
	return [color, red, green, blue, alpha]

def returnMirEdgeColor(mir_acc):
	mir_fatiscan = float( mirna_dict[mir_acc][3] )
	mir_rank = float( mirna_dict[mir_acc][-1] )
	
	if mir_rank < 0:
		red = int(ceil(255*(1 - mir_fatiscan / max_mirEdgeRank)))
		green = int(ceil(255*(1 - mir_fatiscan / max_mirEdgeRank)))
		blue = 255
		alpha = int(ceil(150*(abs(mir_rank) / max_mirRank))) + 50
	else:
		red = 255
		green = int(ceil(255*(1 - mir_fatiscan / max_mirEdgeRank)))
		blue = int(ceil(255*(1 - mir_fatiscan / max_mirEdgeRank)))
		alpha = int(ceil(150*(mir_rank / max_mirRank))) + 50

	color = '#%02x%02x%02x%02x' % tuple([red, green, blue, alpha])
	return [color, red, green, blue, alpha]
	
def returnGeneNodeColor(ucsc_id):
	gene_logp = float( gene_dict[ucsc_id][-2] )
	gene_rank = float( gene_dict[ucsc_id][-1] )
	if gene_logp < 0:
		red = int(ceil(255*(1-gene_logp / min_geneLogp)))
		green = int(ceil(255*(1-gene_logp / min_geneLogp)))
		blue = 255
		alpha = int(ceil(150*(abs(gene_rank) / max_geneRank))) + 50
	else:
		red = 255
		green = int(ceil(255*(1-gene_logp / max_geneLogp)))
		blue = int(ceil(255*(1-gene_logp / max_geneLogp)))
		alpha = int(ceil(150*(gene_rank / max_geneRank))) + 50

	color = '#%02x%02x%02x%02x' % tuple([red, green, blue, alpha])
	if re.search("-", color):
		print >> sys.stderr, ucsc_id + " with gene_rank/max_rank " + str([gene_rank, max_geneRank]) + \
			" and gene_logp/max/min_logp " + str([gene_logp, max_geneLogp, min_geneLogp]) + " produced the invalid hex color " + color
		sys.exit(1)
	return [color, red, green, blue, alpha]

def returnNodeValue(id, dict):
	return abs(float( dict[id][-1] ))

def returnTermNodeColor(term_id):
	category = term_dict[term_id][3]
	alpha = 200
	if category == "KEGG":
		red = 230
		green = 159
		blue = 0
	elif category == "REACTOME":
		red = 0
		green = 158
		blue = 115
	elif category == "BIOCARTA":
		red = 86
		green = 180
		blue = 233
	elif category == "GOBP":
		red = 204
		green = 121
		blue = 167
	
	color = '#%02x%02x%02x%02x' % tuple([red, green, blue, alpha])
	return [color, red, green, blue, alpha]

def returnEdgeColor(rgbList, degree, maxDegree):
	alpha = 255*(log(degree) / log(maxDegree))
	colorList = rgbList + [alpha]
	color = '#%02x%02x%02x%02x' % tuple(colorList)
	return [color] + colorList

def returnLargestDegree(dict):
	maxDegree = 0
	for item in dict.keys():
		degree = dict[item][1] + dict[item][2]
		if degree > maxDegree:
			maxDegree = degree
	return maxDegree

def returnMirNodeHeight(mir_acc):
	degree = mirna_dict[mir_acc][1] + mirna_dict[mir_acc][2]
	value = int(ceil(3*degree / max_mirDegree)) + 1
	return value

def returnGeneNodeHeight(ucsc_id):
	degree = gene_dict[ucsc_id][1] + gene_dict[ucsc_id][2]
	value = int(ceil(3*degree / max_geneDegree)) + 1
	return value

def returnTermNodeHeight(term):
	degree = term_dict[term][1] + term_dict[term][2]
	value = int(ceil(3 * degree / max_termDegree)) + 1
	return value

with open(alias_file, 'r') as alias:
	for line in alias.readlines():
		lineList = line.strip("\n").split("\t")
		name = lineList[2]
		alias = lineList[1]
		if alias != "":
			alias_dict[name] = alias

processDeseq(mirna_file)
processDeseq(transcript_file)

for file in tables:
	processTable(file)

max_mirDegree = returnLargestDegree(mirna_dict)
max_geneDegree = returnLargestDegree(gene_dict)
max_termDegree = returnLargestDegree(term_dict)

with open(node_out, 'w') as output:
	headerList = ["id", "axis", "node_label", "node_value", "node_height",
		"node_color","edge_color","num_Outgoing", "num_Incoming","node_r","node_g","node_b",
		"node_a","edge_r","edge_g","edge_b","edge_a","group","FDR","Rank"]
	header = "\t".join(headerList) + "\n"
	output.write(header)
	axis = "TDP-43-regulated_miRNA"
	for mir_acc in sorted(mirna_dict.keys()):
		node_label, num_out, num_in, edgeRank, fdr, logp, rank = mirna_dict[mir_acc]
		node_height = returnMirNodeHeight(mir_acc)
		node_value = returnNodeValue(mir_acc, mirna_dict)
		if float(logp) < 0:
			group = "upGenesDownMiRNAs"
		else:
			group = "downGenesUpMiRNAs"
		nodeColorList = returnMirNodeColor(mir_acc)
		edgeColorList = returnMirEdgeColor(mir_acc)
		lineList = [mir_acc, axis, node_label, str(node_value), str(node_height), 
			nodeColorList[0], edgeColorList[0], str(num_out), str(num_in)]
		newLine = "\t".join(lineList) + \
			"\t" + "\t".join(str(x) for x in nodeColorList[1:]) + \
			"\t" + "\t".join(str(x) for x in edgeColorList[1:]) + \
			"\t" + "\t".join([group, fdr, str(rank)]) + "\n"
		output.write(newLine)
	
	axis = "Target_Gene"
	for ucsc_id in sorted(gene_dict.keys()):
		node_label, num_out, num_in, fdr, logp, rank = gene_dict[ucsc_id]
		node_height = returnGeneNodeHeight(ucsc_id)
		node_value = returnNodeValue(ucsc_id, gene_dict)
		if float(logp) < 0:
			group = "downGenesUpMiRNAs"
		else:
			group = "upGenesDownMiRNAs"
		nodeColorList = returnGeneNodeColor(ucsc_id)
		#edgeColorList = returnEdgeColor(nodeColorList[1:4], node_value, max_geneDegree)
		edge_color = "#00000000"
		edge_r = edge_g = edge_b = edge_a = ""
		lineList = [ucsc_id, axis, node_label, str(node_value), str(node_height),
			nodeColorList[0], edge_color, str(num_out), str(num_in)]
		newLine = "\t".join(lineList) + \
			"\t" + "\t".join(str(x) for x in nodeColorList[1:]) + \
			"\t" + "\t".join([edge_r, edge_g, edge_b, edge_a, group, fdr, str(rank)]) + "\n"
		output.write(newLine)
		
	axis = "Fatigo"
	for term in sorted(term_dict.keys()):
		node_label, num_out, num_in, category, fdr, rank = term_dict[term]
		node_height = returnTermNodeHeight(term)
		node_value = returnNodeValue(term, term_dict)
		nodeColorList = returnTermNodeColor(term)
		edgeColorList = returnEdgeColor(nodeColorList[1:4], node_value, max_termDegree)
		lineList = [term, axis, node_label, str(node_value), str(node_height),
			nodeColorList[0], edgeColorList[0], str(num_out), str(num_in)]
		newLine = "\t".join(lineList) + \
			"\t" + "\t".join(str(x) for x in nodeColorList[1:]) + \
			"\t" + "\t".join(str(x) for x in edgeColorList[1:]) + \
			"\t" + "\t".join([category, fdr, str(rank)]) + "\n"
		output.write(newLine)

with open(edge_out, 'w') as output:
	#header = "\t".join(["Outgoing Node", "Incoming Node"]) + "\n"
	#output.write(header)
	for target in sorted(edge_dict.keys()):
		nodeList = edge_dict[target]
		for node in nodeList:
			newLine = "\t".join([node, target]) + "\n"
			output.write(newLine)