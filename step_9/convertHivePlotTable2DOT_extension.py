#!/usr/bin/python
import sys, re, itertools, optparse
from collections import OrderedDict

optParser = optparse.OptionParser( 
   
   usage = "python %prog [options] <hiveplot_nodeList> <hiveplot_edgeList> <output_file_prefix>",
   
   description=
      "This scripts converts a table for the HivePlot web tool into a directed graph dot file",
      
   epilog = 
      "Written by Warren McGee (warren-mcgee@fsm.northwestern.edu), Jane Wu Lab" +
      "(Northwestern, Chicago, USA). (c) 2014. Released under the terms of the" +
      "GNU General Public License v3." )
   
if len( sys.argv ) == 1:
   optParser.print_help()
   sys.exit(1)

(opts, args) = optParser.parse_args()

if len( args ) != 3:
   sys.stderr.write( sys.argv[0] + ": Error: Please provide three arguments.\n" )
   sys.stderr.write( "  Call with '-h' to get usage information.\n" )
   sys.exit( 1 )

node_file = args[0]
edge_file = args[1]
out_pre = args[2]
dot_file = out_pre + ".dot"
node_i = out_pre + "_nodeInst.csv"
edge_i = out_pre + "_edgeInst.csv"

#initialize counters
node_dict = OrderedDict()
edge_dict = OrderedDict()
axis_set = set()
value_set = set()
height_set = set()
nodeCol_set = set()
nodeType_set = set()
edgeCol_set = set()
edgeType_set = set()

with open(node_file, 'rU') as input:
	header = input.readline()
	headerList = header.split("\t")
	headerList[-1] = headerList[-1].strip("\n")
	for i in range(0,len(headerList)):
		if headerList[i] == "node_label":
			headerList[i] = "label"
	for line in input:
		lineList = line.strip("\n").split("\t")
		node_id = lineList[0]
		if node_id not in node_dict:
			attr_dict = {}
			for i in range(1, len(lineList)):
				if(re.search(" ",lineList[i])):
					lineList[i] = '"' + lineList[i] + '"'
				#if(re.search("color",headerList[i])):
					#lineList[i] = "#" + lineList[i] + "FF"
				attr_dict[headerList[i]] = lineList[i]
				if headerList[i] == "node_value" and lineList[i] not in value_set:
					value_set.add(lineList[i])
				if headerList[i] == "axis" and lineList[i] not in axis_set:
					axis_set.add(lineList[i])
				if headerList[i] == "node_height" and lineList[i] not in height_set:
					height_set.add(lineList[i])
				if headerList[i] == "node_color" and lineList[i] not in nodeCol_set:
					nodeCol_set.add(lineList[i])
				if headerList[i] == "node_type" and lineList[i] not in nodeType_set:
					nodeType_set.add(lineList[i])
				if headerList[i] == "edge_color" and lineList[i] not in edgeCol_set:
					edgeCol_set.add(lineList[i])
#				if headerList[i] == "group" and lineList[i] not in edgeType_set:
#					edgeType_set.add(lineList[i])
			node_dict[node_id] = attr_dict

with open(edge_file, 'rU') as input:
	for line in input:
		lineList = line.strip("\n").split("\t")
		key = tuple([lineList[0], lineList[1]])
		edge_color = node_dict[lineList[1]]["edge_color"]
		edge_type = node_dict[lineList[1]]["group"]
		attr_dict = {"edge_color":edge_color, "edge_type":edge_type}
		if key not in edge_dict:
			edge_dict[key] = attr_dict
	
with open(dot_file, 'w') as output:
	header = "digraph " + out_pre +  " {\n"
	output.write(header)
	for node in node_dict.keys():
		attr_dict = node_dict[node]
		if(re.search(" ",node)):
			node = '"'+ node + '"'
		attr_str = ' '
		for attr in sorted(attr_dict.keys()):
			attr_str += attr + '=' + attr_dict[attr] + ', '
		attr_str = attr_str[:-2]
		line = node + ' [' + attr_str + ' ];\n'
		output.write('\t'+line)

	for key in edge_dict.keys():
		gene_set, gene = list(key)
		if(re.search(" ",gene_set)):
			gene_set = '"' + gene_set + '"'
		if(re.search(" ",gene)):
			gene = '"' + gene + '"'
		edge = gene_set + ' -> ' + gene
		attr_dict = edge_dict[key]
		attr_str = " "
		for attr in sorted(attr_dict.keys()):
			attr_str += attr + '=' + attr_dict[attr] + ', '
		attr_str = attr_str[:-2]
		line = edge + ' [' + attr_str + '];\n'
		output.write(line)
		#output.write('\t'+edge+'\n')
	output.write('}\n')

with open(node_i, 'w') as output:
	header = "dot.tag,dot.val,hive.tag,hive.val\n"
	output.write(header)
	counter = 3
	for axis in sorted(axis_set):
		line = "axis," + axis + ",axis," + str(counter) + "\n"
		output.write(line)
		counter -= 1
	for color in sorted(nodeCol_set):
		line = "node_color," + color + ",color," + color + "\n"
		output.write(line)
	for value in sorted(value_set):
		line = "node_value," + value + ",radius," + value + "\n"
		output.write(line)
	for height in sorted(height_set):
		line = "node_height," + height + ",size," + height + "\n"
		output.write(line)
	for type in sorted(nodeType_set):
		line = "node_type," + type + ",nodeType," + type + "\n"
		output.write(line)
	
with open(edge_i, 'w') as output:
	header = "dot.tag,dot.val,hive.tag,hive.val\n"
	output.write(header)
	for color in sorted(edgeCol_set):
		line = "edge_color," + color + ",color," + color + "\n"
		output.write(line)
#	for type in sorted(edgeType_set):
#		line = "edge_type," + type + ",lineType," + type + "\n"
#		output.write(line)