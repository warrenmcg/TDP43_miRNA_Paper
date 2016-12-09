#!/opt/loca/bin/python
import sys, re, itertools, optparse
from math import *

optParser = optparse.OptionParser( 
   
   usage = "python %prog [options] <ucscTable1> <>ucscTable2> <ucscTable3> <out_prefix>",
   
   description=
      "This takes multiple kgXref tables from UCSC and combines them into one table, " + 
      "with both old and new entries included. Please list the tables in decreasing " +
      "chronological order (newest table first; oldest table last). Output is an abbreviated " +
      "table, with just UCSC ID, gene name, and description.",
      
   epilog = 
      "Written by Warren McGee (warren-mcgee@fsm.northwestern.edu), " +
      "Jane Wu Lab, Northwestern University, Chicago, USA. (c) 2015. " + 
      "Released under the terms of the GNU General Public License v3." )

(opts, args) = optParser.parse_args()


if len( args ) != 4:
   sys.stderr.write( sys.argv[0] + ": Error: Please provide four arguments.\n" )
   sys.stderr.write( "  Call with '-h' to get usage information.\n" )
   sys.exit( 1 )

kgxref1 = args[0]
kgxref2 = args[1]
kgxref3 = args[2]
out_file = args[3]

ucsc_dict = {}

with open(kgxref1, 'r') as input1:
#	header = input1.readline()
	for line in input1.readlines():
		lineList = line.strip("\n").split("\t")
		kg_id = lineList[0]
		gene_id = lineList[4]
		gene_desc = lineList[7]
		entry = [gene_id, gene_desc]
		assert kg_id not in ucsc_dict, kg_id + " is listed twice in kgXref table!"
		ucsc_dict[kg_id] = entry

with open(kgxref2, 'r') as input2:
#	header = input2.readline()
	for line in input2.readlines():
		lineList = line.strip("\n").split("\t")
		kg_id = lineList[0]
		gene_id = lineList[4]
		gene_desc = lineList[7]
		entry = [gene_id, gene_desc]

		if kg_id not in ucsc_dict:
			ucsc_dict[kg_id] = entry
	
with open(kgxref3, 'r') as input3:
#	header = input3.readline()
	for line in input3.readlines():
		lineList = line.strip("\n").split("\t")
		kg_id = lineList[0]
		gene_id = lineList[4]
		gene_desc = lineList[7]
		entry = [gene_id, gene_desc]

		if kg_id not in ucsc_dict:
			ucsc_dict[kg_id] = entry

with open(out_file, 'w') as output:
	header = "\t".join(["ucscID","geneSymbol","geneDescription"])+"\n"
	output.write(header)
	for kg_id in sorted(ucsc_dict.keys()):
		entry = ucsc_dict[kg_id]
		newLine = "\t".join([kg_id] + entry) + "\n"
		output.write(newLine)
