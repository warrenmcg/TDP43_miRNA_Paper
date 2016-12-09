#!/opt/loca/bin/python
import sys, re, itertools, optparse
from math import *

optParser = optparse.OptionParser( 
   
   usage = "python %prog [options] <bableomics extended Gene List> <TDP43 correlations file> <out_file>",
   
   description=
      "This script takes output from the extract_BabelomicSigGeneLists_TCGAversion.py script" +
      "and extracts the list of miRNAs and their targets that are DE and related to TDP-43. " + 
      "This includes info from UCSC about the transcript gene symbol and description, " + 
      "as well as info about the DESeq rank used for the analysis.",
      
   epilog = 
      "Written by Warren McGee (warren-mcgee@fsm.northwestern.edu), " +
      "Jane Wu Lab, Northwestern University, Chicago, USA. (c) 2014. " + 
      "Released under the terms of the GNU General Public License v3." )

(opts, args) = optParser.parse_args()


if len( args ) != 3:
   sys.stderr.write( sys.argv[0] + ": Error: Please provide three arguments.\n" )
   sys.stderr.write( "  Call with '-h' to get usage information.\n" )
   sys.exit( 1 )

in_file = args[0]
tdp43_file = args[1]
out_file = args[2]

tdp_dict = {}
mir_dict = {}

with open(tdp43_file, 'r') as tdp43:
	for line in tdp43.readlines():
		lineList = line.strip("\n").split("\t")
		name = lineList[0]
		
		assert name not in tdp_dict, name + " is listed twice in the database file."
		tdp_dict[name] = True

with open(in_file, 'r') as input, open(out_file, 'w') as output:
	header = input.readline()
	output.write(header)
	for line in input.readlines():
		lineList = line.split("\t")
		mir_name = lineList[1]
		pNDE = lineList[3]
		if mir_name in tdp_dict and float(pNDE) <= 0.1:
			output.write(line)
