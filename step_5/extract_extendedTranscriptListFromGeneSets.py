#!/opt/loca/bin/python
import sys, re, itertools, optparse
from math import *

optParser = optparse.OptionParser( 
   
   usage = "python %prog [options] <gene sets file> <output prefix>",
   
   description=
      "This script takes a GMT file ready for GSEA and generates an annotation " +
      "file that lists each miRNA that targets a particular transcript (like an " +
      "inverse of the GMT file).",
      
   epilog = 
      "Written by Warren McGee (warren-mcgee@fsm.northwestern.edu), " +
      "Jane Wu Lab, Northwestern University, Chicago, USA. (c) 2014. " + 
      "Released under the terms of the GNU General Public License v3." )

(opts, args) = optParser.parse_args()


if len( args ) != 2:
   sys.stderr.write( sys.argv[0] + ": Error: Please provide two arguments.\n" )
   sys.stderr.write( "  Call with '-h' to get usage information.\n" )
   sys.exit( 1 )

gmt_file = args[0]
out_file = args[1] + "_inverseGMTmirbaseAccExtended.txt"
out2_file = args[1] + "_inverseGMTmirnaNameExtended.txt"

ens_dict = {}
mir_dict = {}

with open(gmt_file, 'r') as input:
	for line in input.readlines():
		lineList = line.strip("\n").split("\t")
		mir_acc = lineList[0]
		mir_name = lineList[1]
		
		mir_dict[mir_acc] = mir_name
		
		ensList = lineList[2:]
		
		for ens_id in ensList:
			if ens_id in ens_dict:
				ens_dict[ens_id].append(mir_acc)
			else:
				ens_dict[ens_id] = [mir_acc]

with open(out_file, 'w') as output:
	for ens_id in sorted(ens_dict.keys()):
		mirList = ens_dict[ens_id]
		for mir in mirList:
			newLine = "\t".join([ens_id, mir]) + "\n"
			output.write(newLine)

with open(out2_file, 'w') as output2:
	for ens_id in sorted(ens_dict.keys()):
		mirList = ens_dict[ens_id]
		for mir in mirList:
			mir_name = mir_dict[mir]
			newLine = "\t".join([ens_id, mir_name]) + "\n"
			output2.write(newLine)
