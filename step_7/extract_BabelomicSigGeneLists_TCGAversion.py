#!/opt/loca/bin/python
import sys, re, itertools, optparse
from math import *

optParser = optparse.OptionParser( 
   
   usage = "python %prog [options] <modified_SPIA_file> <ucscCombinedTable> <DESeqRankFile> <direction> <out_file>",
   
   description=
      "This script takes output from the modified_SPIA R pipeline " +
      "and extracts the list of DE genes that made each miRNA significant. " + 
      "This includes info from UCSC about the transcript gene symbol and description, " + 
      "as well as info about the DESeq rank used for the analysis.",
      
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
ucsc_file = args[1]
rank_file = args[2]

if args[3] == "up":
	up_bool = True
elif args[3] == "down":
	up_bool = False
else:
	sys.stderr.write( sys.argv[0] + ": Error: Fourth argument \"direction\" should be 'up' or 'down'.\n")
	sys.stderr.write( "It was " + args[3] + "\n")
	sys.exit( 1)
	
out_file = args[4]

ucsc_dict = {}
mir_dict = {}

with open(ucsc_file, 'r') as ucsc:
	header = ucsc.readline()
	for line in ucsc.readlines():
		lineList = line.strip("\n").split("\t")
		ucsc_id = lineList[0]
		gene_id = lineList[1]
		desc = lineList[2]
		
		assert ucsc_id not in ucsc_dict, ucsc_id + " is listed twice in the database file."
		ucsc_dict[ucsc_id] = [gene_id, desc]

with open(rank_file, 'r') as rank:
	header = rank.readline()
	for line in rank.readlines():
		lineList = line.strip("\n").split("\t")
		ucsc_id = lineList[0]
		rank = lineList[1]
		ucsc_dict[ucsc_id].append(rank)

with open(in_file, 'r') as input:
	header = input.readline()
	for line in input.readlines():
		lineList = line.split("\t")
		mir_acc = lineList[0]
		mir_name = lineList[1]
		num_targets = lineList[4]
		pNDE = lineList[-5]
		pFatiscan = lineList[-4]
		pGFDR = lineList[-2]
		
		ucsc_list = lineList[7].split(",")
		mir_dict[mir_acc] = [mir_name, num_targets, pNDE, pFatiscan, pGFDR, ucsc_list]
		

with open(out_file, 'w') as output:
	header = "\t".join(["miRBase Accession","Mature miRNA Name","Number Predicted Targets","pNDE","pFatiscan","pGFDR",
		"UCSC ID","DESeq signedLogP","Gene Name","Gene Description"])+"\n"
	output.write(header)
	for mir in sorted(mir_dict.keys()):
		mirList = [mir] + mir_dict[mir][:5]
		ucscList = mir_dict[mir][-1]
		for ucsc_id in ucscList:
			if ucsc_id == "":
				continue
			assert ucsc_id in ucsc_dict, mir + " does not have " + ucsc_id + " as a target in " + in_file
				
			try:
				gene_id, desc, rank = ucsc_dict[ucsc_id]
			except:
				print >> sys.stderr, ucsc_id + " does not have a rank: " + str(ucsc_dict[ucsc_id])
				sys.exit(1)
			if (up_bool and float(rank) > 0) or (not up_bool and float(rank) < 0):
				newLine = "\t".join(mirList + [ucsc_id, rank, gene_id, desc]) + "\n"
				output.write(newLine)
