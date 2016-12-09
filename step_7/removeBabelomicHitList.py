#!/opt/loca/bin/python
import sys, re, itertools, optparse
from math import *

optParser = optparse.OptionParser( 
   
   usage = "python %prog [options] <modified_SPIA_file> <out_file>",
   
   description=
      "This script takes output from the modified_SPIA R pipeline " +
      "and removes the unwieldy and really long list of hits for each miRNA.",
      
   epilog = 
      "Written by Warren McGee (warren-mcgee@fsm.northwestern.edu), " +
      "Jane Wu Lab, Northwestern University, Chicago, USA. (c) 2014. " + 
      "Released under the terms of the GNU General Public License v3." )

(opts, args) = optParser.parse_args()


if len( args ) != 2:
   sys.stderr.write( sys.argv[0] + ": Error: Please provide two arguments.\n" )
   sys.stderr.write( "  Call with '-h' to get usage information.\n" )
   sys.exit( 1 )

in_file = args[0]
out_file = args[1]

with open(in_file, 'r') as input:
	with open(out_file, 'w') as output:
		header = input.readline()
		headerList = header.split("\t")
		badIndex = headerList.index("DEgeneList")
		
		headerList.pop(badIndex)
		newHeader = "\t".join(headerList)
		output.write(newHeader)
		
		for line in input.readlines():
			lineList = line.split("\t")
			
			lineList.pop(badIndex)
			newLine = "\t".join(lineList)
			output.write(newLine)
