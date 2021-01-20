#!/usr/bin/python3

import argparse,re,gzip

############ Input and such ############
desc	='''This program takes a list of loci and outputs the corresponding fasta file.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('i',metavar='infile',help='Path to your loci input file.')
parser.add_argument('c',metavar='catalogue',help='Path to your cstacks tags file (gzipped).')
parser.add_argument('o',metavar='outfile',help='Name of your output file.')
args	=parser.parse_args()
########################################

################ Loci ##################
with open(args.i) as f: myloci={locus.rstrip() for locus in f}
########################################

############## Parsing #################
with gzip.open(args.c,'rt') as f, open(args.o,'w') as out:
	next(f) # Ignoring the header line
	for line in f:
		fields	=line.split('\t')
		ID		=fields[2]
		sequence=fields[9]
		if ID in myloci:
			out.write('>'+ID+'\n'+sequence+'\n')
########################################
print('Done with',args.i)