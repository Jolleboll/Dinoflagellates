#!/usr/bin/python3

import argparse,re,gzip

############ Input and such ############
desc	='''This program takes sequences from a cstacks .tags. output file and returns a corresponding fasta.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('i',metavar='infile',help='Path to your cstacks tags file (gzipped).')
parser.add_argument('o',metavar='outfile',help='Name of your output file.')
args	=parser.parse_args()
########################################

############## Parsing #################
with gzip.open(args.i,'rt') as f, open(args.o,'w') as out:
	next(f) # Ignoring the header line
	fields		=next(f).split('\t') # First/Last entry newline problem, now solved
	ID			=fields[2]
	sequence	=fields[9]
	out.write('>'+ID+'\n'+sequence)
	for line in f:
		fields	=line.split('\t')
		ID		=fields[2]
		sequence=fields[9]
		out.write('\n'+'>'+ID+'\n'+sequence)
	print('Done! (^_^)')
########################################