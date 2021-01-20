#!/usr/bin/python3

import argparse,re,sys,gzip

############ Input and such ############
desc	='''This program takes the output from "parseBayesoutput.py" and from its loci creates a fasta file containing the
corresponding consensus sequences.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('i',metavar='bayesTr',help='Path to input file.')
parser.add_argument('c',metavar='catalog',help='Path to the stacks tags catalog.')
parser.add_argument('-o',metavar='output',help='Output file basename. Default derived from input.')
args	=parser.parse_args()
########################################



########## Parsing BayeScan ############
with open(args.i) as f:
	next(f) # Ignore header. The columns are: locus, probability, log10(PO), qvalue, alpha, fst
	BayeSeqs=dict()
	for line in f:
		# locus,prob,log,qvalue,alpha,fst=line.split()
		locus,*others=line.split()
		BayeSeqs[locus.split('_')[0]]='|'.join(others)
########################################



########## Parsing catalog #############

if args.o: outfile=args.o
# else: outfile='.'.join(args.i.split('.')[:-1])+'.fasta'
else: outfile=args.i.split('_')[0]+'.fasta'

with gzip.open(args.c,'rt') as f, open(outfile,'w') as out:
	next(f) # Ignoring the header line
	for line in f:
		fields	=line.split('\t')
		ID		=fields[2]
		sequence=fields[9]
		if ID in BayeSeqs:
			out.write('>'+ID+' '+BayeSeqs[ID]+'\n'+sequence+'\n')
print('Done! (^_^)')
########################################