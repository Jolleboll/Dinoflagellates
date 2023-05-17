#!/usr/bin/python3

import argparse,re,sys

############ Input and such ############
desc	='''This program takes the fst output file from a BayeScan run and parses it like a champ.
It uses output from a previous stacks/populations run in order to translate loci names (this is the "s" argument).'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('i',metavar='bayesOut',help='Path to input file.')
parser.add_argument('s',metavar='structure',help='Path to the stacks populations output file in structure format.')
parser.add_argument('-o',metavar='output',help='Output file basename. Default derived from input.')
parser.add_argument('-p',metavar='probCut',help='Probability cut-off for the "prob" column in the input file. Default 0.95.',type=float,default=0.95)
args	=parser.parse_args()
########################################



######## Parsing structure file ########
with open(args.s) as f:
	next(f) # Ignore the stacks header
	lociLine=next(f).split() # The header row of the matrix. Not aligned but who cares (not structure).
########################################



########## Parsing/Printing ############

if args.o: outfile=args.o
else: outfile=args.i+'.translated'

with open(args.i) as f, open(outfile,'w') as out:
	next(f) # Ignore header. The columns are: locus, probability, log10(PO), qvalue, alpha, fst
	print("Locus ID","Probability","Log10(PO)","Q-value","Alpha","Fst",sep='\t',file=out) # Better header
	for line in f:
		locus,prob,log,qvalue,alpha,fst=line.split()
		if float(prob) > args.p:
			print(lociLine[int(locus)-1],prob,log,qvalue,alpha,fst,sep='\t',file=out)
########################################