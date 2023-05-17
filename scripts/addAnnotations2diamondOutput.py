#!/usr/bin/python3

import argparse,re,sys

############ Input and such ############
desc	='''
This program uses the SwissProt annotation downloaded from the MMETSP site
and adds it to the diamond output.
Output is input+.annotated'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('d',metavar='diamond',help='Path to your diamond output file.')
parser.add_argument('a',metavar='anno',help='Path to your (swissprot) annotation file.')
args	=parser.parse_args()
########################################



########## Parsing anno file ###########
anno=dict()
with open(args.a) as f:
	next(f) # ignore header
	for line in f:
		fields=line.rstrip().split('\t')
		name,desc=fields[0],fields[8]
		if not name in anno: anno[name]=[desc]
		else: anno[name].append(desc)
########################################

########## Parsing diamond #############
with open(args.d) as f, open(args.d+'.annotated','w') as out:
	print(next(f).rstrip()+'\t'+'Annotation',file=out) # header
	for line in f:
		line=line.rstrip('\n')+'\t'
		mmetsp=line.split()[1].split('___')
		for m in mmetsp:
			if m in anno: line+=' ||| '.join(anno[m])+'___'
		out.write(line[:-3]+'\n')
########################################