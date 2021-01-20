#!/usr/bin/python3

import argparse,re,sys,shutil

############ Input and such ############
desc	='''
This program takes a structure output file from a stacks 1.35 "populations" run
with the --structure argument given (along with the corresponding popmap file),
and preps it for PGDSpider, as well as outputs some useful information.
Output files end up in working directory.
The popmap filename should be something like "popmap_mypopmap.tsv".
The second to last [_.]-delimited field in the filename will be used by default.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('s',metavar='structure',help='Path to your input structure file.')
parser.add_argument('p',metavar='popmap',help='Path to your popmap file used when creating the input structure file.')
parser.add_argument('-o',metavar='output',help='Output files basename. Default derived from popmap filename.')
args	=parser.parse_args()
########################################



############ Pop numbers ###############
# Since the structure format does not support population names as strings, we must abuse the sample names
with open(args.p) as f:
	sample2pop	=dict()
	pop2num		=dict()
	num2pop		=dict()
	num			=1
	for line in f:
		fields=line.split()
		if fields[1] not in pop2num:
			num2pop[num]=fields[1]
			pop2num[fields[1]]=num
			num+=1
		sample2pop[fields[0]]=fields[1]
########################################



##### Parsprinting structure file ######
if args.o: outfilename=args.o
else: outfilename=re.split('[._]',args.p)[-2]

with open(args.s) as f, open(outfilename+'.structure','w') as out:
	next(f) # Ignore the stacks header
	lociLine=next(f).strip() # The header row of the matrix. Not aligned but who cares (not structure).
	out.write(lociLine+'\n')
	lociList=lociLine.split() # Loci, ordered
	lociDips=[0]*len(lociList) # Will contain diploidity frequencies
	for line in f:
		fields	=line.split()
		fields1	=fields[2:] # First two columns are sample and pop, so ignore them
		line	=next(f)
		fields2	=line.split()[2:]
		for i,(f1,f2) in enumerate(zip(fields1,fields2)):
			if not f1==f2:
				fields1[i]='0' # Missing data woop (also applied to those loci to be blacklisted)
				lociDips[i]+=1
				# sys.stdout.write('In '+outfilename+', '+lociList[i]+' at column #'+str(i)+': '+f1+' vs '+f2+'.\n')
		out.write('\t'.join(fields+fields1)+'\n')
########################################



############ Logs and info #############
with open(outfilename+'.information','w') as out:
	out.write("Structure/BayeScan don't use names for the populations, but numbers. Mapping is:\n")
	for key in pop2num:
		out.write(key+"\t"+str(pop2num[key])+"\n")
	out.write("\nLoci and their occurrences of diploidity are listed in the .diploidity file.\n")
with open(outfilename+'.diploidity','w') as out:
	out.write("# Number of samples is "+str(len(sample2pop))+" in case you were interested.\n")
	for locus,freq in zip(lociList,map(str,lociDips)):
		if not freq=='0':
			out.write(str(locus.split('_')[0])+'\t'+freq+'\n')
########################################