#!/usr/bin/python3

import argparse,re,gzip,sys

############ Input and such ############
desc	='''This program takes various stacks-related inputs and converts to BayeScan SNP input.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('p',metavar='popmap',help='Path to your popmap input file that determines what populations you want to consider.')
parser.add_argument('a',metavar='alleles',help='Path to your *catalog.alleles.tsv.gz* (cstacks output) file.')
parser.add_argument('o',metavar='outfile',help='Name of your BayeScan format output file.')
args	=parser.parse_args()
########################################



############ Populations ###############
pops=dict() # pop to samples
popL=list() # pops in a list, index is their number (use enumerate(1,...) so it starts on 1)
with open(args.p) as f:
	for line in f:
		fields	=line.rstrip().split('\t')
		sample	=fields[0]
		pop		=fields[1]
		if pop not in pops:
			pops[pop]=[sample]
			popL.append(pop)
		else:
			pops[pop].append(sample)
print('Populations:',', '.join(popL),file=sys.stderr) # For fun
########################################

############## Loci IDs ################
with open('populations_output/filtered_loci.tsv') as f:
	lociIDs_list=[line.split('\t')[1].rstrip() for line in f]
	lociIDs_set	=set(lociIDs_list)
########################################

############### Alleles ################
alleles=dict()
with gzip.open(args.a,'rt') as f:
	next(f) # ignore header
	for line in f:
		fields	=line.rstrip().split('\t')
		locus	=fields[2]
		allele	=fields[3]
		if locus in lociIDs_set: # ignore loci we're not interested in
			if locus not in alleles:
				alleles[locus]=[allele]
			else:
				alleles[locus].append(allele)
########################################

########## Allele freqs ################
allelefreqs=dict()
# with gzip.open(args.a,'rt') as f:
# 	next(f) # ignore header
# 	for line in f:
# 		fields	=line.rstrip().split('\t')
# 		locus	=fields[2]
# 		allele	=fields[3]
# 		if locus in lociIDs_set: # ignore loci we're not interested in
# 			if locus not in alleles:
# 				alleles[locus]=[allele]
# 			else:
# 				alleles[locus].append(allele)
########################################

########### N "individuals" ############
Ninds=dict() # Given a pop and a locus, how many samples/genes/alleles (haploid...) are represented?
for pop in pops:
	print('Processing population',pop,file=sys.stderr) # For fun
	for sample in pops[pop]:
		processed=set() # To keep track of what loci have been seen already... ugly I know.
		with gzip.open('sstacks_output/testlinks/'+sample+'.matches.tsv.gz','rt') as f:
			next(f) # Skip header
			for line in f:
				locus=line.split('\t')[2]
				if locus in lociIDs_set: # We care only about the loci in the populations program output
					if (pop,locus) not in Ninds:
						Ninds[(pop,locus)]=1
					elif locus not in processed:
						Ninds[(pop,locus)]+=1
					processed.add(locus)
########################################

sys.exit()

############## Printing ################
with open(args.o,'w') as f:
	print('[loci]='+str(len(lociIDs_set))+'\n') # extra newline because
	print('[populations]='+str(len(pops)))
	for ipop,pop in enumerate(popL,1):
		print('\n'+'[pop]='+str(ipop))
		for ilocus,locus in enumerate(lociIDs_list,1):
			print(ilocus,Ninds[(pop,locus)],len(alleles[locus]),*allelefreqs[(pop,locus)])
########################################

'''
### input file start
[loci]=1000

[populations]=10

[pop]=1
	1	40	2	0	40
	2	40	2	34	6
...
[pop]=2
... 
### input file end

#loci id, N individuals, N alleles, N of each allele

loci id, is easily accessible from the output of get_interesting_loci.scr.

N individuals, is now located in the Ninds dictionary's values. Easy game. See the BayeScan manual for info
about this number though.

N alleles, is once the populations program has filtered out unwanted loci, easily accessible with the
bash script in /crex/proj/uppstore2017173/b2014326/Dinopop2020_all/bayeScan. Take this script's output and import
only the wanted loci (two input files -> dictionary in this python script).

N of each allele, is difficult.
'''
