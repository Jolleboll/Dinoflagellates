#!/usr/bin/python3

import argparse,re,sys,gzip,pickle,os

############ Input and such ############
desc	='''
This program takes a popmap file and returns information about how many
loci are shared between populations. It also prints the unique loci to files.
With the -D flag, you only create the sstacks loci dictionary.
It has sample names as keys and sets of loci as values. Its name is sstacks.pkl
and must exist in this script's directory.
Also needed when using -D is the following file at this location:
/crex/proj/uppstore2017173/b2014326/Dinopop2020_all/samples2020.txt
It's just a list of sample names, newline-separated.
It should include all sample names that any popmap includes.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('p',metavar='popmap',help='Path to your popmap file.')
parser.add_argument('s',metavar='sstacks',help='Path to your input sstacks directory.')
parser.add_argument('-D',help='Only create sstacks dictionary.',action='store_true')
args	=parser.parse_args()
########################################



######## Parsing sstacks files #########
if args.D:
	loci=dict() # sets of loci per sample
	with open('/crex/proj/uppstore2017173/b2014326/Dinopop2020_all/samples2020.txt') as f:
		allsamples=[line.rstrip() for line in f]
	for sample in allsamples:
		with gzip.open(args.s.rstrip('/')+'/'+sample+'.matches.tsv.gz','rt') as f:
			next(f) # stacks header
			loci[sample]=set()
			for line in f:
				locus=line.split()[2]
				loci[sample].add(locus)
	with open(sys.path[0]+'/sstacks.pkl','wb') as f:
		pickle.dump(loci,f)
	print('Done dumping sstacks dictionary. Exiting.')
	sys.exit()
########################################

########## Loading sstacks #############
with open(sys.path[0]+'/sstacks.pkl','rb') as f:
	loci=pickle.load(f)
print('\nDone loading sstacks dictionary.\n')
########################################



########## Parsing popmap ##############
pops=dict()
with open(args.p) as f:
	for line in f:
		fields=line.rstrip().split('\t')
		if fields[1] not in pops: pops[fields[1]]={fields[0]}
		else: pops[fields[1]].add(fields[0])
print('\nDone parsing popmap file.\n')
########################################



######### Creating popmap sets #########
popLoci=dict() # pop:set(loci)
for pop in pops:
	popLoci[pop]=set.union(*(loci[sample] for sample in pops[pop]))
########################################



########### Set calculations ###########
popIntersec=set.intersection(*popLoci.values())
# The following is just to check how many loci are present in ALL samples, but the answer is 0 so meh.
# allIntersec=set.intersection(*[loci[sample] for sample in [sample for sublist in pops.values() for sample in sublist]]) # Readability is why we use Python
########################################


########## Printing results ############
popmapname=re.split('[_.]',args.p)[-2]
print('Popmap name:',popmapname)
print('Common loci across all its populations:',len(popIntersec))
print('\nUnique loci for each population...')
blargh=list(popLoci.items())
for i in range(len(blargh)):
	with open(popmapname+'.'+blargh[i][0]+'.uniqueLoci','w') as out:
		smol=[blargh[j][1] for j in range(len(blargh)) if j!=i]
		print('Loci unique to ',blargh[i][0],': ',len(blargh[i][1].difference(*smol)),sep='')
		out.write('\n'.join(blargh[i][1].difference(*smol))+'\n')
# print('Common loci across all samples:',len(allIntersec))
print('\nPairwise intersections...')
for i in range(len(blargh)-1):
	for j in range(i+1,len(blargh)):
		print('Common loci in ',blargh[i][0],' and ',blargh[j][0],': ',len(set.intersection(blargh[i][1],blargh[j][1])),sep='')
########################################