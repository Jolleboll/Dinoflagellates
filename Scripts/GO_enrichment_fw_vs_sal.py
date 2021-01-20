#!/usr/bin/python3

import argparse,sys,math



############ Input and such ############
desc	='''This program calculates ratios between GO term abundances in e.g. the fw and sal unique loci diamond/blast results.
The first file is the nominator.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('i1',metavar='Input file 1',help='The first annotated diamond output file.')
parser.add_argument('i2',metavar='Input file 2',help='The second annotated diamond output file.')
parser.add_argument('i3',metavar='Go file',help='The tab-separated accnr->GO file.')
parser.add_argument('-c',metavar='Cut-off',type=float,help='Include only ratios/genes in this extreme. Default 0.05.',default=0.05)
parser.add_argument('-o1',metavar='Output full',help='Full output.',default='fw_vs_sal_go_full.out')
parser.add_argument('-o2',metavar='Output extr',help='Only top and bottom extremes output.',default='fw_vs_sal_go_extr.out')
parser.add_argument('-o3',metavar='Output excl',help='Only zeroes included.',default='fw_vs_sal_go_excl.out')
args	=parser.parse_args()
########################################



############# Parsing, heh #############
checked_loci=set() # To keep track of whether a locus has been accounted for already
book1,book2=dict(),dict()

t=dict() # To keep track of what transcripts are associated with the accnrs/GO_terms we present in the end

for b,f in zip([book1,book2],[args.i1,args.i2]): # Daym
	with open(f) as f:
		next(f) # Skip!
		for line in f:
			fields=line.split('\t')
			if len(fields)==13: # We disregard hits without annotation
				locus=fields[0] # For exclusion with checked_loci
				trsID=fields[1] # Transcript ID
				fields=fields[-1].split(' ||| ') # Now a list of annotation IDs
				if locus not in checked_loci: # Care only about the first hit/annotation
					accnr=fields[0].split(';')[0].split(':')[1]
					if accnr not in b: b[accnr]=1
					else: b[accnr]+=1
					checked_loci.add(locus)
					if accnr not in t: t[accnr]={trsID}
					else: t[accnr].add(trsID)
########################################
# In other words, multiple distinct loci can have the same top transcript hit,
# but that transcript would be de-duplicated in the transcript book.
# However, if several distinct transcripts have the same top annotation (accnr),
# they will all be associated with this accnr, and then later if a certain GO
# term is associated with several accnrs, the number of transcripts associated with
# that GO term can become quite substantial.
# A few accnrs have a lot of transcripts, up to 155 (a repeat protein with 11 repeated regions!)




########### GO dictionaries ############
acc2go=dict()
with open(args.i3) as f:
	for line in f:
		fields=line.split()
		accnr=fields[0]
		gos=set(fields[1].split(','))
		acc2go[accnr]=gos # Every accnr is now associated with a set of GO terms.

gocount1,gocount2=dict(),dict()
for a,g in zip([book1,book2],[gocount1,gocount2]):
	for accnr in a:
		count=a[accnr]
		if accnr in acc2go:
			for go in acc2go[accnr]:
				if not go in g: g[go]=count
				else: g[go]+=count

				if not go in t: # Remember what transcripts the GO terms came from
					t[go]=t[accnr] # Pray that GO terms are never identical to accession numbers
				else: t[go]|=t[accnr] # Amazing syntax
		else:
			print('Accession number',accnr,'with count',str(count)+', had no associated GO terms.',file=sys.stderr)
########################################
# So, all transcripts whose top annotation was accnr, will be listed as belonging to every GO term associated
# with that same accnr. So, if the final size of t[go] is larger than the number of GO terms amassed,
# the reason must be that the GO belongs to many different accnrs, and probably that several transcripts
# point to the same best accnr. Not sure if that makes biological sense though.

# At this point we have abundances. Now to create ratios...
# But first a silly sorting function because logs are hard
def logsort(tup):
	if		tup[0]=='NA':	return  1000000000+ tup[2]	 /(tup[3]+1)
	elif	tup[0]==0:		return -1000000000+(tup[2]+1)/ tup[3]
	else:					return tup[0]

############# Ratios heh k #############
for gonr in gocount1.keys()|gocount2.keys(): # union of keys to avoid key error later
	if		gonr not in gocount1: gocount1[gonr]=0
	elif	gonr not in gocount2: gocount2[gonr]=0
#	if		gonr not in transcripts1: transcripts1[gonr]={}
#	elif	gonr not in transcripts2: transcripts2[gonr]={}

ratios=[] # empty list fuck ye
for go in gocount1: # Doesn't matter which of the books we loop thru
	if		gocount2[go]==0:	ratios.append(('NA',go,gocount1[go],gocount2[go],len(t[go]),','.join(t[go])))
	else:						ratios.append((gocount1[go]/gocount2[go],go,gocount1[go],gocount2[go],len(t[go]),','.join(t[go])))
ratios.sort(key=logsort) # in-place sorting is a dangerous tool

# Ratios is now a sorted list of tuples for all shared descriptions below the OS taxonomy level.
lower	=math.floor(len(ratios)*args.c)
upper	=len(ratios)-math.floor(len(ratios)*args.c)
Extr	=[tup for i,tup in enumerate(ratios) if (i<lower or i>upper)]

full=[('inf',*tup[1:]) if tup[0]=='NA' else (float('-inf'),*tup[1:]) if tup[0]==0 else (math.log10(tup[0]),*tup[1:]) for tup in ratios]
extr=[('inf',*tup[1:]) if tup[0]=='NA' else (float('-inf'),*tup[1:]) if tup[0]==0 else (math.log10(tup[0]),*tup[1:]) for tup in Extr]
excl=[('inf',*tup[1:]) if tup[0]=='NA' else (float('-inf'),*tup[1:]) for tup in ratios if tup[0] in (0,'NA')]
########################################


################ Output ################
percentage=f'{args.c:3.0%}'
pop1=args.i1.split('.')[-3]
pop2=args.i2.split('.')[-3]
header='# '+'\t'.join(['Log-ratio','GO term',pop1+' counts',pop2+' counts','Number of distinct transcripts','Transcript IDs'])
header1='# Logarithm of GO term "abundance" ratios between {} and {}.'.format(pop1,pop2)
header2='# A ratio of 1.0 for example means 10 times as many were found in the {} blast results.'.format(pop1)
header3='# Step 1: Choose loci respectively unique to {} and {}.'.format(pop1,pop2)
header4='# Step 2: Choose loci that have blast (diamond) hits to the merged transcriptome.'
header5='# Step 3: Choose loci whose transcriptome hits have SwissProt annotations.'
header6='# Step 4: Count the number of times a GO term was ultimately associated with {0} or {1}. Locus -> Best hit -> Best annotation -> One or more GO terms. If all loci had hits, and all hits had annotations, and all annotations had GO terms, the maximum GO count would be the number of loci in either ({0}/{1}) set.'.format(pop1,pop2)
header7='''# Step 5: Divide the GO counts ({}/{}) and log transform where possible.
# The transcript count is the number of distinct transcripts whose best annotation was associated with a particular GO term.
# The reason some transcript counts are much higher than others is probably because the GO term is associated with many different
# proteins, or that the transcriptome for some reason contains many sequences that all point to the same protein (annotation).
# Transcripts present in the transcriptome, but never associated with a locus, are not included.'''.format(pop1,pop2)
with open(args.o1,'w') as o:
	print(header1,header2,header3,header4,header5,header6,header7,header,file=o,sep='\n')
	for ratio in full:
		print(*ratio,sep='\t',file=o)
with open(args.o2,'w') as o:
	header1+=' (only top/bottom '+percentage+' included)'
	print(header1,header2,header3,header4,header5,header6,header7,header,file=o,sep='\n')
	for ratio in extr:
		print(*ratio,sep='\t',file=o)
with open(args.o3,'w') as o:
	header1+=' (only logarithm-incompatible values included)'
	print(header1,header2,header3,header4,header5,header6,header7,header,file=o,sep='\n')
	for ratio in excl:
		print(*ratio,sep='\t',file=o)	
########################################
