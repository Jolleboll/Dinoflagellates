#!/usr/bin/python3

import argparse,sys,math



############ Input and such ############
desc	='''This program calculates ratios between GO term abundances in e.g. the fw and sal unique loci diamond/blast results.
The first file is the nominator.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('i1',metavar='Input file 1',help='The first annotated diamond output file.')
parser.add_argument('i2',metavar='Input file 2',help='The second annotated diamond output file.')
parser.add_argument('i3',metavar='Go file',help='The tab-separated accnr->GO file.')
parser.add_argument('-s',help='Flag for one GO per locus.',action='store_true')
parser.add_argument('-c',metavar='Cut-off',type=float,help='Include only ratios/genes in this extreme. Default 0.05.',default=0.05)
parser.add_argument('-o1',metavar='Output full',help='Full output.',default='fw_vs_sal_go_full.out')
parser.add_argument('-o2',metavar='Output extr',help='Only top and bottom extremes output.',default='fw_vs_sal_go_extr.out')
parser.add_argument('-o3',metavar='Output excl',help='Only zeroes included.',default='fw_vs_sal_go_excl.out')
# parser.add_argument('-hb',metavar='HealingBonus',type=float,help='Bonus healing from talents. Max 1.10.',default=1.0)
# parser.add_argument('-sb',metavar='SpiritBonus',type=float,help='The total healing bonus coeff on Spirit. Max 0.25.',default=0.0)
args	=parser.parse_args()
########################################



############# Parsing, heh #############
checked_loci=set() # To keep track of whether a locus has been accounted for already
book1,book2=dict(),dict()

transcripts1=dict() # To keep track of what transcripts are associated with the accnrs/GO_terms we present in the end
transcripts2=dict()

for b,f,t in zip([book1,book2],[args.i1,args.i2],[transcripts1,transcripts2]): # Daym
	with open(f) as f:
		next(f) # Skip!
		for line in f:
			fields=line.split('\t')
			if len(fields)==13:
				locus=fields[0] # For exclusion with checked_loci
				trsID=fields[1] # Transcript ID
				fields=fields[-1].split(' ||| ')
				if args.s and locus not in checked_loci: # Care only about the first annotation
					accnr=fields[0].split(';')[0].split(':')[1]
					if accnr not in b: b[accnr]=1
					else: b[accnr]+=1
					checked_loci.add(locus)
					if accnr not in t: t[accnr]={trsID}
					else: t[accnr].add(trsID)
				elif not args.s:
					for field in fields:
						accnr=field.split(';')[0].split(':')[1]
						if accnr not in b: b[accnr]=1
						else: b[accnr]+=1
						if accnr not in t: t[accnr]={trsID}
						else: t[accnr].add(trsID)
########################################



########### GO dictionaries ############
acc2go=dict()
with open(args.i3) as f:
	for line in f:
		fields=line.split()
		accnr=fields[0]
		gos=set(fields[1].split(','))
		acc2go[accnr]=gos

gocount1,gocount2=dict(),dict()
for a,g,t in zip([book1,book2],[gocount1,gocount2],[transcripts1,transcripts2]):
	for accnr in a:
		count=a[accnr]
		if accnr in acc2go:
			for go in acc2go[accnr]:
				if not go in g: g[go]=count
				else: g[go]+=count

				if not go in t: # Remember what transcripts the GO terms came from
					t[go]=t[accnr]
				else: t[go]|=t[accnr]
		else:
			print('Accession number',accnr,'with count',str(count)+', had no associated GO terms.',file=sys.stderr)
########################################


# At this point we have abundances. Now to create ratios...
# But first a silly sorting function because logs are hard
def logsort(tup):
	if		tup[0]=='NA':	return  1000000000+ tup[2]	 /(tup[3]+1)
	elif	tup[0]==0:		return -1000000000+(tup[2]+1)/ tup[3]
	else:					return tup[0]

############# Ratios heh k #############
for gonr in gocount1.keys()|gocount2.keys(): # union of keys
	if		gonr not in gocount1: gocount1[gonr]=0
	elif	gonr not in gocount2: gocount2[gonr]=0
	if		gonr not in transcripts1: transcripts1[gonr]={}
	elif	gonr not in transcripts2: transcripts2[gonr]={}

ratios=[] # empty list fuck ye
for go in gocount1: # Doesn't matter which of the books we loop thru
	if		gocount2[go]==0:	ratios.append(('NA',go,gocount1[go],gocount2[go],len(t[go]),','.join(t[go])))
	else:						ratios.append((gocount1[go]/gocount2[go],go,gocount1[go],gocount2[go],len(t[go]),','.join(t[go])))
ratios.sort(key=logsort)

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
yndict={True:'Yes',False:'No'}
header='# '+'\t'.join(['Log-ratio','GO term','Freshwater counts','Saltwater counts','Number of distinct transcripts','Transcript IDs'])

with open(args.o1,'w') as o:
	print(
		'# Logarithm of GO term "abundance" ratios between freshwater and salt water.',
		'# One GO term per locus? '+yndict[args.s]+'.',
		'# A ratio of 1.0 for example means 10 times as many were found in the freshwater blast results.',
		'# Step 1: Choose unique loci for both fw and sal.',
		'# Step 2: Choose loci that had blast hits to the merged transcriptome.',
		'# Step 3: Choose loci whose transcriptome hits have SwissProt annotations.',
		'# Step 4: Count all GO terms associated with the top hit\'s accession number, then divide the former count (fw) by the latter (sal).',
		'# Step 5: Log transform if possible.',header,
	file=o,sep='\n')
	for ratio in full:
		print(*ratio,sep='\t',file=o)
with open(args.o2,'w') as o:
	print(
		'# Logarithm of GO term "abundance" ratios between freshwater and salt water (only top/bottom '+percentage+' included)',
		'# One GO term per locus? '+yndict[args.s]+'.',
		'# A ratio of 1.0 for example means 10 times as many were found in the freshwater blast results.',
		'# Step 1: Choose unique loci for both fw and sal.',
		'# Step 2: Choose loci that had blast hits to the merged transcriptome.',
		'# Step 3: Choose loci whose transcriptome hits have SwissProt annotations.',
		'# Step 4: Count all GO terms associated with the top hit\'s accession number, then divide the former count (fw) by the latter (sal).',
		'# Step 5: Log transform if possible.',header,
	file=o,sep='\n')
	for ratio in extr:
		print(*ratio,sep='\t',file=o)
with open(args.o3,'w') as o:
	print(
		'# Logarithm of GO term "abundance" ratios between freshwater and salt water. Only logarithm-incompatible values included.',
		'# One GO term per locus? '+yndict[args.s]+'.',
		'# A ratio of 1.0 for example means 10 times as many were found in the freshwater blast results.',
		'# Step 1: Choose unique loci for both fw and sal.',
		'# Step 2: Choose loci that had blast hits to the merged transcriptome.',
		'# Step 3: Choose loci whose transcriptome hits have SwissProt annotations.',
		'# Step 4: Count all GO terms associated with the top hit\'s accession number, then divide the former count (fw) by the latter (sal).',
		'# Step 5: Log transform if possible.',header,
			file=o,sep='\n')
	for ratio in excl:
		print(*ratio,sep='\t',file=o)	
########################################