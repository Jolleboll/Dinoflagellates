#!/usr/bin/python3

import argparse,sys,math



############ Input and such ############
desc	='''This program calculates ratios between gene abundances in e.g. the fw and sal unique loci diamond/blast results.
The first file is the nominator.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('i1',metavar='Input file 1',help='The first annotated diamond output file.')
parser.add_argument('i2',metavar='Input file 2',help='The second annotated diamond output file.')
parser.add_argument('-s',help='Flag for one locus per description.',action='store_true')
parser.add_argument('-c',metavar='Cut-off',type=float,help='Include only ratios/genes in this extreme. Default 0.05.',default=0.05)
parser.add_argument('-o1',metavar='Output full',help='Full output.',default='fw_vs_sal_full.out')
parser.add_argument('-o2',metavar='Output extr',help='Only top and bottom extremes output.',default='fw_vs_sal_extr.out')
parser.add_argument('-o3',metavar='Output excl',help='Only zeroes included.',default='fw_vs_sal_excl.out')
# parser.add_argument('-hb',metavar='HealingBonus',type=float,help='Bonus healing from talents. Max 1.10.',default=1.0)
# parser.add_argument('-sb',metavar='SpiritBonus',type=float,help='The total healing bonus coeff on Spirit. Max 0.25.',default=0.0)
args	=parser.parse_args()
########################################



############# Parsing, heh #############
checked_loci=set() # To keep track of whether a locus has been accounted for already
book1,book2=dict(),dict()
for b,f in zip([book1,book2],[args.i1,args.i2]): # Daym
	with open(f) as f:
		next(f) # Skip!
		for line in f:
			fields=line.split('\t')
			if len(fields)==13:
				locus=fields[0] # For exclusion with checked_loci
				fields=fields[-1].split(' ||| ')
				if args.s and locus not in checked_loci: # Care only about the first annotation
					descr=fields[0].split(';')[2].split('=')[1].split(' os ')[0]
					if descr not in b: b[descr]=1
					else: b[descr]+=1
					checked_loci.add(locus)
				else:
					for field in fields:
						descr=field.split(';')[2].split('=')[1].split(' os ')[0]
						if descr not in b: b[descr]=1
						else: b[descr]+=1
########################################

# At this point we have abundances. Now to create ratios...
# But first a silly sorting function because logs are hard
def logsort(tup):
	if		tup[0]=='NA':	return  1000000000+ tup[2]	 /(tup[3]+1)
	elif	tup[0]==0:		return -1000000000+(tup[2]+1)/ tup[3]
	else:					return tup[0]

############# Ratios heh k #############
for descr in book1.keys()|book2.keys(): # union of keys
	if 		descr not in book1:	book1[descr]=0
	elif 	descr not in book2:	book2[descr]=0

ratios=[] # empty list fuck ye
for descr in book1: # Doesn't matter which of the books we loop thru
	if		book2[descr]==0:	ratios.append(('NA',descr,book1[descr],book2[descr]))
	else:						ratios.append((book1[descr]/book2[descr],descr,book1[descr],book2[descr]))
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

with open(args.o1,'w') as o:
	print(
		'# Logarithm of abundance ratios between freshwater and salt water.',
		'# One description per locus? '+yndict[args.s]+'.',
		'# A ratio of 1.0 for example means 10 times as many were found in the freshwater blast results.',
		'# Step 1: Choose unique loci for both fw and sal.',
		'# Step 2: Choose loci that had blast hits to the merged transcriptome.',
		'# Step 3: Choose loci whose transcriptome hits have SwissProt annotations.',
		'# Step 4: Count the annotation description fields and divide the former (fw) by the latter (sal) counts.',
		'# Step 5: Log transform if possible.',
	file=o,sep='\n')
	for ratio in full:
		print(*ratio,sep='\t',file=o)
with open(args.o2,'w') as o:
	print(
		'# Logarithm of abundance ratios between freshwater and salt water (only top/bottom '+percentage+' included)',
		'# One description per locus? '+yndict[args.s]+'.',
		'# A ratio of 1.0 for example means 10 times as many were found in the freshwater blast results.',
		'# Step 1: Choose unique loci for both fw and sal.',
		'# Step 2: Choose loci that had blast hits to the merged transcriptome.',
		'# Step 3: Choose loci whose transcriptome hits have SwissProt annotations.',
		'# Step 4: Count the annotation description fields and divide the former (fw) by the latter (sal) counts.',
		'# Step 5: Log transform if possible.',
	file=o,sep='\n')
	for ratio in extr:
		print(*ratio,sep='\t',file=o)
with open(args.o3,'w') as o:
	print(
		'# Logarithm of abundance ratios between freshwater and salt water. Only logarithm-incompatible values included.',
		'# One description per locus? '+yndict[args.s]+'.',
		'# A ratio of 1.0 for example means 10 times as many were found in the freshwater blast results.',
		'# Step 1: Choose unique loci for both fw and sal.',
		'# Step 2: Choose loci that had blast hits to the merged transcriptome.',
		'# Step 3: Choose loci whose transcriptome hits have SwissProt annotations.',
		'# Step 4: Count the annotation description fields and divide the former (fw) by the latter (sal) counts.',
		'# Step 5: Log transform if possible.',
			file=o,sep='\n')
	for ratio in excl:
		print(*ratio,sep='\t',file=o)	
########################################

#[(book1[descr]/book2[descr],descr) for descr in book1.keys()&book2.keys()],key=lambda x:x[0])