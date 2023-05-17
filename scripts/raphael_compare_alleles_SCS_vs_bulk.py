#!/usr/bin/python3

import argparse,re,sys,gzip,glob

############ Input and such ############
desc	='''
This program compares all the catalog-matched loci for all Raphael\'s SCS samples against
the bulk sample catalog, and looks at what alleles are not present.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('c',metavar='catalog dir',help='Path to catalog\'s directory.')
parser.add_argument('s',metavar='sstacks output dir',help='Path to SCS samples\' sstacks output directory.')
parser.add_argument('n',metavar='stacks sample list',help='Path to file containing all "SCS" samples\' filenames')
parser.add_argument('o',metavar='output file',help='Path to output file')
args	=parser.parse_args()
########################################

'''

Columns in .matches. files:
catalog ID	sample ID (not useful)	Locus ID	Haplotype	Stack depth	Cigar string

Columns in .snps. files:
sample ID (not useful)	Locus ID	Position	Type	Likelihood	Majority nuc	Alternative nuc	...	...


Questions:
2. Are positions aligned? They should be, considering the same restriction enzyme was used.
3. What about snps that occur in the SCS but not in the bulk?
'''


########## Parsing bulk files ##########

# .tags. first
print('Counting catalog loci...',file=sys.stderr)
with gzip.open(args.c+'/'+'catalog.tags.tsv.gz','rt') as f:
	for line in f:
		n_cat_loci=line.split()[1]
n_cat_loci=int(n_cat_loci)

# .alleles. secnd
print('Collecting catalog HZ loci...',file=sys.stderr)
bulk_HZ_loci=dict()
with gzip.open(args.c+'/'+'catalog.alleles.tsv.gz','rt') as f:
	next(f)
	for line in f:
		locus=line.split()[1]
		if locus not in bulk_HZ_loci: bulk_HZ_loci[locus]=1
		else: bulk_HZ_loci[locus]+=1
print('Done.',file=sys.stderr)


# .snps. turd
print('Collecting catalog snps...',file=sys.stderr)
bulk_snps,nbulk_snps=dict(),0
with gzip.open(args.c+'/'+'catalog.snps.tsv.gz','rt') as f:
	next(f)
	for line in f:
		fields	=line.split()
		if fields[3]!='E': continue # Only care about heterozygous ones
		locus	=fields[1]
		position=fields[2]

		nuc1	=fields[5]
		nuc2	=fields[6]
		if not fields[7]=='-':
			nuc3=fields[7]
			if not fields[8]=='-':
				nuc4=fields[8]
			else: nuc4=''
		else:
			nuc3,nuc4='',''
		nucs={nuc for nuc in {nuc1,nuc2,nuc3,nuc4} if nuc}

		if not locus in bulk_snps: bulk_snps[locus]=[(position,nucs)]
		else: bulk_snps[locus].append((position,nucs))
		nbulk_snps+=1
print('Done.',file=sys.stderr)

print('Calculating bulk DP loci...',file=sys.stderr)
n_bulk_HZ_loci	=len(bulk_HZ_loci) # Measures for later use
n_DP_bulk_loci	=sum(bulk_HZ_loci[x]==2 for x in bulk_HZ_loci) # True + True = 2
DP_bulk_loci	={x for x in bulk_HZ_loci if bulk_HZ_loci[x]==2}
print('Done.',file=sys.stderr)
########################################

########## List of paths ###############
with open(args.n) as f:
	samples=[]
	for line in f:
		samples.append(line.rstrip())
########################################

############## Headers #################
with open(args.o,'w') as o:
	print('Printing headers...',file=sys.stderr)
	print('##    "SCS"  =   Single Cell Sample',file=o)
	print('##    "HZ"   =   Heterozygous. A locus is defined HZ if present in the .alleles. file.',file=o)
	print('##    "DP"   =   Diploid. A locus is defined DP if it is HZ and does not have more than two alleles in the .alleles. file.',file=o)
	print('##    "snp"  =   A snp is defined as a single HZ entry in a .snps. file.',file=o)
	print('##    "#"    =   Number of',file=o)
	print('## NOTE that DP_Mrat2 is a lower estimate, since the SCS -> Bulk mapping is unique, but Bulk -> SCS is not.',file=o)
	print('##',file=o)
	print('##',file=o)
	print('## SCS_loci  =   (# SCS loci)',file=o)
	print('##',file=o)
	print('##',file=o)
	print('## HO_vs_HZ  =   (# non-HZ SCS loci that matched HZ catalog loci)',file=o)
	print('## HO_vs_HZ1 =   (# non-HZ SCS loci that matched HZ catalog loci) / (# matched loci)',file=o)
	print('## HO_vs_HZ2 =   (# non-HZ SCS loci that matched HZ catalog loci) / (# HZ bulk loci)',file=o)
	print('## HO_vs_HZ3 =   (# non-HZ SCS loci that matched HZ catalog loci) / (# HZ SCS loci that matched HZ catalog loci)',file=o)
	print('## HO_vs_HZ4 =   (# non-HZ SCS loci that matched HZ catalog loci) / (# HO SCS loci that matched HO catalog loci)',file=o)
	print('##',file=o)
	print('##',file=o)
	print('## SCS_HZ_nomatch    =   (# HZ SCS loci that did not match any catalog locus)',file=o)
	print('## SCS_HZ_nomatch_rat=   (# HZ SCS loci that did not match any catalog locus) / (# HZ SCS loci that matched catalog loci)',file=o)
	print('## SCS_HO_nomatch    =   (# non-HZ SCS loci that did not match any catalog locus)',file=o)
	print('## SCS_HO_nomatch_rat=   (# non-HZ SCS loci that did not match any catalog locus) / (# non-HZ SCS loci that matched catalog loci)',file=o)
	print('##',file=o)
	print('##',file=o)
	print('## SCS_HZ    =   (# HZ SCS loci)',file=o)
	print('## HZ_rat1   =   (# HZ SCS loci)                          / (# SCS loci)',file=o)
	print('## HZ_rat2   =   (# HZ SCS loci)                          / (# HZ bulk loci)',file=o)
	print('## HZ_Mrat1  =   (# HZ SCS loci matching a HZ bulk locus) / (# HZ SCS loci)',file=o)
	print('## HZ_Mrat2  =   (# HZ SCS loci matching a HZ bulk locus) / (# HZ bulk loci)',file=o)
	print('##',file=o)
	print('##',file=o)
	print('## SCS_snp   =   (# SCS snps)',file=o)
	print('## snp_rat   =   (# SCS snps)                     / (# Bulk snps)',file=o)
	print('## snp_Mrat1 =   (# SCS snps matching a bulk snp) / (# SCS snps)',file=o)
	print('## snp_Mrat2 =   (# SCS snps matching a bulk snp) / (# bulk snps)',file=o)
	print('##',file=o)
	print('##',file=o)
	print('## SCS_DP    =   (# DP SCS loci )',file=o)
	print('## DP_rat1   =   (# DP SCS loci )                            / (# SCS loci)',file=o)
	print('## DP_rat2   =   (# DP SCS loci )                            / (# HZ SCS loci)',file=o)
	print('## DP_rat3   =   (# DP SCS loci )                            / (# DP bulk loci)',file=o)
	print('## DP_Mrat1  =   (# DP SCS loci matching a DP bulk locus)    / (# DP SCS loci)',file=o)
	print('## DP_Mrat2  =   (# DP bulk loci "matching" a DP SCS locus)  / (# DP bulk loci)',file=o)
	print('## HO_vs_DP  =   (# non-HZ SCS loci that matched DP catalog loci)',file=o)
	print('## HO_vs_DP1 =   (# non-HZ SCS loci that matched DP catalog loci) / (# matched loci)',file=o)
	print('## HO_vs_DP2 =   (# non-HZ SCS loci that matched DP catalog loci) / (# DP bulk loci)',file=o)
	print('##',file=o)
	print('##',file=o)
	print('## The number of catalog loci is   ',f'{n_cat_loci:,}',file=o)
	print('## The number of HZ catalog loci is',f'{n_bulk_HZ_loci:,}',file=o)
	print('## The number of DP catalog loci is',f'{n_DP_bulk_loci:,}',file=o)
	print('## The number of catalog snps is   ',f'{nbulk_snps:,}',file=o)
	print('##',file=o)
	print('##',file=o)
	print('# Sample','SCS_loci','Matches',
			'HO_vs_HZ','HO_vs_HZ1','HO_vs_HZ2','HO_vs_HZ3','HO_vs_HZ4',
			'SCS_HZ_nomatch','SCS_HZ_nomatch_rat','SCS_HO_nomatch','SCS_HO_nomatch_rat',
			'SCS_HZ','HZ_rat1','HZ_rat2','HZ_Mrat1','HZ_Mrat2',
			'SCS_snp','snp_rat','snp_Mrat1','snp_Mrat2',
			'SCS_DP','DP_rat1','DP_rat2','DP_rat3','DP_Mrat1','DP_Mrat2','HO_vs_DP','HO_vs_DP1','HO_vs_DP2',
			sep='\t',file=o)
print('Done.',file=sys.stderr)
########################################

########## Parsing SCS samples #########
for sample in samples: # Main loop

	basename=sample.split('/')[-1]
	print('Parsing sample',basename+'...',file=sys.stderr)

	print('\tParsing tags...',file=sys.stderr)
	with gzip.open(sample+'.tags.tsv.gz','rt') as f:
		next(f)
		SCS_loci=set()
		for line in f:
			SCS_loci.add(line.split()[1])
	n_SCS_loci=len(SCS_loci)

	print('\tParsing matches...',file=sys.stderr)
	with gzip.open(args.s+'/'+basename+'.matches.tsv.gz','rt') as f:
		sample2bulk=dict()
		next(f)
		for line in f:
			fields		=line.split()
			bulklocus	=fields[0]
			samplelocus	=fields[2]
			sample2bulk[samplelocus]=bulklocus # This way because the other way does not quite feature unique keys (naturally)

	bulk2sample={v:k for k,v in sample2bulk.items()}


	print('\tParsing snps...',file=sys.stderr)
	scs_snps,n_scs_snps=dict(),0	
	with gzip.open(sample+'.snps.tsv.gz','rt') as f:
		next(f)
		for line in f:
			fields	=line.split()
			if fields[3]!='E': continue # Only care about heterozygous ones
			locus	=fields[1]
			position=fields[2]
			nuc1	=fields[5]
			nuc2	=fields[6]

			if not locus in scs_snps: scs_snps[locus]=[(position,{nuc1,nuc2})]
			else: scs_snps[locus].append((position,{nuc1,nuc2})) # SCS snps check!
			n_scs_snps+=1
	

	print('\tParsing alleles...',file=sys.stderr)
	scs_HZ_loci=dict()
	with gzip.open(sample+'.alleles.tsv.gz','rt') as f:
		next(f)
		for line in f:
			locus=line.split()[1]
			if locus not in scs_HZ_loci: scs_HZ_loci[locus]=1
			else: scs_HZ_loci[locus]+=1

	n_DP_SCS_loci	=sum(scs_HZ_loci[x]==2 for x in scs_HZ_loci) # x is a locus
	DP_SCS_loci		={x for x in scs_HZ_loci if scs_HZ_loci[x]==2}


	

	# Calculations
	print('\tCalculating the various measures...',file=sys.stderr)

	SCS_HO_vs_HZ		=sum((locus not in scs_HZ_loci and sample2bulk[locus] in bulk_HZ_loci) for locus in sample2bulk)
	try:	SCS_HO_vs_HZ_rat1	=SCS_HO_vs_HZ/len(sample2bulk)
	except ZeroDivisionError:	SCS_HO_vs_HZ_rat1='na'
	SCS_HO_vs_HZ_rat2	=SCS_HO_vs_HZ/n_bulk_HZ_loci

	SCS_HO				=SCS_loci-scs_HZ_loci.keys()
	SCS_HZ_vs_HZ		=sum(sample2bulk[locus] in bulk_HZ_loci for locus in scs_HZ_loci.keys()&sample2bulk.keys())
	SCS_HO_vs_HO		=sum(sample2bulk[locus] not in bulk_HZ_loci for locus in SCS_HO&sample2bulk.keys())
	try:	SCS_HO_vs_HZ_rat3	=SCS_HO_vs_HZ/SCS_HZ_vs_HZ
	except ZeroDivisionError:	SCS_HO_vs_HZ_rat3='na'
	try:	SCS_HO_vs_HZ_rat4	=SCS_HO_vs_HZ/SCS_HO_vs_HO
	except ZeroDivisionError:	SCS_HO_vs_HZ_rat4='na'

	SCS_HZ_nomatch		=len(scs_HZ_loci.keys()-sample2bulk.keys())
	try:	SCS_HZ_nomatch_rat	=SCS_HZ_nomatch/len(scs_HZ_loci.keys()&sample2bulk.keys())
	except ZeroDivisionError:	SCS_HZ_nomatch_rat='na'
	SCS_HO_nomatch		=len((SCS_loci-scs_HZ_loci.keys())-sample2bulk.keys())
	try:	SCS_HO_nomatch_rat	=SCS_HO_nomatch/len((SCS_loci-scs_HZ_loci.keys())&sample2bulk.keys())
	except ZeroDivisionError:	SCS_HO_nomatch_rat='na'

	n_SCS_HZ_loci		=len(scs_HZ_loci)
	try:	HZ_rat1		=n_SCS_HZ_loci/n_SCS_loci
	except ZeroDivisionError:	HZ_rat1	='na'
	try:	HZ_rat2		=n_SCS_HZ_loci/n_bulk_HZ_loci
	except ZeroDivisionError: HZ_rat2	='na'
	try:	HZ_Mrat1	=SCS_HZ_vs_HZ / n_SCS_HZ_loci # etc
	except ZeroDivisionError: HZ_Mrat1	='na'
	try:	HZ_Mrat2	=SCS_HZ_vs_HZ / n_bulk_HZ_loci
	except ZeroDivisionError: HZ_Mrat2	='na'

	SCS_snps			=n_scs_snps
	snp_rat				=SCS_snps/nbulk_snps
	tmp					=sum( (snp1[0]==snp2[0] and snp1[1]<=snp2[1]) for locus in sample2bulk.keys()&scs_snps.keys() for snp1 in scs_snps[locus] for snp2 in bulk_snps[sample2bulk[locus]])
	try:	snp_Mrat1	=tmp/SCS_snps
	except ZeroDivisionError: snp_Mrat1	='na'
	try:	snp_Mrat2	=tmp/nbulk_snps
	except ZeroDivisionError: snp_Mrat2	='na'

	SCS_DP				=n_DP_SCS_loci
	try:	DP_rat1		=n_DP_SCS_loci/n_SCS_loci
	except ZeroDivisionError:	DP_rat1	='na'
	try:	DP_rat2		=n_DP_SCS_loci/n_SCS_HZ_loci
	except ZeroDivisionError: DP_rat2	='na'
	try:	DP_rat3		=n_DP_SCS_loci/n_DP_bulk_loci
	except ZeroDivisionError:	DP_rat3	='na'
	try:	DP_Mrat1	=sum(sample2bulk[locus] in DP_bulk_loci for locus in DP_SCS_loci&sample2bulk.keys()) / n_DP_SCS_loci
	except ZeroDivisionError:	DP_Mrat1='na'
	try:	DP_Mrat2	=sum(bulk2sample[locus] in DP_SCS_loci for locus in DP_bulk_loci&bulk2sample.keys()) / n_DP_bulk_loci
	except ZeroDivisionError:	DP_Mrat2='na'

	SCS_HO_vs_DP		=sum((locus not in scs_HZ_loci and sample2bulk[locus] in DP_bulk_loci) for locus in sample2bulk)
	try:	SCS_HO_vs_DP_rat1	=SCS_HO_vs_DP/len(sample2bulk)
	except ZeroDivisionError:	SCS_HO_vs_DP_rat1='na'
	SCS_HO_vs_DP_rat2	=SCS_HO_vs_DP/n_DP_bulk_loci

	with open(args.o,'a') as o:
		print(	basename,n_SCS_loci,len(sample2bulk),
				SCS_HO_vs_HZ,SCS_HO_vs_HZ_rat1,SCS_HO_vs_HZ_rat2,SCS_HO_vs_HZ_rat3,SCS_HO_vs_HZ_rat4,
				SCS_HZ_nomatch,SCS_HZ_nomatch_rat,SCS_HO_nomatch,SCS_HO_nomatch_rat,
				n_SCS_HZ_loci,HZ_rat1,HZ_rat2,HZ_Mrat1,HZ_Mrat2,
				SCS_snps,snp_rat,snp_Mrat1,snp_Mrat2,
				SCS_DP,DP_rat1,DP_rat2,DP_rat3,DP_Mrat1,DP_Mrat2,SCS_HO_vs_DP,SCS_HO_vs_DP_rat1,SCS_HO_vs_DP_rat2,
				sep='\t',file=o)
	print('Done with',basename+'.\n',file=sys.stderr)

########################################