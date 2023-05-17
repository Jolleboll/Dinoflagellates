#!/usr/bin/python3

import argparse,re,sys,gzip,glob

############ Input and such ############
desc	='''
This program compares all the catalog-matched loci for all Raphael\'s SCS samples against
the bulk sample catalog, and looks at what alleles are not present.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('n',metavar='stacks file list',help='Path to file containing all stacks files\' filenames')
parser.add_argument('f',metavar='mlem',help='Hurrrr')
parser.add_argument('ba',metavar='bulk alleles',help='Path to bulk sample\' alleles file.')
parser.add_argument('bs',metavar='bulk snps',help='Path to bulk sample\' snps file.')
args	=parser.parse_args()
########################################

########## List of paths ###############
with open(args.n) as f:
	samples=[]
	for line in f:
		if 'snps' in line and 'catalog' not in line and 'std' not in line: # Choose snps files that are not bulk or catalog.
			samples.append('.'.join(line.split('.')[:-3])) # without extensions
########################################

########## Parsing bulk files ##########

# .alleles. first
print('Collecting bulk HZ loci...',file=sys.stderr)
bulk_HZ_loci=dict()
with gzip.open(args.ba,'rt') as f:
	next(f)
	for line in f:
		locus=line.split()[1]
		if locus not in bulk_HZ_loci: bulk_HZ_loci[locus]=1
		else: bulk_HZ_loci[locus]+=1
print('Done.',file=sys.stderr)


# .snps. second
print('Collecting bulk snps...',file=sys.stderr)
bulk_snps,nbulk_snps=dict(),0
with gzip.open(args.bs,'rt') as f:
	next(f)
	for line in f:
		fields	=line.split()
		if fields[3]!='E': continue # Only care about heterozygous ones
		locus	=fields[1]
		position=fields[2]
		nuc1	=fields[5]
		nuc2	=fields[6]
		
		# We ignore third and fourth snps, hmm... perhaps we shouldn't?
		if not locus in bulk_snps: bulk_snps[locus]=[(position,{nuc1,nuc2})]
		else: bulk_snps[locus].append((position,{nuc1,nuc2}))
		nbulk_snps+=1
print('Done.',file=sys.stderr)

print('Calculating bulk DP loci...',file=sys.stderr)
HZ_bulk_loci	=len(bulk_HZ_loci) # Measures for later use
n_DP_bulk_loci	=sum(bulk_HZ_loci[x]==2 for x in bulk_HZ_loci) # True + True = 2
DP_bulk_loci	={x for x in bulk_HZ_loci if bulk_HZ_loci[x]==2}
print('Done.',file=sys.stderr)
########################################

############## Headers #################
print('Printing headers...',file=sys.stderr)
print('##	"SCS"	=	Single Cell Sample')
print('## 	"HZ"	=	Heterozygous. A locus is defined HZ if present in the .alleles. file.')
print('## 	"DP"	=	Diploid. A locus is defined DP if it is HZ and does not have more than two alleles in the .alleles. file.')
print('##	"snp"	=	A snp is defined as a single HZ entry in a .snps. file.')
print('## 	"#"	=	Number of')
print('## NOTE that DP_Mrat2 is a lower estimate, since the SCS -> Bulk mapping is unique, but Bulk -> SCS is not.')
print('##')
print('##')
print('## SCS_HZ	=	(# HZ SCS loci)') # .alleles.
print('## HZ_rat1	=	(# HZ SCS loci)                          / (# SCS loci)') # .alleles.
print('## HZ_rat2	=	(# HZ SCS loci)                          / (# HZ bulk loci)') # .alleles.
print('## HZ_Mrat1	=	(# HZ SCS loci matching a HZ bulk locus) / (# HZ SCS loci)') # .alleles. .matches.
print('## HZ_Mrat2	=	(# HZ SCS loci matching a HZ bulk locus) / (# HZ bulk loci)') # .alleles. .matches.
print('##')
print('##')
print('## SCS_snp	=	(# SCS snps)') # .snps.
print('## snp_rat	=	(# SCS snps)                     / (# Bulk snps)') # .snps.
print('## snp_Mrat1	=	(# SCS snps matching a bulk snp) / (# SCS snps)') # .alleles. .matches.
print('## snp_Mrat2	=	(# SCS snps matching a bulk snp) / (# bulk snps)') # .alleles. .matches.
print('##')
print('##')
print('## SCS_DP	=	(# DP SCS loci )') # .alleles.
print('## DP_rat1	=	(# DP SCS loci )                            / (# SCS loci)') # .alleles.
print('## DP_rat2	=	(# DP SCS loci )                            / (# HZ SCS loci)') # .alleles.
print('## DP_rat3	=	(# DP SCS loci )                            / (# DP bulk loci)') # .alleles.
print('## DP_Mrat1	=	(# DP SCS loci matching a DP bulk locus)    / (# DP SCS loci)') # .alleles. .matches.
print('## DP_Mrat2	=	(# DP bulk loci "matching" a DP SCS locus)  / (# DP bulk loci)') # .alleles. .matches.
print('##')
print('##')
print('## The number of bulk loci is   ',f'{132948:,}')
print('## The number of HZ bulk loci is',f'{HZ_bulk_loci:,}') # Oh shit that syntax hngh
print('## The number of DP bulk loci is',f'{n_DP_bulk_loci:,}')
print('## The number of bulk snps is   ',f'{nbulk_snps:,}')
print('##')
print('##')
print('# Sample',
		'SCS_HZ','HZ_rat1','HZ_rat2','HZ_Mrat1','HZ_Mrat2',
		'SCS_snp','snp_rat','snp_Mrat1','snp_Mrat2',
		'SCS_DP','DP_rat1','DP_rat2','DP_rat3','DP_Mrat1','DP_Mrat2',
		sep='\t')
print('Done.',file=sys.stderr)
########################################

########## Parsing SCS samples #########
with open(args.f) as f:
	for sample,line in zip(samples,f): # Main loop
		basename=sample.split('/')[-1]
		sys.stdout.write(basename+'\t'+line)
########################################







































# #!/usr/bin/python3

# import argparse,re,sys,gzip,glob

# ############ Input and such ############
# desc	='''
# This program compares all the catalog-matched loci for all Raphael\'s SCS samples against
# the bulk sample catalog, and looks at what alleles are not present.'''
# parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
# parser.add_argument('b',metavar='bulk snps',help='Path to bulk sample\' snps file.')
# parser.add_argument('s',metavar='sstacks output dir',help='Path to SCS samples\' sstacks output directory.')
# parser.add_argument('n',metavar='stacks file list',help='Path to file containing all stacks files\' filenames')
# args	=parser.parse_args()
# ########################################

# '''

# OOVVEERRHHAAUULL
# This method doesn't work.
# Ignore the SCS .snps. files. All you need is in the .matches. files and the bulk's .snps. file.
# Get number of bulk loci with snps from the bulk snps file.
# Parse .matches. files for each SCS sample, count number of NOT consensus loci. These matched
# to the bulk catalog snps. Compare the two counts for each SCS sample. Win the game.

# 0. Parse bulk sample's snps. Dictionary locus:alleles
# 1. Parse a SCS sample's snps. Dictionary locus:alleles
# 2. Create dictionary bulk locus:sample locus (order of key:value doesn't matter)
# 3. Compare all matched loci's alleles

# Columns in .matches. files:
# catalog ID	sample ID (not useful)	Locus ID	Haplotype	Stack depth	Cigar string

# Columns in .snps. files:
# sample ID (not useful)	Locus ID	Position	Type	Likelihood	Majority nuc	Alternative nuc	...	...



# Questions:
# 1. What about 3rd and 4th (minority) alleles? Currently ignoring them.
# 2. Are positions aligned? Skip all but 'consensus' loci in .matches. files?
# 3. What about snps that occur in the SCS but not in the bulk?
# '''


# ####### Parsing bulk snps files ########
# bulk_snps=dict()
# with gzip.open(args.b,'rt') as f:
# 	next(f)
# 	for line in f:
# 		fields	=line.split()
# 		if fields[7]!='-': continue # Ignore positions with more than two alleles, to simplify things as a starting point
# 		locus	=fields[1]
# 		position=fields[2]
# 		nuc1	=fields[5]
# 		nuc2	=fields[6]
# 		if not locus in bulk_snps: bulk_snps[locus]=[(position,{nuc1,nuc2})]
# 		else: bulk_snps[locus].append((position,{nuc1,nuc2}))
# ########################################

# ########## List of paths ###############
# with open(args.n) as f:
# 	samples=[]
# 	for line in f:
# 		if 'snps' in line and 'catalog' not in line and 'std' not in line:
# 			samples.append('.'.join(line.split('.')[:-3])) # without extensions
# ########################################

# ########## Parsing samples #############
# for sample in samples: # Loop over SCS sample snps files
# 	basename=sample.split('/')[-1]
# 	scs_snps=dict()
# 	with gzip.open(sample+'.snps.tsv.gz','rt') as f:
# 		next(f)
# 		for line in f:
# 			fields	=line.split()
# 			if fields[5]=='-' or fields[6]=='-': continue # Ignore homozygous entries
# 			locus	=fields[1]
# 			position=fields[2]
# 			nuc1	=fields[5]
# 			nuc2	=fields[6]

# 			if not locus in scs_snps: scs_snps[locus]=[(position,{nuc1,nuc2})]
# 			else: scs_snps[locus].append((position,{nuc1,nuc2})) # SCS snps check!

# 	with gzip.open(args.s+'/'+basename+'.matches.tsv.gz','rt') as f:
# 		sample2bulk=dict()
# 		next(f)
# 		for line in f:
# 			if 'consensus' in line:
# 				fields		=line.split()
# 				bulklocus	=fields[0]
# 				samplelocus	=fields[2]
# 				sample2bulk[samplelocus]=bulklocus # This way because the other way does not quite feature unique keys

# 	# Now to calculate "overlap"
# 	allelematches,total	=0,0
# 	common_loci			=sample2bulk.keys()&scs_snps.keys()
# 	for locus in common_loci:
# 		Salleles=scs_snps[locus]
# 		if sample2bulk[locus] in bulk_snps: # We ignore snps that do not occur in the bulk sample
# 			Balleles=bulk_snps[sample2bulk[locus]]
# 			for snp in Salleles:
# 				total+=1
# 				if snp in Balleles: allelematches+=1

# 	if not total:
# 		print(basename,'No common heterozygous loci.',sep='\t')
# 	else: print(basename,allelematches,total,allelematches/total,sep='\t')
# ########################################