#!/usr/bin/python3

import argparse,sys



############ Input and such ############
desc	='''This program takes the 8 dinoflagellate transcriptomes and merges them into 1 transcriptome.
No duplicate sequences.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
args	=parser.parse_args()
########################################

######## One dict to find them #########
Biggus_Dictus=dict()
########################################

##### One loop to bring them all #######
mmetsps	=['MMETSP0359.cds.fa','MMETSP0360.cds.fa','MMETSP0361.cds.fa','MMETSP0367.cds.fa','MMETSP0368.cds.fa','MMETSP0369.cds.fa','MMETSP0370.cds.fa','MMETSP0371.cds.fa']
prefix	='/crex/proj/uppstore2017173/b2014326/Dino_transcriptomes/'
mmetsps	=[prefix+mmetsp for mmetsp in mmetsps]

for mmetsp in mmetsps:
	with open(mmetsp) as f:
		block,ID=[],next(f).strip('\n>').split()[0] # First line (removed the length= part)
		for line in f: # All other lines
			if line.startswith('>'):
				sequence=''.join(block).upper()
				if sequence not in Biggus_Dictus: Biggus_Dictus[sequence]=ID
				else: Biggus_Dictus[sequence]+='___'+ID
				block,ID=[],line.strip('\n>').split()[0]
			else:
				block.append(line.rstrip())
		sequence=''.join(block).upper()
		if sequence not in Biggus_Dictus: Biggus_Dictus[sequence]=ID # Last one
		else: Biggus_Dictus[sequence]+='___'+ID
########################################

#### and in the darkness bind them #####
with open('dinotranscriptome.fasta','w') as o:
	for seq in Biggus_Dictus:
		o.write('>'+Biggus_Dictus[seq]+'\n'+seq+'\n')
########################################