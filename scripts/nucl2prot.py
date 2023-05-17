#!/usr/bin/python3

import argparse,sys,textwrap

############ Input and such ############
desc	='''This program takes a nucleotide fasta and makes it an aminoacid fasta.
Prints to standard out.'''
parser	=argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('i',metavar='infasta',help='Path to input file.')
args	=parser.parse_args()
########################################

code = {'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
		'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
		'AGA':'R','AGG':'R',
		'AAT':'N','AAC':'N','GAT':'D','GAC':'D',
		'TGT':'C','TGC':'C',
		'CAA':'Q','CAG':'Q','GAA':'E','GAG':'E',
		'GGT':'G','GGC':'G','GGA':'G','GGG':'G',
		'CAT':'H','CAC':'H',
		'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
		'TTA':'L','TTG':'L','CTT':'L','CTC':'L',
		'CTA':'L','CTG':'L',
		'AAA':'K','AAG':'K','TTT':'F','TTC':'F',
		'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
		'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
		'AGT':'S','AGC':'S',
		'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
		'TGG':'W','TAT':'Y','TAC':'Y',
		'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
		'TAA':'*','TGA':'*','TAG':'*'}

########## Parsing nucleos #############
with open(args.i) as f:
	for line in f:
		if line.startswith('>'):
			sys.stdout.write(line)
		else:
			block=[]
			for codon in [line[i:i+3] for i in range(0,len(line),3)]:
				if codon in code: block.append(code[codon])
			sys.stdout.write(''.join(block)+'\n')
########################################