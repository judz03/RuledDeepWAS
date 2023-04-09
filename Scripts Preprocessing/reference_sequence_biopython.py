# -*- coding: utf-8 -*-
"""
Spyder Editor

This script is used to experiment with the Biopython library and the 
GCA_000001405.29_GRCh38.p14_genomic.fna file
"""
#%% Import the libraries
from Bio import SeqIO

#%% Import the reference sequence
databases_path = '/home/mosotelo/Project Files/Experimentos/Databases/'
reference_genome_path = databases_path + 'GCA_000001405.29_GRCh38.p14_genomic.fna'
sequence = SeqIO.index(reference_genome_path,"fasta")

# %%
# Load the GRCh38 reference genome sequence in FASTA format
#genome_file = "GRCh38.fasta"
genome_file = reference_genome_path
genome = SeqIO.index(genome_file, "fasta") # This generates a dictionary-like object

# Get the keys and values from 'genome' creating a new dictionary
genome_dict = dict(genome)

'''
# Get the reference genome sequence for chromosome 1
chromosome = "chr1"
start = 0
end = 1000
ref_seq = genome[chromosome][start:end].seq

# Print the reference genome sequence
print(ref_seq)
'''