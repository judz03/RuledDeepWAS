# -*- coding: utf-8 -*-
"""
Spyder Editor

This script is used to experiment with the Biopython library and the 
GCA_000001405.28_GRCh38.p13_genomic.fna file
"""
#%% Import the libraries
from Bio import SeqIO
# To be able to modify sequences
from Bio.Seq import MutableSeq

#%% Import the reference sequence
databases_path = '/home/msr/Documents/Experimentos/Databases/'
reference_genome_path = databases_path + 'GCA_000001405.28/GCA_000001405.28_GRCh38.p13_genomic.fna'
reference_genome = SeqIO.index(reference_genome_path,"fasta")
#%% We will work with the chromosome 21 because it's the shortest
chromosome_21_record = reference_genome["CM000683.2"] # The key is the GenBank identifier
chromosome_21_sequence = chromosome_21_record.seq
# Convert it into a mutable sequence
mutable_chromosome_21 = MutableSeq(chromosome_21_sequence)
# Import the ensembl database through a personalized function and extract the chromosome 21 variants
# In this line we import the database
chromosome_21_variants = ensembl_df.loc[ensembl_df['seqid']=='21']
# Consider that the loci counts start from 1 and Python starts counting from 0
