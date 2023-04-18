# -*- coding: utf-8 -*-
"""
Spyder Editor

This script is used to experiment with the Biopython library and the 
GCA_000001405.28_GRCh38.p13_genomic.fna file
"""
#%% Import the libraries
from Bio import SeqIO
# To be able to modify sequences
from Bio.Seq import MutableSeq, Seq
# Import SeqRecord
from Bio.SeqRecord import SeqRecord
# To import the databases directly
from database_io import load_db
# Import pandas
import pandas as pd
from utils import replace_seq_info

#%% Import the reference sequence
databases_path = '/home/msr/Documents/Experimentos/Databases/'
reference_genome_path = databases_path + 'GCA_000001405.28/GCA_000001405.28_GRCh38.p13_genomic.fna'
reference_genome = SeqIO.index(reference_genome_path,"fasta")

ensembl_db = load_db(databases_path+'Ensembl_database_ready.csv')
#%% We will work with the chromosome 21 because it's the shortest
chromosome_21_record = reference_genome["CM000683.2"] # The key is the GenBank identifier
chromosome_21_sequence = chromosome_21_record.seq
# Convert it into a mutable sequence
mutable_chromosome_21 = MutableSeq(chromosome_21_sequence)
# Import the ensembl database through a personalized function and extract the chromosome 21 variants
# In this line we import the database
chromosome_21_variants = ensembl_db.loc[ensembl_db['seqid']=='21']
# Consider that the loci counts start from 1 and Python starts counting from 0

# %% We are going to start only with th SNVs
chromosome_21_snvs = chromosome_21_variants.loc[chromosome_21_variants['type']=='SNV']
#test = chromosome_21_snvs.iloc[1] # The first element (0) of the dataframe has two alleles
# Extract all the SNVs with only 1 variant allele
chromosome_21_single_allele_snv = chromosome_21_snvs.loc[chromosome_21_snvs['Variant_seq'].str.len()==1]
# %% Replace the variant information in the reference sequence
variant_sequences_data = {'sequence':[],
                          'id':[],
                          'chromosome':[],
                          'start':[],
                          'end':[],
                          'Dbxref':[],
                          'clinical_significance':[]}

for index in range(len(chromosome_21_single_allele_snv)):
    variant_sequence_info = replace_seq_info(mutable_chromosome_21,chromosome_21_single_allele_snv.iloc[index])
    variant_sequences_data['sequence'].append(variant_sequence_info.seq)
    variant_sequences_data['id'].append(variant_sequence_info.id)
    variant_sequences_data['chromosome'].append(variant_sequence_info.annotations['Chromosome'])
    variant_sequences_data['start'].append(variant_sequence_info.annotations['Start of mutation'])
    variant_sequences_data['end'].append(variant_sequence_info.annotations['End of mutation'])
    variant_sequences_data['Dbxref'].append(variant_sequence_info.dbxrefs)
    variant_sequences_data['clinical_significance'].append(variant_sequence_info.annotations['Clinical_significance'])
    
# %%
