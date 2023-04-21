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
from utils import replace_seq_info, export_variants_df

#%% Import the reference sequence
databases_path = '/home/msr/Documents/Experimentos/Databases/'
reference_genome_path = databases_path + 'GCA_000001405.28/GCA_000001405.28_GRCh38.p13_genomic.fna'
reference_genome = SeqIO.index(reference_genome_path,"fasta")

ensembl_db = load_db(databases_path+'Ensembl_database_ready.csv')

# %% Extract the SeqRecord of the chromosomes into a dict and then put them into a dataframe
# DO NOT RUN UNLESS ENOUGH RAM
reference_genome_keys = list(reference_genome.keys())
chromosome_keys = []
for key in reference_genome_keys:
    if 'CM' in key:
        chromosome_keys.append(key)

chromosomes_dict = {}

for key in chromosome_keys:
    chromosomes_dict[key] = reference_genome[key] # Extract the correspondent SeqRecords of each chromosome


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

# To extract all the SeqRecords in a dictionary and save them in a df


chromosome_21_modified_df = export_variants_df(chromosome_21_sequence,chromosome_21_single_allele_snv)
