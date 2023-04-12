#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 21:27:25 2023

@author: msr
This script is meant to load the pre-processed data we have been working on.
Ensembl
dbSNP
GWAS Catalog
"""
# %% Import necessary modules
import pandas as pd
import re
import matplotlib.pyplot as plt

# %% Load databses' data frames (just 100 lines to save memory)
#ensembl_df = pd.read_csv('/mnt/DATA/Databases/Ensembl/Ensembl_database_ready.csv', dtype = str)
data_path = '/home/msr/Documents/Experimentos/Databases/' #Change depending on the PC the script is being run
ensembl_db_path = data_path+'Ensembl_database_ready.csv'
#ensembl_df = pd.read_csv(ensembl_db_path,
#                         dtype = str,
#                         nrows=100)
ensembl_df = pd.read_csv(ensembl_db_path,
                         dtype = str)

'''
dbsnp_df = pd.read_csv('/mnt/DATA/Databases/dbSNP DATA/dbsnp_b154_million_samples.csv',
                       dtype = str,
                       nrows=100)
#gwas_catalog_df = pd.read_csv('/mnt/DATA/Databases/GWAS Catalog DATA/gwas_catalog_v1.0.2-associations_e108_r2022-12-21.tsv',
                              delimiter='\t',
                              usecols=['SNP_ID_CURRENT','SNPS','STRONGEST SNP-RISK ALLELE','CHR_ID','CHR_POS','REGION','CONTEXT','INTERGENIC','DISEASE/TRAIT'],
                              dtype = str,
                              nrows = 100)
'''
# %% Load the reference genome
# read in the reference genome string
with open("/mnt/DATA/Databases/GCA_000001405.29_GRCh38.p14_genomic.fna") as f:
    ref_genome = f.read()
# %% Replace variant information in the reference genome
# define the variant information
chromosome = "chr1"
position = 12345
ref_base = "A"
alt_base = "C"

# determine the position and length of the region to replace
start = position - 1
end = position

# replace the nucleotides in the reference genome string
ref_genome = ref_genome[:start] + alt_base + ref_genome[end:]

# write out the modified reference genome string
with open("modified_reference_genome.fa", "w") as f:
    f.write(ref_genome)