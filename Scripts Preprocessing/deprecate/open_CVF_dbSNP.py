#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 09:46:44 2023

@author: msr
This is a script to load and print a range of specific lines from the dbSNP 
database. It can't be fully loaded due to its size (118GB)
"""
import pandas as pd
import vcf
# %% Load the dbSNP database VCF file

# BIG FILE, HIGH RAM USAGE! (30GB)
#with open('/mnt/DATA/Databases/dbSNP DATA/GCF_000001405.38.vcf', 'r') as f:
#    dbsnp_vcf = vcf.Reader(f)
#    dbsnp_data = []
#    for record in dbsnp_vcf:
#        dbsnp_data.append([record.CHROM, record.POS, record.ID, record.REF, record.ALT])
#    dbsnp_df = pd.DataFrame(data, columns=['CHROM','POS','ID','REF','ALT'])

# Open the VCF file line by line
with open("/mnt/DATA/Databases/dbSNP DATA/GCF_000001405.38.vcf", "r") as file:
    # Read the first 1000 lines
    lines = [next(file) for i in range(10000000)]

# Print the first 1000 lines
#for line in lines:
#    print(line)
# %% Split and organize the data
dbsnp_dc = {}
dbsnp_pragmas = lines[:36] # This is hardcoded as well, the pragmas start with ## and end with \n
dbsnp_headers = lines[37] # This line is important because I can iterate trough the dictionary with its values
dbsnp_plain_data = lines[38:]

dbsnp_headers = dbsnp_headers.split('\t')

if dbsnp_headers[-1] != 'INFO':
    dbsnp_headers[-1]='INFO' #It has a space at the end of the word
    
for header in dbsnp_headers: # This loop creates the dictionary keys based on the headers extracted from the file
    dbsnp_dc[header]=[]

dbsnp_data = []
for variant in dbsnp_plain_data:
    dbsnp_data.append(variant.split('\t'))
    
for variant in range(len(dbsnp_data)): # This can be better, it's harcoded, not generalized
    dbsnp_dc['#CHROM'].append(dbsnp_data[variant][0])
    dbsnp_dc['POS'].append(dbsnp_data[variant][1])
    dbsnp_dc['ID'].append(dbsnp_data[variant][2])
    dbsnp_dc['REF'].append(dbsnp_data[variant][3])
    dbsnp_dc['ALT'].append(dbsnp_data[variant][4])
    dbsnp_dc['QUAL'].append(dbsnp_data[variant][5])
    dbsnp_dc['FILTER'].append(dbsnp_data[variant][6])
    dbsnp_dc['INFO'].append(dbsnp_data[variant][7])
# %% Create DataFrame with all the previous information
dbsnp_df = pd.DataFrame(data=dbsnp_dc)