#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 18:29:15 2023

@author: msr
This script is used for exploring and integrating the data from the dbsnp and 
ensembl once they have been converted into dataframes and exported in csv format.
Refer to open_GVF.py and open_CVF_dbSNP.py files
"""
import pandas as pd
import re #Let's try to import the pragmas using regular expressions
# %% Import dbSNP and Ensembl dataframes

## Let's consider include both file's pragmas just for reference
# Just loading 1000 rows for memory saving and efficiency
ensembl_db = pd.read_csv('/mnt/DATA/Databases/Ensembl/ensembl_data.csv',dtype=str,nrows=1000)
dbsnp_db = pd.read_csv('/mnt/DATA/Databases/dbSNP DATA/10milldbSNP.csv',nrows=1000)
# %% Split the information contained in the INFO column of the dbsnp data
dbsnp_info = dbsnp_db['INFO']
dbsnp_list_info=[]
for variant in dbsnp_info:
    dbsnp_list_info.append(variant.replace('\n','').split(';'))
    
# %% Extract all the kinds of variants from dbsnp_data

# IN THIS ONE THE PROCESS IS CORRECT BECAUSE I'M EXTRACTING ONLY ONE FEATURE
dbsnp_kinds_of_variants = []
# We will search the type of variant using regular expressions
type_of_variant_re = 'VC=(\w+)'
pattern = re.compile(type_of_variant_re)

for variant in dbsnp_list_info:
    for attribute in variant:
        match = pattern.search(attribute)
        if match:
            dbsnp_kinds_of_variants.append(match.group(1)) # Indicating the 1 captures everything after the '=' symbol
            # for this re in specific
            pass # Once the pattern has been found go to next variant

# Insert the new 'types of variants' into the original dataframe and save it in csv
#dbsnp_db.insert(2,'VC',dbsnp_kinds_of_variants)
# %% Organize ensembl attributes
# Separate the attributes using ';' as delimeter per variant and save the resulting list in another list
ensembl_variant_attributes = []
ensembl_attributes = ensembl_db['attributes']
for variant in ensembl_attributes:
    ensembl_variant_attributes.append(variant.split(';'))
# %% Extract all the keys contained in the attributes
# HERE, I NEED TO EXTRACT ALL THE ATTRIBUTES JUST TO KNOW HOW MANY THERE ARE
attribute_pattern = '(\w+)=((\w+)(,\w*)*)'
pattern = re.compile(attribute_pattern)
ensembl_attribute_values = []
ensembl_attribute_keys = []
for variant in ensembl_variant_attributes:
    current_value_list = []
    current_key_list = []
    for attribute in variant:
        match = pattern.search(attribute)
        if match:
            current_value_list.append(match.group(2))# Indicating the 2 captures everything after the '=' symbol
            current_key_list.append(match.group(1))
            # for this re in specific
    ensembl_attribute_values.append(current_value_list)
    ensembl_attribute_keys.append(current_key_list) # Here, I save the attribute keys contained per variant
        
# %% To check and count the attributes in ensembl database
keys_dic ={}
# This for loop counts and adds all the attributes contained in the whole Ensembl databse
for c_list in ensembl_attribute_keys:
    for c_attribute in c_list:
        if c_attribute in keys_dic:
            keys_dic[c_attribute]=keys_dic[c_attribute]+1
        else:
            keys_dic[c_attribute]=1
# %% Now, let's add all the attributes per variant in a DataFrame

