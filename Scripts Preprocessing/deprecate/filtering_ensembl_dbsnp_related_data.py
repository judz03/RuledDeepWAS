#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 20:36:58 2023

@author: msr
This script is meant to filter the Ensembl database.
The objective is to filter the dbSNP-related variants with 7 attributes, since
they represent a great percentage of the database.
Many of the code here was extracted from another script I wrote, "organize_data.py"

We have to extract the data from the attributes contained in the variants from the Ensembl database.
"""
# %% Load necessary modules
import pandas as pd
import re
# %% Load the Ensembl database dataframe
ensembl_db = pd.read_csv('/mnt/DATA/Databases/Ensembl/ensembl_data.csv',dtype=str)
# %% Quit all the variants that don't come from dbSNP
ensembl_dbsnp_related = ensembl_db.loc[ensembl_db['source']=='dbSNP']
# %% Organize ensembl attributes
# Separate the attributes using ';' as delimeter per variant and save the resulting list in another list
ensembl_splitted_attributes = []
ensembl_attributes = ensembl_dbsnp_related['attributes']
for variant in ensembl_attributes:
    ensembl_splitted_attributes.append(variant.split(';'))
# %% Extract the attributes keys and count them
keys_pattern = '(\w+)='
pattern = re.compile(keys_pattern)
# %% Extract the attributes now to rearrange them in columns taking the filtered data
'''
# Separate the attributes using ';' as delimeter per variant and save the resulting list in another list
filtered_splitted_attributes = []
filtered_attributes = filtered_data['attributes']
for variant in filtered_attributes:
    filtered_splitted_attributes.append(variant.split(';'))

# HERE, I NEED TO EXTRACT ALL THE ATTRIBUTES JUST TO KNOW HOW MANY THERE ARE
attribute_pattern = '(\w+)=((\w+)(,\w*)*)'
pattern = re.compile(attribute_pattern)
ensembl_attribute_values = []
ensembl_attribute_keys = []
for variant in filtered_splitted_attributes:
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
'''