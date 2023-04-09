#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 19:02:43 2023

@author: msr

This script is used to extract the attributes contained in each variant from Ensebl
We are going to extract the following attributes:
    - Dbxref
    - Reference_seq
    - Variant_seq
    - clinical_signifcance

I have a test set with 100 samples named test.csv in the following directory:
    /mnt/DATA/Databases/Ensembl/test.csv
"""
# %% Load the necessary modules
import pandas as pd
import re
# %% Import the test.csv with 100 samples
test_df = pd.read_csv('/mnt/DATA/Databases/Ensembl/test.csv') # Consider add the parameter dtype=str
# %% Check all the existing attributes within the strings in all variants
def get_attributes(database_dataframe: pd.DataFrame) -> dict:
    attributes = database_dataframe['attributes'] # extract only the attributes column from the complete dataframe
    attribute_regex = '(\w+)=((\w+)(,\w*)*)' # All the attributes have this pattern
    pattern = re.compile(attribute_regex)
    keys_dict = {} # This dictionary will save all the existing attributes and the number of times they repeat in the db
    
    # Iterate trough the attributes DF and extract the keys contained in the string
    for attribute in attributes:
        matches = re.findall(pattern, attribute) # matches is a list of tuples containing the groups found by the regex
        for group in matches:
            key = group[0]
            if key not in keys_dict:
                keys_dict[key] = 1
            else:
                keys_dict[key]+=1
    return keys_dict
# %% Get the attributes that we are looking for
# Save the found keys in a list
keys_dict = get_attributes(test_df)
keys_list = list(keys_dict.keys())