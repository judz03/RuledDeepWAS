#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 29 19:58:14 2023

@author: msr
This script is to plot some of the characteristics about the data of the
ensembl and dbsnp databases.
"""

import pandas as pd
import matplotlib.pyplot as plt

# function to add value labels
def addlabels(x,y):
    for i in range(len(x)):
        plt.text(i,y[i],y[i])
# %% Import dbSNP and Ensembl dataframes

## Let's consider include both file's pragmas just for reference
ensembl_db = pd.read_csv('/mnt/DATA/Databases/Ensembl/ensembl_data.csv')
dbsnp_db = pd.read_csv('/mnt/DATA/Databases/dbSNP DATA/dbsnp_b154_million_samples.csv')
# %% Check the types of variants in the database and how many there are
variant_count_dbsnp = dbsnp_db['VC'].value_counts()
        
plt.bar(variant_count_dbsnp.index, variant_count_dbsnp.values)
addlabels(variant_count_dbsnp.index, variant_count_dbsnp.values)

#%% Check the most common variation's location
variant_pos_dbsnp = dbsnp_db['POS'].value_counts().sort_index(0)
plt.stem(variant_pos_dbsnp.index,variant_pos_dbsnp.values,markerfmt='o')
plt.xlabel('Posición')
plt.ylabel('Número de instancias')
plt.show()
# %% Now let's check the ensembl data
ensembl_variation_types = ensembl_db['type'].value_counts()
plt.barh(ensembl_variation_types.index,ensembl_variation_types.values)
#addlabels(ensembl_variation_types.index,ensembl_variation_types.values)
plt.ylabel('Tipo de variación')
plt.xlabel('Número de instancias')
plt.title('Tipos y número de variaciones contempladas en Ensembl')
plt.show()
# %% Ensembl's variation positions
ensembl_variation_pos = ensembl_db['start'].value_counts().sort_index(0)
plt.stem(ensembl_variation_pos.index,ensembl_variation_pos.values)
plt.xlabel('Posición')
plt.ylabel('Número de instancias')
plt.title('Número y posición de las variantes (Ensembl)')
plt.show()