# %%
# -*- coding: utf-8 -*-
# Import libraries
"""
This is an script divided in cells which main objective is to import, load and explore the data from the dbSNP, Ensembl
and GWAS Catalog databases.
"""
import numpy as np
import pandas as pd
import vcf #PyVCF
import gffutils
# %% Load Ensembl data from plain text (GVF File)
## Import and read the file
with open('/mnt/DATA/Databases/Ensembl/homo_sapiens_phenotype_associated.gvf',mode = 'r') as ensembl_file:
    ensembl_content=ensembl_file.read()
# Recortar archivo para mejor manejo
plain_pragmas = ensembl_content[:1165] #Hay que generalizar, los pragmas son líneas que empiezan con ## y terminan con \n
plain_variants = ensembl_content[1166:]
ensembl_plain_variants_info_split = plain_variants.split('\n')
ensembl_info_dict = {'seqid':[],
                     'source':[],
                     'type':[],
                     'start':[],
                     'end':[],
                     'attributes':[]}
for i in range(len(ensembl_plain_variants_info_split)-1):
    # The -1 is because the last element in the ensembl_plain_variants_info_split variable is empty
    current_list = ensembl_plain_variants_info_split[i].split('\t')
    ensembl_info_dict['seqid'].append(current_list[0])
    ensembl_info_dict['source'].append(current_list[1])
    ensembl_info_dict['type'].append(current_list[2])
    ensembl_info_dict['start'].append(current_list[3])
    ensembl_info_dict['end'].append(current_list[4])
    ensembl_info_dict['attributes'].append(current_list[8])

ensembl_df = pd.DataFrame(data=ensembl_info_dict)
# Save in dataframe in a csv. Normally commented
#ensembl_df.to_csv('/mnt/DATA/Databases/Ensembl/ensembl_data.csv', index=False)

# Provided by Chat GPT # It is a good reference, but I need to organize the data differently
#db = gffutils.create_db("/mnt/DATA/Databases/Ensembl/homo_sapiens_phenotype_associated.gvf", dbfn='/mnt/DATA/Databases/Ensembl/EnsemblData.db', force=True)
#for feature in db.features_of_type('variant'):
#    print(feature)

# Check the contents of the GVF File
#fn = gffutils.example_filename('/mnt/DATA/Databases/Ensembl/homo_sapiens_phenotype_associated.gvf')
#print(open(fn).read())

# %% Load the GWAS Catalog Variants-Associations File (specific columns)
# Import the file as a DataFrame
gwas_catalog_df = pd.read_csv('/mnt/DATA/Databases/GWAS Catalog DATA/gwas_catalog_v1.0.2-associations_e108_r2022-12-21.tsv',
                              delimiter='\t',
                              usecols=['SNP_ID_CURRENT','SNPS','STRONGEST SNP-RISK ALLELE','CHR_ID','CHR_POS','REGION','CONTEXT','INTERGENIC','DISEASE/TRAIT'])

# %% Load the GWAS Catalog Ancestries File
# Import the file as a DataFrame
# Uncomment the line below later
#gwas_catalog_ancestry_data_DF = pd.read_csv('/mnt/DATA/Databases/GWAS Catalog DATA/gwas_catalog-ancestry_r2023-01-14.tsv',delimiter='\t')

# %% Load the ClinVar databse VCF file
# Just to check how the VCF looks
#clinvar_vcf_file = open("/mnt/DATA/Databases/ClinVar/clinvar_20230121.vcf", mode='r')
#clinvar_plain_text = clinvar_vcf_file.read()
#clinvar_vcf_file.close()
#linvar_splitted = clinvar_plain_text.split('\n')

# Provided by ChatGPT
with open('/mnt/DATA/Databases/ClinVar/clinvar_20230121.vcf', 'r') as f:
    vcf_reader = vcf.Reader(f)
    data = []
    for record in vcf_reader:
        data.append([record.CHROM, record.POS, record.ID, record.REF, record.ALT])
    clinvar_df = pd.DataFrame(data, columns=['CHROM','POS','ID','REF','ALT'])