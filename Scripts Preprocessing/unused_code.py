#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 15 15:53:08 2023

@author: msr
I will put here the code that I didn't use in case I need to refer to it in the future
"""

# Extract the indexes where the length of attributes is equal to 7
indexes_list = []
for index in range(len(ensembl_splitted_attributes)):
    if len(ensembl_splitted_attributes[index]) == 7:
        indexes_list.append(index)
    else:
        pass
