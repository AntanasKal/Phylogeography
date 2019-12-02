#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:39:10 2019

@author: Yuxuan
"""

import pandas as pd
from pandas import read_excel
import random

##import the DNA sequence file
df = pd.read_fwf('/WNV/WNV.fasta', sep = ';', header=None)
n = len(df)  
dna = pd.DataFrame(columns=['traits','DNA'], index = range(int(n/2)))
for i in range(int(n/2)):
    dna['traits'][i] = df[0][2*i][1:]
    dna['DNA'][i] = df[0][2*i + 1]

#generate the samples in the east
my_sheet = 'Sheet1'
file_name1 = '/WNV/WNV_lat_long_east.xlsx' # name of your excel file
east = read_excel(file_name1, sheet_name = my_sheet)

east_merged = pd.merge(east, dna, on='traits', how='inner')

dna_east = pd.DataFrame(columns=[0], index = range(2*len(east_merged)))
for i in range(int(n/2)):
    dna_east[0][2*i] = ">" + east_merged['traits'][i]
    dna_east[0][2*i + 1] = east_merged['DNA'][i]

#generate samples in the west
my_sheet = 'Sheet1'
file_name2 = '/WNV/WNV_lat_long_west.xlsx' # name of your excel file
west = read_excel(file_name2, sheet_name = my_sheet)

west_merged = pd.merge(west, dna, on='traits', how='inner')

dna_west = pd.DataFrame(columns=[0], index = range(2*len(west_merged)))
for i in range(int(n/2)):
    dna_west[0][2*i] = ">" + west_merged['traits'][i]
    dna_west[0][2*i + 1] = west_merged['DNA'][i]
    
##add 20 samples randomly from west into east and output the combined sample
eastrandom20fromwest = pd.DataFrame(columns=east_merged.columns, index = range(0))
for i in range(len(east_merged)):
    eastrandom20fromwest.loc[i] = east_merged.iloc[i]
    
t = len(west_merged)
for i in random.sample(range(len(east_merged)), 20):
    eastrandom20fromwest.loc[t] = west_merged.iloc[i]
    t +=1
    
dna_new = pd.DataFrame(columns=[0], index = range(2*len(eastrandom20fromwest)))
for i in range(len(east_merged)):
    dna_new[0][2*i] = ">" + eastrandom20fromwest['traits'][i]
    dna_new[0][2*i + 1] = eastrandom20fromwest['DNA'][i]    
    
for i in range(len(west_merged),len(eastrandom20fromwest)):
    dna_new[0][2*i] = ">" + eastrandom20fromwest['traits'][i]
    dna_new[0][2*i + 1] = len(eastrandom20fromwest['DNA'][0])*'N' 
    
    
##add 20 samples randomly from east into the west and output the combined sample   
westrandom20fromeast = pd.DataFrame(columns=west_merged.columns, index = range(0))
for i in range(len(west_merged)):
    westrandom20fromeast.loc[i] = west_merged.iloc[i]
    
t = len(east_merged)
for i in random.sample(range(len(west_merged)), 20):
    westrandom20fromeast.loc[t] = east_merged.iloc[i]
    t +=1
    
dna_new = pd.DataFrame(columns=[0], index = range(2*len(westrandom20fromeast)))
for i in range(len(east_merged)):
    dna_new[0][2*i] = ">" + westrandom20fromeast['traits'][i]
    dna_new[0][2*i + 1] = westrandom20fromeast['DNA'][i]    
    
for i in range(len(east_merged),len(westrandom20fromeast)):
    dna_new[0][2*i] = ">" + westrandom20fromeast['traits'][i]
    dna_new[0][2*i + 1] = len(westrandom20fromeast['DNA'][0])*'N'     
      
