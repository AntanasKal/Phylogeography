#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 14:51:09 2019

@author: Yuxuan
"""


import pandas as pd
from pandas import read_excel
import matplotlib.pyplot as plt
import random

##import the DNA sequence file
df = pd.read_fwf('/YFV/YFV_sequences.fasta', sep = ';', header=None)
n = len(df)  
dna = pd.DataFrame(columns=['traits','DNA'], index = range(int(n/2)))
for i in range(int(n/2)):
    dna['traits'][i] = df[0][2*i][1:]
    dna['DNA'][i] = df[0][2*i + 1]

##draw a picture of separating the points into north and south
my_sheet = 'Sheet1'
file_name3 = '/YFV/YFV_coordinates.xlsx' # name of your excel file
pic = read_excel(file_name3, sheet_name = my_sheet)    
   

x = list()
y = list()
for i in range(len(pic['lat'])):
    y.append(pic['lat'][i])
    x.append(pic['long'][i])
plt.scatter(x, y, alpha=0.6)  
plt.plot([-48,-40],[-19,-19] , color='red')
plt.show()  
    
#generate the north files
my_sheet = 'Sheet1'
file_name4 = '/YFV/YFV_coordinates_north.xlsx' # name of your excel file
north = read_excel(file_name4, sheet_name = my_sheet)

north_merged = pd.merge(north, dna, on='traits', how='inner')

dna_north = pd.DataFrame(columns=[0], index = range(2*len(north_merged)))
for i in range(int(n/2)):
    dna_north[0][2*i] = ">" + north_merged['traits'][i]
    dna_north[0][2*i + 1] = north_merged['DNA'][i]

#generate the south files
my_sheet = 'Sheet1'
file_name5 = '/YFV/YFV_coordinates_south.xlsx' # name of your excel file
south = read_excel(file_name5, sheet_name = my_sheet)

south_merged = pd.merge(south, dna, on='traits', how='inner')

dna_south = pd.DataFrame(columns=[0], index = range(2*len(south_merged)))
for i in range(int(n/2)):
    dna_south[0][2*i] = ">" + south_merged['traits'][i]
    dna_south[0][2*i + 1] = south_merged['DNA'][i]    
    
    
##add 10 samples randomly from south into north
northrandom10fromsouth = pd.DataFrame(columns=north_merged.columns, index = range(0))
for i in range(len(north_merged)):
    northrandom10fromsouth.loc[i] = north_merged.iloc[i]
    
t = len(north_merged)
for i in random.sample(range(len(south_merged)), 10):
    northrandom10fromsouth.loc[t] = south_merged.iloc[i]
    t +=1
    
dna_10s = pd.DataFrame(columns=[0], index = range(2*len(northrandom10fromsouth)))
for i in range(len(north_merged)):
    dna_10s[0][2*i] = ">" + northrandom10fromsouth['traits'][i]
    dna_10s[0][2*i + 1] = northrandom10fromsouth['DNA'][i]    
    
for i in range(len(north_merged),len(northrandom10fromsouth)):
    dna_10s[0][2*i] = ">" + northrandom10fromsouth['traits'][i]
    dna_10s[0][2*i + 1] = len(northrandom10fromsouth['DNA'][0])*'N' 
    
    
##add 10 samples randomly from north  into south
southrandom10fromnorth = pd.DataFrame(columns=south_merged.columns, index = range(0))
for i in range(len(south_merged)):
    southrandom10fromnorth.loc[i] = south_merged.iloc[i]
    
t = len(south_merged)
for i in random.sample(range(len(north_merged)), 10):
    southrandom10fromnorth.loc[t] = north_merged.iloc[i]
    t +=1
    
dna_10n = pd.DataFrame(columns=[0], index = range(2*len(southrandom10fromnorth)))
for i in range(len(south_merged)):
    dna_10n[0][2*i] = ">" + southrandom10fromnorth['traits'][i]
    dna_10n[0][2*i + 1] = southrandom10fromnorth['DNA'][i]    
    
for i in range(len(south_merged),len(southrandom10fromnorth)):
    dna_10n[0][2*i] = ">" + southrandom10fromnorth['traits'][i]
    dna_10n[0][2*i + 1] = len(southrandom10fromnorth['DNA'][0])*'N' 
