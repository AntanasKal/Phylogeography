#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 11:01:37 2019

@author: Yuxuan
"""
import pandas as pd
import re

##read the YFV_cauchy_RRW_north.trees
df = pd.read_fwf('/output/YFV_output/YFV_cauchy_RRW_north.trees', sep = ' ', header=None)

##add 800080000 to each tree state number so that the tree can be added to YFV_cauchy_RRW_north2.trees 
for j in range(len(df[1])):
    data = df[1][j].split() #split string into a list
    data[0] = int(re.sub('STATE_', '', data[0]))  
    data[0] +=800080000
    data[0] = str(data[0])
    data[0] = 'STATE_' + data[0]

    char = data[0]
    for i in range(1 ,len(data)):
        char = char + ' ' + data[i]
    
    df[1][j] = char

##This will create a file which can be added to YFV_cauchy_RRW_north2.trees   
file = open("output/YFV_output/new", "w")
for i in range(len(df)):
    file.write(df[0][i] + ' ' + df[1][i] +'\n')

file.close()
