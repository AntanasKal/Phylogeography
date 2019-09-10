# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 11:18:17 2019

@author: Antanas
"""

import argparse
import dendropy
import numpy as np

import os


parser = argparse.ArgumentParser(description='Run simulations')
#parser.add_argument('-file', action="store", type=int, dest="file_name", default="output/c_beast/sampled1/beast_input/beast0.xml", help='File name')
parser.add_argument('-sample_index', action="store", type=int, dest="sample_index", default = 1, help='index of sampling scenario for analysis (default: 1)')
parser.add_argument('-index', action="store", type=int, dest="index", default = 0, help='index of sample (default: 0)')

args = parser.parse_args()

index = args.index
sample_index = args.sample_index

os.system('beast -overwrite -seed 123456795 "output/c_beast/sampled'+str(sample_index)+'/beast_input/beast'+str(index)+'.xml"')

sample_taxon_labels = []

for line in open("output/c_beast/sampled"+str(sample_index)+"/root_data/labels"+str(index)+".txt","r"):
    if len(line)>1:
        sample_taxon_labels.append(line[:-1])
        print(line[:-1])
    
    
initialize_dict = False
taxon_dict = {}
for line in open("output/c_beast/sampled"+str(sample_index)+"/beast_output/beast"+str(index)+".trees.txt"):
    
    if initialize_dict==True and line.startswith("		;"):
        break           
    
    if initialize_dict==True:            
        if line[-2]==',':
            [a, b]=line[2:-2].split()
        else:
            [a, b]=line[2:-1].split()            
        
        taxon_dict.update({b : a}) 
    
    if line.startswith("	Translate"):
        initialize_dict = True 
translated_sample_taxon_labels = []
for label in sample_taxon_labels:
    translated_sample_taxon_labels.append(taxon_dict.get(label))        
        
        
file = open("output/c_beast/sampled"+str(sample_index)+"/root_data/observed_roots"+str(index)+".txt", "w")    
    
for line in open("output/c_beast/sampled"+str(sample_index)+"/beast_output/beast"+str(index)+".trees.txt"):
    
    if line.startswith("tree"):
        start_index = 0
        while True:
            if line[start_index:start_index+4]=="[&R]":
                break
            start_index=start_index+1
        single_tree=dendropy.Tree.get(data=line[start_index:], schema="newick", extract_comment_metadata=True)
        for leaf in single_tree.leaf_node_iter():
            print("label: "+leaf.taxon.label)
        mrca = single_tree.mrca(taxon_labels=translated_sample_taxon_labels)
        file.write(mrca.annotations.require_value("location")[0]+"\t"+mrca.annotations.require_value("location")[1]+'\n')
file.close()