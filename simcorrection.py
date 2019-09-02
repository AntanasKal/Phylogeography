# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 16:47:06 2019

@author: Antanas
"""

import argparse
import dendropy
import numpy as np

import sampling
import beastxmlwriter
import treegenerator

import os


parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-N', action="store", type=int, dest="num_trees", default=5, help='number of simulations (default 5)')
parser.add_argument('-jobi', action="store", type=int, dest="job_index", default=1, help='job index')

args = parser.parse_args()

num_trees = args.num_trees
job_index = args.job_index



num_tips = 50
sample_size = 5
other_sample_size = 5
sigma = 1
seq_len = 100 # 10000

#num_simulations = 1




for i in range(num_trees*(job_index-1), num_trees*job_index):
    
    tree = treegenerator.generate_yule_tree(num_tips)    
    tree = treegenerator.simulate_brownian(tree, sigma, 2) 
    tree= treegenerator.simulate_brownian(tree, sigma, dimension=2) 
    #t=calculate_time_to_tips(t)
    
    
    
    sampled_t = sampling.sample_biased_extreme(tree, sample_ratio = sample_size/num_tips, dimension = 2)    
    d = dendropy.model.discrete.hky85_chars(kappa=3, mutation_rate=0.1, seq_len=seq_len, tree_model=sampled_t, retain_sequences_on_tree=False)    
    beastxmlwriter.write_BEAST_xml_corrected(tree, sampled_t, d, i=i,  mcmc=1000000, log_every=1000, beast_input_string ="output/c_beast/beast_input/beast", beast_output_string="output/c_beast/beast_output/beast", other_sample_size=other_sample_size, seq_len=seq_len)
    
    beastxmlwriter.write_BEAST_xml(sampled_t, i=i, dimension=2, mcmc=1000000, log_every=1000, beast_input_string ="output/c_beast/beast_input/nbeast", beast_output_string="output/c_beast/beast_output/nbeast")
    
    
    os.system('beast -overwrite -seed 123456795 "output/c_beast/beast_input/beast'+str(i)+'.xml"')
    os.system('beast -overwrite -seed 123456795 "output/c_beast/beast_input/nbeast'+str(i)+'.xml"')
    
    file = open("output/actual_root"+str(i)+".txt", "w")
    file.write(str(sampled_t.seed_node.X)+'\n')
    file.write(str(sampled_t.seed_node.Y)+'\n')
    file.close()
    
    file = open("output/observed_roots"+str(i)+".txt", "w")
    f = open("output/c_beast/beast_output/beast"+str(i)+".trees.txt", "r")
    for line in f:
        if line.startswith("tree"):
            single_t = dendropy.Tree.get_from_string(line, extract_comment_metadata=True, suppress_internal_node_taxa=False, schema='newick')
            file.write(single_t.seed_node.annotations.require_value("X")+"\t"+single_t.seed_node.annotations.require_value("Y")+'\n')
    file.close()
    
    
    