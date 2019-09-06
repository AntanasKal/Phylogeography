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
import time
import os


parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-N', action="store", type=int, dest="num_trees", default=5, help='number of simulations (default 5)')
parser.add_argument('-jobi', action="store", type=int, dest="job_index", default=1, help='job index')

args = parser.parse_args()

num_trees = args.num_trees
job_index = args.job_index



num_tips = 500
sample_size = 50
other_sample_size = 50
sigma = 1
seq_len = 10000 # 10000
mcmc = 1e7


#num_simulations = 1




for i in range(num_trees*(job_index-1), num_trees*job_index):
    
    tree = treegenerator.generate_yule_tree(num_tips)    
    tree = treegenerator.simulate_brownian(tree, sigma, 2) 
    tree= treegenerator.simulate_brownian(tree, sigma, dimension=2) 
    #t=calculate_time_to_tips(t)
    
    
    
    sampled_t = sampling.sample_biased_extreme(tree, sample_ratio = sample_size/num_tips, dimension = 2)    
    d = dendropy.model.discrete.hky85_chars(kappa=3, mutation_rate=0.05, seq_len=seq_len, tree_model=sampled_t, retain_sequences_on_tree=False)    
    beastxmlwriter.write_BEAST_xml_corrected(tree, sampled_t, d, i=i,  mcmc=10000, log_every=1000, beast_input_string ="output/c_beast/beast_input/beast", beast_output_string="output/c_beast/beast_output/beast", other_sample_size=other_sample_size, seq_len=seq_len)
    
    beastxmlwriter.write_BEAST_xml(sampled_t, i=i, dimension=2, mcmc=100, log_every=10, beast_input_string ="output/c_beast/beast_input/nbeast", beast_output_string="output/c_beast/beast_output/nbeast")
    
    
    os.system('beast -overwrite -seed 123456795 "output/c_beast/beast_input/beast'+str(i)+'.xml"')
    os.system('beast -overwrite -seed 123456795 "output/c_beast/beast_input/nbeast'+str(i)+'.xml"')
    
    file = open("output/root_data/actual_root"+str(i)+".txt", "w")
    file.write(str(sampled_t.seed_node.X)+'\n')
    file.write(str(sampled_t.seed_node.Y)+'\n')
    file.close()
    
    
    time.sleep(2)
    
    sample_taxon_labels = []
    
    for leaf in sampled_t.leaf_node_iter():
        sample_taxon_labels.append(leaf.taxon.label)
    
    initialize_dict = False
    taxon_dict = {}
    for line in open("output/c_beast/beast_output/beast"+str(i)+".trees.txt"):
        
        if initialize_dict==True and line.startswith("		;"):
            break           
        
        if initialize_dict==True:
#            print(line[2:-2])
            
            if line[-2]==',':
                [a, b]=line[2:-2].split()
            else:
                [a, b]=line[2:-1].split()            
#            print(a)
#            print(b)
            
            taxon_dict.update({b : a})
            #break   
        
        if line.startswith("	Translate"):
            #print("adsfasdf")
            initialize_dict = True 
    
    print(taxon_dict)
    translated_sample_taxon_labels = []
    for label in sample_taxon_labels:
        translated_sample_taxon_labels.append(taxon_dict.get(label))
        print(taxon_dict.get(label))
    
      
        
    file = open("output/root_data/observed_roots"+str(i)+".txt", "w")    
    
    for line in open("output/c_beast/beast_output/beast"+str(i)+".trees.txt"):
        
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
        
#    treelist = dendropy.TreeList.get(path="output/c_beast/beast_output/beast"+str(i)+".trees.txt", extract_comment_metadata=True, schema="nexus")
#    for single_tree in treelist:
#        
#        mrca = single_tree.mrca(taxon_labels=sample_taxon_labels)
#        file.write(mrca.annotations.require_value("location")[0]+"\t"+mrca.annotations.require_value("location")[1]+'\n')
#        
#    
    
    
    file = open("output/root_data/old_observed_roots"+str(i)+".txt", "w")    
    treelist = dendropy.TreeList.get(path="output/c_beast/beast_output/nbeast"+str(i)+".trees.txt", extract_comment_metadata=True, schema="nexus")
    for single_tree in treelist:
        file.write(single_tree.seed_node.annotations.require_value("location")[0]+"\t"+single_tree.seed_node.annotations.require_value("location")[1]+'\n')
    file.close() 