# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:37:34 2019

@author: Antanas
"""

#the code that wirtes BEAST xml files is in different file, "beastxmlwriter.py"
import math 
import argparse
import dendropy
import random
import numpy as np


import sampling
import beastxmlwriter
import phyrexxmlwriter
import treegenerator

import os





parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-dims', dest='dimension', type=int, default=1, help='number of dimensions (1 or 2) for which the random walk is generated (default: 1)')
parser.add_argument('-N', action="store", type=int, dest="num_trees", default=5, help='number of simulations (default 5)')
parser.add_argument('-treetype', action="store", dest="tree_type", default='nuc', help='type of tree generated \n "nuc" - nonultrametric coalescent \n "uc" - ultrametric coalescent \n "bd" - birth-death tree')
parser.add_argument('-mcmc', action="store", type=int, dest="mcmc", default=5000, help='number of simulations (default 5000)')
parser.add_argument('-logevery', action="store", type=int, dest="log_every", default=10, help='number of simulations (default 10)')
parser.add_argument('-ntips', action="store", type=int, dest="ntips", default=100, help='number of tips for ultrametric coalescent tree(default 100)')
parser.add_argument('-ntipspp', action="store", type=int, dest="num_tips_per_period", default=20, help='number of tips for per period for nonultrametric coalescent tree (default 20)')
parser.add_argument('-nps', action="store", type=int, dest="num_periods", default=10, help='number of periods for non-ultrametric coalescent tree (default 10)')
parser.add_argument('-lambda', action="store", type=float, dest="lamb", default=1, help='rate of coalescent for coalescent trees (default 1)')
parser.add_argument('-pl', action="store", type=float, dest="period_length", default=1, help='period length for nonultrametric coalescent tree (default 1)')
parser.add_argument('-sigma', action="store", type=float, dest="sigma", default=1, help='standart deviation for brownian diffusion simulation for a branch of length 1 (default 1)')
parser.add_argument('-jobi', action="store", type=int, dest="job_index", default=1, help='job index')
parser.add_argument('--linux', dest='linux', action='store_const', const=True, default=False, help='is the program run on Linux? (default: False)')
parser.add_argument('--annotate', dest='run_tree_annotator', action='store_const', const=True, default=False, help='run tree annotator (default: False)')

sample_size = 50
other_sample_size = 50
seq_len = 10000



args = parser.parse_args()
dimension = args.dimension
num_trees = args.num_trees
mcmc = args.mcmc
log_every = args.log_every
sigma = args.sigma
num_tips = args.ntips
num_tips_per_period = args.num_tips_per_period
num_periods=args.num_periods
lamb=args.lamb
period_length=args.period_length
job_index = args.job_index
linux = args.linux
run_tree_annotator = args.run_tree_annotator


#####
sample_ratio=0.05
num_sampling = 4
corr_beast_mcmc=1e7
generate_corrected_files=True
#####

generate_sample_files=True
run_sample_analysis=True

if not os.path.exists("output"):
    os.makedirs("output")
if not os.path.exists("output/beast"):
    os.makedirs("output/beast")
if not os.path.exists("output/c_beast"):
    os.makedirs("output/c_beast")  
if not os.path.exists("output/phyrex"):
    os.makedirs("output/phyrex")   
    
if not os.path.exists("output/beast/no_sampling"):
    os.makedirs("output/beast/no_sampling")    
if not os.path.exists("output/beast/no_sampling/beast_input"):
    os.makedirs("output/beast/no_sampling/beast_input")
if not os.path.exists("output/beast/no_sampling/beast_output"):
    os.makedirs("output/beast/no_sampling/beast_output")
if not os.path.exists("output/beast/no_sampling/generated_trees"):
    os.makedirs("output/beast/no_sampling/generated_trees")
if not os.path.exists("output/beast/no_sampling/annotated_trees"):
    os.makedirs("output/beast/no_sampling/annotated_trees")
    


for output_index in range(1, num_sampling+1):
    if not os.path.exists("output/beast/sampled"+str(output_index)):
        os.makedirs("output/beast/sampled"+str(output_index))   
    if not os.path.exists("output/beast/sampled"+str(output_index)+"/beast_input"):
        os.makedirs("output/beast/sampled"+str(output_index)+"/beast_input")
    if not os.path.exists("output/beast/sampled"+str(output_index)+"/beast_output"):
        os.makedirs("output/beast/sampled"+str(output_index)+"/beast_output")
    if not os.path.exists("output/beast/sampled"+str(output_index)+"/generated_trees"):
        os.makedirs("output/beast/sampled"+str(output_index)+"/generated_trees")
    if not os.path.exists("output/beast/sampled"+str(output_index)+"/annotated_trees"):
        os.makedirs("output/beast/sampled"+str(output_index)+"/annotated_trees")
    if not os.path.exists("output/beast/sampled"+str(output_index)+"/root_data"):
        os.makedirs("output/beast/sampled"+str(output_index)+"/root_data")    
        
    if not os.path.exists("output/c_beast/sampled"+str(output_index)):
        os.makedirs("output/c_beast/sampled"+str(output_index))   
    if not os.path.exists("output/c_beast/sampled"+str(output_index)+"/beast_input"):
        os.makedirs("output/c_beast/sampled"+str(output_index)+"/beast_input")
    if not os.path.exists("output/c_beast/sampled"+str(output_index)+"/beast_output"):
        os.makedirs("output/c_beast/sampled"+str(output_index)+"/beast_output")
    if not os.path.exists("output/c_beast/sampled"+str(output_index)+"/generated_trees"):
        os.makedirs("output/c_beast/sampled"+str(output_index)+"/generated_trees")
    if not os.path.exists("output/c_beast/sampled"+str(output_index)+"/annotated_trees"):
        os.makedirs("output/c_beast/sampled"+str(output_index)+"/annotated_trees")
    if not os.path.exists("output/c_beast/sampled"+str(output_index)+"/root_data"):
        os.makedirs("output/c_beast/sampled"+str(output_index)+"/root_data") 

    
    if not os.path.exists("output/phyrex/sampled"+str(output_index)):
        os.makedirs("output/phyrex/sampled"+str(output_index))        
    if not os.path.exists("output/phyrex/sampled"+str(output_index)+"/phyrex_output"):
        os.makedirs("output/phyrex/sampled"+str(output_index)+"/phyrex_output")
    if not os.path.exists("output/phyrex/sampled"+str(output_index)+"/phyrex_input"):
        os.makedirs("output/phyrex/sampled"+str(output_index)+"/phyrex_input")

for i in range(num_trees*(job_index-1), num_trees*job_index):
    
    #print("index is" + str(i))
    
    t=dendropy.Tree()
    
    
    if args.tree_type == "nuc":
        t = treegenerator.generate_nonultrametric_coalescent_tree(num_tips_per_period, num_periods, period_length, lamb)
    elif args.tree_type == "uc":
        t = treegenerator.generate_ultrametric_coalescent_tree(num_tips, lamb)
    elif args.tree_type == "bd":
        t = treegenerator.generate_birthdeath_tree(1, 0.1, num_tips)
    elif args.tree_type == "yule":
        t = treegenerator.generate_yule_tree(num_tips)
    elif args.tree_type == "star":
        t = treegenerator.generate_star_tree()
    else:
        print(args.tree_type+" is invalid tree type")
        break
                
    t= treegenerator.simulate_brownian(t, sigma, dimension) 
    t= treegenerator.calculate_time_to_tips(t)
    
    max_coordinate = 0
    
    for node in t.preorder_node_iter():
        max_coordinate = max(max_coordinate, abs(node.X))
        if dimension == 2:
            max_coordinate = max(max_coordinate, abs(node.Y))
     
    beastxmlwriter.write_BEAST_xml(t, i, dimension, mcmc, log_every, beast_input_string="output/beast/no_sampling/beast_input/beast", beast_output_string="output/beast/no_sampling/beast_output/beast")    
    
    if generate_sample_files:    
        
        for output_index in range(1, num_sampling+1):
            if output_index ==1:
                sampled_t=sampling.sample_unbiased(t, dimension, sample_ratio=sample_ratio)
            elif output_index==2:
                sampled_t=sampling.sample_biased_most_central(t, dimension, sample_ratio=sample_ratio)
            elif output_index==3:
                sampled_t=sampling.sample_biased_diagonal(t, dimension, sample_ratio=sample_ratio)            
            elif output_index==4:
                sampled_t=sampling.sample_biased_extreme(t, dimension, sample_ratio=sample_ratio)
            sampled_t=treegenerator.calculate_time_to_tips(sampled_t)
            beastxmlwriter.write_BEAST_xml(sampled_t, i, dimension, mcmc, log_every, "output/beast/sampled"+str(output_index)+"/beast_input/beast", beast_output_string="output/beast/sampled"+str(output_index)+"/beast_output/beast")
            for node in sampled_t.preorder_node_iter():
                node.annotations.add_bound_attribute("time")
                node.annotations.add_bound_attribute("X")
                node.annotations.add_bound_attribute("time_to_tips")
                if dimension==2:
                    node.annotations.add_bound_attribute("Y")
                    
                    
            ######
            if generate_corrected_files:
                d = dendropy.model.discrete.hky85_chars(kappa=3, mutation_rate=0.05, seq_len=seq_len, tree_model=sampled_t, retain_sequences_on_tree=False)    
                beastxmlwriter.write_BEAST_xml_corrected(t, sampled_t, d, i=i,  mcmc=corr_beast_mcmc, log_every=1000, beast_input_string ="output/c_beast/sampled"+str(output_index)+"/beast_input/beast", beast_output_string="output/c_beast/sampled"+str(output_index)+"/beast_output/beast", other_sample_size=other_sample_size, seq_len=seq_len)
            ######
            
            
            sampled_t.write(path="output/beast/sampled"+str(output_index)+"/generated_trees/tree"+str(i)+".txt", schema="nexus", suppress_internal_taxon_labels=True)
            if dimension ==2:
                phyrexxmlwriter.write_phyrex_input(sampled_t, i, input_string="output/phyrex/sampled"+str(output_index)+"/phyrex_input/" , output_string="output/phyrex/sampled"+str(output_index)+"/phyrex_output/", bound=2*max_coordinate) 
            

            file = open("output/beast/sampled"+str(output_index)+"/root_data/actual_root"+str(i)+".txt", "w")
            file.write(str(sampled_t.seed_node.X)+'\n')
            file.write(str(sampled_t.seed_node.Y)+'\n')
            file.close()
                
            
            
            
            if generate_corrected_files:
                file = open("output/c_beast/sampled"+str(output_index)+"/root_data/actual_root"+str(i)+".txt", "w")
                file.write(str(sampled_t.seed_node.X)+'\n')
                file.write(str(sampled_t.seed_node.Y)+'\n')
                file.close()
                
                file = open("output/c_beast/sampled"+str(output_index)+"/root_data/labels"+str(i)+".txt", "w")
                for leaf in sampled_t.leaf_node_iter():
                    file.write(leaf.taxon.label+'\n')
                file.close()
        
    for node in t.preorder_node_iter():
        node.annotations.add_bound_attribute("time")
        node.annotations.add_bound_attribute("X")
        node.annotations.add_bound_attribute("time_to_tips")
        if dimension==2:
            node.annotations.add_bound_attribute("Y")
    t.write(path="output/beast/no_sampling/generated_trees/tree"+str(i)+".txt", schema="nexus", suppress_internal_taxon_labels=True)
    
    run_analysis= False

    burnin=int(mcmc/10)
    if run_analysis:
        if linux:
            os.system('beast -overwrite -seed 123456795 "output/beast/no_sampling/beast_input/beast'+str(i)+'.xml"')
        else:
            os.system('cmd /c java -jar beast.jar -overwrite -seed 123456795 "output/beast/no_sampling/beast_input/beast'+str(i)+'.xml"')
            
            
        if run_tree_annotator:
            os.system('treeannotator -burnin '+str(burnin)+' "output/beast/no_sampling/beast_output/beast'+str(i)+'.trees.txt" "output/beast/no_sampling/annotated_trees/beast'+str(i)+'.tree.txt"')
        
    
    if run_sample_analysis:
        for output_index in range(1, num_sampling+1):
            if linux:
                os.system('beast -overwrite -seed 123456795 "output/beast/sampled'+str(output_index)+'/beast_input/beast'+str(i)+'.xml"')
            else:
                os.system('cmd /c java -jar beast.jar -overwrite -seed 123456795 "output/beast/sampled'+str(output_index)+'/beast_input/beast'+str(i)+'.xml"')
            if run_tree_annotator:
                #os.system('cmd /c ""C:/Users/Antanas/Desktop/BEAST v1.10.4/bin/treeannotator" -burnin '+str(burnin)+' "C:/Users/Antanas/Phylogeny Simulation/output/sampled_beast_output'+str(output_index)+'/sampled_beast'+str(i)+'.trees.txt" "C:/Users/Antanas/Phylogeny Simulation/output/annotated_sampled_trees'+str(output_index)+'/sampled_beast'+str(i)+'.tree.txt""')
                os.system('treeannotator -burnin '+str(burnin)+' "output/beast/sampled'+str(output_index)+'beast_output/beast'+str(i)+'.trees.txt" "output/beast/sampled'+str(output_index)+'/annotated_trees/beast'+str(i)+'.tree.txt"')
               
            file = open("output/beast/sampled"+str(output_index)+"/root_data/observed_roots"+str(i)+".txt", "w")   
            for line in open("output/beast/sampled"+str(output_index)+"/beast_output/beast"+str(i)+".trees.txt"):
                
                if line.startswith("tree"):
                    start_index = 0
                    while True:
                        if line[start_index:start_index+4]=="[&R]":
                            break
                        start_index=start_index+1
                    single_tree=dendropy.Tree.get(data=line[start_index:], schema="newick", extract_comment_metadata=True)
                    file.write(single_tree.seed_node.annotations.require_value("location")[0]+"\t"+single_tree.seed_node.annotations.require_value("location")[1]+'\n')
            file.close()
        
        
        #def analyze_tree_list(tree, i, dimension, mcmc): 
##   command to execute treeannotator       
#    burnin=int(mcmc/10)
#    os.system('cmd /c ""C:/Users/Antanas/Desktop/BEAST v1.10.4/bin/treeannotator" -burnin '+str(burnin)+' "C:/Users/Antanas/Phylogeny Simulation/output/beast_output/beast'+str(i)+'.trees.txt" "C:/Users/Antanas/Phylogeny Simulation/annotated_trees/beast'+str(i)+'.tree.txt""')
#    for node in tree.preorder_node_iter():
#        node.averageX=0
#        if dimension == 2:
#            node.averageY=0
#        
#    annotated_tree=dendropy.Tree.get(path='output/annotated_trees/beast'+str(i)+'.tree.txt', taxon_namespace=tree.taxon_namespace, extract_comment_metadata=True, suppress_internal_node_taxa=False, schema="nexus")     
#           
#    nodes_tree = [nd for nd in tree.preorder_node_iter()]
#    nodes_annotated_tree = [nd for nd in annotated_tree.preorder_node_iter()]
#    for i, n in enumerate(nodes_tree):            
#        if dimension==2:
#            nodes_tree[i].averageX=float(nodes_annotated_tree[i].annotations.require_value("location1"))
#            nodes_tree[i].averageY=float(nodes_annotated_tree[i].annotations.require_value("location2"))
#        else:
#            nodes_tree[i].averageX=float(nodes_annotated_tree[i].annotations.require_value("X"))
#            
#    return tree                  

