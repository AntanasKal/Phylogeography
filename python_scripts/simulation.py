# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:37:34 2019

@author: Antanas
"""

import argparse
import dendropy
#from dendropy import treecalc
import os
from Bio import Phylo
from cStringIO import StringIO
import itertools
import math
import numpy as np

import sampling
import beastxmlwriter
import phyrexxmlwriter
import treegenerator


#argparses file to extract arguments from command line
parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-dims', dest='dimension', type=int, default=2, help='number of dimensions (1 or 2) for which the random walk is generated (default: 2)')
parser.add_argument('-N', action="store", type=int, dest="num_trees", default=5, help='number of simulations (default 5)')
parser.add_argument('-treetype', action="store", dest="tree_type", default='yule', help='type of tree generated \n "nuc" - nonultrametric coalescent \n "uc" - ultrametric coalescent \n "bd" - birth-death tree \n "yule" - yule tree')
parser.add_argument('-mcmc', action="store", type=int, dest="mcmc", default=5000, help='number of mcmc steps (default 5000)')
parser.add_argument('-logevery', action="store", type=int, dest="log_every", default=10, help='log every these many steps (default 10)')
parser.add_argument('-ntips', action="store", type=int, dest="ntips", default=100, help='number of tips for ultrametric coalescent tree (default 100)')
parser.add_argument('-ntipspp', action="store", type=int, dest="num_tips_per_period", default=20, help='number of tips for per period for nonultrametric coalescent tree (default 20)')
parser.add_argument('-nps', action="store", type=int, dest="num_periods", default=10, help='number of periods for non-ultrametric coalescent tree (default 10)')
parser.add_argument('-lambda', action="store", type=float, dest="lamb", default=1, help='rate of coalescent for coalescent trees (default 1)')
parser.add_argument('-pl', action="store", type=float, dest="period_length", default=1, help='period length for nonultrametric coalescent tree (default 1)')
parser.add_argument('-sigma', action="store", type=float, dest="sigma", default=1, help='standart deviation for brownian diffusion simulation for a branch of length 1 (default 1)')
parser.add_argument('-jobi', action="store", type=int, dest="job_index", default=1, help='job index')
parser.add_argument('--linux', dest='linux', action='store_const', const=True, default=False, help='is the program run on Linux? (default: False)')
parser.add_argument('--annotate', dest='run_tree_annotator', action='store_const', const=True, default=False, help='run tree annotator (default: False)')
parser.add_argument('-sample_ratio', dest='sample_ratio', action="store", type=float, default=0.05, help='ratio of samples (default: 0.05)')
parser.add_argument('-seq_len', action="store", type=int, dest="seq_len", default=10000, help='alignment length for corrected beast (default 10000)')
parser.add_argument('--c_beast', dest='c_beast', action='store_const', const=True, default=False, help='should the files be generated for corrected BEAST? (default: False)')
parser.add_argument('--re_run', dest='re_run', action='store_const', const=True, default=False, help='is this a re-run, so we should not generate new trees but only run beast again on the previous simulations that did not complete? (default: False)')
parser.add_argument('--corr', dest='corr', action='store_const', const=True, default=False, help='just calculate genetic-geographic distance correlations.')

parser.add_argument('-br', action="store", type=float, dest="br", default=1, help='birth rate for birth-death or Yule trees (default 1)')
parser.add_argument('-dr', action="store", type=float, dest="dr", default=1, help='death rate for birth-death trees (default 0.1)')

args = parser.parse_args()
dimension = args.dimension
num_trees = args.num_trees
mcmc = args.mcmc
log_every = args.log_every
sigma = args.sigma
num_tips = args.ntips
num_tips_per_period = args.num_tips_per_period
num_periods=args.num_periods

br = args.br
dr = args.dr
lamb=args.lamb
period_length=args.period_length
job_index = args.job_index


linux = args.linux

run_tree_annotator = args.run_tree_annotator

#ratio of samples
sample_ratio=args.sample_ratio


#variables for corrected BEAST
#other_sample_size is the number of additional trees nodes that are added
other_sample_size = 50

#seq_len is the length of generated alignment
seq_len = args.seq_len

#mcmc chain length for corrected BEAST
corr_beast_mcmc=int(1e7)



#number of sampling scenarios (set to 4).
# scenario 1 - unbiased sampling
# scenario 2 - central sampling
# scenario 3 - diagonal sampling
# scenario 4 - one-sided sampling
num_sampling = 4


##Bellow are several variables on whether to run particular inference or not.
##useful if we are interested in running some particular parts of inference, but not every part

#Should I only run beast again for the cases that did not run before?
reRun_beast_only=args.re_run

#Should I just calculate correlations?
correlations_only=args.corr

#boolean variable for whether to generate files for sampled scenarios
generate_sample_files=True

if reRun_beast_only:
	generate_sample_files=False

#boolean variable for whether to generate files for corrected beast
generate_corrected_files=args.c_beast


#boolean variable for whether to run BEAST on sampled scenarios
run_sample_analysis=True

#boolean variable for whether to run BEAST on full sampling
run_analysis= False


#code to create needed folders for output
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

corrs=[[],[],[],[]]
#loop is useful if we want to split a number of simulations into several jobs
#the program runs num_trees of simulations
for i in range(num_trees*(job_index), num_trees*(job_index+1)):
    print("new tree")
    
    #generating a tree according to the parameteres given from the command line
    if args.tree_type == "nuc":
        t = treegenerator.generate_nonultrametric_coalescent_tree(num_tips_per_period, num_periods, period_length, lamb)
    elif args.tree_type == "uc":
        t = treegenerator.generate_ultrametric_coalescent_tree(num_tips, lamb)
    elif args.tree_type == "bd":
        t = treegenerator.generate_birthdeath_tree(num_tips, br, dr)
    elif args.tree_type == "yule":
        t = treegenerator.generate_yule_tree(num_tips, br)
    elif args.tree_type == "star":
        t = treegenerator.generate_star_tree()
    else:
        print(args.tree_type+" is invalid tree type")
        break
    
    #simulating brownian motion along entire tree
    t= treegenerator.simulate_brownian(t, sigma, dimension) 
    
    #calculating max coordinate to give the bounds of habitat for phyrex
    max_coordinate = 0    
    for node in t.preorder_node_iter():
        max_coordinate = max(max_coordinate, abs(node.X))
        if dimension == 2:
            max_coordinate = max(max_coordinate, abs(node.Y))
    
    #line of code for writing xml file for full sampling (this is not needed necessary)
    beastxmlwriter.write_BEAST_xml(t, i, dimension, mcmc, log_every, beast_input_string="output/beast/no_sampling/beast_input/beast", beast_output_string="output/beast/no_sampling/beast_output/beast")    
    
    
    
    
    if generate_sample_files or (not os.path.exists("output/beast/sampled1/beast_input/beast"+str(i)+".xml")):    
        
        for output_index in range(1, num_sampling+1):
            if output_index ==1:
                sampled_t=sampling.sample_unbiased(t, dimension, sample_ratio=sample_ratio)
            elif output_index==2:
                sampled_t=sampling.sample_biased_most_central(t, dimension, sample_ratio=sample_ratio)
            elif output_index==3:
                sampled_t=sampling.sample_biased_diagonal(t, dimension, sample_ratio=sample_ratio)            
            elif output_index==4:
                sampled_t=sampling.sample_biased_extreme(t, dimension, sample_ratio=sample_ratio)
            
            if correlations_only:
            	gen_dists=[]
            	geo_dists=[]
            	newickT=sampled_t.as_string(schema="newick",suppress_rooting=True)
            	tree = Phylo.read(StringIO(newickT), "newick")
            	d = {}
            	for x, y in itertools.combinations(tree.get_terminals(), 2):
            		v = tree.distance(x, y)
            		d[x.name] = d.get(x.name, {})
            		d[x.name][y.name] = v
            		d[y.name] = d.get(y.name, {})
            		d[y.name][x.name] = v
            	for x in tree.get_terminals():
            		d[x.name][x.name] = 0
            	for leaf in sampled_t.leaf_node_iter():
            		for leaf2 in sampled_t.leaf_node_iter():
            			if leaf.taxon.label>leaf2.taxon.label:
            				gen_dists.append(d[leaf.taxon.label][leaf2.taxon.label])
            				geo_dists.append(math.sqrt(((leaf.X-leaf2.X)**2) + ((leaf.Y-leaf2.Y)**2)))
            	print(np.corrcoef(gen_dists,geo_dists)[0][1])
            	corrs[output_index-1].append(np.corrcoef(gen_dists,geo_dists)[0][1])
            	continue
            
            beastxmlwriter.write_BEAST_xml(sampled_t, i, dimension, mcmc, log_every, "output/beast/sampled"+str(output_index)+"/beast_input/beast", beast_output_string="output/beast/sampled"+str(output_index)+"/beast_output/beast")
            for node in sampled_t.preorder_node_iter():
                node.annotations.add_bound_attribute("time")
                node.annotations.add_bound_attribute("X")
                if dimension==2:
                    node.annotations.add_bound_attribute("Y")
                    
                    
            #
            #mutation rate for corrected beast 
            mutation_rate =0.01
            
            if generate_corrected_files:
                #using dendropy implementation for generating alignments
                d = dendropy.model.discrete.hky85_chars(kappa=3, mutation_rate=mutation_rate, seq_len=seq_len, tree_model=sampled_t, retain_sequences_on_tree=False)    
                beastxmlwriter.write_BEAST_xml_corrected(t, sampled_t, d, i=i,  mcmc=corr_beast_mcmc, log_every=1000, log_every_tree=10000, beast_input_string ="output/c_beast/sampled"+str(output_index)+"/beast_input/beast", beast_output_string="output/c_beast/sampled"+str(output_index)+"/beast_output/beast", other_sample_size=other_sample_size, seq_len=seq_len)
            ######
            
            #writing the full tree with annotations containig locations and times of the nodes
            #(currently not necessary but might be useful to check how, for instance internal nodes are inferred)
            sampled_t.write(path="output/beast/sampled"+str(output_index)+"/generated_trees/tree"+str(i)+".txt", schema="nexus", suppress_internal_taxon_labels=True)
            if dimension ==2:
                phyrexxmlwriter.write_phyrex_input(sampled_t, i, input_string="output/phyrex/sampled"+str(output_index)+"/phyrex_input/" , output_string="output/phyrex/sampled"+str(output_index)+"/phyrex_output/", bound=2*max_coordinate) 
            
            #writing the actual root data in convenient form
            file = open("output/beast/sampled"+str(output_index)+"/root_data/actual_root"+str(i)+".txt", "w")
            file.write(str(sampled_t.seed_node.X)+'\n')
            file.write(str(sampled_t.seed_node.Y)+'\n')
            file.close()   
            
            
            #for corrected BEAST we write the real root locations
            #and write the labels of the taxa in the subsampled tree
            #in order to keep track when inference is launched for corrected BEAST which nodes are actual in the tree
            if generate_corrected_files:
                file = open("output/c_beast/sampled"+str(output_index)+"/root_data/actual_root"+str(i)+".txt", "w")
                file.write(str(sampled_t.seed_node.X)+'\n')
                file.write(str(sampled_t.seed_node.Y)+'\n')
                file.close()
                
                file = open("output/c_beast/sampled"+str(output_index)+"/root_data/labels"+str(i)+".txt", "w")
                for leaf in sampled_t.leaf_node_iter():
                    file.write(leaf.taxon.label+'\n')
                file.close()
        
        if correlations_only:
        	continue
                
        #code for writing the full tree in nexus format (this is not necessary now, but might be useful)
        for node in t.preorder_node_iter():
        	node.annotations.add_bound_attribute("time")
        	node.annotations.add_bound_attribute("X")
        	if dimension==2:
        		node.annotations.add_bound_attribute("Y")
        t.write(path="output/beast/no_sampling/generated_trees/tree"+str(i)+".txt", schema="nexus", suppress_internal_taxon_labels=True)   
    
    
    
    
    
    
    #now that all input files have been generated, we run BEAST inference
    burnin=int(mcmc/10)
    if run_sample_analysis:
        #we cycle over the 4 sampling scenarios
        for output_index in range(1, num_sampling+1):
            if linux:
            	#NICOLA: changed here so that only analyses that have not already done are now executed.
            	#NICOLA: also I have added "ulimit -c unlimited" so that maybe the core dumping problem is solved?
            	if (not reRun_beast_only) or (not os.path.exists("output/beast/sampled"+str(output_index)+"/beast_output/beast"+str(i)+".log.txt")):
            		#os.system('ulimit -c unlimited; beast -overwrite -seed 123456795 "output/c_beast/sampled'+str(sample_index)+'/beast_input/beast'+str(index)+'.xml"')
            		os.system('ulimit -c unlimited; beast -overwrite -seed 123456795 "output/beast/sampled'+str(output_index)+'/beast_input/beast'+str(i)+'.xml"')
                	
            else:
                os.system('cmd /c java -jar beast.jar -overwrite -seed 123456795 "output/beast/sampled'+str(output_index)+'/beast_input/beast'+str(i)+'.xml"')
                
            #code for running tree annotator (unneeded)
            if run_tree_annotator:                
                os.system('treeannotator -burnin '+str(burnin)+' "output/beast/sampled'+str(output_index)+'beast_output/beast'+str(i)+'.trees.txt" "output/beast/sampled'+str(output_index)+'/annotated_trees/beast'+str(i)+'.tree.txt"')            
            
            
            #print root positions from .trees.txt file to more convenient format
            #only the locations of the roots are important from the .trees.txt file
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
            

        
        
    if run_analysis:
        if linux:
            os.system('beast -overwrite -seed 123456795 "output/beast/no_sampling/beast_input/beast'+str(i)+'.xml"')
        else:
            os.system('cmd /c java -jar beast.jar -overwrite -seed 123456795 "output/beast/no_sampling/beast_input/beast'+str(i)+'.xml"')           
        if run_tree_annotator:
            os.system('treeannotator -burnin '+str(burnin)+' "output/beast/no_sampling/beast_output/beast'+str(i)+'.trees.txt" "output/beast/no_sampling/annotated_trees/beast'+str(i)+'.tree.txt"')
        
if correlations_only:
	for output_index in range(num_sampling):
		print("correlation: "+str(np.average(corrs[output_index])))
exit()