# -*- coding: utf-8 -*-
"""
@author: Nicola De Maio
"""

#Run extra simulations to investigate intermediate levels of bias and extreme levels of corrected samples.

import argparse
import dendropy
#from dendropy import treecalc
import os
#from Bio import Phylo
#from cStringIO import StringIO
#import itertools
import math
import numpy as np

import sampling
import beastxmlwriter
import phyrexxmlwriter
import treegenerator


#argparses file to extract arguments from command line
parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-dims', dest='dimension', type=int, default=2, help='number of dimensions (1 or 2) for which the random walk is generated (default: 2)')
parser.add_argument('-N', action="store", type=int, dest="num_trees", default=1, help='number of simulations (default 5)')
parser.add_argument('-mcmc', action="store", type=int, dest="mcmc", default=10000, help='number of mcmc steps (default 5000)')
parser.add_argument('-logevery', action="store", type=int, dest="log_every", default=10, help='log every these many steps (default 10)')
parser.add_argument('-lambda', action="store", type=float, dest="lamb", default=1, help='rate of coalescent for coalescent trees (default 1)')
parser.add_argument('-sigma', action="store", type=float, dest="sigma", default=1, help='standart deviation for brownian diffusion simulation for a branch of length 1 (default 1)')
parser.add_argument('-jobi', action="store", type=int, dest="job_index", default=1, help='job index')
parser.add_argument('-sample_ratio', dest='sample_ratio', action="store", type=float, default=0.05, help='ratio of samples (default: 0.05)')
parser.add_argument('-seq_len', action="store", type=int, dest="seq_len", default=10000, help='alignment length for corrected beast (default 10000)')
parser.add_argument('--c_beast', dest='c_beast', action='store_const', const=True, default=False, help='should the files be generated for corrected BEAST? (default: False)')
parser.add_argument('--re_run', dest='re_run', action='store_const', const=True, default=False, help='is this a re-run, so we should not generate new trees but only run beast again on the previous simulations that did not complete? (default: False)')
parser.add_argument('--corr', dest='corr', action='store_const', const=True, default=False, help='just calculate genetic-geographic distance correlations.')
parser.add_argument('-br', action="store", type=float, dest="br", default=1, help='birth rate for birth-death or Yule trees (default 1)')
parser.add_argument('-extraSamples', action='store_const', const=True, default=False, help='run simulations of extra samples?')
parser.add_argument('-gradient', action="store", type=int, default=1, help='intensity of the gradient for simulations')
parser.add_argument('-biasIntensity', action='store_const', const=True, default=False, help='run simulations of extra samples?')
parser.add_argument('-gradientNum', action="store", type=int, default=5, help='total number of intensity of the gradient for simulations')

args = parser.parse_args()
dimension = args.dimension
num_trees = args.num_trees
mcmc = args.mcmc
log_every = args.log_every
sigma = args.sigma
extraSamples=args.extraSamples
biasIntensity=args.biasIntensity
gradient=args.gradient

br = args.br
lamb=args.lamb
job_index = args.job_index

#ratio of samples
sample_ratio=args.sample_ratio

#variables for corrected BEAST
#other_sample_size is the number of additional trees nodes that are added
other_sample_size = 50

#seq_len is the length of generated alignment
seq_len = args.seq_len

#mcmc chain length for corrected BEAST
corr_beast_mcmc=int(1e7)


#Levels of intermediate extreme sampling
num_sampling = args.gradientNum

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


for output_index in range(1, num_sampling+1):
    if not os.path.exists("output/beast/sampled_gradient"+str(output_index)):
        os.makedirs("output/beast/sampled_gradient"+str(output_index))   
    if not os.path.exists("output/beast/sampled_gradient"+str(output_index)+"/beast_input"):
        os.makedirs("output/beast/sampled_gradient"+str(output_index)+"/beast_input")
    if not os.path.exists("output/beast/sampled_gradient"+str(output_index)+"/beast_output"):
        os.makedirs("output/beast/sampled_gradient"+str(output_index)+"/beast_output")
    if not os.path.exists("output/beast/sampled_gradient"+str(output_index)+"/generated_trees"):
        os.makedirs("output/beast/sampled_gradient"+str(output_index)+"/generated_trees")
    if not os.path.exists("output/beast/sampled_gradient"+str(output_index)+"/annotated_trees"):
        os.makedirs("output/beast/sampled_gradient"+str(output_index)+"/annotated_trees")
    if not os.path.exists("output/beast/sampled_gradient"+str(output_index)+"/root_data"):
        os.makedirs("output/beast/sampled_gradient"+str(output_index)+"/root_data")    
        
    if not os.path.exists("output/c_beast/sampled_gradient"+str(output_index)):
        os.makedirs("output/c_beast/sampled_gradient"+str(output_index))   
    if not os.path.exists("output/c_beast/sampled_gradient"+str(output_index)+"/beast_input"):
        os.makedirs("output/c_beast/sampled_gradient"+str(output_index)+"/beast_input")
    if not os.path.exists("output/c_beast/sampled_gradient"+str(output_index)+"/beast_output"):
        os.makedirs("output/c_beast/sampled_gradient"+str(output_index)+"/beast_output")
    if not os.path.exists("output/c_beast/sampled_gradient"+str(output_index)+"/generated_trees"):
        os.makedirs("output/c_beast/sampled_gradient"+str(output_index)+"/generated_trees")
    if not os.path.exists("output/c_beast/sampled_gradient"+str(output_index)+"/annotated_trees"):
        os.makedirs("output/c_beast/sampled_gradient"+str(output_index)+"/annotated_trees")
    if not os.path.exists("output/c_beast/sampled_gradient"+str(output_index)+"/root_data"):
        os.makedirs("output/c_beast/sampled_gradient"+str(output_index)+"/root_data") 


#loop is useful if we want to split a number of simulations into several jobs
#the program runs num_trees of simulations
for i in range(num_trees*(job_index), num_trees*(job_index+1)):
	print("new tree")
	#t = treegenerator.generate_yule_tree(num_tips, br)

	output_index=gradient
	if extraSamples:
		nTips2=100
		t2=treegenerator.generate_yule_tree(nTips2, br)
		t2= treegenerator.simulate_brownian(t2, sigma, dimension) 
		step2=1.0/(num_sampling)
		mutation_rate =0.01
		other_sample_size=int(step2*(output_index-1)*nTips2)
		print("Tree generated")
		print(other_sample_size)
		print(step2)
		print("")
		sampled_t2=sampling.sample_biased_side(t2, dimension=2, sample_ratio=step2, ratioExtreme=1.0)
		print("Tree sampled")
		d = dendropy.model.discrete.hky85_chars(kappa=3, mutation_rate=mutation_rate, seq_len=seq_len, tree_model=sampled_t2, retain_sequences_on_tree=False)    
		beastxmlwriter.write_BEAST_xml_corrected(t2, sampled_t2, d, i=i,  mcmc=corr_beast_mcmc*2, log_every=2000, log_every_tree=20000, beast_input_string ="output/c_beast/sampled_gradient"+str(output_index)+"/beast_input/beast", beast_output_string="output/c_beast/sampled_gradient"+str(output_index)+"/beast_output/beast", other_sample_size=other_sample_size, seq_len=seq_len)
		sampled_t2.write(path="output/c_beast/sampled_gradient"+str(output_index)+"/generated_trees/tree"+str(i)+".txt", schema="nexus", suppress_internal_taxon_labels=True)
		file = open("output/c_beast/sampled_gradient"+str(output_index)+"/root_data/actual_root"+str(i)+".txt", "w")
		file.write(str(sampled_t2.seed_node.X)+'\n')
		file.write(str(sampled_t2.seed_node.Y)+'\n')
		file.close()
		file = open("output/c_beast/sampled_gradient"+str(output_index)+"/root_data/labels"+str(i)+".txt", "w")
		for leaf in sampled_t2.leaf_node_iter():
			file.write(leaf.taxon.label+'\n')
		file.close()
		os.system('ulimit -c unlimited; beast -overwrite -seed 1 "output/c_beast/sampled_gradient'+str(output_index)+'/beast_input/beast'+str(i)+'.xml"')
		file = open("output/c_beast/sampled_gradient"+str(output_index)+"/root_data/observed_roots"+str(i)+".txt", "w")   
		for line in open("output/c_beast/sampled_gradient"+str(output_index)+"/beast_output/beast"+str(i)+".trees.txt"):                
			if line.startswith("tree"):
				start_index = 0
				while True:
					if line[start_index:start_index+4]=="[&R]":
						break
					start_index=start_index+1
				single_tree=dendropy.Tree.get(data=line[start_index:], schema="newick", extract_comment_metadata=True)
				file.write(single_tree.seed_node.annotations.require_value("location")[0]+"\t"+single_tree.seed_node.annotations.require_value("location")[1]+'\n')
		file.close()
	
	if biasIntensity:
		nTips1=10000
		t1=treegenerator.generate_yule_tree(nTips1, br)
		t1= treegenerator.simulate_brownian(t1, sigma, dimension) 
		step=1.0/(num_sampling-1)
		print("Tree generated")
		print(step*(output_index-1))
		print(step)
		print("")
		sampled_t=sampling.sample_biased_side(t1, dimension=2, sample_ratio=0.01, ratioExtreme=step*(output_index-1))
		beastxmlwriter.write_BEAST_xml(sampled_t, i, dimension, mcmc, log_every, "output/beast/sampled_gradient"+str(output_index)+"/beast_input/beast", beast_output_string="output/beast/sampled_gradient"+str(output_index)+"/beast_output/beast")
		for node in sampled_t.preorder_node_iter():
			node.annotations.add_bound_attribute("time")
			node.annotations.add_bound_attribute("X")
			if dimension==2:
				node.annotations.add_bound_attribute("Y")
		sampled_t.write(path="output/beast/sampled_gradient"+str(output_index)+"/generated_trees/tree"+str(i)+".txt", schema="nexus", suppress_internal_taxon_labels=True)
		file = open("output/beast/sampled_gradient"+str(output_index)+"/root_data/actual_root"+str(i)+".txt", "w")
		file.write(str(sampled_t.seed_node.X)+'\n')
		file.write(str(sampled_t.seed_node.Y)+'\n')
		print("Actual root location")
		print(sampled_t.seed_node.X)
		print(sampled_t.seed_node.Y)
		file.close()
		burnin=int(mcmc/10)
		os.system('ulimit -c unlimited; beast -overwrite -seed 1 "output/beast/sampled_gradient'+str(output_index)+'/beast_input/beast'+str(i)+'.xml"')
		file = open("output/beast/sampled_gradient"+str(output_index)+"/root_data/observed_roots"+str(i)+".txt", "w") 
		Xs=[]
		Ys=[]
		for line in open("output/beast/sampled_gradient"+str(output_index)+"/beast_output/beast"+str(i)+".trees.txt"):                
			if line.startswith("tree"):
				start_index = 0
				while True:
					if line[start_index:start_index+4]=="[&R]":
						break
					start_index=start_index+1
				single_tree=dendropy.Tree.get(data=line[start_index:], schema="newick", extract_comment_metadata=True)
				Xroot=single_tree.seed_node.annotations.require_value("location")[0]
				Xs.append(float(Xroot))
				Yroot=single_tree.seed_node.annotations.require_value("location")[1]
				Ys.append(float(Yroot))
				file.write(Xroot+"\t"+Yroot+'\n')
		file.close()
		print(np.mean(Xs))
		print(np.std(Xs))
		print(np.mean(Ys))
		print(np.std(Ys))

exit()