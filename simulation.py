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

#this function generates borwnian motion for a given tree
def simulate_brownian(t, sigma, dimension):
    #t is the tree
    #sigma is the standard deviation of the brownian motion
    #dimension is the number of dimensions in which we generate the random walk
    
    positions = {}
    
    for node in t.preorder_node_iter():   
        
        if node.parent_node is None:
            node.X = 0
            node.displacementX = 0
            if dimension==2:
                node.Y = float (0)
                node.displacementY = float(0)
        else:
 #           node.displacementX = random.gauss(mu=0, sigma=sigma*math.sqrt(node.edge_length))            
            node.displacementX = random.gauss(0, sigma*math.sqrt(node.edge_length))
            #node.displacementX = np.random.normal()*math.sqrt(node.edge_length)
            
            
            node.X = node.parent_node.X+node.displacementX      
            #node.X = random.gauss(node.parent_node.X, sigma*math.sqrt(node.edge_length))   
            if dimension==2:
                node.displacementY = random.gauss(0, sigma*math.sqrt(node.edge_length))
                node.Y = node.parent_node.Y+node.displacementY   
        positions.update({node.taxon.label: node.X})
    return t

    
def analyze_tree_list(tree, i, dimension, mcmc): 
#   command to execute treeannotator       
    burnin=int(mcmc/10)
    os.system('cmd /c ""C:/Users/Antanas/Desktop/BEAST v1.10.4/bin/treeannotator" -burnin '+str(burnin)+' "C:/Users/Antanas/Phylogeny Simulation/output/beast_output/beast'+str(i)+'.trees.txt" "C:/Users/Antanas/Phylogeny Simulation/annotated_trees/beast'+str(i)+'.tree.txt""')
    for node in tree.preorder_node_iter():
        node.averageX=0
        if dimension == 2:
            node.averageY=0
        
    annotated_tree=dendropy.Tree.get(path='output/annotated_trees/beast'+str(i)+'.tree.txt', taxon_namespace=tree.taxon_namespace, extract_comment_metadata=True, suppress_internal_node_taxa=False, schema="nexus")
#    treelist = dendropy.TreeList.get(path="beast_output\\beast"+str(i)+'.trees.txt', schema="nexus")
    
#    print(annotated_tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True))
#    print(tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True))       
           
    nodes_tree = [nd for nd in tree.preorder_node_iter()]
    nodes_annotated_tree = [nd for nd in annotated_tree.preorder_node_iter()]
    for i, n in enumerate(nodes_tree):            
        if dimension==2:
            nodes_tree[i].averageX=float(nodes_annotated_tree[i].annotations.require_value("location1"))
            nodes_tree[i].averageY=float(nodes_annotated_tree[i].annotations.require_value("location2"))
        else:
            nodes_tree[i].averageX=float(nodes_annotated_tree[i].annotations.require_value("X"))
            
            
#        if nodes_annotated_tree[i].taxon is not None:
#            if nodes_tree[i].taxon is not None:
#                print("%s : %s : %s : %s " % (nodes_tree[i].taxon.label, nodes_annotated_tree[i].taxon.label, nodes_tree[i].X, nodes_tree[i].averageX))
#            else:
#                print("None : None : %s : %s " % (nodes_tree[i].X, nodes_tree[i].averageX))
#                              
#        else:
#            if nodes_tree[i].taxon is not None:
#                print("%s : None : %s : %s " % (nodes_tree[i].taxon.label, nodes_tree[i].X, nodes_tree[i].averageX))
#            else:
#                print("None : None : %s : %s " % (nodes_tree[i].X, nodes_tree[i].averageX))
    
#    n=0
#    average_diffusion_rate = 0
#    
#    for i, t in enumerate(treelist):
#        if i > 100:
#            diffusion_rate = 0
#            tree_length = 0
#            for node in t.preorder_node_iter():
#                if node.parent_node is not None:
#                    tree_length=tree_length + node.edge_length
#                    if dimension ==2:
#                        displacementX=float(node.annotations.require_value("location")[0])-float(node.parent_node.annotations.require_value("location")[0])
#                        displacementY=float(node.annotations.require_value("location")[1])-float(node.parent_node.annotations.require_value("location")[1])
#                        diffusion_rate = diffusion_rate + math.sqrt(displacementX*displacementX+displacementY*displacementY)
#                    else:
#                        displacementX=float(node.annotations.require_value("X"))-float(node.parent_node.annotations.require_value("X"))
#                        diffusion_rate=diffusion_rate + abs(displacementX)
#            diffusion_rate= diffusion_rate/tree_length
#            average_diffusion_rate=average_diffusion_rate+diffusion_rate
#            n=n+1
#    average_diffusion_rate = average_diffusion_rate/n
#    print("Average diffusion rate: "+str(average_diffusion_rate)+"\n")
    return tree                  

def calculate_time_to_tips(tree):
    for leaf in tree.leaf_node_iter():
        leaf.time_to_tips = 0
    for node in tree.postorder_node_iter():
        if not hasattr(node, "time_to_tips"):
            node.time_to_tips=min(set([child.time_to_tips+child.edge_length for child in node.child_nodes()]))  
    return tree

def sample_biased_extreme(tree, dimension):
    for node in tree.preorder_node_iter():
        node.leave = False
        
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips/5)
    X_positions = np.zeros(num_tips)
    index = 0
    for leaf in tree.leaf_node_iter():
        X_positions[index]=leaf.X
        index=index+1
        
        
    X_positions=X_positions[X_positions.argsort()]
    cutoff=X_positions[-sample_size]

    

    index = 0
    for leaf in tree.leaf_node_iter():
        if leaf.X >= cutoff:
            leaf.leave = True
        index = index+1    
            
    labels = []        
    for node in tree.postorder_node_iter(): 
        for child in node.child_node_iter():
            if child.leave:
                node.leave =True
#        print(node.taxon.label+ " " +str(node.leave))
        
        if not node.leave:
            labels.append(node.taxon.label)
            
            
    t1 = tree.extract_tree_without_taxa_labels(labels=labels)

    #seed_position_X = tree.find_node_with_taxon_label(t1.seed_node.taxon.label).X
    #if dimension ==2:
        #seed_position_Y = tree.find_node_with_taxon_label(t1.seed_node.taxon.label).Y
    for node in t1.preorder_node_iter():
        node.X=tree.find_node_with_taxon_label(node.taxon.label).X#-seed_position_X
        node.time=tree.find_node_with_taxon_label(node.taxon.label).time
        if dimension ==2:
            node.Y=tree.find_node_with_taxon_label(node.taxon.label).Y#-seed_position_Y
    return t1


def sample_biased_halfline(tree, dimesnion):
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
        if leaf.X >=0:
            leaf.leave = True
            
            
    labels = []        
    for node in tree.postorder_node_iter(): 
        for child in node.child_node_iter():
            if child.leave:
                node.leave =True
        
        if not node.leave:
            labels.append(node.taxon.label)

    t1 = tree.extract_tree_without_taxa_labels(labels=labels)
    #seed_position = tree.find_node_with_taxon_label(t1.seed_node.taxon.label).X
    for node in t1.preorder_node_iter():
        node.X=tree.find_node_with_taxon_label(node.taxon.label).X#-seed_position
        if dimension==2:
            node.Y=tree.find_node_with_taxon_label(node.taxon.label).Y
        node.time=tree.find_node_with_taxon_label(node.taxon.label).time
    return t1
          

def sample_biased_quadrant(tree, dimension):
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
#        print(leaf.taxon.label+ " " +str(leaf.X))
        if leaf.X >=0 and leaf.Y >= 0:
            leaf.leave = True
            
            
    labels = []        
    for node in tree.postorder_node_iter(): 
        for child in node.child_node_iter():
            if child.leave:
                node.leave =True
        
        if not node.leave:
            labels.append(node.taxon.label)
                
    t1 = tree.extract_tree_without_taxa_labels(labels=labels)

#    print(t1.taxon_namespace)
#    print(tree.taxon_namespace)
    #seed_position_X = tree.find_node_with_taxon_label(t1.seed_node.taxon.label).X
    #seed_position_Y = tree.find_node_with_taxon_label(t1.seed_node.taxon.label).Y
    for node in t1.preorder_node_iter():
        node.X=tree.find_node_with_taxon_label(node.taxon.label).X#-seed_position_X
        node.Y=tree.find_node_with_taxon_label(node.taxon.label).Y#-seed_position_Y
        node.time=tree.find_node_with_taxon_label(node.taxon.label).time
    return t1
          

def sample_biased_radius_central(tree, dimension):
    radius = 1
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
        if dimension==2:            
            if math.sqrt(leaf.X*leaf.X+leaf.Y*leaf.Y) <= radius:
                leaf.leave = True
        else:
            if abs(leaf.X) <= radius:
                leaf.leave = True            
            
    labels = []        
    for node in tree.postorder_node_iter(): 
        for child in node.child_node_iter():
            if child.leave:
                node.leave =True
        
        if not node.leave:
            labels.append(node.taxon.label)
                
    t1 = tree.extract_tree_without_taxa_labels(labels=labels)


    for node in t1.preorder_node_iter():
        node.X=tree.find_node_with_taxon_label(node.taxon.label).X#-seed_position_X
        if dimension==2:
            node.Y=tree.find_node_with_taxon_label(node.taxon.label).Y
        node.time=tree.find_node_with_taxon_label(node.taxon.label).time
    return t1

def sample_biased_most_central(tree, dimension):
    for node in tree.preorder_node_iter():
        node.leave = False        
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips/5)
    distances = np.zeros(num_tips)
    index = 0
    for leaf in tree.leaf_node_iter():
        if dimension ==2:
            leaf.distance=math.sqrt(leaf.X*leaf.X+leaf.Y*leaf.Y)
            distances[index]=leaf.distance
            index = index+1
        else:
            leaf.distance=abs(leaf.X)
            distances[index]=leaf.distance
            index=index+1
    distances=distances[distances.argsort()]
    cutoff=distances[sample_size-1]
    
    for leaf in tree.leaf_node_iter():
        if leaf.distance <= cutoff:
            leaf.leave = True  

    labels = []        
    for node in tree.postorder_node_iter(): 
        for child in node.child_node_iter():
            if child.leave:
                node.leave =True        
        if not node.leave:
            labels.append(node.taxon.label)            
            
    t1 = tree.extract_tree_without_taxa_labels(labels=labels)
    for node in t1.preorder_node_iter():
        node.X=tree.find_node_with_taxon_label(node.taxon.label).X#-seed_position_X
        node.time=tree.find_node_with_taxon_label(node.taxon.label).time        
        if dimension ==2:
            node.Y=tree.find_node_with_taxon_label(node.taxon.label).Y#-seed_position_Y
    return t1



def sample_unbiased(tree, dimension):
    
    for node in tree.preorder_node_iter():
        node.leave = False    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips/5)
    sample_set = random.sample(range(num_tips), sample_size)
    index = 0
    for leaf in tree.leaf_node_iter():
        if index in sample_set:
            leaf.leave = True
        index = index+1    
            
    labels = []        
    for node in tree.postorder_node_iter(): 
        for child in node.child_node_iter():
            if child.leave:
                node.leave =True
        
        if not node.leave:
            labels.append(node.taxon.label)
                
    t1 = tree.extract_tree_without_taxa_labels(labels=labels)
    for node in t1.preorder_node_iter():
        node.X=tree.find_node_with_taxon_label(node.taxon.label).X
        node.time=tree.find_node_with_taxon_label(node.taxon.label).time
        if dimension ==2:
            node.Y=tree.find_node_with_taxon_label(node.taxon.label).Y
    return t1

def sample_biased_diagonal(tree, dimension):
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
        if abs(leaf.X-leaf.Y) <= 1:
            leaf.leave = True
                        
    labels = []        
    for node in tree.postorder_node_iter(): 
        for child in node.child_node_iter():
            if child.leave:
                node.leave =True
        
        if not node.leave:
            labels.append(node.taxon.label)
                
    t1 = tree.extract_tree_without_taxa_labels(labels=labels)


    for node in t1.preorder_node_iter():
        node.X=tree.find_node_with_taxon_label(node.taxon.label).X#-seed_position_X
        node.Y=tree.find_node_with_taxon_label(node.taxon.label).Y#-seed_position_Y
        node.time=tree.find_node_with_taxon_label(node.taxon.label).time
    return t1


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
job_index = args.linux
linux = args.linux
print(linux)
run_tree_annotator = args.run_tree_annotator
#import os
if not os.path.exists("output"):
    os.makedirs("output")
if not os.path.exists("output/beast_input"):
    os.makedirs("output/beast_input")
if not os.path.exists("output/beast_output"):
    os.makedirs("output/beast_output")
if not os.path.exists("output/generated_trees"):
    os.makedirs("output/generated_tree")
if not os.path.exists("output/annotated_trees"):
    os.makedirs("output/annotated_trees")
if not os.path.exists("output/phyrex_input"):
    os.makedirs("output/phyrex_input")
if not os.path.exists("output/phyrex_output"):
    os.makedirs("output/phyrex_output")

num_sampling = 4

for i in range(0, num_sampling):
    if not os.path.exists("output/sampled_beast_input"+str(i+1)):
        os.makedirs("output/sampled_beast_input"+str(i+1))
    if not os.path.exists("output/sampled_beast_output"+str(i+1)):
        os.makedirs("output/sampled_beast_output"+str(i+1))
    if not os.path.exists("output/generated_sampled_trees"+str(i+1)):
        os.makedirs("output/generated_sampled_trees"+str(i+1))
    if not os.path.exists("output/annotated_sampled_trees"+str(i+1)):
        os.makedirs("output/annotated_sampled_trees"+str(i+1))

for i in range(num_trees*(job_index-1), num_trees*job_index):
    #random.seed = 1357+i
    print("index is" + str(i))
    t=dendropy.Tree()
    
    if args.tree_type == "nuc":
        t = treegenerator.generate_nonultrametric_coalescent_tree(num_tips_per_period, num_periods, period_length, lamb)
    elif args.tree_type == "uc":
        t = treegenerator.generate_ultrametric_coalescent_tree(num_tips, lamb)
    elif args.tree_type == "bd":
        t = treegenerator.generate_birthdeath_tree(1, 0.1, num_tips)
        #print(t.as_ascii_plot(plot_metric='length', show_internal_node_labels=True))
    elif args.tree_type == "star":
        t = treegenerator.generate_star_tree()
    else:
        print(args.tree_type+" is invalid tree type")
        break
                
    t= simulate_brownian(t, sigma, dimension) 
    t=calculate_time_to_tips(t)
     
    beastxmlwriter.write_BEAST_xml(t, i, dimension, mcmc, log_every)   
    phyrexxmlwriter.write_phyrex_input(t, i)     
    run_sample_analysis=False
    
    if run_sample_analysis:    
        
        for output_index in range(1, num_sampling+1):
            if output_index ==1:
                sampled_t=sampling.sample_unbiased(t, dimension, sample_ratio=0.05)
            elif output_index==2:
                sampled_t=sampling.sample_biased_vertical(t, dimension, sample_ratio=0.05)
            elif output_index==3:
                sampled_t=sampling.sample_biased_most_central(t, dimension, sample_ratio=0.05)            
            elif output_index==4:
                sampled_t=sampling.sample_biased_extreme(t, dimension, sample_ratio=0.05)
            sampled_t=calculate_time_to_tips(sampled_t)
            beastxmlwriter.write_BEAST_xml(sampled_t, i, dimension, mcmc, log_every, "output/sampled_beast_input"+str(output_index)+"/sampled_beast", beast_output_string="output/sampled_beast_output"+str(output_index)+"/sampled_beast")
            for node in sampled_t.preorder_node_iter():
                node.annotations.add_bound_attribute("time")
                node.annotations.add_bound_attribute("X")
                node.annotations.add_bound_attribute("time_to_tips")
                if dimension==2:
                    node.annotations.add_bound_attribute("Y")
            sampled_t.write(path="output/generated_sampled_trees"+str(output_index)+"/sampled_tree"+str(i)+".txt", schema="nexus", suppress_internal_taxon_labels=True)

 
        
    for node in t.preorder_node_iter():
        node.annotations.add_bound_attribute("time")
        node.annotations.add_bound_attribute("X")
        node.annotations.add_bound_attribute("time_to_tips")
        if dimension==2:
            node.annotations.add_bound_attribute("Y")
    t.write(path="output/generated_trees/tree"+str(i)+".txt", schema="nexus", suppress_internal_taxon_labels=True)
    
    run_analysis= True

    burnin=int(mcmc/10)
    if run_analysis:
        if linux:
            os.system('beast -overwrite -seed 123456795 "output/beast_input/beast'+str(i)+'.xml"')
        else:
            os.system('cmd /c java -jar beast.jar -overwrite -seed 123456795 "output/beast_input/beast'+str(i)+'.xml"')
        if run_tree_annotator:
            os.system('treeannotator -burnin '+str(burnin)+' "output/beast_output/beast'+str(i)+'.trees.txt" "output/annotated_trees/beast'+str(i)+'.tree.txt"')
        
        
    if run_sample_analysis:
        for output_index in range(1, num_sampling+1):
            if linux:
                os.system('beast -overwrite -seed 123456795 "output/sampled_beast_input'+str(output_index)+'/sampled_beast'+str(i)+'.xml"')
            else:
                os.system('cmd /c java -jar beast.jar -overwrite -seed 123456795 "output/sampled_beast_input'+str(output_index)+'/sampled_beast'+str(i)+'.xml"')
            if run_tree_annotator:
                #os.system('cmd /c ""C:/Users/Antanas/Desktop/BEAST v1.10.4/bin/treeannotator" -burnin '+str(burnin)+' "C:/Users/Antanas/Phylogeny Simulation/output/sampled_beast_output'+str(output_index)+'/sampled_beast'+str(i)+'.trees.txt" "C:/Users/Antanas/Phylogeny Simulation/output/annotated_sampled_trees'+str(output_index)+'/sampled_beast'+str(i)+'.tree.txt""')
                os.system('treeannotator -burnin '+str(burnin)+' "output/sampled_beast_output'+str(output_index)+'/sampled_beast'+str(i)+'.trees.txt" "output/annotated_sampled_trees'+str(output_index)+'/sampled_beast'+str(i)+'.tree.txt"')
           

