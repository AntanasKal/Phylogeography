# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:37:34 2019

@author: Antanas
"""

#the code that wirtes BEAST xml files is in different file, "beastxmlwriter.py"
import beastxmlwriter
import treegenerator

import math 
import argparse
import dendropy
import random
import numpy as np


import os

#this function generates borwnian motion for a given tree
def simulate_brownian(t, sigma, dimension):
    #t is the tree
    #sigma is the standard deviation of the brownian motion
    #dimension is the number of dimensions in which we generate the random walk
    #random.seed = 134679
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
    return t

    
def analyze_tree_list(tree, i): 
#   command to execute treeannotator       
    os.system('cmd /c ""C:\\Users\\Antanas\\Desktop\\BEAST v1.10.4\\bin\\treeannotator" -burnin 1000 "C:\\Users\\Antanas\\Phylogeny Simulation\\beast_output\\beast'+str(i)+'.trees.txt" "C:\\Users\\Antanas\\Phylogeny Simulation\\annotated_trees\\beast'+str(i)+'.tree.txt""')
    for node in tree.preorder_node_iter():
        node.averageX=0
        
    annotated_tree=dendropy.Tree.get(path='annotated_trees\\beast'+str(i)+'.tree.txt', taxon_namespace=tree.taxon_namespace, extract_comment_metadata=True, suppress_internal_node_taxa=False, schema="nexus")
    
    print(annotated_tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True))
    print(tree.as_ascii_plot(plot_metric='length', show_internal_node_labels=True))              
    nodes_tree = [nd for nd in tree.preorder_node_iter()]
    nodes_annotated_tree = [nd for nd in annotated_tree.preorder_node_iter()]
    for i, n in enumerate(nodes_tree):            
        nodes_tree[i].averageX=float(nodes_annotated_tree[i].annotations.require_value("X"))
        if nodes_annotated_tree[i].taxon is not None:
            if nodes_tree[i].taxon is not None:
                print("%s : %s : %s : %s " % (nodes_tree[i].taxon.label, nodes_annotated_tree[i].taxon.label, nodes_tree[i].X, nodes_tree[i].averageX))
            else:
                print("None : None : %s : %s " % (nodes_tree[i].X, nodes_tree[i].averageX))
                              
        else:
            if nodes_tree[i].taxon is not None:
                print("%s : None : %s : %s " % (nodes_tree[i].taxon.label, nodes_tree[i].X, nodes_tree[i].averageX))
            else:
                print("None : None : %s : %s " % (nodes_tree[i].X, nodes_tree[i].averageX))                  

parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-dims', dest='dimension', type=int, default=1, help='number of dimensions (1 or 2) for which the random walk is generated (default: 1)')
parser.add_argument('-N', action="store", type=int, dest="num_trees", default=5, help='number of simulations (default 5)')
parser.add_argument('-treetype', action="store", dest="tree_type", default='nuc', help='type of tree generated \n "nuc" - nonultrametric coalescent \n "uc" - ultrametric coalescent \n "bd" - birth-death tree')



args = parser.parse_args()
dimension = args.dimension
num_trees = args.num_trees

sigma = 1


for i in range(num_trees):
    #random.seed = 1357+i
    t=dendropy.Tree()
    if args.tree_type == "nuc":
        t = treegenerator.generate_coalescent_nonultrametric_tree()
    elif args.tree_type == "uc":
        t = treegenerator.generate_ultrametric_coalescent_tree(50)
    elif args.tree_type == "bd":
        t = treegenerator.generate_birthdeath_tree(1, 0.5, 30)
    elif args.tree_type == "star":
        t = treegenerator.generate_star_tree2()
    else:
        print(args.tree_type+" is invalid tree type")
        break
                
#    if args.tree_type != "star":    
    t= simulate_brownian(t, sigma, dimension)    
    
    #simulate hky (currently not needed)
    d=dendropy.model.discrete.hky85_chars(kappa=3, mutation_rate=0.01, seq_len=1000,tree_model=t, retain_sequences_on_tree=False)
 #   d= dendropy.DnaCharacterMatrix()
    beastxmlwriter.write_BEAST_xml(t, d, i, dimension)

    d.write(path='output1/fasta_output'+str(i)+'.txt', schema="fasta")
    file1 = open("output6/dates_output"+str(i)+".txt","w") 
    file2 = open("output7/position_output"+str(i)+".txt","w") 
    file2.write("\tX\tY\n")
    for tax in d:
        file1.write(tax.label+"\t"+str(t.find_node_for_taxon(tax).time)+"\n")
        if dimension==2:
            file2.write(tax.label+"\t"+str(t.find_node_for_taxon(tax).X)+"\t"+str(t.find_node_for_taxon(tax).Y)+"\n")   
        else:
            file2.write(tax.label+"\t"+str(t.find_node_for_taxon(tax).X)+"\n")
    file1.close() 
    file2.close()    
    
    
    os.system('cmd /c java -jar beast.jar -overwrite -seed 123456795 "beast_input\\beast'+str(i)+'.xml"')
    
#    analyze_tree_list(t, i)
    
#    for node in t.leaf_node_iter():
#        if node.parent_node is not None:
#            print("%s : %s : %s" % (node.parent_node.taxon.label, node.taxon.label, node.edge_length))
            
    mean=0
    n=0
    mle1 = 0
    for node in t.leaf_node_iter():
        n=n+1
        mean=mean+node.X
        mle1=mle1+node.X*node.X
    mle1=mle1/n
    mean=mean/n
    print("%s : %s : %s" %(n, mean, mle1))
    mle2=0
    for node in t.leaf_node_iter():            
        mle2=mle2+(node.X-mean)*(node.X-mean)
    mle2=mle2/n
    print("%s : %s : %s" %(n, mean, mle2))      
    

    