# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 15:37:34 2019

@author: Antanas
"""

#the code that wirtes BEAST xml files is in different file, "beastxmlwriter.py"
import beastxmlwriter

import math 
import argparse
import dendropy
import random
from dendropy.simulate import treesim

#this function generates borwnian motion for a given tree
def simulate_brownian(t, sigma, dimension):
    #t is the tree
    #sigma is the standard deviation of the brownian motion
    #dimension is the number of dimensions in which we generate the random walk
    for node in t.preorder_node_iter():
        if node.parent_node is None:
            node.X = float(0)
            node.displacementX = float(0)
            if dimension==2:
                node.Y = float (0)
                node.displacementY = float(0)
        else:
            node.displacementX = random.gauss(0, sigma*math.sqrt(node.edge.length))
#            node.displacementX = random.gauss(0, sigma*node.edge.length)
            node.X = node.parent_node.X+node.displacementX            
            if dimension==2:
                node.displacementY = random.gauss(0, sigma*math.sqrt(node.edge.length))
                node.Y = node.parent_node.Y+node.displacementY            
    return t

#this function calculates the time from the seed node to the present (seed node has time 0)
def calculate_times(t):
    for node in t.preorder_node_iter():
        if node.parent_node is None:
            node.time = 0
        else:
            node.time = node.parent_node.time+node.edge.length            
    return t 


def generate_birthdeath_tree(br, dr, num_extinct):
    t = treesim.birth_death_tree(birth_rate=br, death_rate=dr, num_extinct_tips=num_extinct, is_retain_extinct_tips=True, is_add_extinct_attr=True)
    #t.print_plot()    
    
    index = 0
    namespace = [];
    
    for node in t.preorder_node_iter():
        index=index+1
        namespace.append("T"+str(index))
    
    #name all nodes instead of just leaves
    taxon_namespace = dendropy.TaxonNamespace(namespace)
    t.taxon_namespace=taxon_namespace
    index=0
    for node in t.preorder_node_iter():
        index=index+1
        node.taxon=t.taxon_namespace.get_taxon("T"+str(index))
    
    t =prune_nodes(t)
    
    #distance to root
    t=calculate_times(t)
        
    return t

#this function checks which of the nodes is an extinct leaf or an ancestor of an extinct leaf
#and only leaves these nodes in the tree
def prune_nodes(t):
    for leaf in t.leaf_node_iter():    
        if hasattr(leaf, 'is_extinct'):
            leaf.extinct_ancestor = True
        else:
            leaf.extinct_ancestor = False
        
    for node in t.postorder_node_iter(): 
        if not hasattr(node, 'extinct_ancestor'):
            child_extinct = False
            for child in node.child_node_iter():
                if child.extinct_ancestor:
                    child_extinct =True
            node.extinct_ancestor = child_extinct
    labels = set([taxon.label for taxon in t.taxon_namespace
        if not t.find_node_for_taxon(taxon).extinct_ancestor])
    t1 = t.extract_tree_without_taxa_labels(labels=labels)
    return t1


#function to generate coalescent trees (this is the ultrametric case)
def generate_ultrametric_coalescent_tree(num_tips):
    names = []
    lamb = 1
    
    #if there are N tips, there must be 2N-1 nodes
    for i in range(2*num_tips-1):
        names.append("T"+str(i))
    
    taxon_namespace = dendropy.TaxonNamespace(names)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    time_from_present = 0
    
    #current_nodes is a list of nodes that are currently not merged
    current_nodes = []
    for i in range(num_tips):
        node = dendropy.Node(taxon=taxon_namespace.get_taxon("T"+str(i)))
        current_nodes.append(node)
        node.age = 0 
    
    
    #if there are N leaves, there must be N-1 merges
    for merges in range(num_tips-1):
        #calculating time to the next coalescent
        time_to_coalescent=random.expovariate(lamb*len(current_nodes)*(len(current_nodes)-1)/2)
        
        time_from_present=time_from_present+time_to_coalescent
        
        #choosing 2 indices of nodes that will be merged  randomly
        merging_branches = random.sample(range(len(current_nodes)),2)
        node = dendropy.Node(taxon=taxon_namespace.get_taxon("T"+str(merges+num_tips)))
        
        #if it is the last merge, instead of creting a new node, we set the node of the merge to be the seed node
        if merges == num_tips-2:
            node=tree.seed_node
            node.taxon=taxon_namespace.get_taxon("T"+str(merges+num_tips))
        node.age = time_from_present
        current_nodes[merging_branches[0]].edge.length=time_from_present-current_nodes[merging_branches[0]].age
        current_nodes[merging_branches[1]].edge.length=time_from_present-current_nodes[merging_branches[1]].age
        node.set_child_nodes([current_nodes[merging_branches[0]], current_nodes[merging_branches[1]]])
        
        #deleting the nodes that have been merging from the list of nodes
        current_nodes.pop(max(merging_branches))
        current_nodes.pop(min(merging_branches))
        current_nodes.append(node)
        
    tree=calculate_times(tree)
    return tree




#function to generate coalescent nonultrametric trees
def generate_coalescent_nonultrametric_tree():
    #a number num_tips_per_period is added every period_length of time units for num_periods times
    lamb=1
    period_length=1
    num_tips_per_period = 5
    num_periods = 3
    num_tips = num_tips_per_period*num_periods
    names = []
    
    #if there are N tips, there must be 2N-1 nodes
    for i in range(2*num_tips-1):
        names.append("T"+str(i))

    taxon_namespace = dendropy.TaxonNamespace(names)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    time_from_present = 0
    current_nodes = []
    
    #index variable simply keeps track of the number of nodes added just to choose the right name
    index = 0
    
    #simulating coalescent during every period
    for current_period in range(num_periods):
        time_from_present=current_period*period_length
        
        #adding new tips during this period
        for i in range(num_tips_per_period):
            node = dendropy.Node(taxon=taxon_namespace.get_taxon("T"+str(index)))
            current_nodes.append(node)
            index= index+1
            node.age = time_from_present
        
        
        current_num_tips = len(current_nodes)
        
        #iterating over the merges (if currently there are N unmerged nodes, there can be at most N-1 merges)
        for merges in range(current_num_tips-1):
            time_to_coalescent=random.expovariate(lamb*len(current_nodes)*(len(current_nodes)-1)/2)
            time_from_present=time_from_present+time_to_coalescent
            
            #if current time exceeds the end of current period, we stop the simulation of the period
            #and start simulating the next period
            if current_period < num_periods-1 and time_from_present > (current_period+1)*period_length:
                break
            else:
                merging_branches = random.sample(range(len(current_nodes)),2)
                
                if merges == current_num_tips-2 and current_period==num_periods-1:
                    node=tree.seed_node
                    node.taxon=taxon_namespace.get_taxon("T"+str(index))
                else:
                    node = dendropy.Node(taxon=taxon_namespace.get_taxon("T"+str(index)))
                index=index+1
                    
                node.age = time_from_present
                current_nodes[merging_branches[0]].edge.length=time_from_present-current_nodes[merging_branches[0]].age
                current_nodes[merging_branches[1]].edge.length=time_from_present-current_nodes[merging_branches[1]].age
                node.set_child_nodes([current_nodes[merging_branches[0]], current_nodes[merging_branches[1]]])
        
                current_nodes.pop(max(merging_branches))
                current_nodes.pop(min(merging_branches))
                current_nodes.append(node)
#     print(tree.as_ascii_plot(show_internal_node_labels=True, plot_metric='length'))
    tree=calculate_times(tree)
    return tree
    

parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-dims', dest='dimension', type=int, default=1, help='number of dimensions (1 or 2) for which the random walk is generated (default: 1)')
parser.add_argument('-N', action="store", type=int, dest="num_trees", default=5, help='number of simulations (default 5)')
parser.add_argument('-treetype', action="store", dest="tree_type", default='nuc', help='type of tree generated \n "nuc" - nonultrametric coalescent \n "uc" - ultrametric coalescent \n "bd" - birth-death tree')



args = parser.parse_args()
dimension = args.dimension
num_trees = args.num_trees



for i in range(num_trees):
    t=dendropy.Tree()
    if args.tree_type == "nuc":
        t = generate_coalescent_nonultrametric_tree()
    elif args.tree_type == "uc":
        t = generate_ultrametric_coalescent_tree()
    elif args.tree_type == "bd":
        t = generate_birthdeath_tree(1, 0.5, 20)
    else:
        print(args.tree_type+" is invalid tree type")
        break
                
        
    t= simulate_brownian(t, 1, dimension)
    
    #simulate hky (currently not needed)
    d=dendropy.model.discrete.hky85_chars(kappa=3, mutation_rate=0.01, seq_len=1000,tree_model=t, retain_sequences_on_tree=False)
 #   d= dendropy.DnaCharacterMatrix()
    beastxmlwriter.write_BEAST_xml(t, d, i, dimension)
    
    import os
    os.system('cmd /c java -jar beast.jar -overwrite -seed 1234 "beast_input\\beast'+str(i)+'.xml"')