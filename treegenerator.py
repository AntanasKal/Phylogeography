# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 14:29:19 2019

@author: Antanas
"""


import math 
import dendropy
import random
from dendropy.simulate import treesim
import numpy as np


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
        namespace.append("s"+str(index))
    
    #name all nodes instead of just leaves
    taxon_namespace = dendropy.TaxonNamespace(namespace)
    t.taxon_namespace=taxon_namespace
    index=0
    for node in t.preorder_node_iter():
        index=index+1
        node.taxon=t.taxon_namespace.get_taxon("s"+str(index))
    
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
def generate_ultrametric_coalescent_tree(num_tips, lamb):
    names = []
    #lamb = 1
    
    #if there are N tips, there must be 2N-1 nodes
    for i in range(2*num_tips-1):
        names.append("s"+str(i))
    
    taxon_namespace = dendropy.TaxonNamespace(names)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    time_from_present = 0
    
    #current_nodes is a list of nodes that are currently not merged
    current_nodes = []
    for i in range(num_tips):
        node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(i)))
        current_nodes.append(node)
        node.age = 0 
    
    
    #if there are N leaves, there must be N-1 merges
    for merges in range(num_tips-1):
        #calculating time to the next coalescent
        time_to_coalescent=random.expovariate(lamb*len(current_nodes)*(len(current_nodes)-1)/2)
        
        time_from_present=time_from_present+time_to_coalescent
        
        #choosing 2 indices of nodes that will be merged  randomly
        merging_branches = random.sample(range(len(current_nodes)),2)
        node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(merges+num_tips)))
        
        #if it is the last merge, instead of creting a new node, we set the node of the merge to be the seed node
        if merges == num_tips-2:
            node=tree.seed_node
            node.taxon=taxon_namespace.get_taxon("s"+str(merges+num_tips))
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


def generate_birth_death_custom(br, dr, num_extinct):
    names = []
    #lamb = 1
    
    successful = False
    while not successful:
        #generation_done=False
        alive_nodes=[]
        tree = dendropy.Tree()
        allive_nodes.append(tree.seed_node)
        #while not generation_done:
            
    
    
    
    #if there are N tips, there must be 2N-1 nodes
    for i in range(2*num_tips-1):
        names.append("s"+str(i))
    
    taxon_namespace = dendropy.TaxonNamespace(names)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    time_from_present = 0
    
    #current_nodes is a list of nodes that are currently not merged
    current_nodes = []
    for i in range(num_tips):
        node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(i)))
        current_nodes.append(node)
        node.age = 0 
    
    
    
    #if there are N leaves, there must be N-1 merges
    for merges in range(num_tips-1):
        #calculating time to the next coalescent
        time_to_coalescent=random.expovariate(lamb*len(current_nodes)*(len(current_nodes)-1)/2)
        
        time_from_present=time_from_present+time_to_coalescent
        
        #choosing 2 indices of nodes that will be merged  randomly
        merging_branches = random.sample(range(len(current_nodes)),2)
        node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(merges+num_tips)))
        
        #if it is the last merge, instead of creting a new node, we set the node of the merge to be the seed node
        if merges == num_tips-2:
            node=tree.seed_node
            node.taxon=taxon_namespace.get_taxon("s"+str(merges+num_tips))
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


def generate_star_tree2():
    num_tips = 10
    branch_length = 1
    names = []
    
    for i in range(num_tips+1):
        names.append("s"+str(i))
    
    taxon_namespace = dendropy.TaxonNamespace(names)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    
    index = 0
    for i in range(num_tips+1):
        if index == 0:
            tree.seed_node.taxon=taxon_namespace.get_taxon("s"+str(0))
            
            tree.seed_node.X = 0
            tree.seed_node.time = 0
        else:
            node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(index)))
            node.edge_length = branch_length
            node.X = random.gauss(0,1)
            node.time = branch_length
            tree.seed_node.add_child(node)
    return tree

def generate_star_tree():
    num_tips = 100
    branch_length = 1
    fake_step = 0.000000001
    names = [] 
      
    for i in range(num_tips):
        names.append("s"+str(i))   
        
    taxon_namespace = dendropy.TaxonNamespace(names)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    current_seed=dendropy.Node()
    current_seed.edge_length = fake_step
    i=0
    while i <num_tips:
        if i ==0:
            node1 = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(i)))
            node1.edge_length = branch_length
            node1.X = random.gauss(0,branch_length)
            current_seed.add_child(node1)
            i=i+1
            
            node2 = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(i)))
            node2.edge_length = branch_length
            node2.X = random.gauss(0,branch_length)
            current_seed.add_child(node2)
            current_seed.X=0
            i=i+1
        elif i == num_tips-1:
            node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(i)))
            node.edge_length = branch_length
            node.X = random.gauss(0,branch_length)
            i=i+1
            tree.seed_node.X=0
            tree.seed_node.add_child(node)
            tree.seed_node.add_child(current_seed)
        else:
            current_seed2=dendropy.Node()
            current_seed2.edge_length = fake_step
            
            node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(i)))
            node.edge_length = branch_length
            node.X = random.gauss(0,branch_length)
            i=i+1
            
            current_seed2.add_child(node)
            current_seed2.add_child(current_seed)
            current_seed= current_seed2
            current_seed.X=0
    for node in tree.internal_nodes():
        node.taxon = taxon_namespace.get_taxon("s"+str(i))
        i=i+1
    tree=calculate_times(tree)      
    return tree


#function to generate coalescent nonultrametric trees
def generate_nonultrametric_coalescent_tree(num_tips_per_period, num_periods, period_length, lamb):
    
    #a number num_tips_per_period is added every period_length of time units for num_periods times
    num_tips = num_tips_per_period*num_periods
    names = []
    
    #if there are N tips, there must be 2N-1 nodes
    for i in range(2*num_tips-1):
        names.append("s"+str(i))

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
            node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(index)))
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
                    node.taxon=taxon_namespace.get_taxon("s"+str(index))
                else:
                    node = dendropy.Node(taxon=taxon_namespace.get_taxon("s"+str(index)))
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