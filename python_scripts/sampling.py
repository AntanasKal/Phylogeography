# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 10:05:24 2019

@author: Antanas

Modified by Nicola De Maio
"""

import numpy as np
import math
import random

#Code for sampling


def delete_tips(tree, dimension=2):
    #the function takes a tree and returns a subtree containing all the nodes
    #that have a descendant node with "leave" attribute set to node.leave=True
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

def sample_biased_extreme(tree, dimension=2, sample_ratio=0.1):
    #function returns a subtree with the all of the leaves that have the largest x coordinate
    for node in tree.preorder_node_iter():
        node.leave = False       
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips*sample_ratio)
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
        
    return delete_tips(tree, dimension)
    
#NICOLA: this is the new function to sample partly extreme tips and partly random tips.
def sample_biased_side(tree, dimension=2, sample_ratio=0.1, ratioExtreme=0.5):
    #function returns a subtree with partly randomly sampled leaves and partly leaves that have the largest x coordinate
    for node in tree.preorder_node_iter():
        node.leave = False       
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size_extreme = int(num_tips*sample_ratio*ratioExtreme)
    X_positions = np.zeros(num_tips)
    index = 0
    for leaf in tree.leaf_node_iter():
        X_positions[index]=leaf.X
        index=index+1        
        
    X_positions=X_positions[X_positions.argsort()]
    cutoff=X_positions[-sample_size_extreme]  
    if sample_size_extreme==0:
    	cutoff=float("inf")

    #index = 0
    leaves=[]
    for leaf in tree.leaf_node_iter():
        if leaf.X >= cutoff:
            leaf.leave = True
        else:
        	leaves.append(leaf)
        #index = index+1  
        
    sample_size =   int(num_tips*sample_ratio*(1.0-ratioExtreme)) 
    
    sample_set = random.sample(range(num_tips-sample_size_extreme), sample_size)
    print("values")
    print(len(leaves))
    print(num_tips)
    print(sample_size_extreme)
    print(sample_size)
    print(sample_set)
    print(len(sample_set))
    for i in range(len(sample_set)):
    	#print(sample_set[i])
    	leaves[sample_set[i]].leave = True 
        
    return delete_tips(tree, dimension)


def sample_biased_halfline(tree, dimension=2, sample_ratio=1):
    #function returns the subtree with all leaves with positive x coordinate
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
        if leaf.X >=0:
            leaf.leave = True
            
    return delete_tips(tree, dimension)
          

def sample_biased_quadrant(tree, dimension=2, sample_ratio=0.1):
    #function returns the subtree with leaves in the 1st quadrant
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
#        print(leaf.taxon.label+ " " +str(leaf.X))
        if leaf.X >=0 and leaf.Y >= 0:
            leaf.leave = True
            
    return delete_tips(tree, dimension)

def sample_biased_radius_central(tree, dimension=2, radius=1, sample_ratio=1):
    #function returns the subtree with all the leaves inside a disc of given radius    
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
        if dimension==2:            
            if math.sqrt(leaf.X*leaf.X+leaf.Y*leaf.Y) <= radius:
                leaf.leave = True
        else:
            if abs(leaf.X) <= radius:
                leaf.leave = True        
                
    return delete_tips(tree, dimension)

def sample_biased_most_central(tree, dimension=2, sample_ratio=0.1):
    #function returns the subtree with the most central leaves
    for node in tree.preorder_node_iter():
        node.leave = False        
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips*sample_ratio)
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

    return delete_tips(tree, dimension)



def sample_unbiased(tree, dimension=2, sample_ratio=0.1):
    #function returns a subtree with uniformly randomly sampled leaves
    
    for node in tree.preorder_node_iter():
        node.leave = False    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips*sample_ratio)
    sample_set = random.sample(range(num_tips), sample_size)
    index = 0
    for leaf in tree.leaf_node_iter():
        if index in sample_set:
            leaf.leave = True
        index = index+1    
        
    return delete_tips(tree, dimension)

def sample_biased_diagonal_strip(tree, dimension=2, strip_width=1, sample_ratio=1):
    #function returns the subtree with leaves within given distance from the diagonal
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
        if abs(leaf.X-leaf.Y) <= strip_width/2:
            leaf.leave = True
            
    return delete_tips(tree, dimension)

def sample_biased_diagonal(tree, dimension=2, sample_ratio=0.1):
    #function returns the subtree with leaves closest to the diagonal line
    for node in tree.preorder_node_iter():
        node.leave = False        
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips*sample_ratio)
    distances = np.zeros(num_tips)
    index = 0
    for leaf in tree.leaf_node_iter():

        leaf.distance=abs(leaf.X-leaf.Y)/math.sqrt(2)
        distances[index]=leaf.distance
        index = index+1

    distances=distances[distances.argsort()]
    cutoff=distances[sample_size-1]
    
    for leaf in tree.leaf_node_iter():
        if leaf.distance <= cutoff:
            leaf.leave = True  
            
    return delete_tips(tree, dimension)


def sample_biased_horizontal(tree, dimension=2, sample_ratio=0.1):
    #function returns the subtree with leaves closest to x axis    
    for node in tree.preorder_node_iter():
        node.leave = False        
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips*sample_ratio)
    distances = np.zeros(num_tips)
    index = 0
    for leaf in tree.leaf_node_iter():

        leaf.distance=abs(leaf.Y)
        distances[index]=leaf.distance
        index = index+1
#        else:
#            leaf.distance=abs(leaf.X)
#            distances[index]=leaf.distance
#            index=index+1
    distances=distances[distances.argsort()]
    cutoff=distances[sample_size-1]
    
    for leaf in tree.leaf_node_iter():
        if leaf.distance <= cutoff:
            leaf.leave = True  
            
    return delete_tips(tree, dimension)


def sample_biased_vertical(tree, dimension=2, sample_ratio=0.1):
    #function returns the subtree with leaves closest to y axis
    for node in tree.preorder_node_iter():
        node.leave = False        
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips*sample_ratio)
    distances = np.zeros(num_tips)
    index = 0
    for leaf in tree.leaf_node_iter():

        leaf.distance=abs(leaf.X)
        distances[index]=leaf.distance
        index = index+1

    distances=distances[distances.argsort()]
    cutoff=distances[sample_size-1]
    
    for leaf in tree.leaf_node_iter():
        if leaf.distance <= cutoff:
            leaf.leave = True  
            
    return delete_tips(tree, dimension)


def sample_biased_2_points(tree, dimension=2, sample_ratio=0.1):
    #function the subtree with leaves closest to one of 2 points
    point_1_X=2
    point_1_Y=0
    
    point_2_X=-2
    point_2_Y=0
    
    for node in tree.preorder_node_iter():
        node.leave = False        
    
    num_tips = 0
    for leaf in tree.leaf_node_iter():
        num_tips = num_tips +1
    
    sample_size = int(num_tips*sample_ratio)
    distances = np.zeros(num_tips)
    index = 0
    for leaf in tree.leaf_node_iter():
        distance_1=math.sqrt((leaf.X-point_1_X)*(leaf.X-point_1_X)+(leaf.Y-point_1_Y)*(leaf.Y-point_1_Y))
        distance_2=math.sqrt((leaf.X-point_2_X)*(leaf.X-point_2_X)+(leaf.Y-point_2_Y)*(leaf.Y-point_2_Y))
        leaf.distance=min(distance_1, distance_2)
        distances[index]=leaf.distance
        index = index+1
    distances=distances[distances.argsort()]
    cutoff=distances[sample_size-1]
    
    for leaf in tree.leaf_node_iter():
        if leaf.distance <= cutoff:
            leaf.leave = True  
            
    return delete_tips(tree, dimension)
