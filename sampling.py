# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 10:05:24 2019

@author: Antanas
"""

import numpy as np
import math
import random

def delete_tips(tree, dimension):
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

def sample_biased_extreme(tree, dimension, sample_ratio=0.1):
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


def sample_biased_halfline(tree, dimension, sample_ratio=1):
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
        if leaf.X >=0:
            leaf.leave = True
            
    return delete_tips(tree, dimension)
          

def sample_biased_quadrant(tree, dimension, sample_ratio=0.1):
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
#        print(leaf.taxon.label+ " " +str(leaf.X))
        if leaf.X >=0 and leaf.Y >= 0:
            leaf.leave = True
            
    return delete_tips(tree, dimension)

def sample_biased_radius_central(tree, dimension, radius=1, sample_ratio=1):
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

def sample_biased_most_central(tree, dimension, sample_ratio=0.1):
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



def sample_unbiased(tree, dimension, sample_ratio=0.1):
    
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

def sample_biased_diagonal_strip(tree, dimension, strip_width=1, sample_ratio=1):
    for node in tree.preorder_node_iter():
        node.leave = False
    for leaf in tree.leaf_node_iter():
        if abs(leaf.X-leaf.Y) <= strip_width/2:
            leaf.leave = True
            
    return delete_tips(tree, dimension)

def sample_biased_diagonal(tree, dimension, sample_ratio=0.1):
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


def sample_biased_horizontal(tree, dimension, sample_ratio=0.1):
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


def sample_biased_vertical(tree, dimension, sample_ratio=0.1):
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


def sample_biased_2_points(tree, dimension, sample_ratio=0.1):
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
