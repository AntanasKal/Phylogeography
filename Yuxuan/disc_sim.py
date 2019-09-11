#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  1 14:56:35 2019

@author: Yuxuan
"""

import ercs
import discsim
import time
import math
from random import random
import numpy as np
 
start = time.clock() ##record the time of running the programme
rad = 10
''' rad is the diameter of the torus we are simulating on. 
    This defines the size of the 1D or 2D space that lineages
    can move around in.'''

sim = discsim.Simulator(rad)

a = [None]

##sample n points uniformly in a unit square and centering at (rad/2,rad/2)
n = 22
x = np.zeros(n) 
y = np.zeros(n) 
for i in range(n):
    x[i] = random() + rad/2 - 0.5
    y[i] = random() + rad/2 - 0.5 
    a.append((x[i], y[i]))
sim.sample = a

##According to the Latex file 'Diffusion_rate' provided, 
##we are going to simulate the event class under disc model
R = 0.01    
lambda1 = 1e11
mu = rad**2/(lambda1*(1 - (128/(45*math.pi))**2)*math.pi*R**4)

sim.event_classes = [ercs.DiscEventClass(r = R, u = mu, rate = lambda1)]
''' All individuals within distance R of the centre of an event
    have probability mu of dying in the event 
    and parents are thrown down uniformly within this disc.'''

sim.run()

##get the tree and the root of location of the root
pi, tau = sim.get_history()
al = sim.get_population()[0][0]

print("root location  = ", al)

### If we do not want oriented tree format, just delete the next two lines
print("pi  = ", pi)
print("tau = ", tau)

##print the time used
elapsed = (time.clock() - start)
print("Time used:",elapsed)

##transfer the tree format into newick format
##The python script newick.py can be downloaded at https://github.com/tyjo/newick.py
import os
os.system('python newick.py -b -p "%s" -t "%s"' %(pi, tau))

