print("test")
import ercs
import discsim
import time
import math
from random import random
import numpy as np
print("test2") 
R = 10
''' R is the diameter of the torus we are simulating on. 
    This defines the size of the 1D or 2D space that lineages
    can move around in.'''
sim = discsim.Simulator(R)
a = [None]
n = 100
x = np.zeros(n) 
y = np.zeros(n) 
for i in range(n):
    x[i] = random() + 4.5
    y[i] = random() + 4.5
    a.append((x[i], y[i]))
sim.sample = a
rad = 0.5    
lambda1 = 1e3
mu = 1/(lambda1*R**2*(1 - (128/(45*math.pi))**2)*math.pi*rad**4)
print("test3")
sim.event_classes = [ercs.DiscEventClass(r = rad, u = mu, rate = lambda1)]
''' All individuals within distance r of the centre of an event
    have probability u of dying in the event 
    and parents are thrown down uniformly within this disc.'''
sim.run()
pi, tau = sim.get_history()
al = sim.get_population()[0][0]
print("test4")
file = open("output102.txt", "w")
file.write("ancestor location  = "+str( al)+"\n")

### If we do not want oriented tree format, just delete the next two lines
file.write("pi  = "+str(pi)+"\n")
file.write("tau = "+str(tau)+"\n")
##print the time used

file.close()

#import os
#os.system('python newick.py -b -p "%s" -t "%s"' %(pi, tau))
