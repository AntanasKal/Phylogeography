import ercs
import discsim
import numpy as np
import random
import math
import argparse

parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-jobi', action="store", type=int, dest="job_index", default=1, help='job index')
parser.add_argument('-N', action="store", type=int, dest="num_simulations", default=1, help='number of simulations (default 1)')

args = parser.parse_args()

job_index = args.job_index
num_simulations = args.num_simulations


file = open("output/discs"+str(job_index)+".txt", "w")


for index in range((job_index-1)*num_simulations, job_index*num_simulations):
	L = 50
	''' R is the diameter of the torus we are simulating on. 
	    This defines the size of the 1D or 2D space that lineages
	    can move around in.'''
	sim = discsim.Simulator(L)
	a = [None]
	n = 2
	#x = np.zeros(n) 
	#y = np.zeros(n) 
	#for i in range(n):
	#    x[i] = 2*random.random() + 4
	#    y[i] = 2*random.random() + 4
	#    a.append((x[i], y[i]))
	a=[None, (20,25), (30, 25)]
	sim.sample = a
	print(a)
	rad = 0.1    
	#lambda1 = 5000000
	#mu = R**2/(lambda1*4*math.pi*rad**4)
	#mu = R**2/(lambda1*rad**4)
	#mu= R**2/(lambda1*rad**4*math.pi*(1-(128/(45*math.pi))**2))
	mu=0.9
	lambda1=L**2/(mu*rad**4*math.pi*(1-(128/(45*math.pi))**2))
	print(index)
	print("lambda = "+str(lambda1))
	sim.event_classes = [ercs.DiscEventClass(r = rad, u = mu, rate = lambda1)]
	''' All individuals within distance r of the centre of an event
	    have probability u of dying in the event 
	    and parents are thrown down uniformly within this disc.'''
	sim.run()
	pi, tau = sim.get_history()
	al = sim.get_population()[0][0]
	print(sim.get_population())
	file.write(str(tau[0][3])+"\n")
	print("mu = "+str(mu))
	print("rad = "+str(rad))


	print(str(al))
	print(str(pi))
	print(str(tau))

file.close()


