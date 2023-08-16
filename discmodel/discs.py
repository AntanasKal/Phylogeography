import ercs
import discsim
import numpy as np
import random
import math
import argparse
import newick
import discs_phyrexxmlwriter
import beastxmlwriter
import dendropy
import treegenerator
import os


#tree by this is not completely ultrametric
#there are some small on the branch ages
#this functions make so branches longer if needed
def ultrametrize(tree):
	#print(tree.as_ascii_plot())
	max_time =0
	for leaf in tree.leaf_node_iter():
		print(leaf.time)
		print(max_time)
		max_time = max(leaf.time, max_time)
	for leaf in tree.leaf_node_iter():
		leaf.edge_length = leaf.edge_length +max_time-leaf.time
		leaf.time = max_time
	return tree

parser = argparse.ArgumentParser(description='Run simulations')
parser.add_argument('-jobi', action="store", type=int, dest="job_index", default=1, help='job index')
parser.add_argument('-N', action="store", type=int, dest="num_simulations", default=1, help='number of simulations (default 1)')
parser.add_argument('--re_run', dest='re_run', action='store_const', const=True, default=False, help='is this a re-run, so we should not generate new trees but only run beast again on the previous simulations that did not complete? (default: False)')

args = parser.parse_args()

job_index = args.job_index
num_simulations = args.num_simulations
#Should I only run beast again for the cases that did not run before?
reRun_beast_only=args.re_run





if not os.path.exists("output"):
	os.makedirs("output")
if not os.path.exists("output/beast"):
	os.makedirs("output/beast")
if not os.path.exists("output/beast/LV"):
	os.makedirs("output/beast/LV")    
if not os.path.exists("output/beast/LV/beast_input"):
	os.makedirs("output/beast/LV/beast_input")
if not os.path.exists("output/beast/LV/beast_output"):
	os.makedirs("output/beast/LV/beast_output")

if not os.path.exists("output/phyrex"):
	os.makedirs("output/phyrex")

if not os.path.exists("output/phyrex/LV"):
	os.makedirs("output/phyrex/LV")        
if not os.path.exists("output/phyrex/LV/phyrex_output"):
	os.makedirs("output/phyrex/LV/phyrex_output")
if not os.path.exists("output/phyrex/LV/phyrex_input"):
	os.makedirs("output/phyrex/LV/phyrex_input")
if not os.path.exists("output/LV/root_data"):
	os.makedirs("output/LV/root_data")

for index in range(job_index*num_simulations, (job_index+1)*num_simulations):
	if (not reRun_beast_only) or (not os.path.exists("output/phyrex/LV/phyrex_input/phyrex"+str(index)+".xml")):
		L = 100
		''' R is the diameter of the torus we are simulating on. 
			This defines the size of the 1D or 2D space that lineages
			can move around in.'''
		sim = discsim.Simulator(L)
	
		#array for sample locations
		a = [None]

		#number of samples
		n = 100
		x = np.zeros(n) 
		y = np.zeros(n) 
		for i in range(n):
			if index < 100:
				x[i] = (random.uniform(25, 75))%L
				y[i] = (random.uniform(25, 75))%L
				a.append((x[i], y[i]))
			else:
				x[i] = (random.uniform(45, 55))%L
				y[i] = (random.uniform(45, 55))%L
				a.append((x[i], y[i]))

		sim.sample = a

		rad = 0.1    

		mu=0.1
		lambda1=2*L**2/(mu*rad**4*math.pi)

		sim.event_classes = [ercs.DiscEventClass(r = rad, u = mu, rate = lambda1)]
		''' All individuals within distance r of the centre of an event
			have probability u of dying in the event 
			and parents are thrown down uniformly within this disc.'''
	
	
		sim.run()
	
	
		pi, tau = sim.get_history()
		al = sim.get_population()[0][0]


		file = open("output/LV/root_data/actual_root"+str(index)+".txt", "w")
		file.write(str(al[0])+"\n"+str(al[1]))
		file.close()

		tree_newick = newick.convert_to_newick(pi, tau, True)
		tree = dendropy.Tree.get(data=tree_newick, schema="newick")

		tree = treegenerator.calculate_times(tree)
		tree = ultrametrize(tree)
		for leaf_index in range(1,n+1):
			leaf_label ="s"+str(leaf_index)
			if leaf_index < 10:
				leaf_label = "s000"+str(leaf_index)
			elif leaf_index < 100:
				leaf_label = "s00"+str(leaf_index)
			elif leaf_index < 1000:
				leaf_label = "s0"+str(leaf_index)	
			node = tree.find_node_with_taxon_label(leaf_label)
			node.X = a[leaf_index][0]
			node.Y = a[leaf_index][1]
			node.annotations.add_bound_attribute("X")
			node.annotations.add_bound_attribute("Y")	
			node.annotations.add_bound_attribute("time")	
		beastxmlwriter.write_BEAST_xml(tree, index, dimension=2, mcmc=10000, log_every=10, beast_input_string="output/beast/LV/beast_input/beast", beast_output_string="output/beast/LV/beast_output/beast")
		discs_phyrexxmlwriter.write_phyrex_input(tree, index, input_string="output/phyrex/LV/phyrex_input/" , output_string="output/phyrex/LV/phyrex_output/", bound=L)
		
		os.system('ulimit -c unlimited; beast -overwrite -seed 123456795 "output/beast/LV/beast_input/beast'+str(index)+'.xml"')

	#line to run BEAST
	#NICOLA: added the "ulimit -c unlimited;" again to try to avoid core dumping
	elif (not reRun_beast_only) or (not os.path.exists("output/beast/LV/beast_output/beast"+str(index)+".log.txt")):
		os.system('ulimit -c unlimited; beast -overwrite -seed 123456795 "output/beast/LV/beast_input/beast'+str(index)+'.xml"')
		
