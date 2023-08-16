# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 11:18:17 2019

@author: Antanas
"""

import argparse
import dendropy
import numpy as np

import os

#NICOLA: I created this variant to run beast on genetic sequences but without the extra samples, for better comparison.

parser = argparse.ArgumentParser(description='Run simulations')
#parser.add_argument('-file', action="store", type=int, dest="file_name", default="output/c_beast/sampled1/beast_input/beast0.xml", help='File name')
#parser.add_argument('-sample_index', action="store", type=int, dest="sample_index", default = 1, help='index of sampling scenario for analysis (default: 1)')
#parser.add_argument('-index', action="store", type=int, dest="index", default = 0, help='index of sample (default: 0)')
parser.add_argument('--re_run', dest='re_run', action='store_const', const=True, default=False, help='is this a re-run, so we should not generate new trees but only run beast again on the previous simulations that did not complete? (default: False)')


args = parser.parse_args()

#index = args.index
#sample_index = args.sample_index

#Should I only run beast again for the cases that did not run before?
reRun_beast_only=args.re_run

for s in range(4):
	sample_index=s+1
	for index in range(100):

		#NICOLA: changed here so that only analyses that have not already done are now executed.
		#NICOLA: also I have added "ulimit -c unlimited" so that maybe the core dumping problem is solved?
		if (not reRun_beast_only) or ((not os.path.exists("output/unc_beast/sampled"+str(sample_index)+"/beast_output/beast"+str(index)+".log.txt")) or (os.path.getsize("output/unc_beast/sampled"+str(sample_index)+"/beast_output/beast"+str(index)+".log.txt")<1000)):
			if (not os.path.exists("output/unc_beast/sampled"+str(sample_index)+"/beast_input/beast"+str(index)+".xml")):
				file=open("output/c_beast/sampled"+str(sample_index)+"/beast_input/beast"+str(index)+".xml")
				fileO=open("output/unc_beast/sampled"+str(sample_index)+"/beast_input/beast"+str(index)+".xml","w")
				line=file.readline()
				fileO.write(line)
				line=file.readline()
				fileO.write(line)
				line=file.readline()
				fileO.write(line)
				line=file.readline()
				lines={}
				for i in range(100):
					#print(line)
					#print(line.split("\""))
					name=line.split("\"")[1]
					#print(name)
					lines[name]=[]
					for j in range(12):
						lines[name].append(line)
						line=file.readline()
				linesToAdd=[]
				linesToAdd.append(line)
				line=file.readline()
				linesToAdd.append(line)
				line=file.readline()
				lines2={}
				counts={}
				for i in range(100):
					oldline=line
					line=file.readline()
					name=line.split("\"")[1]
					lines2[name]=[]
					lines2[name].append(oldline)
					lines2[name].append(line)
					line=file.readline()
					counts[name]=-1
					counts[name]=line.count('A')
					lines2[name].append(line)
					line=file.readline()
					lines2[name].append(line)
					line=file.readline()
				#print(counts)
				for n in lines.keys():
					if counts[n]>0:
						for l in lines[n]:
							fileO.write(l)
				for l in linesToAdd:
					fileO.write(l)
				for n in lines.keys():
					if counts[n]>0:
						for l in lines2[n]:
							fileO.write(l)
		
				while line!="":
					fileO.write(line)
					line=file.readline().replace("/c_beast","/unc_beast")
					#.replace("10000000","10000")
				#exit()
				file.close()
				fileO.close()
		
			os.system('bsub -g /unc_beast -o console_output/unc_beast_out'+str(sample_index)+'_'+str(index)+'.txt -e console_output/unc_beast_err'+str(sample_index)+'_'+str(index)+'.txt beast -overwrite -seed 123456795 "output/unc_beast/sampled'+str(sample_index)+'/beast_input/beast'+str(index)+'.xml"')
	
			#NICOLA: just take root location
			#print root positions from .trees.txt file to more convenient format
			#only the locations of the roots are important from the .trees.txt file
			file = open("output/unc_beast/sampled"+str(sample_index)+"/root_data/observed_roots"+str(index)+".txt", "w")   
			for line in open("output/unc_beast/sampled"+str(sample_index)+"/beast_output/beast"+str(index)+".trees.txt"):                
				if line.startswith("tree"):
					start_index = 0
					while True:
						if line[start_index:start_index+4]=="[&R]":
							break
						start_index=start_index+1
					single_tree=dendropy.Tree.get(data=line[start_index:], schema="newick", extract_comment_metadata=True)
					file.write(single_tree.seed_node.annotations.require_value("location")[0]+"\t"+single_tree.seed_node.annotations.require_value("location")[1]+'\n')
			file.close()
	
	