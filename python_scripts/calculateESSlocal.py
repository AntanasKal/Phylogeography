# -*- coding: utf-8 -*-
"""
Created on Sep 30 2020

@author: Nicola De Maio
"""

import argparse
#import dendropy
import os
import numpy as np
from scipy.stats import norm
import pymc3 as pm
import tensorflow as tf
import tensorflow_probability as tfp
tfd = tfp.distributions

# target = tfd.MultivariateNormalDiag(scale_diag=[1., 2.])
# # Get 1000 states from one chain.
# states = tfp.mcmc.sample_chain(
#     num_burnin_steps=20,
#     num_results=100,
#     current_state=tf.constant([0., 0.]),
#     trace_fn=None,
#     kernel=tfp.mcmc.HamiltonianMonteCarlo(
#       target_log_prob_fn=target.log_prob,
#       step_size=0.05,
#       num_leapfrog_steps=2))
#print(states)
#print(tfp.mcmc.effective_sample_size(states))
#print(tfp.mcmc.effective_sample_size(states)[0])

N=100

analyses=["No Bias BMP, beast","Central Sampling BMP, beast", "Diagonal Sampling BMP, beast", "One sided sampling BMP, beast", "No Bias BMP, beast with extra samples","Central Sampling BMP, beast with extra samples", "Diagonal Sampling BMP, beast with extra samples", "One sided sampling BMP, beast with extra samples", "Broad sampling LFV, beast", "Narrow sampling LFV, beast", "No Bias BMP, phyrex","Central Sampling BMP, phyrex", "Diagonal Sampling BMP, phyrex", "One sided sampling BMP, phyrex", "Broad sampling LFV, phyrex", "Narrow sampling LFV, phyrex", "No Bias BMP, beast without extra samples","Central Sampling BMP, beast without extra samples", "Diagonal Sampling BMP, beast without extra samples", "One sided sampling BMP, beast without extra samples"]
foldersBEAST=["output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/", "output/beast/LV/","output/beast/LV/","output/phyrex/sampled1/","output/phyrex/sampled2/","output/phyrex/sampled3/","output/phyrex/sampled4/", "output/phyrex/LV/","output/phyrex/LV/","output/unc_beast/sampled1/","output/unc_beast/sampled2/","output/unc_beast/sampled3/","output/unc_beast/sampled4/"]
foldersRoots=["output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/","output/LV/","output/LV/","output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/LV/","output/LV/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/"]
toSkip=[0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,100,0,0,0,0]



#code for computing ESS from https://github.com/jwalton3141/jwalton3141.github.io/blob/master/assets/posts/ESS/rwmh.py
def my_ESS(x):
    """ Compute the effective sample size of estimand of interest. Vectorised implementation. """
    m_chains, n_iters = x.shape
    variogram = lambda t: ((x[:, t:] - x[:, :(n_iters - t)])**2).sum() / (m_chains * (n_iters - t))
    post_var = my_gelman_rubin(x)
    t = 1
    rho = np.ones(n_iters)
    negative_autocorr = False
    # Iterate until the sum of consecutive estimates of autocorrelation is negative
    while not negative_autocorr and (t < n_iters):
        rho[t] = 1 - variogram(t) / (2 * post_var)
        if not t % 2:
            negative_autocorr = sum(rho[t-1:t+1]) < 0
        t += 1
    return int(m_chains*n_iters / (1 + 2*rho[1:t].sum()))

def my_gelman_rubin(x):
    """ Estimate the marginal posterior variance. Vectorised implementation. """
    m_chains, n_iters = x.shape
    # Calculate between-chain variance
    B_over_n = ((np.mean(x, axis=1) - np.mean(x))**2).sum() / (m_chains - 1)
    # Calculate within-chain variances
    W = ((x - x.mean(axis=1, keepdims=True))**2).sum() / (m_chains*(n_iters - 1))
    # (over) estimate of variance
    s2 = W * (n_iters - 1) / n_iters + B_over_n
    return s2
		
path="/Users/demaio/Desktop/phylogeography_2019_summer_project/check_LFV_convergence/"
#Phyrex
analyses=["No Bias BMP, phyrex", "Central Sampling BMP, phyrex", "Diagonal Sampling BMP, phyrex", "One sided sampling BMP, phyrex", "Broad sampling LFV, phyrex", "Narrow sampling LFV, phyrex"]
foldersBEAST=["NoBias_phyrex/","other_phyrex2/","other_phyrex3/","other_phyrex4/", "FLV_phyrex/","FLV_phyrex/"]
#foldersRoots=["output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/","output/LV/","output/LV/","output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/LV/","output/LV/"]
toSkip=[0,0,0,0,0,100]
lessThan200=0
lessThan100=0
total=0
for i in range(len(analyses)):
	print("\n")
	print(analyses[i])
	for j in range(N):
		skip=False
		#print(j)
		sigmas=[]
		sigmas2=[]
		sigmas3=[]
		sigmas4=[]
		lnPs=[]
		rootLats=[]
		rootLons=[]
		#file=open(foldersBEAST[i]+"phyrex_output/out_new2_phyrex_stats_"+str(j)+".txt")
		file=open(path+foldersBEAST[i]+"out_new2_phyrex_stats_"+str(j+toSkip[i])+"_extract.txt")
		line=file.readline()
		linelist=line.split()
		if line!="" and line!="\n":
			rootLoni=linelist.index("rootLon")
			rootLati=linelist.index("rootLat")
			lnPi=linelist.index("lnP")
			sigsqi=linelist.index("sigSq")
			sigsqobsi=linelist.index("realsigsqroot")
			#sigsqobs2i=linelist.index("realsigsqtips")
			sigsqobs2i=linelist.index("realsigsqtipsbis")
			sigsqobs3i=linelist.index("realsigsqtipster")
			line=file.readline()
			#count=0
			#fileO.write("rootLon\t"+"rootLat\t"+"lnP\t"+"sigSq\t"+"realsigsqroot\t"+"realsigsqtipsbis\t"+"realsigsqtipster\n")
			while line!="" and line!="\n":
				#count+=1
				#if count%thinPhyrex==0:
				linelist=line.split()
				sigsq=float(linelist[sigsqi])
				sigsqobs=float(linelist[sigsqobsi])
				sigsqobs3=float(linelist[sigsqobs2i])
				sigsqobs4=float(linelist[sigsqobs3i])
				lnP=float(linelist[lnPi])
				rootLat=float(linelist[rootLati])
				rootLon=float(linelist[rootLoni])
				sigmas.append(sigsq)
				sigmas2.append(sigsqobs)
				sigmas3.append(sigsqobs3)
				sigmas4.append(sigsqobs4)
				lnPs.append(lnP)
				rootLats.append(rootLat)
				rootLons.append(rootLon)
					#fileO.write(linelist[rootLoni]+"\t"+linelist[rootLati]+"\t"+linelist[lnPi]+"\t"+linelist[sigsqi]+"\t"+linelist[sigsqobsi]+"\t"+linelist[sigsqobs2i]+"\t"+linelist[sigsqobs3i]+"\n")
				line=file.readline()
		
			file.close()
			values=[sigmas,sigmas2,sigmas3,sigmas4,lnPs,rootLats,rootLons]
		
			#print(j)
			ESSs=[]
			ESSs2=[]
			found=False
			found2=False
			for v in range(len(values)):
				array1=np.array(values[v])
				lenAr=len(array1)
				#print(lenAr)
				array1.shape=(lenAr,1)
				#print(my_ESS(array1))
				ESS=tfp.mcmc.effective_sample_size(array1)
				if ESS[0].numpy()<200:
					found=True
				if ESS[0].numpy()<100:
					found2=True
				#ESS2=tfp.mcmc.effective_sample_size(array1[500:])
				#print(tfp.mcmc.effective_sample_size(array1))
				ESSs.append(ESS.numpy()[0])
				#ESSs2.append(ESS2.numpy()[0])
				#print(pm.effective_n(values[v]))
		else:
				found=True
				found2=True
				ESSs=[0.0,0.0,0.0,0.0,0.0,0.0,0.0]
		file=open(path+foldersBEAST[i]+"out_new2_phyrex_stats_"+str(j+toSkip[i])+"_ESS.txt","w")
		if found:
			lessThan200+=1
			print("\n")
			print(str(j))
			print(ESSs)
			#print(ESSs2)
		if found2:
			lessThan100+=1
			file.write("0")
		else:
			file.write("1")
		file.close()
		total+=1
		
print(total)
print(lessThan200)
print(lessThan100)
exit()
	

	
	