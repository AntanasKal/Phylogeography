# -*- coding: utf-8 -*-
"""
Created on Oct 4 2019

@author: Nicola De Maio
"""

import argparse
#import dendropy
import os
import numpy

N=100

#Check the new version of the phyrex output and create files accordingly

rootFiles=True
sigmaFiles=True
BEAST=True
ChangePhyrex=False

analyses=["No Bias BMP, beast","Central Sampling BMP, beast", "Diagonal Sampling BMP, beast", "One sided sampling BMP, beast", "No Bias BMP, beast with extra samples","Central Sampling BMP, beast with extra samples", "Diagonal Sampling BMP, beast with extra samples", "One sided sampling BMP, beast with extra samples", "Broad sampling LFV, beast", "Narrow sampling LFV, beast", "No Bias BMP, phyrex","Central Sampling BMP, phyrex", "Diagonal Sampling BMP, phyrex", "One sided sampling BMP, phyrex", "Broad sampling LFV, phyrex", "Narrow sampling LFV, phyrex", "No Bias BMP, beast without extra samples","Central Sampling BMP, beast without extra samples", "Diagonal Sampling BMP, beast without extra samples", "One sided sampling BMP, beast without extra samples"]
foldersBEAST=["output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/", "output/beast/LV/","output/beast/LV/","output/phyrex/sampled1/","output/phyrex/sampled2/","output/phyrex/sampled3/","output/phyrex/sampled4/", "output/phyrex/LV/","output/phyrex/LV/","output/unc_beast/sampled1/","output/unc_beast/sampled2/","output/unc_beast/sampled3/","output/unc_beast/sampled4/"]
foldersRoots=["output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/","output/LV/","output/LV/","output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/LV/","output/LV/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/"]
toSkip=[0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,100,0,0,0,0]

if rootFiles:
	for i in range(len(analyses)):
		print("\n")
		print(analyses[i])
		fileX=open("plots/"+analyses[i].replace(" ","_")+"_rootX.txt","w")
		fileY=open("plots/"+analyses[i].replace(" ","_")+"_rootY.txt","w")
		print("plots/"+analyses[i].replace(" ","_")+"_rootX.txt")
		for j in range(N):
			#print(j)
			file=open(foldersRoots[i]+"root_data/actual_root"+str(j+toSkip[i])+".txt")
			line=file.readline()
			xLoc=float(line.replace("\n",""))
			line=file.readline()
			yLoc=float(line.replace("\n",""))
			file.close()
			trueRoots=[xLoc,yLoc]
		
			estimatedXRoots=[]
			estimatedYRoots=[]
			if i<8 or i>15:
				file=open(foldersBEAST[i]+"root_data/observed_roots"+str(j)+".txt")
				#print(foldersBEAST[i]+"root_data/observed_roots"+str(j)+".txt")
				line=file.readline()
				while line!="" and line!="\n":
					linelist=line.split()
					estimatedXRoot=float(linelist[0])-trueRoots[0]
					estimatedYRoot=float(linelist[1])-trueRoots[1]
					estimatedXRoots.append(estimatedXRoot)
					estimatedYRoots.append(estimatedYRoot)
					fileX.write(str(estimatedXRoot))
					fileY.write(str(estimatedYRoot))
					fileX.write("\t")
					fileY.write("\t")
					line=file.readline()
				fileX.write("\n")
				fileY.write("\n")
				file.close()
			elif i==8 or i==9:
				file=open(foldersBEAST[i]+"beast_output/beast"+str(j+toSkip[i])+".trees.txt")
				#print(foldersBEAST[i]+"beast_output/beast"+str(j+toSkip[i])+".trees.txt")
				line=file.readline()
				linelist=line.split()
				while len(linelist)<2 or linelist[0]!="tree":
					line=file.readline()
					linelist=line.split()
					if line=="":
						break
				while line!="" and line!="\n" and len(linelist)>1:
					#print(line)
					locs=linelist[-1].split("{")[-1].replace("}];","")
					estimatedXRoot=float(locs.split(",")[0])-trueRoots[0]
					estimatedYRoot=float(locs.split(",")[1])-trueRoots[1]
					estimatedXRoots.append(estimatedXRoot)
					estimatedYRoots.append(estimatedYRoot)
					fileX.write(str(estimatedXRoot))
					fileY.write(str(estimatedYRoot))
					fileX.write("\t")
					fileY.write("\t")
					line=file.readline()
					linelist=line.split()
				fileX.write("\n")
				fileY.write("\n")
				file.close()
			if i>9 and i<=15:
				file=open(foldersBEAST[i]+"phyrex_output/out_phyrex_stats_"+str(j+toSkip[i])+".txt")
				#print(foldersBEAST[i]+"phyrex_output/out_phyrex_stats_"+str(j+toSkip[i])+".txt")
				line=file.readline()
				linelist=line.split()
				while len(linelist)<2 or linelist[0]!="sample":
					line=file.readline()
					linelist=line.split()
					if line=="":
						break
				line=file.readline()
				linelist=line.split()
				while line!="" and line!="\n":
					estimatedXRoot=float(linelist[19].replace("f",""))-trueRoots[0]
					estimatedYRoot=float(linelist[20])-trueRoots[1]
					estimatedXRoots.append(estimatedXRoot)
					estimatedYRoots.append(estimatedYRoot)
					fileX.write(str(estimatedXRoot))
					fileY.write(str(estimatedYRoot))
					fileX.write("\t")
					fileY.write("\t")
					line=file.readline()
					linelist=line.split()
				fileX.write("\n")
				fileY.write("\n")
				file.close()




			#print(len(estimatedXRoots))
			#if len(estimatedXRoots)>100:
			#	print(str(numpy.percentile(estimatedXRoots, 2.5))+" "+str(numpy.percentile(estimatedXRoots, 50))+" "+str(numpy.percentile(estimatedXRoots, 97.5)))
			#	print(str(numpy.percentile(estimatedYRoots, 2.5))+" "+str(numpy.percentile(estimatedYRoots, 50))+" "+str(numpy.percentile(estimatedYRoots, 97.5)))
			#else:
			#	print("Less than 100 samples!")
		fileX.close()
		fileY.close()


if sigmaFiles:
	if BEAST:
		analyses=["No Bias BMP, beast","Central Sampling BMP, beast", "Diagonal Sampling BMP, beast", "One sided sampling BMP, beast", "No Bias BMP, beast with extra samples","Central Sampling BMP, beast with extra samples", "Diagonal Sampling BMP, beast with extra samples", "One sided sampling BMP, beast with extra samples", "No Bias BMP, beast without extra samples","Central Sampling BMP, beast without extra samples", "Diagonal Sampling BMP, beast without extra samples", "One sided sampling BMP, beast without extra samples", "Broad sampling LFV, beast", "Narrow sampling LFV, beast"]
		foldersBEAST=["output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/","output/unc_beast/sampled1/","output/unc_beast/sampled2/","output/unc_beast/sampled3/","output/unc_beast/sampled4/", "output/beast/LV/","output/beast/LV/"]
		#foldersRoots=["output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/","output/LV/","output/LV/"]
		truth=[1.0,0.0,1.0]
		toSkip=[0,0,0,0,0,0,0,0,0,0,0,0,0,100]
		for i in range(len(analyses)):
			print("\n")
			print(analyses[i])
			fileX=open("plots/"+analyses[i].replace(" ","_")+"_sigmaX.txt","w")
			fileCorr=open("plots/"+analyses[i].replace(" ","_")+"_Corr.txt","w")
			fileY=open("plots/"+analyses[i].replace(" ","_")+"_sigmaY.txt","w")
			#fileSigma=open("plots/"+analyses[i].replace(" ","_")+"_Sigma.txt","w")
			#fileXoverY=open("plots/"+analyses[i].replace(" ","_")+"_XoverY.txt","w")
			for j in range(N):
				#print(j)
		
				estimatedSigmaX=[]
				estimatedSigmaY=[]
				estimatedCov=[]
				estimatedSigma=[]
				estimatedCorr=[]
				treeLs=[]
				XoverYs=[]
				if i<12:
					file=open(foldersBEAST[i]+"beast_output/beast"+str(j)+".log.txt")
				elif i==12 or i==13:
					file=open(foldersBEAST[i]+"beast_output/beast"+str(j+toSkip[i])+".log.txt")
				for k in range(4):
					line=file.readline()
				linelist=line.split()
				treeLi=linelist.index("treeLength")
				corri=linelist.index("correlation")
				xi=linelist.index("location.varCovar.location.precision.col11")
				yi=linelist.index("location.varCovar.location.precision.col22")
				#covi=linelist.index("location.varCovar.location.precision.col12")
				#diffi=linelist.index("location.diffusionRate")
				line=file.readline()
				while line!="" and line!="\n":
					linelist=line.split()
					estimatedX=float(linelist[xi])/float(linelist[treeLi])
					estimatedY=float(linelist[yi])/float(linelist[treeLi])
					#estimatedC=float(linelist[covi])/float(linelist[treeLi])
					estimatedC2=float(linelist[corri])
					#estimatedD=float(linelist[diffi])
					treeL=float(linelist[treeLi])
					estimatedSigmaX.append(estimatedX)
					estimatedSigmaY.append(estimatedY)
					#estimatedCov.append(estimatedC)
					estimatedCorr.append(estimatedC2)
					#estimatedSigma.append(estimatedD)
					treeLs.append(treeL)
					fileX.write(str(estimatedX))
					fileY.write(str(estimatedY))
					fileCorr.write(str(estimatedC2))
					#fileSigma.write(str(estimatedD))
					#fileXoverY.write(str(estimatedX/estimatedY))
					#XoverYs.append(estimatedX/estimatedY)
					fileX.write("\t")
					fileY.write("\t")
					fileCorr.write("\t")
					#fileSigma.write("\t")
					#fileXoverY.write("\t")
					line=file.readline()
				fileX.write("\n")
				fileY.write("\n")
				fileCorr.write("\n")
				#fileSigma.write("\n")
				#fileXoverY.write("\n")
				file.close()

				#print(len(estimatedSigmaX))
				#if len(estimatedSigmaX)>100:
				#	print(str(numpy.percentile(estimatedSigmaX, 2.5))+" "+str(numpy.percentile(estimatedSigmaX, 50))+" "+str(numpy.percentile(estimatedSigmaX, 97.5)))
				#	print(str(numpy.percentile(estimatedSigmaY, 2.5))+" "+str(numpy.percentile(estimatedSigmaY, 50))+" "+str(numpy.percentile(estimatedSigmaY, 97.5)))
					#print(str(numpy.percentile(estimatedCov, 2.5))+" "+str(numpy.percentile(estimatedCov, 50))+" "+str(numpy.percentile(estimatedCov, 97.5)))
				#	print(str(numpy.percentile(estimatedCorr, 2.5))+" "+str(numpy.percentile(estimatedCorr, 50))+" "+str(numpy.percentile(estimatedCorr, 97.5)))
					#print(str(numpy.percentile(estimatedSigma, 2.5))+" "+str(numpy.percentile(estimatedSigma, 50))+" "+str(numpy.percentile(estimatedSigma, 97.5)))
					#print(str(numpy.percentile(XoverYs, 2.5))+" "+str(numpy.percentile(XoverYs, 50))+" "+str(numpy.percentile(XoverYs, 97.5)))
					#print(str(numpy.percentile(treeLs, 2.5))+" "+str(numpy.percentile(treeLs, 50))+" "+str(numpy.percentile(treeLs, 97.5)))
				#else:
				#	print("Less than 100 samples!")
			fileX.close()
			fileY.close()
			fileCorr.close()
			#fileSigma.close()
			#fileXoverY.close()
		
		
	
	#Phyrex
	analyses=["No Bias BMP, phyrex", "Central Sampling BMP, phyrex", "Diagonal Sampling BMP, phyrex", "One sided sampling BMP, phyrex", "Broad sampling LFV, phyrex", "Narrow sampling LFV, phyrex"]
	foldersBEAST=["output/phyrex/sampled1/","output/phyrex/sampled2/","output/phyrex/sampled3/","output/phyrex/sampled4/", "output/phyrex/LV/","output/phyrex/LV/"]
	#foldersRoots=["output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/c_beast/sampled1/","output/c_beast/sampled2/","output/c_beast/sampled3/","output/c_beast/sampled4/","output/LV/","output/LV/","output/beast/sampled1/","output/beast/sampled2/","output/beast/sampled3/","output/beast/sampled4/","output/LV/","output/LV/"]
	toSkip=[0,0,0,0,0,100]
	for i in range(len(analyses)):
		print("\n")
		print(analyses[i])
		fileSigma=open("plots/"+analyses[i].replace(" ","_")+"_SigmaFormula.txt","w")
		fileSigma2=open("plots/"+analyses[i].replace(" ","_")+"_SigmaObservedRoot.txt","w")
		fileSigma3=open("plots/"+analyses[i].replace(" ","_")+"_SigmaObservedTips.txt","w")
		for j in range(N):
			skip=False
			#print(j)
			sigmas=[]
			sigmas2=[]
			sigmas3=[]
			file=open(foldersBEAST[i]+"phyrex_output/out_phyrex_stats_"+str(j)+".txt")
			line=file.readline()
			linelist=line.split()
			while len(linelist)<3 or linelist[0]!="sample":
				line=file.readline()
				linelist=line.split()
				if line=="":
					skip=True
					break
			if skip:
				continue
			sigsqi=linelist.index("sigsq")
			sigsqobsi=linelist.index("realsigsqroot")
			sigsqobs2i=linelist.index("realsigsqtipsbis")
			line=file.readline()
			while line!="" and line!="\n":
				linelist=line.split()
				sigsq=float(linelist[sigsqi])
				sigsqobs=float(linelist[sigsqobsi])
				sigsqobs3=float(linelist[sigsqobs2i])
				sigmas.append(sigsq)
				sigmas2.append(sigsqobs)
				sigmas3.append(sigsqobs3)
				fileSigma.write(str(sigsq))
				fileSigma2.write(str(sigsqobs))
				fileSigma3.write(str(sigsqobs3))
				fileSigma.write("\t")
				fileSigma2.write("\t")
				fileSigma3.write("\t")
				line=file.readline()
			fileSigma.write("\n")
			fileSigma2.write("\n")
			fileSigma3.write("\n")
			file.close()

			#print(len(sigmas))
			#if len(sigmas)>100:
			#	print(str(numpy.percentile(sigmas, 2.5))+" "+str(numpy.percentile(sigmas, 50))+" "+str(numpy.percentile(sigmas, 97.5)))
			#	print(str(numpy.percentile(sigmas2, 2.5))+" "+str(numpy.percentile(sigmas2, 50))+" "+str(numpy.percentile(sigmas2, 97.5)))
			#	print(str(numpy.percentile(sigmas3, 2.5))+" "+str(numpy.percentile(sigmas3, 50))+" "+str(numpy.percentile(sigmas3, 97.5)))
				#print(str(numpy.percentile(treeLs, 2.5))+" "+str(numpy.percentile(treeLs, 50))+" "+str(numpy.percentile(treeLs, 97.5)))
			#else:
			#	print("Less than 100 samples!")
		fileSigma.close()
		fileSigma2.close()
		fileSigma3.close()
		
if ChangePhyrex:
	analyses=["Broad sampling LFV, phyrex", "Narrow sampling LFV, phyrex"]
	foldersBEAST=["output/phyrex/LV/","output/phyrex/LV/"]
	toSkip=[0,100]
	for i in range(len(analyses)):
		for j in range(N):
			file=open(foldersBEAST[i]+"phyrex_input/phyrex"+str(j+toSkip[i])+".xml")
			fileNew=open(foldersBEAST[i]+"phyrex_input/phyrex"+str(j+toSkip[i])+"_new.xml","w")
			text=file.read()
			#print(line)
			#print(line.replace('mcmc.sample.every="1000"','mcmc.sample.every="100"'))
			textNew=text.replace('mcmc.sample.every="1000"','mcmc.sample.every="100"').replace('mcmc.burnin="10000"','mcmc.burnin="1000"')
			file.close()
			fileNew.write(textNew)
			fileNew.close()
	
	
	
	
	