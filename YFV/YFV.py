import sys
import math
import argparse
import os


file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_coordinates.txt")
line=file.readline()
samples={}
line=file.readline()
while line!="" and line!="\n":
	linelist=line.split()
	samples[linelist[0]]=[float(linelist[1]),float(linelist[2])]
	line=file.readline()
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_sequences.fasta")
line=file.readline()
#line=file.readline()
while line!="" and line!="\n":
	linelist=line.split()
	#print(linelist[0])
	line=file.readline()
	linelist2=line.split()
	#print(samples.keys())
	#print(linelist[0].replace(">","") in samples.keys())
	#print(linelist[0].replace(">",""))
	samples[linelist[0].replace(">","")].append(linelist2[0])
	line=file.readline()
file.close()

#print(samples)

samList=[]



for s in samples.keys():
	newList=list(samples[s])
	newList.append(s)
	#print(s)
	#print(samples[s])
	
	#print(newList)
	samList.append(newList)
	#print(samples[s][0])

#print(samList[0])
samList.sort(key=lambda x: x[0])
#print(samList[0])
southList=[]
northList=[]
downList=[]
boundary=-19
for s in samList:
	downList.append(s)
	if s[0]<-19:
		southList.append(s)
	else:
		northList.append(s)
		
file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_sequences_downSamp.fasta","w")
for s in downList:
	file.write(">"+s[3]+"\n"+s[2]+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_coordinates_downSamp.txt","w")
file.write("traits	lat	long\n")
for s in downList:
	file.write(s[3]+"\t"+str(s[0])+"\t"+str(s[1])+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_sequences_downSamp2.fasta","w")
for s in downList:
	file.write(">"+s[3]+"\n"+s[2]+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_coordinates_downSamp2.txt","w")
file.write("traits	lat	long\n")
for s in downList:
	file.write(s[3]+"\t"+str(s[0])+"\t"+str(s[1])+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_sequences_south.fasta","w")
for s in southList:
	file.write(">"+s[3]+"\n"+s[2]+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_coordinates_south.txt","w")
file.write("traits	lat	long\n")
for s in southList:
	file.write(s[3]+"\t"+str(s[0])+"\t"+str(s[1])+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_sequences_north.fasta","w")
for s in northList:
	file.write(">"+s[3]+"\n"+s[2]+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_coordinates_north.txt","w")
file.write("traits	lat	long\n")
for s in northList:
	file.write(s[3]+"\t"+str(s[0])+"\t"+str(s[1])+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_sequences_north+S.fasta","w")
for s in northList:
	file.write(">"+s[3]+"\n"+s[2]+"\n")
for s in southList:
	file.write(">"+s[3]+"\n"+("N"*len(s[2]))+"\n")
file.close()

file=open("/Users/demaio/Desktop/phylogeography_2019_summer_project/Yuxuan/YFV_input/YFV_sequences_south+N.fasta","w")
for s in southList:
	file.write(">"+s[3]+"\n"+s[2]+"\n")
for s in northList:
	file.write(">"+s[3]+"\n"+("N"*len(s[2]))+"\n")
file.close()
	
exit()



