
import argparse
import os


parser = argparse.ArgumentParser(description='Run simulations')
#parser.add_argument('-file', action="store", type=int, dest="file_name", default="output/c_beast/sampled1/beast_input/beast0.xml", help='File name')
#parser.add_argument('-sample_index', action="store", type=int, dest="sample_index", default = 1, help='index of sampling scenario for analysis (default: 1)')
#parser.add_argument('-index', action="store", type=int, dest="index", default = 0, help='index of sample (default: 0)')
parser.add_argument('--re_run', dest='re_run', action='store_const', const=False, default=True, help='is this a re-run, so we should not generate new trees but only run beast again on the previous simulations that did not complete? (default: False)')


args = parser.parse_args()

#index = args.index
#sample_index = args.sample_index

#Should I only run beast again for the cases that did not run before?
reRun_beast_only=args.re_run

#NICOLA: changed here so that only analyses that have not already done are now executed.
#NICOLA: also I have added "ulimit -c unlimited" so that maybe the core dumping problem is solved?
#if sample_index==4:
for index in range(200):
	if (not reRun_beast_only) or ((not os.path.exists("console_output/phyrex_out_"+str(index)+".txt")) or (os.path.getsize("console_output/phyrex_out_"+str(index)+".txt")<10000)):
		if os.path.exists("console_output/phyrex_out_"+str(index)+".txt"):
			os.system("rm console_output/phyrex_out_"+str(index)+".txt")
		os.system('bsub -g /LV_phyrex -J LV'+str(index)+' -o console_output/phyrex_out_'+str(index)+'.txt -e console_output/phyrex_err_'+str(index)+'.txt /nfs/research1/goldman/demaio/phylogeography/phyml/src/phyrex --xml=output/phyrex/LV/phyrex_input/phyrex'+str(index)+'_new.xml')

for index in range(100):
	for sam in range(4):
		if (not reRun_beast_only) or ((not os.path.exists("console_output/phyrex_out"+str(sam+1)+"_"+str(index)+".txt")) or (os.path.getsize("console_output/phyrex_out"+str(sam+1)+"_"+str(index)+".txt")<10000)):
			if os.path.exists("console_output/phyrex_out"+str(sam+1)+"_"+str(index)+".txt"):
				os.system("rm console_output/phyrex_out"+str(sam+1)+"_"+str(index)+".txt")
			os.system('bsub -g /yule/phyrex -J '+str(sam+1)+'s'+str(index)+' -o console_output/phyrex_out'+str(sam+1)+'_'+str(index)+'.txt -e console_output/phyrex_err'+str(sam+1)+'_'+str(index)+'.txt /nfs/research1/goldman/demaio/phylogeography/phyml/src/phyrex --xml=output/phyrex/sampled'+str(sam+1)+'/phyrex_input/phyrex'+str(index)+'.xml')


#else:
#	if (not reRun_beast_only) or ((not os.path.exists("console_output/phyrex_out"+str(sample_index+1)+"_"+str(index)+".txt")) or (os.path.getsize("console_output/phyrex_out"+str(sample_index+1)+"_"+str(index)+".txt")<1000)):
#		os.system('/nfs/research1/goldman/demaio/phylogeography/phyml/src/phyrex --xml=output/phyrex/sampled'+str(sample_index+1)+'/phyrex_input/phyrex'+str(index)+'.xml')
	
exit()
	