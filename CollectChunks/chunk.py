import ROOT as R
import os
import subprocess as sp
import shlex
import sys
from samples_Janik import samplelist_Janik
from samples_Markus import samplelist_Markus
import string
import argparse
#from runUtils import checkProxy as cP
import datetime

def main(): 
	#proxy_cmd = 'voms-proxy-init --voms cms' 
	#print "executing: {0}".format(proxy_cmd)
	#os.system(proxy_cmd)

	checkDirExist(directory='./proxy')
	print "copying proxy"
	os.system("cp /tmp/x509up_u3522 proxy/")
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', dest='sample_name', help='Enter the name of the sample to be merged', type=str, metavar = 'SAMPLE_NAME', choices=['SingleElectron','SingleMuon','Tau','VBF','GluGlu','DY','EWK','ST_t','W','ZZ','TTT'],default='')
	parser.add_argument('-t', dest='submit', help='Where to submit', type=str, metavar = 'SUBMIT', choices=['local','hephybatch'],default='hephybatch')
	parser.add_argument('-o', dest='output_directory', help='output directory', type=str, metavar = 'OUTDIR', default='create')
	parser.add_argument('-f', dest='which_file_list', help='which file list', type=str, metavar = 'FILELIST', choices=['Janik','Markus'], default='Janik')
	args = parser.parse_args()
			    
	
	submit = args.submit

	output_choice = args.output_directory
	sample_name = args.sample_name

	if args.which_file_list == 'Janik' :
		samplelist = samplelist_Janik
	else :
		samplelist = samplelist_Markus

	if output_choice == 'create' :
		day_stamp = str(datetime.datetime.today()).split()[0]
		outdir = "/afs/hephy.at/data/higgs02/haddJanik/{0}/".format(day_stamp)
	else :
		if '/afs/' in output_choice :
			outdir = output_choice
		else :
			outdir = "/afs/hephy.at/data/higgs02/haddJanik/{0}/".format(output_choice)

	print 'output directory: {0}'.format(outdir)
	checkDirExist(outdir)

	sources = samplelist.keys()
	sources.sort()
	counter = 0
	for source in sources:
		splsource = source.split("/")
		name = "_".join(splsource[1:3]).replace("RunIIFall17MiniAODv2-","").replace("-", "_")
		if not(sample_name in name) :
			continue

		file_list =  getfiles(samplelist[source][0] )
		n = len(file_list)
		done = (n == samplelist[source][1] )
		
		# check whether sample has finished running
		if done :
			print "sample {0} has finished running".format(name)
			#print "id {0}".format(counter)
			output_dir = outdir+name+'/'
			checkDirExist(output_dir)
			# if the output directory is empty it means that sample has finished
			# and can be merged. If it is not empty it has been already merged.
			if not os.listdir(output_dir) :
    			#print "Directory is empty"
				if sample_name in name :
					print 'going to merge the sample {0}'.format(name)
					merge(outdir=output_dir,files=file_list,submit_run=submit)
					
			else :
				"sample already merged - see: {0}".format(output_dir)
		counter += 1

def checkDirExist (directory) :
	if not os.path.exists(directory) :
    		os.makedirs(directory)

def getfiles(source):

    proc = sp.Popen(shlex.split("dpns-ls -R {0}".format(source) ), stdout=sp.PIPE, stderr=sp.PIPE, shell=False)
    (out, err) = proc.communicate()
    folder = {}
    current = ""
    for line in out.splitlines():
    	if ":" in line:
    		current = line.replace(":","")
    		folder[current] = []
    	if current and ".root" in line:
    		folder[current].append(line)

    filelist = []
    for fold in folder.keys():
    	if not folder[fold]:
    		folder.pop(fold, None)
    	else:
    		for f in folder[fold]:
    			filelist.append(  "/".join(["root://hephyse.oeaw.ac.at", fold, f ])  )

    return filelist



def merge(outdir,files,submit_run):
	test = True

	nevents = 0
	merge_pack = []
	pack = []
	
	if "31Mar2018" in outdir :
		threshold = 1500000
	else :
		threshold = 500000

	for f in files:
		ch = R.TChain("Events")
		ch.Add(f)
		nevents += ch.GetEntries()
		if nevents > threshold:
			nevents = 0
			merge_pack.append(pack)
			pack = []
		pack.append(f)
	
	if nevents < 200000 and len(merge_pack) != 0 : # fix the case where merge_pack is empty
		merge_pack[-1] += pack
	else:
		merge_pack.append(pack)

	for i, pack in enumerate(merge_pack):
	
		if submit_run == 'local' :
			hadd_cmd = "haddnano.py "+outdir+"nano_{0}.root {1}".format(i, " ".join(pack) )
			print i
			print hadd_cmd
			print '-'*100
			os.system(hadd_cmd)
		if submit_run == 'hephybatch' :
			with open("submit_hephybatch.sh".format('batch') ) as FSO:
            			templ = string.Template( FSO.read() )
			runscript = templ.substitute(cmd = outdir+"nano_{0}.root {1}".format(i, " ".join(pack) ))
			rundir_name = "./rundir/"+outdir.split("/")[-2]
			checkDirExist(rundir_name+'/')

			os.chdir(rundir_name)
        		with open("submit_"+str(i)+".sh","w") as FSO:
                		FSO.write( runscript )
            
			os.chmod("submit_"+str(i)+".sh", 0777)
	
       			#print("sbatch submit_"+str(i)+".sh") 
			os.system( "sbatch submit_"+str(i)+".sh" )
			os.chdir("../..")
			if test :
				break

if __name__ == '__main__':
	main()
