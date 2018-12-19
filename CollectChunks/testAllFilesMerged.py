import os
import subprocess as sp
import shlex
import sys
from samples_Janik import samplelist_Janik
from samples_Markus import samplelist_Markus
import string
import numpy as np
import argparse
from chunk import checkDirExist


parser = argparse.ArgumentParser()
parser.add_argument('-n', dest='sample_name', help='Enter the name of the sample to be tested for the correct merging', type=str, metavar = 'SAMPLE_NAME', choices=['SingleElectron','SingleMuon','Tau','VBF','GluGlu','DY','EWK','ST_t','W','ZZ','TTT','Wplus','Wminus','ZH'],default='')
parser.add_argument('-f', dest='which_file_list', help='which file list', type=str, metavar = 'FILELIST', choices=['Janik','Markus'], default='Janik')
parser.add_argument('-d', dest='directory', help='directory', type=str, metavar = 'DIR')

args = parser.parse_args()

sample_name = args.sample_name
if args.which_file_list == 'Janik' :
	samplelist = samplelist_Janik
else :
	samplelist = samplelist_Markus
directory = args.directory



sources = samplelist.keys()
sources.sort()
counter = 0
for source in sources:
	splsource = source.split("/")
	name = "_".join(splsource[1:3]).replace("RunIIFall17MiniAODv2-","").replace("-", "_")
	if not(sample_name in name) :
		continue

	with open("fullMerge.sh".format('batch') ) as FSO:
		templ = string.Template( FSO.read() )
		
	runscript = templ.substitute(sample_name = sample_name, outputdir=directory, samplelist=args.which_file_list)
	
	rundir_name = "./out/outstamp_{0}".format(directory.split('/')[-2])
	checkDirExist(rundir_name+'/')

	os.chdir(rundir_name)
	with open("submit_{0}.sh".format(name),"w") as FSO:
		FSO.write( runscript )
        
	os.chmod("submit_{0}.sh".format(name), 0777)
	print("sbatch submit_{0}.sh".format(name)) 
	os.system( "sbatch submit_{0}.sh".format(name) )
	os.chdir("../../")