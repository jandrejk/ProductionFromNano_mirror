import ROOT as R
import os
import subprocess as sp
import shlex
import sys
from samples import samplelist
import string
import numpy as np
def main():
	for index in np.arange(0,48,1):
		#print index
		with open("fullMerge.sh".format('batch') ) as FSO:
            		templ = string.Template( FSO.read() )
		runscript = templ.substitute(a = index)
		
		rundir_name = "./out"
		checkDirExist(rundir_name+'/')

		os.chdir(rundir_name)
        	with open("submit_"+str(index)+".sh","w") as FSO:
                	FSO.write( runscript )
            
		os.chmod("submit_"+str(index)+".sh", 0777)

        	#print("sbatch submit_"+str(index)+".sh") 
		os.system( "sbatch submit_"+str(index)+".sh" )
		os.chdir("../")

def checkDirExist (directory) :
	if not os.path.exists(directory) :
    		os.makedirs(directory)

		

if __name__ == '__main__':
	main()
