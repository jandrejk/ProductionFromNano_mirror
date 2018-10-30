import sys, os
from glob import glob
import shlex
import string
import json
from runUtils import checkProxy, checkTokens, useToken, getSystem, getHeplxPublicFolder
import run as r

# if the output directory does not exist, create it
if not os.path.exists("./MCpuHistograms") :
    os.makedirs("./MCpuHistograms")
if not os.path.exists("./rundir") :
    os.makedirs("./rundir")
rundir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/rundir"

base_dir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/"

if not r.checkProxy(): sys.exit()

file_list = [f for f in glob(base_dir+'/samples/mc/*/*')]


job_counter = 0
for file_name in file_list :
    for sample in open(file_name,'r') :
	os.chdir(base_dir)
	with open("template.sh".format('batch') ) as FSO:
            templ = string.Template( FSO.read() )
	runscript = templ.substitute(sample_name = sample)
	os.chdir('./rundir')
        with open("submit_"+str(job_counter)+".sh","w") as FSO:
            FSO.write( runscript )
            
	os.chmod("submit_"+str(job_counter)+".sh", 0777)

        #print("sbatch submit_"+str(job_counter)+".sh") 
        os.system( "sbatch submit_"+str(job_counter)+".sh" )
        job_counter += 1

