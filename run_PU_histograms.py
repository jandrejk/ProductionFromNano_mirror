import sys, os
import shutil
import threading
from glob import glob
from time import sleep, time
import subprocess as sp
import multiprocessing as mp
import shlex
import argparse
import string
import json
from runUtils import checkProxy, checkTokens, useToken, getSystem, getHeplxPublicFolder


# if the output directory does not exist, create it
if not os.path.exists("./MCpuHistograms") :
    os.makedirs("./MCpuHistograms")
if not os.path.exists("./rundir") :
    os.makedirs("./rundir")
    rundir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/rundir"
    shutil.copytree("proxy", "/".join([rundir,"proxy"]))

base_dir = "/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/"

#if not checkProxy(): sys.exit()
#print checkProxy()

file_list = [f for f in glob(base_dir+'/samples/mc/*/*')]

# use only small sample for the beginning
#file_list = file_list[:1]

test = True
job_counter = 0
for file_name in file_list :
    for sample in open(file_name,'r') :
	os.chdir(base_dir)

	if test :
	    #test = False
            
	    with open("template.sh".format('batch') ) as FSO:
                templ = string.Template( FSO.read() )
	    #print templ

	    runscript = templ.substitute(sample_name = sample)
	    os.chdir('./rundir')
            with open("submit_"+str(job_counter)+".sh","w") as FSO:
                FSO.write( runscript )
            
            #print("sbatch submit_"+str(job_counter)+".sh") 
	    os.chmod("submit_"+str(job_counter)+".sh", 0777)

            os.system( "sbatch submit_"+str(job_counter)+".sh" )

        job_counter += 1


def checkProxy():
    proxy_path = glob("/tmp/x509*_u{0}".format( os.getuid() ) )
    if len(proxy_path) == 1 and os.path.exists( proxy_path[0] ):

        p = sp.Popen( shlex.split("voms-proxy-info --timeleft"), stdout=sp.PIPE, stderr=sp.PIPE )
        (out,err) =  p.communicate()


        if err:
            print err
            return False
        if out:
            if int(out) > 0:
                if not os.path.exists("proxy"):
                    os.mkdir("proxy")
                shutil.copyfile(proxy_path[0], "proxy/x509_proxy")
                return True
            else:
                print "Proxy not valid! Get a new one with 'voms-proxy-init --voms cms'"
                return False

    
    print "No proxy found! Get a new one with 'voms-proxy-init --voms cms'" 
    return False






