import sys, os
from glob import glob
import string
import argparse


def main() :
    print 'start searching for missing files...'
    parser = argparse.ArgumentParser()
    parser.add_argument('-r',  dest='resub', metavar='resub',  help='Used to resubmit production of missing histograms', default = False)

    args = parser.parse_args()
    resubmit = args.resub
    base_dir = "{CMSSW_BASE}/src/WawTools/NanoAODTools/".format(**os.environ)

    file_list = [f for f in glob(base_dir+'/samples/mc/*/*')]
    root_files = []

    MC_hist_list = [mc[112:] for mc in glob(base_dir+'/MCpuHistograms/*')]
    counter = 0
    loopcounter = 0
    for file_name in file_list :
        for sample in open(file_name,'r') :
            sample = sample.replace('/','_')[31:].rstrip()
            root_files.append(sample)
	    os.chdir(base_dir)
	    if not any(sample in x for x in MC_hist_list) :
	        print 'missing file \n', sample, '\n'
	        print 'look at cat ./rundir/submit_'+str(loopcounter)+'.sh'
	        if resubmit == 'yes':
		    print 'resubmitting file'
		    resubmit_func(counter=loopcounter)
		counter += 1
		print '-'*30
 	    loopcounter += 1    

    if counter == 0 :
        print 'all Pileup histograms are there!' 
    else :
        print counter, 'files are missing' 


def resubmit_func (counter) :
    submitbatch_cmd = ("sbatch submit_"+str(counter)+".sh") 
    #print submitbatch_cmd
    os.chdir("./rundir")
    os.system(submitbatch_cmd)
    os.chdir("../")	

if __name__ == '__main__':
    main()
