import ROOT as R
import os
import subprocess as sp
import shlex
import sys
from samples import samplelist
import string
def main(a):
	s = True
	outputdir = '/afs/hephy.at/data/higgs02/haddJanik/2018-11-23/'
	sources = samplelist.keys()
	sources.sort()
	ListOfIncompleteMerges = []
	b = a + 1
	#print a
	#print b
	counter = 0
	for source in sources[a:b] :
		
		splsource = source.split("/")
		name = "_".join(splsource[1:3]).replace("RunIIFall17MiniAODv2-","").replace("-", "_")
		
		print "Sample name: {0}".format(name)
		#print counter 
		counter += 1
		file_list =  getfiles(samplelist[source][0] )
		n = len(file_list)	
		done = (n == samplelist[source][1] )
        	n = str(n)
        	#print  n + " "*(6 - len(n)) + (not done)*"    " + done*"DONE",  name 
        	
		
		if done :
			try :
				outputFile_list = [outputdir+name+'/'+f for f in  os.listdir(outputdir+name+'/')]
	
				outxt = open(name+'.txt','w')
				nEventsBeforeMerging = GetNumberofEvents(files=file_list)
				print 'No. of events before merging: {0}'.format(nEventsBeforeMerging)			
				nEventsAfterMerging = GetNumberofEvents(files=outputFile_list)
				print 'No. of events after merging: {0}'.format(nEventsAfterMerging)
				if nEventsBeforeMerging - nEventsAfterMerging != 0 :
					ListOfIncompleteMerges.append(name)
			
			except OSError :
				print 'No merged files for this sample'
				s=False
				break

		
		print '-'*100
			
	if done and s:		
		print '_'*50
		print 'Summary:'
		if len(ListOfIncompleteMerges) == 0 :
			print "All Events are merged correctly"
			outxt.write('merging successful /t')
			outxt.write('No. of events before merging: {0} /t '.format(nEventsBeforeMerging))
			outxt.write('No. of events after merging: {0} /t '.format(nEventsAfterMerging))
		else :
			print "In the following samples not all events were merged: "
			outxt.write('merging failed')
			for incomp in ListOfIncompleteMerges :
				print incomp
		print '_'*50
		
def GetNumberofEvents (files) :
	nevents = 0
	for f in files:
                ch = R.TChain("Events")
                ch.Add(f)
                nevents += ch.GetEntries()
	return nevents


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


if __name__ == '__main__':
	main(a = int(sys.argv[1]))

