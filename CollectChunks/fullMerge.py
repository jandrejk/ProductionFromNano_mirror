import ROOT as R
import os
import subprocess as sp
import shlex
import sys
import string
from samples_Janik import samplelist_Janik
from samples_Markus import samplelist_Markus

def main(sample_name,outputdir,whatsamplelist):
	if 'Janik' in whatsamplelist :
		samplelist = samplelist_Janik
	else :
		samplelist = samplelist_Markus
	s = True
	
	sources = samplelist.keys()
	sources.sort()
	ListOfIncompleteMerges = []
	
	for source in sources:
		
		splsource = source.split("/")
		name = "_".join(splsource[1:3]).replace("RunIIFall17MiniAODv2-","").replace("-", "_")
		
		
		if not(sample_name in name) :
			continue	
		
		print "Sample name: {0}".format(name)
		file_list =  getfiles(samplelist[source][0] )
		n = len(file_list)	
		done = (n == samplelist[source][1] )

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
			outxt.write('merging successful	')
			outxt.write('No. of events before merging: {0}	'.format(nEventsBeforeMerging))
			outxt.write('No. of events after merging: {0}	'.format(nEventsAfterMerging))
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
	main(sample_name = sys.argv[1], outputdir=sys.argv[2],whatsamplelist=sys.argv[3] )

