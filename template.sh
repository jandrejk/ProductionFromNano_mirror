#! /bin/sh

#export X509_USER_PROXY=./proxy/x509_proxy
export X509_USER_PROXY=/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/proxy/x509up_u3522
echo "---------------------"
echo "Grid certificate"
voms-proxy-info --all
echo "---------------------"

eval `scramv1 runtime -sh`

python ../produce_MCpuHistograms.py ${sample_name} 
