#! /bin/sh

#export X509_USER_PROXY=./proxy/x509_proxy
eval `scramv1 runtime -sh`
export X509_USER_PROXY=$$CMSSW_BASE/src/WawTools/NanoAODTools/proxy/x509_proxy
echo "---------------------"
echo "Grid certificate"
voms-proxy-info --all
echo "---------------------"

eval `scramv1 runtime -sh`

python ../produce_MCpuHistograms.py ${sample_name} 
