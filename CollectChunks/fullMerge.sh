#! /bin/sh

export X509_USER_PROXY=/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/CollectChunks/proxy/x509up_u3522
echo "---------------------"
echo "Grid certificate"
voms-proxy-info --all
echo "---------------------"

eval `scramv1 runtime -sh`

python ../testAllFilesMerged.py ${a}  
