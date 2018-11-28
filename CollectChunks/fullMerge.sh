#! /bin/sh
eval 'scramv1 runtime -sh'
export X509_USER_PROXY=$$CMSSW_BASE/src/WawTools/NanoAODTools/CollectChunks/proxy/x509up_u3522
echo "---------------------"
echo "Grid certificate"
voms-proxy-info --all
echo "---------------------"

eval `scramv1 runtime -sh`

python ../../fullMerge.py ${sample_name} ${outputdir} ${samplelist}
