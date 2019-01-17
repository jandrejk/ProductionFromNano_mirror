#!/usr/bin/env bash

date=yyyy_mm_dd
name=name
forMarkus=/afs/hephy.at/user/j/jandrejkovic/public/forMarkus/datacards
dc_dir=/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/Datacards

out_dir=${forMarkus}/${date}/${name}/Vienna
datacards_dir=/afs/hephy.at/work/j/jandrejkovic/MC_histWriting/CMSSW_9_4_4_fromNano/src/WawTools/NanoAODTools/Datacards/gitlab_datacards/2019_01_17/SM-ML-2017

mkdir -p ${out_dir}

python ${dc_dir}/readShift.py -c et -y 2016 -o ${out_dir} -in1 ${datacards_dir}/Vienna -in2 ${datacards_dir}/KIT |& tee -a ${out_dir}/log.log 

# additional options:
# --all --data_obs --copy