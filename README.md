# ProductionFromNano

Tools for running WAW ntuple production from NanoAOD by CMS

Port of [tools by A. Kalinowski, M. Bluj et al. based on KLUB/LLR trees](https://github.com/akalinow/Production.git) to CMS NanoAOD

---

* NanoEventsSkeleton.{h,C}: interfacte to NanoAOD ntuples produced from NanoAOD with MakeClass tool, it needs be updated after each modification of NanoAOD format.
* HTauTauTreeFromNanoBase.{h,C}: base class to translate to WAW format
* HMuTauhTreeFromNano.{h,C}: specialization for the mu+tau channel
* HTauhTauhTreeFromNano.{h,C}: specialization for the di-tau channel
* HTTEvent.{h,cxx}: definition of WAW analysis classes
* PropertyEnum.h, TriggerEnum.h: definition of enums, (re)generated by the tool
* AnalysisEnums.h, SelectionBitsEnum.h: definition of enums
* convertNano.py: script to run conversion
* Missing: production tools, need be taken modified from old repo

---

Installation recipe for CMSSW_9_2_4
```
scram project -n CMSSW_9_4_2_fromNano CMSSW CMSSW_9_4_2
cd CMSSW_9_4_2_fromNano/src/
cmsenv
# NanoAOD and tools 
git cms-addpkg PhysicsTools/NanoAOD #not mandatory, but it initializes git for CMSSW which is already useful
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools # not used for now, but can be in future, e.g. JES?
# SVFit
git clone https://github.com/svfit/ClassicSVfit.git TauAnalysis/ClassicSVfit 
git clone https://github.com/svfit/SVfitTF.git TauAnalysis/SVfitTF
# MET recoil corrections
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections
# WAW production tools from NanoAOD
git clone https://github.com/mbluj/ProductionFromNano.git WawTools/NanoAODTools
# compile
scram b -j 4
```


---
Release notes:
* 16.01.2018, M.Bluj, initial version for 80X (2016) inputs with CMSSW_9_4_2
