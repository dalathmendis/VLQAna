#!/bin/bash                                                                                                                                                                                                                           

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src
eval `scram runtime -sh`

cd ${_CONDOR_SCRATCH_DIR}
echo ${_CONDOR_SCRATCH_DIR}
let "sample=${1}+1"
cp /uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/tZtZ1000MU1skim/SingleMu4/SingleMu4_${sample}.py .
cp /uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/*.txt .
cp /uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/*.root .
cp /uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/inputFiles_cfi.py .
cp /uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/*.csv .
cp /uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt .
cmsRun SingleMu4_${sample}.py
rm btag-eff-subjet.root
rm scalefactors_v4.root
rm PU*
rm Run*
rm *.csv
rm *.txt
rm dataset*
rm os2lana*
rm scale*
rm Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
xrdcp *.root root://cmseos.fnal.gov//store/group/lpcbprime/noreplica/dmendis/0530elnew/Muons2/SingleMu4
#xrdcp *.root root://cmseos.fnal.gov//store/user/dmendis/ntuples2016/reminiaod/batch2/Muons/SingleMu4
rm SingleMu4_${sample}.py
rm *.root
ls
echo "DONE!"
