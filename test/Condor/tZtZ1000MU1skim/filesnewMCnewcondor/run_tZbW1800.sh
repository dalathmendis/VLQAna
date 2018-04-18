#!/bin/bash                                                                                                                                                                                                                           
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc530
eval `scramv1 project CMSSW CMSSW_8_0_20`
cd CMSSW_8_0_20/src/
eval `scram runtime -sh`
echo "CMSSW: "$CMSSW_BASE

cd ${_CONDOR_SCRATCH_DIR}
echo ${_CONDOR_SCRATCH_DIR}
let "sample=${1}+1"
xrdcp root://cmseos.fnal.gov//store/user/dmendis/toCondor/CMSSW_8_0_20.tar.gz .
echo "executing ..."
echo "tar -xf CMSSW_8_0_20.tar.gz"
tar -xf CMSSW_8_0_20.tar.gz
rm CMSSW_8_0_20.tar.gz
cd CMSSW_8_0_20/src/Analysis/VLQAna
mkdir test
cd test
cp ../data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt .
xrdcp root://cmseos.fnal.gov//store/user/dmendis/toCondor/muon/tZbW1800/tZbW1800_${sample}.py .
echo "cmsRun"
cmsRun tZbW1800_${sample}.py

rm Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
echo "xrdcp *.root root://cmseos.fnal.gov//store/group/lpcbprime/noreplica/dmendis/0530elnew/Muons2/tZbW1800"

xrdcp *.root root://cmseos.fnal.gov//store/group/lpcbprime/noreplica/dmendis/0530elnew/Muons2/tZbW1800
rm tZbW1800_${sample}.py
rm *.root
rm *.py
cd ../../../../../../..
rm -rf CMSSW_8_0_20
ls
echo "DONE!"
