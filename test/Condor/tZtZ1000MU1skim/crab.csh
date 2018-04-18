#source /uscmst1/prod/grid/gLite_SL6.csh
#cd /uscms_data/d2/skhalil/BPrimeBoost/CMSSW_5_3_3/src/
#setenv SCRAM_ARCH slc5_amd64_gcc462
#cd /uscms/home/skhalil/SHyFT/CMSSW_3_8_6/src
#setenv SCRAM_ARCH slc5_ia32_gcc434
#cmsenv
#source /uscmst1/prod/grid/CRAB_2_7_1/crab.csh
#source /uscmst1/prod/grid/CRAB_2_7_2_p1/crab.csh
#source /uscmst1/prod/grid/CRAB_2_7_5/crab.csh
#source /uscmst1/prod/grid/CRAB_2_8_4/crab.csh
#source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.csh
source /cvmfs/cms.cern.ch/crab3/crab.csh
grid-proxy-init -debug -verify
voms-proxy-init -voms cms --valid 168:00
