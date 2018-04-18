#!/bin/python                                                                                                                        
import sys
import subprocess
import os

#path    = '/eos/uscms/store/group/lpcbprime/noreplica/dmendis/0530elnew/Electrons/'
path1   = '/uscms_data/d3/dmendis/80xB2Gana6/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/tZtZ1000MU1skim/'

Samples = [
       "dy_ht100-200",
       "dy_ht200-400",
       "dy_ht400-600",
       "dy_ht600-800",
       "dy_ht800-1200",
       "dy_ht1200-2500", 
       "dy_ht2500-Inf",
       "dy_pt100-250",
       "dy_pt250-400",
       "dy_pt400-650",
       "dy_pt650-Inf",
      "dy_pt100-250a",
       "dy_pt250-400a",
       "dy_pt400-650a",
       "dy_pt650-Infa",
       "dy_pt100-250b",
       "dy_pt250-400b",

       "tZtZ800",
       'tZtZ1000',
       'tZtZ1200',
       'tZtZ900',
       'tZtZ1100',
       'tZtZ1300',
       'tZtZ1400',
       'tZtZ1500',
       'tZtZ1600',
       'tZtZ1700',
       'tZtZ1800',
       'tZtH800',
       'tZtH1000',
       'tZtH1200',
       'tZtH900',
       'tZtH1100',
       'tZtH1300',
       'tZtH1400',
       'tZtH1500',
       'tZtH1600',
       'tZtH1700',
       'tZtH1800',
       "tZbW800",
       'tZbW1000',
       'tZbW1200',
       'tZbW900',
       'tZbW1100',
       'tZbW1300',
       'tZbW1400',
       'tZbW1500',
       'tZbW1600',
       'tZbW1700',
       'tZbW1800',
       'DoubleMu1',
       'DoubleMu2',
       'DoubleMu3',
       'WW1',
       'WZ1',
       'ZZ1',
       'ZZ',
       'SingleMu1',
       'SingleMu1p1',
       'SingleMu2',
       'SingleMu3',
       'SingleMu4',
       'SingleMu5',
       'SingleMu6',
       #'SingleMu6p1',
       'SingleMu7',
       'SingleMu7p1',
       'SingleMu7p2',
       'ttbar',
       'ttbar1',
       ]  #'WW','WZto2','WZto3','ZZto2','ZZto4',"Del_Prompt"]  
       
os.chdir(path1)

for isample in Samples:
       print "deleting "+str(isample)
       os.chdir(str(isample))
       subprocess.call( 'rm *' , shell=True, executable='/bin/tcsh')
       os.chdir('../')
       print "------------------------------>" + "done " + str(isample)
