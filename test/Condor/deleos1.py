#!/bin/python                                                                                                                        
import sys
import subprocess
import os

path    = '/store/group/lpcbprime/noreplica/dmendis/0530elnew/Muons2/'

Samples = [
      # "dy_ht100-200",
      # "dy_ht200-400",
      # "dy_ht400-600",
      # "dy_ht600-800",
      # "dy_ht800-1200",
      # "dy_ht1200-2500", 
      # "dy_ht2500-Inf",
       "dy_Inc",
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
       'tZtH800',
       'tZtH1000',
       'tZtH1200',
       "tZbW800",
       'tZbW1000',
       'tZbW1200',
       #'DoubleMu1',
       #'DoubleMu2',
       #'DoubleMu3',
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
       "tHtH800",
       "tHtH1000",
       "tHtH1200",
       "tHtH900",
       "tHtH1100",
       "tHtH1300",
       "tHtH1400",
       "tHtH1500",
       "tHtH1600",
       "tHtH1700",
       "tHtH1800",

       "bWbW800",
       "bWbW1000",
       "bWbW1200",
       "bWbW900",
       "bWbW1100",
       "bWbW1300",
       "bWbW1400",
       "bWbW1500",
       "bWbW1600",
       "bWbW1700",
       "bWbW1800",

       "tHbW800",
       "tHbW1000",
       "tHbW1200",
       "tHbW900",
       "tHbW1100",
       "tHbW1300",
       "tHbW1400",
       "tHbW1500",
       "tHbW1600",
       "tHbW1700",
       "tHbW1800",




       'ttbar',
       'ttbar1',
       'WW1',
       'WZ1',
       'ZZ1',
       'ZZ',
       'tant',
       'tt',
       'twant',
       'twt',
       's',
       'SingleMu1',
       'SingleMu1p1',                                                                                                                                                
      'SingleMu2',                                                                                                                                                    
       'SingleMu4',   
       'SingleMu5',
       'SingleMu6',
       'SingleMu7',
       'SingleMu7p1',
       'SingleMu7p2',
       'SingleMu3',


       ]  #'WW','WZto2','WZto3','ZZto2','ZZto4',"Del_Prompt"]  


for isample in Samples:
       print "deleting "+str(isample)
       command = 'eosrm -r '
       subprocess.call( command+path+str(isample) , shell=True, executable='/bin/tcsh')
       print "------------------------------>" + "done " + str(isample)

       command1 = 'eosmkdir '
       subprocess.call( command1+path+str(isample) , shell=True, executable='/bin/tcsh')
       
       print "------------------------------>" + "new directory " + str(isample) + " crated"
