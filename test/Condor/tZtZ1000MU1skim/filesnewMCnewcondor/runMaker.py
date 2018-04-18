import sys, math, os, re
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('pyInput', 'muon',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "path to read config"
)
#options.register('output', '/store/user/dmendis/ntuples2016/muons', 
options.register('output', '/store/group/lpcbprime/noreplica/dmendis/0530elnew/Muons2',  
                VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "path to output"
)
options.register('test', '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test', 
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "path to test directory"
)
options.parseArguments()

pyInput = options.pyInput
output = options.output
test = options.test

toMake = [
    'ttbar',
    'ttbar1',
    'dy_pt100-250', 
    'dy_pt250-400', 
    'dy_pt400-650',
    'dy_pt650-Inf',
    'dy_pt100-250a',
    'dy_pt250-400a',
    'dy_pt400-650a',
    'dy_pt650-Infa',
    'dy_pt100-250b',
    'dy_pt250-400b',
    'tZtZ800',
    'tZtZ900',
    'tZtZ1000',
    'tZtZ1100',
    'tZtZ1200',
    'tZtZ1300',
    'tZtZ1400',
    'tZtZ1500',
    'tZtZ1600',
    'tZtZ1700',
    'tZtZ1800',

    'tZtH800',
    'tZtH900',
    'tZtH1000',
    'tZtH1100',
    'tZtH1200',
    'tZtH1300',
    'tZtH1400',
    'tZtH1500',
    'tZtH1600',
    'tZtH1700',
    'tZtH1800',

    'tZbW800',
    'tZbW900',
    'tZbW1000',
    'tZbW1100',
    'tZbW1200',
    'tZbW1300',
    'tZbW1400',
    'tZbW1500',
    'tZbW1600',
    'tZbW1700',
    'tZbW1800',

    'tHtH800',
    'tHtH900',
    'tHtH1000',
    'tHtH1100',
    'tHtH1200',
    'tHtH1300',
    'tHtH1400',
    'tHtH1500',
    'tHtH1600',
    'tHtH1700',
    'tHtH1800',


    'bWbW800',
    'bWbW900',
    'bWbW1000',
    'bWbW1100',
    'bWbW1200',
    'bWbW1300',
    'bWbW1400',
    'bWbW1500',
    'bWbW1600',
    'bWbW1700',
    'bWbW1800',

 'tHbW800',
    'tHbW900',
    'tHbW1000',
    'tHbW1100',
    'tHbW1200',
    'tHbW1300',
    'tHbW1400',
    'tHbW1500',
    'tHbW1600',
    'tHbW1700',
    'tHbW1800',

    'bZbZ800',
    'bZbZ900',
    'bZbZ1000',
    'bZbZ1100',
    'bZbZ1200',
    'bZbZ1300',
    'bZbZ1400',
    'bZbZ1500',
    'bZbZ1600',
    'bZbZ1700',
    'bZbZ1800',


    'bZbH800',
    'bZbH900',
    'bZbH1000',
    'bZbH1100',
    'bZbH1200',
    'bZbH1300',
    'bZbH1400',
    'bZbH1500',
    'bZbH1600',
    'bZbH1700',
    'bZbH1800',

    'bZtW800',
    'bZtW900',
    'bZtW1000',
    'bZtW1100',
    'bZtW1200',
    'bZtW1300',
    'bZtW1400',
    'bZtW1500',
    'bZtW1600',
    'bZtW1700',
    'bZtW1800',


    'bHbH800',
    'bHbH900',
    'bHbH1000',
    'bHbH1100',
    'bHbH1200',
    'bHbH1300',
    'bHbH1400',
    'bHbH1500',
    'bHbH1600',
    'bHbH1700',
    'bHbH1800',


    'tWtW800',
    'tWtW900',
    'tWtW1000',
    'tWtW1100',
    'tWtW1200',
    'tWtW1300',
    'tWtW1400',
    'tWtW1500',
    'tWtW1600',
    'tWtW1700',
    'tWtW1800',


    'bHtW800',
    'bHtW900',
    'bHtW1000',
    'bHtW1100',
    'bHtW1200',
    'bHtW1300',
    'bHtW1400',
    'bHtW1500',
    'bHtW1600',
    'bHtW1700',
    'bHtW1800',




    'WW1',
    'WZ1',
    'ZZ1', 
    'ZZ',

    'tant',
    'tt',
    'twant',
    'twt',
    's',


   # 'SingleMu1',
   # 'SingleMu1p1',
   # 'SingleMu2',
   # 'SingleMu3',
   # 'SingleMu4',
   # 'SingleMu5',
   # 'SingleMu6',
   # 'SingleMu7',
   # 'SingleMu7p1',
   # 'SingleMu7p2',
    ]

for n in toMake:
    inputFile = open('runDummy.py')
    outputFile = open('run_'+str(n)+'.sh', 'w')
    for line in inputFile:
        line = line.replace('SAMPLE', n)
        line = line.replace('PYTHON', pyInput)
        line = line.replace('OUTPUT', output)
        line = line.replace('TEST', test)
        outputFile.writelines(line)
    inputFile.close()
    outputFile.close()
