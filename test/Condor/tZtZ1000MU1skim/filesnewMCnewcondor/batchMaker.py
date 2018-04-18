import sys, math, os, re
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('Path', '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/tZtZ1000MU1skim/Condorout',
                VarParsing.multiplicity.singleton,
                VarParsing.varType.string,
                "path to store"
)
options.parseArguments()

PATH = options.Path
print PATH

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
    inputFile = open('batchDummy.py')
    outputFile = open('batch_'+str(n)+'.jdl', 'w')
    if n == 'ttbar':
        QUEUE = '984'
        EXE = 'run_'+n+'.sh'
    
    if n == 'ttbar1':
        QUEUE = '1192'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt100-250':
        QUEUE = '34'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt250-400':
        QUEUE = '11'
        EXE = 'run_'+n+'.sh'
    
    if n == 'dy_pt400-650':
        QUEUE = '18'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt650-Inf':
        QUEUE = '15'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt100-250a':
        QUEUE = '34'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt250-400a':
        QUEUE = '9'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt400-650a':
        QUEUE = '7'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt650-Infa':
        QUEUE = '7'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt100-250b':
        QUEUE = '1056'
        EXE = 'run_'+n+'.sh'

    if n == 'dy_pt250-400b':
        QUEUE = '379'
        EXE = 'run_'+n+'.sh'


    if n == 'tZtZ800':
        QUEUE = '24'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH800':
        QUEUE = '24'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW800':
        QUEUE = '24'
        EXE = 'run_'+n+'.sh'

    if n == 'tHtH800':
        QUEUE = '24'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW800' or n == 'tHbW800':
        QUEUE = '24'
        EXE = 'run_'+n+'.sh'


    if n == 'tZtZ900':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH900':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW900':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tHtH900':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW900' or n == 'tHbW900':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtZ1000':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1000':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1000':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tHtH1000':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1000' or n == 'tHbW1000':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'


    if n == 'tZtZ1100':
        QUEUE = '19'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1100':
        QUEUE = '19'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1100':
        QUEUE = '19'
        EXE = 'run_'+n+'.sh'

    if n == 'tHtH1100':
        QUEUE = '19'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1100' or n == 'tHbW1100':
        QUEUE = '19'
        EXE = 'run_'+n+'.sh'


    if n == 'tZtZ1200':
        QUEUE = '36'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1200':
        QUEUE = '36'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1200':
        QUEUE = '36'
        EXE = 'run_'+n+'.sh'


    if n == 'tHtH1200':
        QUEUE = '36'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1200' or n == 'tHbW1200':
        QUEUE = '36'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtZ1300':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1300':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1300':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tHtH1300':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1300' or n == 'tHbW1300':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'


    if n == 'tZtZ1400':
        QUEUE = '33'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1400':
        QUEUE = '33'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1400':
        QUEUE = '33'
        EXE = 'run_'+n+'.sh'


    if n == 'tHtH1400':
        QUEUE = '33'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1400' or n == 'tHbW1400':
        QUEUE = '33'
        EXE = 'run_'+n+'.sh'
    if n == 'tZtZ1500':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1500':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1500':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tHtH1500':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1500' or n == 'tHbW1500':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtZ1600':
        QUEUE = '35'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1600':
        QUEUE = '35'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1600':
        QUEUE = '35'
        EXE = 'run_'+n+'.sh'

    if n == 'tHtH1600':
        QUEUE = '35'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1600' or n == 'tHbW1600':
        QUEUE = '35'
        EXE = 'run_'+n+'.sh'


    if n == 'tZtZ1700':
        QUEUE = '29'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1700':
        QUEUE = '29'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1700':
        QUEUE = '29'
        EXE = 'run_'+n+'.sh'



    if n == 'tHtH1700':
        QUEUE = '29'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1700' or n == 'tHbW1700':
        QUEUE = '29'
        EXE = 'run_'+n+'.sh'
    



    if n == 'tZtZ1800':
        QUEUE = '22'
        EXE = 'run_'+n+'.sh'

    if n == 'tZtH1800':
        QUEUE = '22'
        EXE = 'run_'+n+'.sh'

    if n == 'tZbW1800':
        QUEUE = '22'
        EXE = 'run_'+n+'.sh'
        
    if n == 'tHtH1800':
        QUEUE = '22'
        EXE = 'run_'+n+'.sh'

    if n == 'bWbW1800' or n == 'tHbW1800':
        QUEUE = '22'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ800' or n == 'bZbH800' or n == 'bZtW800'or n == 'bHbH800'or n == 'tWtW800'or n == 'bHtW800':
        QUEUE = '26'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ900' or n == 'bZbH900' or n == 'bZtW900'or n == 'bHbH900'or n == 'tWtW900'or n == 'bHtW900':
        QUEUE = '17'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ1000' or n == 'bZbH1000' or n == 'bZtW1000'or n == 'bHbH1000'or n == 'tWtW1000'or n == 'bHtW1000':
        QUEUE = '19'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ1100' or n == 'bZbH1100' or n == 'bZtW1100'or n == 'bHbH1100'or n == 'tWtW1100'or n == 'bHtW1100':
        QUEUE = '28'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ1200' or n == 'bZbH1200' or n == 'bZtW1200'or n == 'bHbH1200'or n == 'tWtW1200'or n == 'bHtW1200':
        QUEUE = '24'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ1300' or n == 'bZbH1300' or n == 'bZtW1300'or n == 'bHbH1300'or n == 'tWtW1300'or n == 'bHtW1300':
        QUEUE = '22'
        EXE = 'run_'+n+'.sh'
    
    if n == 'bZbZ1400' or n == 'bZbH1400' or n == 'bZtW1400'or n == 'bHbH1400'or n == 'tWtW1400'or n == 'bHtW1400':
        QUEUE = '27'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ1500' or n == 'bZbH1500' or n == 'bZtW1500'or n == 'bHbH1500'or n == 'tWtW1500'or n == 'bHtW1500':
        QUEUE = '18'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ1600' or n == 'bZbH1600' or n == 'bZtW1600'or n == 'bHbH1600'or n == 'tWtW1600'or n == 'bHtW1600':
        QUEUE = '15'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ1700' or n == 'bZbH1700' or n == 'bZtW1700'or n == 'bHbH1700'or n == 'tWtW1700'or n == 'bHtW1700':
        QUEUE = '23'
        EXE = 'run_'+n+'.sh'

    if n == 'bZbZ1800' or n == 'bZbH1800' or n == 'bZtW1800'or n == 'bHbH1800'or n == 'tWtW1800'or n == 'bHtW1800':
        QUEUE = '20'
        EXE = 'run_'+n+'.sh'

    if n == 'WW1':
        QUEUE = '114'
        EXE = 'run_'+n+'.sh'
    
    if n == 'WZ1':
        QUEUE = '128'
        EXE = 'run_'+n+'.sh'

    if n == 'ZZ1':
        QUEUE = '66'
        EXE = 'run_'+n+'.sh'

    if n == 'ZZ':
        QUEUE = '11'
        EXE = 'run_'+n+'.sh'

    if n == 'tant':
        QUEUE = '168'
        EXE = 'run_'+n+'.sh'

    if n == 'tt':
        QUEUE = '122'
        EXE = 'run_'+n+'.sh'
   
    if n == 'twant':
        QUEUE = '17'
        EXE = 'run_'+n+'.sh'

    if n == 'twt':
        QUEUE = '14'
        EXE = 'run_'+n+'.sh'

    if n == 's':
        QUEUE = '21'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu1':
        QUEUE = '57'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu1p1':
        QUEUE = '1277'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu2':
        QUEUE = '713'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu3':
        QUEUE = '1060'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu4':
        QUEUE = '886'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu5':
        QUEUE = '768'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu6':
        QUEUE = '1430'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu7':
        QUEUE = '14'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu7p1':
        QUEUE = '1928'
        EXE = 'run_'+n+'.sh'

    if n == 'SingleMu7p2':
        QUEUE = '52'
        EXE = 'run_'+n+'.sh'

    for line in inputFile:
        line = line.replace('queue', QUEUE)
        line = line.replace('path', PATH+'/'+n)
        line = line.replace('exe', EXE)
        outputFile.writelines(line)
    inputFile.close()
    outputFile.close()
