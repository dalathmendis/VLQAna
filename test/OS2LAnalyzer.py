#! /usr/bin/env python
from __future__ import print_function
import os
import glob
import math
import copy
import subprocess
from ROOT import gROOT,std,ROOT,TFile,TH1F,TH2F,TStopwatch,TLorentzVector,TMath
#import ROOT
gROOT.Macro("~/rootlogon.C")
import sys
from DataFormats.FWLite import Events, Handle

from optparse import OptionParser
parser = OptionParser()

parser.add_option('--files', metavar='F', type='string', action='store',
                  default = 'OS2LAnaEvts_mu.root',
                  dest='files',
                  help='Input files')

parser.add_option("--runMuons", action='store_true',
                  default=True,
                  dest="runMuons",
                  help="Electrons(0), Muons(1)")

parser.add_option('--runMC', metavar='I', type='int', action='store',
                  default=1,
                  dest='runMC',
                  help='0 (all bkg MC), 1 (TPrime), 2 (BPrime), 3 (Data)')

parser.add_option('--pickChannel', metavar='S', type='string', action='store',
                  default='ZtZt',
                  dest='pickChannel',
                  help='VLQ channel to be filtered: ZtZt, WbZt, HtZt')

# Parse and get arguments
(options, args) = parser.parse_args()
print('options', options)

runbkgMC = 0; runTPrime = 0; runBPrime = 0; runData = 0;

if options.runMC == 0:
    runbkgMC = 1
elif options.runMC == 1:
    runTPrime = 1
elif options.runMC == 2:
    runBPrime = 1
elif options.runMC == 3:
    runData = 1

print ('runbkgMC: {0:2d}, runTPrime : {1:2d}, runBPrime : {2:2d}, runData : {3:2d}'.format (runbkgMC, runTPrime, runBPrime, runData))

print ('Getting files from this dir: ' + options.files)

if options.files:
    files = glob.glob( options.files )
else: files =[]

# Get the FWLite "Events"
events = Events (files)

#get the handles and labels
if runTPrime:
    TTtoHtHt_H = Handle("unsigned"); TTtoHtHt_L = ("GenInfo", "TTtoHtHt")
    TTtoHtZt_H = Handle("unsigned"); TTtoHtZt_L = ("GenInfo", "TTtoHtZt")
    TTtoWbHt_H = Handle("unsigned"); TTtoWbHt_L = ("GenInfo", "TTtoWbHt")
    TTtoWbWb_H = Handle("unsigned"); TTtoWbWb_L = ("GenInfo", "TTtoWbWb")
    TTtoWbZt_H = Handle("unsigned"); TTtoWbZt_L = ("GenInfo", "TTtoWbZt")
    TTtoZtZt_H = Handle("unsigned"); TTtoZtZt_L = ("GenInfo", "TTtoZtZt")
    
# Keep some timing information
nEventsAnalyzed = 0
timer = TStopwatch()
timer.Start()

# loop over events
i = 0
for event in events:
    i = i + 1
    if i % 1000 == 0 :
        print("EVENT ", i)
    nEventsAnalyzed = nEventsAnalyzed + 1

    # Get the objects
    if runTPrime:
        event.getByLabel(TTtoHtHt_L, TTtoHtHt_H); TTtoHtHt = TTtoHtHt_H.product()[0]
        event.getByLabel(TTtoHtZt_L, TTtoHtZt_H); TTtoHtZt = TTtoHtZt_H.product()[0]
        event.getByLabel(TTtoWbHt_L, TTtoWbHt_H); TTtoWbHt = TTtoWbHt_H.product()[0]
        event.getByLabel(TTtoWbWb_L, TTtoWbWb_H); TTtoWbWb = TTtoWbWb_H.product()[0]
        event.getByLabel(TTtoWbZt_L, TTtoWbZt_H); TTtoWbZt = TTtoWbZt_H.product()[0]
        event.getByLabel(TTtoZtZt_L, TTtoZtZt_H); TTtoZtZt = TTtoZtZt_H.product()[0]
        if (TTtoHtHt or TTtoWbHt or TTtoWbWb): continue
        #print ("TPrime event with one Z  = {0:1d}".format(TTtoZtZt+TTtoHtZt+TTtoWbZt))
        if (options.pickChannel == 'ZtZt' and (TTtoHtZt+TTtoWbZt == 1)): continue
        print ("ZtZt: {0:2d}, HtZt : {1:2d}, WbZt : {2:2d}".format (TTtoZtZt, TTtoHtZt, TTtoWbZt))

# Done processing the events!
# Stop our timer
timer.Stop()

# Print out our timing information
rtime = timer.RealTime(); # Real time (or "wall time")
ctime = timer.CpuTime(); # CPU time
print("Analyzed events: {0:6d}".format(nEventsAnalyzed))
print("RealTime={0:6.2f} seconds, CpuTime={1:6.2f} seconds".format(rtime,ctime))
print("{0:4.2f} events / RealTime second .".format( nEventsAnalyzed/rtime))
print("{0:4.2f} events / CpuTime second .".format( nEventsAnalyzed/ctime))
subprocess.call( ["ps aux | grep skhalil | cat > memory.txt", ""], shell=True )
