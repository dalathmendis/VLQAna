#!/usr/bin/env python
import sys
import os
import subprocess
from array import array
from ROOT import TH1D,TH2D,TFile,TMath,TCanvas,THStack,TLegend,TPave,TLine,TLatex
from ROOT import gROOT,gStyle,gPad,gStyle
from ROOT import Double,kBlue,kRed,kOrange,kMagenta,kYellow,kCyan,kGreen,kGray,kBlack,kTRUE

gROOT.Macro("~/rootlogon.C")
gStyle.SetOptStat(0)

# ===============
# options
# ===============
from optparse import OptionParser
parser = OptionParser()
parser.add_option('--Lumi', metavar='D', type='float', action='store',
                  default= 35900.,# 36800.,37900.,35737., 35710., 27165.,   12900., 2263.,    2200., 538.134,
                  dest='Lumi',
                  help='Data Luminosity in pb-1')

parser.add_option('--globalSF', metavar='SF', type='float',
                  default = 1.0,
                  dest='globalSF',
                  help='Global trigger SF (%default default)')

parser.add_option('--plotDir', metavar='P', type='string', action='store',
                  #default='newlogplots/mu/nobug/dr2',
                  #default='data/MU/catsig2',
                  #default='eltriggers/Ele32_sig2',
                  default='dyinc/EL/test3pre',
                  dest='plotDir',
                  help='output directory of plots')

parser.add_option('--skimType', metavar='S', type='string', action='store',
                  default='CR_Zelel',
                  dest='skimType',
                  help='Skim type: CR_Zelel, CR_Zmumu, SR_Zelel, SR_Zmumu')

parser.add_option('--processDir', metavar='pD', type='string', action='store',
                  default='ana/pre',
                  dest='processDir',
                  help='directory to read histograms from: ana/pre, ana/cnt, ana/sig')

parser.add_option('--var', metavar='T', type='string', action='store',
                  default='ht_zsel',#cutflow, st
                  dest='var',
                  help='variable to plot')

parser.add_option('--sys', metavar='T', type='string', action='store',
                  default='nominal',
                  dest='sys',
                  help='nominal, BTagSFup, BTagSFdown, ScaleUp, ScaleDown, MatchUp, MatchDown')

parser.add_option('--verbose',action='store_true',
                  default=True,
                  dest='verbose',
                  help='verbose switch')

parser.add_option('--rebin', metavar='T', type='int', action='store',
                  default='1',
                  dest='rebin',
                  help='rebin the histograms')



(options,args) = parser.parse_args()
# ==========end: options =============
var = options.var
lumi = options.Lumi
gSF  = options.globalSF
rebinS = options.rebin
pDir = options.processDir
plotDir = options.plotDir
skimType = options.skimType
verbose = options.verbose

print "var is = ", var
if 'elel' in skimType: title = 'e^{#pm}e^{#mp}+jets'
elif 'mumu' in skimType: title = '#mu^{#pm}#mu^{#mp}+jets'
else: title = ''
# ==========add the input ============

execfile("input.py")

# === prepare the input labels and legends ===========
dataLabel     = 'Data_'
topLabel      = 'Top_'
dyLabel       = 'DY_'
dibosonLabel  = 'Diboson_'
wjLabel       = 'WJets_'
sTLabel       = 'ST_'
vvLabel       = 'VV_'

dataLeg       = 'Data'
topLeg        = 't#bar{t}'
dyLeg         = 'Drell-Yan'
dibosonLeg    = 'Diboson'
wjLeg         = 'W+Jets'
sTLeg         = 'Single top'
vvLeg         = 'Diboson'
# === create structure ============
data     = [
            #[f_Data_Oct2015, 1., 1., 1.], # this corresponds to .53 fb-1
             [f_Data1, 1., 1., 1.],  #addign this should give 2.2 fb-1  
             #[f_Data1p1, 1., 1., 1.], 
             [f_Data2, 1., 1., 1.],  #addign this should give 2.2 fb-1             
             [f_Data3, 1., 1., 1.],  #addign this should give 2.2 fb-1                
             [f_Data4, 1., 1., 1.],  #addign this should give 2.2 fb-1                                                                                                                            
             [f_Data5, 1., 1., 1.],  #addign this should give 2.2 fb-1                                                                                                                            
             [f_Data6, 1., 1., 1.],  #addign this should give 2.2 fb-1 
             [f_Data7, 1., 1., 1.],  #addign this should give 35.733 fb-1  
             [f_Data7p1, 1., 1., 1.],
             #[f_Data7p2, 1., 1., 1.],

             ]

top      = [[f_ttbar,         Top_xs,            Top_num,            lumi]]

dy       = [[f_DYmcnlo,       DY_xs,             DYmcnlo_num,        lumi]]

#dy       = [[f_DYmad,       DY_xs,             DYmad_num,        lumi]] 

#dy       = [
 #           [f_DY100to200,    DY100to200_xs,      DY100to200_num,    lumi],
  #          [f_DY200to400,    DY200to400_xs,      DY200to400_num,    lumi],
   #         [f_DY400to600,    DY400to600_xs,      DY400to600_num,    lumi],
    #        [f_DY600to800,    DY600to800_xs,      DY600to800_num,    lumi], 
     #       [f_DY800to1200,   DY800to1200_xs,     DY800to1200_num,    lumi],
      #      [f_DY1200to2500,  DY1200to2500_xs,    DY1200to2500_num,   lumi],
       #     [f_DY2500toInf,   DY2500toInf_xs,     DY2500toInf_num,    lumi],    
        #    ]

#dy       = [
 #           [f_DY100to250,    DY100to250_xs,      DY100to250_num,    lumi],
 #           [f_DY250to400,    DY250to400_xs,      DY250to400_num,    lumi],
 #           [f_DY400to650,    DY400to650_xs,      DY400to650_num,    lumi],
 #           [f_DY650toInf,    DY650toInf_xs,      DY650toInf_num,    lumi],
           # [f_DY800to1200,   DY800to1200_xs,     DY800to1200_num,    lumi],
           # [f_DY1200to2500,  DY1200to2500_xs,    DY1200to2500_num,   lumi],
           # [f_DY2500toInf,   DY2500toInf_xs,     DY2500toInf_num,    lumi],
 #           ]

diboson  =[                                                                                                                                             
             [      f_WW,             WW_xs,              WW_num,    lumi],                                                                            
             [      f_WZ,             WZ_xs,              WZ_num,    lumi],
             [      f_ZZ,             ZZ_xs,              ZZ_num,    lumi],

             ]

#diboson  =[
 #            [      f_WW,             WW_xs,              WW_num,    lumi],
  #           [      f_WZto2,           WZto2_xs,          WZto2_num,    lumi],
   #          [      f_WZto3,           WZto3_xs,          WZto3_num,    lumi],
    #         [      f_ZZto2,           ZZto2_xs,          ZZto2_num,    lumi],
     #        [      f_ZZto4,           ZZto4_xs,          ZZto4_num,    lumi],   
      #       ] 
#wjets    = [
#            [f_WJ100to200,    WJ100to200_xs,      WJ100to200_num,    lumi],
#            [f_WJ200to400,    WJ200to400_xs,      WJ200to400_num,    lumi],
#            [f_WJ400to600,    WJ400to600_xs,      WJ400to600_num,    lumi],
#            [f_WJ600to800,    WJ600to800_xs,      WJ600to800_num,    lumi],     
#            [f_WJ800to1200,   WJ800to1200_xs,     WJ800to1200_num,   lumi],
#            [f_WJ1200to2500,  WJ1200to2500_xs,    WJ1200to2500_num,  lumi],
#            [f_WJ2500toInf,   WJ2500toInf_xs,     WJ2500toInf_num,   lumi], 
#            ]

#st       = [
 #           [f_tant,     tant_xs,       tant_num,     lumi],
  #          [f_tt,     tt_xs,       tt_num,     lumi],    
   #         [f_twant,     twant_xs,       twant_num,     lumi],
    #        [f_twt,     twt_xs,       twt_num,     lumi],        
     #       [f_s,     s_xs,       s_num,     lumi],    
      #      ]

#vv       = [
#            [f_ZZTo2L2Nu,     ZZTo2L2Nu_xs,       ZZTo2L2Nu_num,    lumi],
#            [f_WZTo2L2Q,      WZTo2L2Q_xs,        WZTo2L2Q_num,     lumi],
#            [f_WWTo2L2Nu,     WWTo2L2Nu_xs,       WWTo2L2Nu_num,    lumi],
#           ]
#tZtZ_700 = [[f_BpBp_tZtZ_700, TpTp700_xs,         TpTp700_num,       lumi]]
#tZbW_700 = [[f_BpTp_tZbW_700, TpTp700_xs,         TpTp700_num,       lumi]]
#tZtH_700 = [[f_BpTp_tZtH_700, TpTp700_xs,         TpTp700_num,       lumi]]
#tZtZ_800 = [[f_TpTp_tZtZ_800, TpTp800_xs,         TpTp800tZtZ_num,       lumi]]
#tZbW_800 = [[f_TpTp_tZbW_800, TpTp800_xs,         TpTp800tZbW_num,       lumi]]
#tZtH_800 = [[f_TpTp_tZtH_800, TpTp800_xs,         TpTp800tZtH_num,       lumi]]
tZtZ_1000 = [[f_TpTp_tZtZ_1000, TpTp1000_xs,         TpTp1000tZtZ_num,       lumi]]
tZbW_1000 = [[f_TpTp_tZbW_1000, TpTp1000_xs,         TpTp1000tZbW_num,       lumi]]
tZtH_1000 = [[f_TpTp_tZtH_1000, TpTp1000_xs,         TpTp1000tZtH_num,       lumi]]

#tZ_1000 = [[f_TpTp_tZ_1000, Tp1000_xs,         TpTp1000tZ_num,       lumi]]

#tZtZ_1200 = [[f_TpTp_tZtZ_1200, TpTp1200_xs,         TpTp1200tZtZ_num,       lumi]]
#tZbW_1200 = [[f_TpTp_tZbW_1200, TpTp1200_xs,         TpTp1200tZbW_num,       lumi]]
#tZtH_1200 = [[f_TpTp_tZtH_1200, TpTp1200_xs,         TpTp1200tZtH_num,       lumi]]


#tZtZ_1400 = [[f_TpTp_tZtZ_1400, TpTp1400_xs,         TpTp1400tZtZ_num,       lumi]]
#tZbW_1400 = [[f_TpTp_tZbW_1400, TpTp1400_xs,         TpTp1400tZbW_num,       lumi]]
#tZtH_1400 = [[f_TpTp_tZtH_1400, TpTp1400_xs,         TpTp1400tZtH_num,       lumi]]


#bZbZ_700 = [[f_BpBp_bZbZ_700, BpBp700_xs,         BpBp700_num,       lumi]]
#bZtW_700 = [[f_BpBp_bZtW_700, BpBp700_xs,         BpBp700_num,       lumi]]
#bZbH_700 = [[f_BpBp_bZbH_700, BpBp700_xs,         BpBp700_num,       lumi]]
#bZbZ_800 = [[f_BpBp_bZbZ_800, BpBp800_xs,         BpBp800_num,       lumi]]
#bZtW_800 = [[f_BpBp_bZtW_800, BpBp800_xs,         BpBp800_num,       lumi]]
#bZbH_800 = [[f_BpBp_bZbH_800, BpBp800_xs,         BpBp800_num,       lumi]]
#bZbZ_1000 = [[f_BpBp_bZbZ_1000, BpBp1000_xs,         BpBp1000_num,       lumi]]
#bZtW_1000 = [[f_BpBp_bZtW_1000, BpBp1000_xs,         BpBp1000_num,       lumi]]
#bZbH_1000 = [[f_BpBp_bZbH_1000, BpBp1000_xs,         BpBp1000_num,       lumi]]


#bZbZ_1200 = [[f_BpBp_bZbZ_1200, BpBp1200_xs,         BpBp1200_num,       lumi]]
#bZtW_1200 = [[f_BpBp_bZtW_1200, BpBp1200_xs,         BpBp1200_num,       lumi]]
#bZbH_1200 = [[f_BpBp_bZbH_1200, BpBp1200_xs,         BpBp1200_num,       lumi]]


if len(data)>0:
    h_data     = getHisto(dataLabel,       dataLeg,        pDir, var,  data,     kBlack,     verbose)



h_top      = getHisto(topLabel,        topLeg,         pDir, var,  top,      8,          verbose)
h_dy       = getHisto(dyLabel,         dyLeg,          pDir, var,  dy,       90,         verbose)
h_diboson  = getHisto(dibosonLabel,  dibosonLeg,       pDir, var,  diboson,   kBlue,         verbose)
#h_wjets    = getHisto(wjLabel,         wjLeg,          pDir, var,  wjets,    kBlue,      verbose)
#h_st       = getHisto(sTLabel,         sTLeg,          pDir, var,  st,       kCyan,      verbose)
#h_vv       = getHisto(vvLabel,         vvLeg,          pDir, var,  vv,       kRed,       verbose)

#h_tZtZ_700 = getHisto('TT_tZtZ_M700_', 'TT_tZtZ_M700', pDir, var,  tZtZ_700, kRed,    verbose)
#h_tZbW_700 = getHisto('TT_tZbW_M700_', 'TT_tZbW_M700', pDir, var,  tZbW_700, kCyan+4,  verbose)
#h_tZtH_700 = getHisto('TT_tZtH_M700_', 'TT_tZtH_M700', pDir, var,  tZtH_700, kBlue+4, verbose)
#h_tZtZ_800 = getHisto('TT_tZtZ_M800_', 'TT_tZtZ_M800', pDir, var,  tZtZ_800, kGreen+4,    verbose)
#h_tZbW_800 = getHisto('TT_tZbW_M800_', 'TT_tZbW_M800', pDir, var,  tZbW_800, kGreen+3,  verbose)
#h_tZtH_800 = getHisto('TT_tZtH_M800_', 'TT_tZtH_M800', pDir, var,  tZtH_800, kGreen+2, verbose)
h_tZtZ_1000 = getHisto('TT_tZtZ_M1000_', 'TT_tZtZ_M1000', pDir, var,  tZtZ_1000, kBlue-9,    verbose)
h_tZbW_1000 = getHisto('TT_tZbW_M1000_', 'TT_tZbW_M1000', pDir, var,  tZbW_1000, kOrange-9,  verbose)
h_tZtH_1000 = getHisto('TT_tZtH_M1000_', 'TT_tZtH_M1000', pDir, var,  tZtH_1000, kMagenta+1, verbose)

#h_tZ_1000 = getHisto('TT_tZ_M1000_', 'TT_tZ_M1000', pDir, var,  tZ_1000,  kRed,    verbose)

#h_tZtZ_1200 = getHisto('TT_tZtZ_M1200_', 'TT_tZtZ_M1200', pDir, var,  tZtZ_1200, kBlue+4,    verbose)
#h_tZbW_1200 = getHisto('TT_tZbW_M1200_', 'TT_tZbW_M1200', pDir, var,  tZbW_1200, kBlue+3,    verbose)
#h_tZtH_1200 = getHisto('TT_tZtH_M1200_', 'TT_tZtH_M1200', pDir, var,  tZtH_1200, kBlue+2,    verbose)

#h_tZtZ_1400 = getHisto('TT_tZtZ_M1400_', 'TT_tZtZ_M1400', pDir, var,  tZtZ_1400, kBlue+4,    verbose)
#h_tZbW_1400 = getHisto('TT_tZbW_M1400_', 'TT_tZbW_M1400', pDir, var,  tZbW_1400, kBlue+3,    verbose)
#h_tZtH_1400 = getHisto('TT_tZtH_M1400_', 'TT_tZtH_M1400', pDir, var,  tZtH_1400, kBlue+2,    verbose)



#h_bZbZ_700 = getHisto('BB_bZbZ_M700_', 'BB_bZbZ_M700', pDir, var,  bZbZ_700, kRed,    verbose)
#h_bZtW_700 = getHisto('BB_bZtW_M700_', 'BB_bZtW_M700', pDir, var,  tZbW_700, kCyan+4,  verbose)
#h_bZbH_700 = getHisto('BB_bZbH_M700_', 'BB_bZbH_M700', pDir, var,  bZbH_700, kBlue+4, verbose)
#h_bZbZ_800 = getHisto('BB_bZbZ_M800_', 'BB_bZbZ_M800', pDir, var,  bZbZ_800, kRed+4,    verbose)
#h_bZtW_800 = getHisto('BB_bZtW_M800_', 'BB_bZtW_M800', pDir, var,  bZtW_800, kRed+3,  verbose)
#h_bZbH_800 = getHisto('BB_bZbH_M800_', 'BB_bZbH_M800', pDir, var,  bZbH_800, kRed+2, verbose)
#h_bZbZ_1000 = getHisto('BB_bZbZ_M1000_', 'BB_bZbZ_M1000', pDir, var,  bZbZ_1000, kMagenta+4,    verbose)
#h_bZtW_1000 = getHisto('BB_bZtW_M1000_', 'BB_bZtW_M1000', pDir, var,  bZtW_1000, kOrange+11,    verbose)
#h_bZbH_1000 = getHisto('BB_bZbH_M1000_', 'BB_bZbH_M1000', pDir, var,  bZbH_1000, kCyan+1,    verbose)
#h_bZbZ_1200 = getHisto('BB_bZbZ_M1200_', 'BB_bZbZ_M1200', pDir, var,  bZbZ_1200, kMagenta+4,    verbose)
#h_bZtW_1200 = getHisto('BB_bZtW_M1200_', 'BB_bZtW_M1200', pDir, var,  bZtW_1200, kMagenta+3,    verbose)
#h_bZbH_1200 = getHisto('BB_bZbH_M1200_', 'BB_bZbH_M1200', pDir, var,  bZbH_1200, kMagenta+2,    verbose)


#h_data.GetXaxis().SetRangeUser(0, 200)
#h_dy.GetXaxis().SetRangeUser(0, 200)
#h_top.GetXaxis().SetRangeUser(0, 200)
#h_diboson.GetXaxis().SetRangeUser(0, 200)
#h_tZtZ_1000.GetXaxis().SetRangeUser(0, 200)
#h_tZbW_1000.GetXaxis().SetRangeUser(0, 200)
#h_tZtH_1000.GetXaxis().SetRangeUser(0, 200)

templates = []
#templates.append(h_data)
templates.append(h_dy)
templates.append(h_top)
templates.append(h_diboson)
#templates.append(h_vv)
#templates.append(h_st)
#templates.append(h_wjets)
#templates.append(h_tZtZ_700)
#templates.append(h_tZbW_700)
#templates.append(h_tZtH_700)

#templates.append(h_tZtZ_800)
#templates.append(h_tZbW_800)
#templates.append(h_tZtH_800)
templates.append(h_tZtZ_1000)
templates.append(h_tZbW_1000)
templates.append(h_tZtH_1000)

#templates.append(h_dy)
#templates.append(h_tZ_1000)
#templates.append(h_tZtZ_1200)
#templates.append(h_tZbW_1200)
#templates.append(h_tZtH_1200)

#templates.append(h_tZtZ_1400)
#templates.append(h_tZbW_1400)
#templates.append(h_tZtH_1400)

#templates.append(h_bZbZ_700)
#templates.append(h_tZbW_700)
#templates.append(h_bZbH_700)
#templates.append(h_bZbZ_800)
#templates.append(h_bZtW_800)
#templates.append(h_bZbH_800)
#templates.append(h_bZbZ_1000)
#templates.append(h_bZtW_1000)
#templates.append(h_bZbH_1000)
#templates.append(h_bZbZ_1200)
#templates.append(h_bZtW_1200)
#templates.append(h_bZbH_1200)


'''
for x in range (1,17):
    print '{0:<5} & {1:<5.1f} & {2:<5.1f} & {3:<5.1f} & {4:<5.1f} & {5:<5.1f} & {6:<5.1f}  \\\\ '.format('Category'+str(x),templates[0].GetBinContent(x),templates[1].GetBinContent(x),templates[2].GetBinContent(x),templates[3].GetBinContent(x),templates[4].GetBinContent(x),templates[5].GetBinContent(x))
#f = TFile(plotDir+"/"+skimType+"/"+var+".root", "RECREATE")
#for ihist in templates:
#    ihist.Write()
'''
#f.Close()
'''
f = TFile(plotDir+"/"+skimType+"/"+var+".root", "RECREATE")                                                                                             

#f = TFile(plotDir+"/"+"/"+var+".root", "RECREATE")
    
if h_data in templates:
    
    h_data.SetName(var+"__DATA")
    h_data.Write("")
                                                 

if h_dy in templates:
        
    h_dy.SetName(var+"__DY")
    h_dy.Write("")
        
if h_diboson in templates:
            
    h_diboson.SetName(var+"__VV")
    h_diboson.Write("")
            
if h_top in templates:                                                                                                                                                                        
        
    h_top.SetName(var+"__Top")                                                                                                                                                               
    h_top.Write("")    


if h_tZtZ_800 in templates:

    h_tZtZ_800.SetName(var+"__TT_tZtZ_M800")
    h_tZtZ_800.Write("")

if h_tZtH_800 in templates:

    h_tZtH_800.SetName(var+"__TT_tZtH_M800")
    h_tZtH_800.Write("")

if h_tZbW_800 in templates:

    h_tZbW_800.SetName(var+"__TT_tZbW_M800")
    h_tZbW_800.Write("")


#single Tprime
#if h_tZ_1000 in templates:
 #   h_tZ_1000.SetName(var+"__TT_tZ_M1000")
  #  h_tZ_1000.Write("")


if h_tZtZ_1000 in templates:

    h_tZtZ_1000.SetName(var+"__TT_tZtZ_M1000")
    h_tZtZ_1000.Write("")

if h_tZtH_1000 in templates:

    h_tZtH_1000.SetName(var+"__TT_tZtH_M1000")
    h_tZtH_1000.Write("")

if h_tZbW_1000 in templates:

    h_tZbW_1000.SetName(var+"__TT_tZbW_M1000")
    h_tZbW_1000.Write("")


if h_tZtZ_1200 in templates:
        
    h_tZtZ_1200.SetName(var+"__TT_tZtZ_M1200")
    h_tZtZ_1200.Write("")

if h_tZtH_1200 in templates:

    h_tZtH_1200.SetName(var+"__TT_tZtH_M1200")
    h_tZtH_1200.Write("")

if h_tZbW_1200 in templates:
    
    h_tZbW_1200.SetName(var+"__TT_tZbW_M1200")
    h_tZbW_1200.Write("")
     

f.Close()
'''

#get background uncertainty
h_bkg = h_top.Clone()
h_bkg.Reset()
h_bkg.SetName("total bkg")
h_bkg.Add(h_diboson)
h_bkg.Add(h_dy)
h_bkg.Add(h_top)

#h_bkg.Add(h_vv)
#h_bkg.Add(h_wjets)
#h_bkg.Add(h_st)

#histo properties
nBins = h_bkg.GetNbinsX()
bMin = h_bkg.GetBinLowEdge(1)
bMax = h_bkg.GetBinLowEdge(nBins+1)
bin1 = h_bkg.GetXaxis().FindBin(bMin)
bin2 = h_bkg.GetXaxis().FindBin(bMax)

for ibin in range(0,nBins+1):    
    iTop     = h_top.GetBinContent(ibin)
    iDY      = h_dy.GetBinContent(ibin)
    iDiboson = h_diboson.GetBinContent(ibin)
 #   iWJ      = h_wjets.GetBinContent(ibin)
  #  iVV      = h_vv.GetBinContent(ibin)
    # stat error
    stat_err = (h_bkg.GetBinError(ibin))**2 
    # add approximate systematic uncertainty to each bin
    lumi_err = 0.026**2
    btag_err = 0.017**2
    ID_err   = 0.03**2
    JES_err  = 0.05*0.05
    dy_err   = (0.15*iDY)**2
    top_err  = (0.15*iTop)**2
    diboson_err= (0.2*iDiboson)**2
   # st_err   = (0.3*iTop)**2
   # wjet_err = (0.1*iWJ)**2
   # vv_err   = (0.3*iVV)**2
    scale_err = ((iDY+iTop+iDiboson)*0.3)**2
    #new_err = stat_err + lumi_err + btag_err + ID_err + JES_err + dy_err + top_err   + diboson_err# +st_err  + wjet_err +st_err + vv_err
    new_err = stat_err
    if h_bkg.GetBinError(ibin) != 0: h_bkg.SetBinError(ibin, TMath.Sqrt(new_err))

h_bkg.SetMarkerSize(0)
h_bkg.SetLineWidth(2)
h_bkg.SetFillColor(14)
h_bkg.SetLineColor(0)
h_bkg.SetFillStyle(3244)

#histogram to print the total background with stat uncertainty
h_tot = h_top.Clone()
h_tot.Reset()
h_tot.SetName("Total_"+h_tot.GetName().split('_',1)[1])
h_tot.Add(h_top)
h_tot.Add(h_dy)
h_tot.Add(h_diboson)
#h_tot.Add(h_wjets)
#h_tot.Add(h_st)

#if len(data)>0:
 #   overUnderFlow(h_data)
#overUnderFlow(h_bkg)
#overUnderFlow(h_tot)
print h_tot.GetName().split('_',1)[1]

## =========Drawing==============
#integralError = Double(5)
# print the latex table:
print '\\begin{tabular}{|c|c| }'
print '\hline'
print 'Sample     & Events  \\\\ '
print '\hline'
count = 0
for ihist in templates :
    integralError = Double(0.05)
    
    #if var != 'cutflow' or var != 'cutflow4':overUnderFlow(ihist)
    count = count+1
    if count == 10:
        print '\hline'
    if count == 14:
        print '\hline'
    ihist.IntegralAndError(bin1,bin2,integralError)
    if 'TT' in ihist.GetName() or 'BB' in ihist.GetName():
       # print "hist name", ihist.GetName()
       # print '{0:<5} & {1:<5.1f} $\pm$ {2:<5.1f} \\\\ '.format(ihist.GetName().split('_')[1], ihist.Integral(bin1,bin2), integralError)
        if len(ihist.GetName())==15:
            print '{0:<5} & {1:<5.1f} $\pm$ {2:<5.1f} \\\\ '.format(ihist.GetName()[0:12], ihist.Integral(bin1,bin2), integralError)
        else:
            print '{0:<5} & {1:<5.1f} $\pm$ {2:<5.1f} \\\\ '.format(ihist.GetName()[0:13], ihist.Integral(bin1,bin2), integralError)  
    else:      
        print '{0:<5} & {1:<5.1f} $\pm$ {2:<5.1f} \\\\ '.format(ihist.GetName().split('_')[0], ihist.Integral(bin1,bin2), integralError)

print '\hline'

print "Tot Bkg",  h_tot.IntegralAndError(bin1,bin2,integralError),"$\pm$", integralError
print '{0:<5} & {1:<5.0f} \\\\ '.format('Tot Bkg', h_tot.Integral(bin1,bin2), integralError)
#print '{0:<5} & {1:<5.0f} \\\\ '.format('Tot Bkg', h_tot.IntegralAndError(bin1,bin2,integralError), integralError)
print '\hline'
if len(data)>0:
    print '{0:<5} & {1:<5.0f} \\\\ '.format(h_data.GetName().split('_')[0], h_data.Integral())
print '\end{tabular}'
print 'bkg : ', h_bkg.Integral(ibin,bin2), 'tot : ', h_tot.Integral(ibin,bin2)



for h in templates :
    bmin = h.GetXaxis().FindBin(120.0)
    bmax = h.GetXaxis().FindBin(bin2)
    h.IntegralAndError(bmin,bmax,integralError)
    print  '{0:<5} & {1:<5.1f} $\pm$ {2:<5.1f}\\\\'.format(h.GetName(), h.Integral(bmin,bmax), integralError)
  #  values.append(h.Integral(bmin,bmax))                                                                                                                            
   # errors.append(integralError)                                                                                                                                    
print "total bkg in the range ",h_bkg.IntegralAndError(bmin,bmax,integralError),"$\pm$",integralError
if len(data)>0:
    print '{0:<5} & {1:<5.0f} \\\\ '.format(h_data.GetName().split('_')[0], h_data.IntegralAndError(bmin,bmax,integralError)) 


hs = THStack("","")

for ihist in reversed(templates[0:3]):# 5 was here
    hs.Add(ihist)
    print 'histo added', ihist.GetName()

# Canvas
c1 = TCanvas('c1', 'c1', 800, 600)

if len(data)>0:
    c1.Divide(1,2)
scale = (1.0 - 0.3)/0.35

# prepare top pad for original plot
pad = c1.cd(1)
if len(data)>0:
    pad.SetPad(0, 0.3, 1, 1)
    pad.SetTopMargin(0.1) 
else:
    pad.SetPad(0, 0.0, 1, 1)
#pad.SetTopMargin(0.1)

#if len(data) == 0:
 #   pad.SetBottomMargin(0.0)
if len(data) > 0:   
    pad.SetBottomMargin(0.005)
    pad.SetTickx(1)
t = pad.GetTopMargin()
 

# prepare the 2nd pad
pad = c1.cd(2)
if len(data)>0:
    pad.SetPad(0, 0.0, 1, 0.3)
#else:
 #   pad.SetPad(0, 0.0, 1, 0.0)
    pad.SetTopMargin(0.06)
    pad.SetBottomMargin(0.4)
    pad.SetTickx(1)
    pad.SetTicky(1)
c1.cd(1)

#hs.SetMaximum(hs.GetMaximum()*5)
#hs.SetMaximum(hs.GetMaximum()+5.0) 
#hs.SetMaximum(hs.GetMaximum("nostack"))

hs.SetMinimum(0.1)#
gPad.SetLogy()
hs.Draw("Hist")
#h_dy.Draw("l same")
#if var == "massak4jet1_pre" or var == "massak4jet2_pre" or "massak4jet1_cnt" or var == "massak4jet2_cnt" or "massak4jet1" or var == "massak4jet2"  :
 #   hs.GetXaxis().SetRangeUser(0, 200)
  #  print "var 333 ", var , "hist range" ,  hs.GetXaxis().GetXmin(),",",  hs.GetXaxis().GetXmax()

#print "var 444 ", var , "hist range" ,  hs.GetXaxis().GetXmin(),",",  hs.GetXaxis().GetXmax() 
#if var == "cutflow":
#hs.GetXaxis().SetRangeUser(3,10)
#hs.GetXaxis().SetRangeUser(6,10) 
#hs.GetXaxis().SetRangeUser(0,200)
#hs.GetXaxis().SetRangeUser(600,1000)  
hs.GetXaxis().Draw()

h_bkg.Draw("e2 same")
if len(data)>0:
    h_data.Draw("same")
maximum = 0
for ihist in reversed(templates[3:13]): #5:8
    print 'overlaying, ', ihist.GetName() 
    ihist.Draw("ehist same")
# put this part to change cnavas size according to higest vlaues of all histograms including signal
    if  ihist.GetBinContent(ihist.GetMaximumBin()) > maximum:
        maximum = ihist.GetBinContent(ihist.GetMaximumBin())
    if maximum > hs.GetMaximum(""):
        hs.SetMaximum(maximum*50)
    else:
        hs.SetMaximum(hs.GetMaximum()*50)
if len(data)>0:
    if hs.GetMaximum("") < h_data.GetBinContent(h_data.GetMaximumBin()):
        hs.SetMaximum(h_data.GetMaximum()*50)

#hs.SetMaximum(h_data.GetMaximum()*10450)

xTitle= h_top.GetXaxis().GetTitle() # unmakrk
yTitle= h_top.GetYaxis().GetTitle()
setTitle(hs, xTitle)
#hs.GetXaxis().SetTitle(xTitle)
 

gPad.RedrawAxis()

ll = TLatex()
ll.SetNDC(kTRUE)
ll.SetTextSize(0.05)
if len(data)>0:
    ll.DrawLatex(0.78,0.92, "35.9 fb^{-1} (13 TeV)");#35.9,2.2 0r 0.53
else:
    ll.DrawLatex(0.74,0.92, "35.9 fb^{-1} (13 TeV)");#2.2 0r 0.53 

cms = TLatex()
cms.SetNDC(kTRUE)
cms.SetTextFont(61)
cms.SetTextSize(0.08)
cms.DrawLatex(0.12, 1-t+0.2*t,"CMS")

sel = TLatex()
sel.SetNDC(kTRUE)
sel.SetTextSize(0.065)

chan = TLatex()
chan.SetNDC(kTRUE)
chan.SetTextSize(0.065)
chan.DrawLatex(0.50, 0.76, title)

prel = TLatex()
prel.SetNDC(kTRUE)
prel.SetTextFont(52)
prel.SetTextSize(0.75*t*0.76)
if len(data)>0:
    prel.DrawLatex(0.22,0.92,"Preliminary")
else:
    prel.DrawLatex(0.24,0.92,"Preliminary")

leg.Draw()
gPad.RedrawAxis()




c1.cd(2)
# add the systematic band
if len(data)>0:
    h_ratio = h_data.Clone()
    h_ratio_bkg = h_bkg.Clone()
    h_ratio_bkg.SetDirectory(0)
    h_ratio.SetDirectory(0)
    h_ratio.Divide(h_data, h_tot)
    h_ratio_bkg.Divide(h_bkg, h_tot)
    
    for ibin in range(1, nBins+1):
        if h_bkg.GetBinContent(ibin) == 0: h_ratio_bkg.SetBinContent(ibin,1)
        
        prepareRatio(h_ratio, h_ratio_bkg, scale, xTitle)
        
        line = TLine(bMin, 1, bMax, 1)
        #line = TLine(75, 1, 105, 1)     
        #line.SetLineColor(kBlack)
       # if var == "massak4jet1_pre" or var == "massak4jet2_pre" or "massak4jet1_cnt" or var == "massak4jet2_cnt" or "massak4jet1" or var == "massak4jet2"  :        
        #h_ratio.GetXaxis().SetRangeUser(0, 200);
        #h_ratio_bkg.GetXaxis().SetRangeUser(0, 200);
        #line = TLine(0, 1, 200, 1) 
        #if var == "cutflow":
        #h_ratio.GetXaxis().SetRangeUser(3,10)
        #h_ratio_bkg.GetXaxis().SetRangeUser(3,10)
        #line = TLine(2.5, 1, bMax, 1)
        
        #h_ratio.GetXaxis().SetRangeUser(6,10)                                                                       
        #h_ratio_bkg.GetXaxis().SetRangeUser(6,10)                                                                  
        #line = TLine(5.5, 1, bMax, 1) 

        #h_ratio.GetXaxis().SetRangeUser(600, 1000);                                                                                                 
        #h_ratio_bkg.GetXaxis().SetRangeUser(600, 1000);                                                                                      
        #line = TLine(600, 1, 1000, 1)               

        h_ratio.Draw("")
        h_ratio_bkg.Draw("e2same")
        h_ratio.Draw("same")
        line.Draw()
        #h_ratio.GetXaxis().SetBinLabel(6,"Pre Selection")
        #leg.AddEntry(h_ratio_bkg, 'stat unc.', 'f')
        gPad.RedrawAxis()
leg.AddEntry(h_ratio_bkg, 'stat unc.', 'f')
#create a directory if it doesn't exist
m_1 = 'mkdir '+plotDir
m_2 = 'mkdir '+plotDir+"/"+skimType
if not os.path.isdir(plotDir):
    subprocess.call( [m_1], shell=True )
if not os.path.isdir(plotDir+"/"+skimType):
    subprocess.call( [m_2], shell=True )    
    
#c1.SaveAs(plotDir+"/"+skimType+"/"+var+"_.pdf")
c1.SaveAs(plotDir+"/"+"/"+var+"_.pdf")
#if len(data)==0:
 #   c1.cd(1).SaveAs(plotDir+"/"+skimType+"/"+var+str(1)+"_.pdf")

#c1.SaveAs(plotDir+"/"+skimType+"/"+var+"_.gif")

#raw_input("hold on")
