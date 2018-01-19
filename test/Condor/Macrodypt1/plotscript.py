#! /usr/bin/env python
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
                  default= 35700., #2263.,    2200., 538.134,
                  dest='Lumi',
                  help='Data Luminosity in pb-1')

parser.add_option('--globalSF', metavar='SF', type='float',
                  default = 1.0,
                  dest='globalSF',
                  help='Global trigger SF (%default default)')

parser.add_option('--plotDir', metavar='P', type='string', action='store',
                  default='cutflowel',
                  dest='plotDir',
                  help='output directory of plots')

parser.add_option('--skimType', metavar='S', type='string', action='store',
                  default='CR_Zelel',
                  dest='skimType',
                  help='Skim type: CR_Zelel, CR_Zmumu, SR_Zelel, SR_Zmumu')

parser.add_option('--processDir', metavar='pD', type='string', action='store',
                  default='ana',
                  dest='processDir',
                  help='directory to read histograms from: ana, anabcUp, anabcDown, analightUp, analightDown,anaJecUp,anaJecDown, anaJerUp,anaJerDown, anaPileupUp,anaPileupDown')

parser.add_option('--processtype', metavar='pT', type='string', action='store',
                  default='',
                  dest='processtype',
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
pType = options.processtype
plotDir = options.plotDir
skimType = options.skimType
verbose = options.verbose

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
             [f_Data2, 1., 1., 1.],  #addign this should give 2.2 fb-1             
             [f_Data3, 1., 1., 1.],  #addign this should give 2.2 fb-1                
             [f_Data4, 1., 1., 1.],  #addign this should give 2.2 fb-1   
             [f_Data5, 1., 1., 1.],  #addign this should give 2.2 fb-1   
             [f_Data6, 1., 1., 1.],  #addign this should give 2.2 fb-1   
             [f_Data7, 1., 1., 1.],  #addign this should give 2.2 fb-1   
             ]

top      = [[f_ttbar,         Top_xs,            Top_num,            lumi]]

#dy       = [[f_DYmcnlo,       DY_xs,             DYmcnlo_num,        lumi]]

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

dy       = [
    [f_DY100to250,    DY100to250_xs,      DY100to250_num,    lumi],
    [f_DY250to400,    DY250to400_xs,      DY250to400_num,    lumi],
    [f_DY400to650,    DY400to650_xs,      DY400to650_num,    lumi],
    [f_DY650toInf,    DY650toInf_xs,      DY650toInf_num,    lumi],
    # [f_DY800to1200,   DY800to1200_xs,     DY800to1200_num,    lumi],                                                                                                       
    # [f_DY1200to2500,  DY1200to2500_xs,    DY1200to2500_num,   lumi],                                                                                                       
    # [f_DY2500toInf,   DY2500toInf_xs,     DY2500toInf_num,    lumi],                                                                                                       
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
#            [f_ST_tW_top,     ST_tW_top_xs,       ST_tW_top_num,     lumi],
#            [f_ST_tW_antitop, ST_tW_antitop_xs,   ST_tW_antitop_num, lumi],
#            [f_ST_t,          ST_t_xs,            ST_t_num,          lumi],
#            [f_ST_t_ex1,      ST_t_xs,            ST_t_ex1_num,      lumi],
#            [f_ST_s,          ST_s_xs,            ST_s_num,          lumi], 
#           ]

#vv       = [
#            [f_ZZTo2L2Nu,     ZZTo2L2Nu_xs,       ZZTo2L2Nu_num,    lumi],
#            [f_WZTo2L2Q,      WZTo2L2Q_xs,        WZTo2L2Q_num,     lumi],
#            [f_WWTo2L2Nu,     WWTo2L2Nu_xs,       WWTo2L2Nu_num,    lumi],
#           ]
#tZtZ_700 = [[f_BpBp_tZtZ_700, TpTp700_xs,         TpTp700_num,       lumi]]
#tZbW_700 = [[f_BpTp_tZbW_700, TpTp700_xs,         TpTp700_num,       lumi]]
#tZtH_700 = [[f_BpTp_tZtH_700, TpTp700_xs,         TpTp700_num,       lumi]]
tZtZ_800 = [[f_TpTp_tZtZ_800, TpTp800_xs,         TpTp800tZtZ_num,       lumi]]
tZbW_800 = [[f_TpTp_tZbW_800, TpTp800_xs,         TpTp800tZbW_num,       lumi]]
tZtH_800 = [[f_TpTp_tZtH_800, TpTp800_xs,         TpTp800tZtH_num,       lumi]]

tZtZ_1000 = [[f_TpTp_tZtZ_1000, TpTp1000_xs,         TpTp1000tZtZ_num,       lumi]]

tZbW_1000 = [[f_TpTp_tZbW_1000, TpTp1000_xs,         TpTp1000tZbW_num,       lumi]]
tZtH_1000 = [[f_TpTp_tZtH_1000, TpTp1000_xs,         TpTp1000tZtH_num,       lumi]]

tZtZ_1200 = [[f_TpTp_tZtZ_1200, TpTp1200_xs,         TpTp1200tZtZ_num,       lumi]]
tZbW_1200 = [[f_TpTp_tZbW_1200, TpTp1200_xs,         TpTp1200tZbW_num,       lumi]]
tZtH_1200 = [[f_TpTp_tZtH_1200, TpTp1200_xs,         TpTp1200tZtH_num,       lumi]]

tZtZ_1100 = [[f_TpTp_tZtZ_1100, TpTp1100_xs,         TpTp1100tZtZ_num,       lumi]]
tZbW_1100 = [[f_TpTp_tZbW_1100, TpTp1100_xs,         TpTp1100tZbW_num,       lumi]]
tZtH_1100 = [[f_TpTp_tZtH_1100, TpTp1100_xs,         TpTp1100tZtH_num,       lumi]]

tZtZ_900 = [[f_TpTp_tZtZ_900, TpTp900_xs,         TpTp900tZtZ_num,       lumi]]
tZbW_900 = [[f_TpTp_tZbW_900, TpTp900_xs,         TpTp900tZbW_num,       lumi]]
tZtH_900 = [[f_TpTp_tZtH_900, TpTp900_xs,         TpTp900tZtH_num,       lumi]]

tZtZ_1300 = [[f_TpTp_tZtZ_1300, TpTp1300_xs,         TpTp1300tZtZ_num,       lumi]]
tZbW_1300 = [[f_TpTp_tZbW_1300, TpTp1300_xs,         TpTp1300tZbW_num,       lumi]]
tZtH_1300 = [[f_TpTp_tZtH_1300, TpTp1300_xs,         TpTp1300tZtH_num,       lumi]]

tZtZ_1400 = [[f_TpTp_tZtZ_1400, TpTp1400_xs,         TpTp1400tZtZ_num,       lumi]]
tZbW_1400 = [[f_TpTp_tZbW_1400, TpTp1400_xs,         TpTp1400tZbW_num,       lumi]]
tZtH_1400 = [[f_TpTp_tZtH_1400, TpTp1400_xs,         TpTp1400tZtH_num,       lumi]]

tZtZ_1500 = [[f_TpTp_tZtZ_1500, TpTp1500_xs,         TpTp1500tZtZ_num,       lumi]]
tZbW_1500 = [[f_TpTp_tZbW_1500, TpTp1500_xs,         TpTp1500tZbW_num,       lumi]]
tZtH_1500 = [[f_TpTp_tZtH_1500, TpTp1500_xs,         TpTp1500tZtH_num,       lumi]]

tZtZ_1600 = [[f_TpTp_tZtZ_1600, TpTp1600_xs,         TpTp1600tZtZ_num,       lumi]]
tZbW_1600 = [[f_TpTp_tZbW_1600, TpTp1600_xs,         TpTp1600tZbW_num,       lumi]]
tZtH_1600 = [[f_TpTp_tZtH_1600, TpTp1600_xs,         TpTp1600tZtH_num,       lumi]]

tZtZ_1700 = [[f_TpTp_tZtZ_1700, TpTp1700_xs,         TpTp1700tZtZ_num,       lumi]]
tZbW_1700 = [[f_TpTp_tZbW_1700, TpTp1700_xs,         TpTp1700tZbW_num,       lumi]]
tZtH_1700 = [[f_TpTp_tZtH_1700, TpTp1700_xs,         TpTp1700tZtH_num,       lumi]]

tZtZ_1800 = [[f_TpTp_tZtZ_1800, TpTp1800_xs,         TpTp1800tZtZ_num,       lumi]]
tZbW_1800 = [[f_TpTp_tZbW_1800, TpTp1800_xs,         TpTp1800tZbW_num,       lumi]]
tZtH_1800 = [[f_TpTp_tZtH_1800, TpTp1800_xs,         TpTp1800tZtH_num,       lumi]]


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
    h_data     = getHisto(dataLabel,       dataLeg,        pDir+'/'+pType, var,  data,     kBlack,     verbose)



h_top      = getHisto(topLabel,        topLeg,         pDir+'/'+pType, var,  top,      8,          verbose)
h_dy       = getHisto(dyLabel,         dyLeg,          pDir+'/'+pType, var,  dy,       90,         verbose)
h_diboson  = getHisto(dibosonLabel,  dibosonLeg,       pDir+'/'+pType, var,  diboson,   kBlue,         verbose)
#h_wjets    = getHisto(wjLabel,         wjLeg,          pDir+'/'+pType, var,  wjets,    kBlue,      verbose)
#h_st       = getHisto(sTLabel,         sTLeg,          pDir+'/'+pType, var,  st,       kCyan,      verbose)
#h_vv       = getHisto(vvLabel,         vvLeg,          pDir+'/'+pType, var,  vv,       kRed,       verbose)

#h_tZtZ_700 = getHisto('TT_tZtZ_M700_', 'TT_tZtZ_M700', pDir+'/'+pType, var,  tZtZ_700, kRed,    verbose)
#h_tZbW_700 = getHisto('TT_tZbW_M700_', 'TT_tZbW_M700', pDir+'/'+pType, var,  tZbW_700, kCyan+4,  verbose)
#h_tZtH_700 = getHisto('TT_tZtH_M700_', 'TT_tZtH_M700', pDir+'/'+pType, var,  tZtH_700, kBlue+4, verbose)
h_tZtZ_800 = getHisto('TT_tZtZ_M800_', 'TT_tZtZ_M800', pDir+'/'+pType, var,  tZtZ_800, kGreen+4,    verbose)
h_tZbW_800 = getHisto('TT_tZbW_M800_', 'TT_tZbW_M800', pDir+'/'+pType, var,  tZbW_800, kGreen+3,  verbose)
h_tZtH_800 = getHisto('TT_tZtH_M800_', 'TT_tZtH_M800', pDir+'/'+pType, var,  tZtH_800, kGreen+2, verbose)
h_tZtZ_1000 = getHisto('TT_tZtZ_M1000_', 'TT_tZtZ_M1000', pDir+'/'+pType, var,  tZtZ_1000, kBlue-9,    verbose)
h_tZbW_1000 = getHisto('TT_tZbW_M1000_', 'TT_tZbW_M1000', pDir+'/'+pType, var,  tZbW_1000, kOrange-9,  verbose)
h_tZtH_1000 = getHisto('TT_tZtH_M1000_', 'TT_tZtH_M1000', pDir+'/'+pType, var,  tZtH_1000, kMagenta+1, verbose)
h_tZtZ_1200 = getHisto('TT_tZtZ_M1200_', 'TT_tZtZ_M1200', pDir+'/'+pType, var,  tZtZ_1200, kBlue+4,    verbose)
h_tZbW_1200 = getHisto('TT_tZbW_M1200_', 'TT_tZbW_M1200', pDir+'/'+pType, var,  tZbW_1200, kBlue+3,    verbose)
h_tZtH_1200 = getHisto('TT_tZtH_M1200_', 'TT_tZtH_M1200', pDir+'/'+pType, var,  tZtH_1200, kBlue+2,    verbose)


h_tZtZ_900 = getHisto('TT_tZtZ_M900_', 'TT_tZtZ_M900', pDir+'/'+pType, var,  tZtZ_900, kGreen+4,    verbose)
h_tZbW_900 = getHisto('TT_tZbW_M900_', 'TT_tZbW_M900', pDir+'/'+pType, var,  tZbW_900, kGreen+3,  verbose)
h_tZtH_900 = getHisto('TT_tZtH_M900_', 'TT_tZtH_M900', pDir+'/'+pType, var,  tZtH_900, kGreen+2, verbose)
'''
h_tZtZ_1100 = getHisto('TT_tZtZ_M1100_', 'TT_tZtZ_M1100', pDir+'/'+pType, var,  tZtZ_1100, kGreen+4,    verbose)
h_tZbW_1100 = getHisto('TT_tZbW_M1100_', 'TT_tZbW_M1100', pDir+'/'+pType, var,  tZbW_1100, kGreen+3,  verbose)
h_tZtH_1100 = getHisto('TT_tZtH_M1100_', 'TT_tZtH_M1100', pDir+'/'+pType, var,  tZtH_1100, kGreen+2, verbose)

h_tZtZ_1300 = getHisto('TT_tZtZ_M1300_', 'TT_tZtZ_M1300', pDir+'/'+pType, var,  tZtZ_1300, kGreen+4,    verbose)
h_tZbW_1300 = getHisto('TT_tZbW_M1300_', 'TT_tZbW_M1300', pDir+'/'+pType, var,  tZbW_1300, kGreen+3,  verbose)
h_tZtH_1300 = getHisto('TT_tZtH_M1300_', 'TT_tZtH_M1300', pDir+'/'+pType, var,  tZtH_1300, kGreen+2, verbose)

h_tZtZ_1400 = getHisto('TT_tZtZ_M1400_', 'TT_tZtZ_M1400', pDir+'/'+pType, var,  tZtZ_1400, kGreen+4,    verbose)
h_tZbW_1400 = getHisto('TT_tZbW_M1400_', 'TT_tZbW_M1400', pDir+'/'+pType, var,  tZbW_1400, kGreen+3,  verbose)
h_tZtH_1400 = getHisto('TT_tZtH_M1400_', 'TT_tZtH_M1400', pDir+'/'+pType, var,  tZtH_1400, kGreen+2, verbose)
'''
h_tZtZ_1500 = getHisto('TT_tZtZ_M1500_', 'TT_tZtZ_M1500', pDir+'/'+pType, var,  tZtZ_1500, kGreen+4,    verbose)
h_tZbW_1500 = getHisto('TT_tZbW_M1500_', 'TT_tZbW_M1500', pDir+'/'+pType, var,  tZbW_1500, kGreen+3,  verbose)
h_tZtH_1500 = getHisto('TT_tZtH_M1500_', 'TT_tZtH_M1500', pDir+'/'+pType, var,  tZtH_1500, kGreen+2, verbose)
'''
h_tZtZ_1600 = getHisto('TT_tZtZ_M1600_', 'TT_tZtZ_M1600', pDir+'/'+pType, var,  tZtZ_1600, kGreen+4,    verbose)
h_tZbW_1600 = getHisto('TT_tZbW_M1600_', 'TT_tZbW_M1600', pDir+'/'+pType, var,  tZbW_1600, kGreen+3,  verbose)
h_tZtH_1600 = getHisto('TT_tZtH_M1600_', 'TT_tZtH_M1600', pDir+'/'+pType, var,  tZtH_1600, kGreen+2, verbose)

h_tZtZ_1700 = getHisto('TT_tZtZ_M1700_', 'TT_tZtZ_M1700', pDir+'/'+pType, var,  tZtZ_1700, kGreen+4,    verbose)
h_tZbW_1700 = getHisto('TT_tZbW_M1700_', 'TT_tZbW_M1700', pDir+'/'+pType, var,  tZbW_1700, kGreen+3,  verbose)
h_tZtH_1700 = getHisto('TT_tZtH_M1700_', 'TT_tZtH_M1700', pDir+'/'+pType, var,  tZtH_1700, kGreen+2, verbose)

h_tZtZ_1800 = getHisto('TT_tZtZ_M1800_', 'TT_tZtZ_M1800', pDir+'/'+pType, var,  tZtZ_1800, kGreen+4,    verbose)
h_tZbW_1800 = getHisto('TT_tZbW_M1800_', 'TT_tZbW_M1800', pDir+'/'+pType, var,  tZbW_1800, kGreen+3,  verbose)
h_tZtH_1800 = getHisto('TT_tZtH_M1800_', 'TT_tZtH_M1800', pDir+'/'+pType, var,  tZtH_1800, kGreen+2, verbose)
'''
#h_bZbZ_700 = getHisto('BB_bZbZ_M700_', 'BB_bZbZ_M700', pDir+'/'+pType, var,  bZbZ_700, kRed,    verbose)
#h_bZtW_700 = getHisto('BB_bZtW_M700_', 'BB_bZtW_M700', pDir+'/'+pType, var,  tZbW_700, kCyan+4,  verbose)
#h_bZbH_700 = getHisto('BB_bZbH_M700_', 'BB_bZbH_M700', pDir+'/'+pType, var,  bZbH_700, kBlue+4, verbose)
#h_bZbZ_800 = getHisto('BB_bZbZ_M800_', 'BB_bZbZ_M800', pDir+'/'+pType, var,  bZbZ_800, kRed+4,    verbose)
#h_bZtW_800 = getHisto('BB_bZtW_M800_', 'BB_bZtW_M800', pDir+'/'+pType, var,  bZtW_800, kRed+3,  verbose)
#h_bZbH_800 = getHisto('BB_bZbH_M800_', 'BB_bZbH_M800', pDir+'/'+pType, var,  bZbH_800, kRed+2, verbose)
#h_bZbZ_1000 = getHisto('BB_bZbZ_M1000_', 'BB_bZbZ_M1000', pDir+'/'+pType, var,  bZbZ_1000, kMagenta+4,    verbose)
#h_bZtW_1000 = getHisto('BB_bZtW_M1000_', 'BB_bZtW_M1000', pDir+'/'+pType, var,  bZtW_1000, kOrange+11,    verbose)
#h_bZbH_1000 = getHisto('BB_bZbH_M1000_', 'BB_bZbH_M1000', pDir+'/'+pType, var,  bZbH_1000, kCyan+1,    verbose)
#h_bZbZ_1200 = getHisto('BB_bZbZ_M1200_', 'BB_bZbZ_M1200', pDir+'/'+pType, var,  bZbZ_1200, kMagenta+4,    verbose)
#h_bZtW_1200 = getHisto('BB_bZtW_M1200_', 'BB_bZtW_M1200', pDir+'/'+pType, var,  bZtW_1200, kMagenta+3,    verbose)
#h_bZbH_1200 = getHisto('BB_bZbH_M1200_', 'BB_bZbH_M1200', pDir+'/'+pType, var,  bZbH_1200, kMagenta+2,    verbose)




templates = []
templates.append(h_data)
templates.append(h_dy)
templates.append(h_top)
#templates.append(h_diboson)
#templates.append(h_vv)
#templates.append(h_st)
#templates.append(h_wjets)
#templates.append(h_tZtZ_700)
#templates.append(h_tZbW_700)
#templates.append(h_tZtH_700)
templates.append(h_tZtZ_800)
templates.append(h_tZbW_800)
templates.append(h_tZtH_800)
templates.append(h_tZtZ_1000)
templates.append(h_tZbW_1000)
templates.append(h_tZtH_1000)
templates.append(h_tZtZ_1200)
templates.append(h_tZbW_1200)
templates.append(h_tZtH_1200)

templates.append(h_tZtZ_900)
templates.append(h_tZbW_900)
templates.append(h_tZtH_900)
'''
templates.append(h_tZtZ_1100)
templates.append(h_tZbW_1100)
templates.append(h_tZtH_1100)

templates.append(h_tZtZ_1300)
templates.append(h_tZbW_1300)
templates.append(h_tZtH_1300)

templates.append(h_tZtZ_1400)
templates.append(h_tZbW_1400)
templates.append(h_tZtH_1400)
'''
templates.append(h_tZtZ_1500)
templates.append(h_tZbW_1500)
templates.append(h_tZtH_1500)
'''
templates.append(h_tZtZ_1600)
templates.append(h_tZbW_1600)
templates.append(h_tZtH_1600)

templates.append(h_tZtZ_1700)
templates.append(h_tZbW_1700)
templates.append(h_tZtH_1700)

templates.append(h_tZtZ_1800)
templates.append(h_tZbW_1800)
templates.append(h_tZtH_1800)
'''

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

#f = TFile(plotDir+"/"+skimType+"/"+var+".root", "RECREATE")
#for ihist in templates:
#    ihist.Write()

#f.Close()


m_1 = 'mkdir '+plotDir
#m_2 = 'mkdir '+plotDir+"/"+skimType                                                                                                                                              
m_2 = 'mkdir '+plotDir+"/"+pDir
if not os.path.isdir(plotDir):
    subprocess.call( [m_1], shell=True )
#if not os.path.isdir(plotDir+"/"+skimType):
 #   subprocess.call( [m_2], shell=True )

#if not os.path.isdir(plotDir+"/"+pDir):                                                                                                                                    
 #   subprocess.call( [m_2], shell=True )   


print "variable" , var
#print " pdir 0 nad 1 ", pDir.split('/',1)[0], pDir.split('/',1)[1]
print "set of variable ", var[0:2]
#f = TFile(plotDir+"/"+skimType+"/"+var+".root", "RECREATE")                                                                                             

#if len(var)>6:
 #   name = str(var[6:14])
#else:
    #name = str(var[3:6])

name = str(var[0:3])
print "name " , name
#f = TFile(plotDir+"/"+pDir+"/"+"mumu"+name+".root", "RECREATE")
print os.path.isfile(plotDir+"/"+"nominal.root")

if os.path.isfile(plotDir+"/"+"nominal.root") == False : f = TFile(plotDir+"/"+"nominal.root", "RECREATE")
else : f = TFile(plotDir+"/"+"nominal.root", "update")  

if pDir == 'anabcUp': sub = "__bc__plus"
if pDir== 'anabcDown': sub = "__bc__minus"
if pDir== 'analightUp': sub = "__light__plus"
if pDir== 'analightDown': sub = "__light__minus"
if pDir== 'anaJecUp': sub = "__jec__plus"
if pDir== 'anaJecDown': sub = "__jec__minus"
if pDir== 'anaJerUp': sub = "__jer__plus"
if pDir== 'anaJerDown': sub = "__jer__minus"
if pDir== 'anaPileupUp': sub = "__pileup__plus"
if pDir== 'anaPileupDown': sub = "__pileup__minus"
if pDir== 'ana': sub = ""



if h_data in templates:
    
    h_data.SetName("ee"+name+"__Data"+sub)
    h_data.Write("")
                                  

if h_dy in templates:
    h_dy.SetName("ee"+name+"__DY"+sub)
    h_dy.Write("")
        
#if h_diboson in templates:
            
 #   h_diboson.SetName("ee"+name+"__VV+sub")
  #  h_diboson.Write("")
            
if h_top in templates:                                                                                                                                                                        
        
    h_top.SetName("ee"+name+"__Top"+sub)                                                                                                                                                               
    h_top.Write("")    


if h_tZtZ_800 in templates:

    h_tZtZ_800.SetName("ee"+name+"__TT_tZtZ_M800"+sub)
    h_tZtZ_800.Write("")

if h_tZtH_800 in templates:

    h_tZtH_800.SetName("ee"+name+"__TT_tZtH_M800"+sub)
    h_tZtH_800.Write("")

if h_tZbW_800 in templates:

    h_tZbW_800.SetName("ee"+name+"__TT_tZbW_M800"+sub)
    h_tZbW_800.Write("")

if h_tZtZ_1000 in templates:

    h_tZtZ_1000.SetName("ee"+name+"__TT_tZtZ_M1000"+sub)
    h_tZtZ_1000.Write("")

if h_tZtH_1000 in templates:

    h_tZtH_1000.SetName("ee"+name+"__TT_tZtH_M1000"+sub)
    h_tZtH_1000.Write("")

if h_tZbW_1000 in templates:

    h_tZbW_1000.SetName("ee"+name+"__TT_tZbW_M1000"+sub)
    h_tZbW_1000.Write("")

if h_tZtZ_1200 in templates:
        
    h_tZtZ_1200.SetName("ee"+name+"__TT_tZtZ_M1200"+sub)
    h_tZtZ_1200.Write("")

if h_tZtH_1200 in templates:

    h_tZtH_1200.SetName("ee"+name+"__TT_tZtH_M1200"+sub)
    h_tZtH_1200.Write("")

if h_tZbW_1200 in templates:
    
    h_tZbW_1200.SetName("ee"+name+"__TT_tZbW_M1200"+sub)
    h_tZbW_1200.Write("")
        

if h_tZtZ_900 in templates:

    h_tZtZ_900.SetName("ee"+name+"__TT_tZtZ_M900"+sub)
    h_tZtZ_900.Write("")

if h_tZtH_900 in templates:

    h_tZtH_900.SetName("ee"+name+"__TT_tZtH_M900"+sub)
    h_tZtH_900.Write("")

if h_tZbW_900 in templates:

    h_tZbW_900.SetName("ee"+name+"__TT_tZbW_M900"+sub)
    h_tZbW_900.Write("")
'''

if h_tZtZ_1100 in templates:

    h_tZtZ_1100.SetName("ee"+name+"__TT_tZtZ_M1100"+sub)
    h_tZtZ_1100.Write("")

if h_tZtH_1100 in templates:

    h_tZtH_1100.SetName("ee"+name+"__TT_tZtH_M1100"+sub)
    h_tZtH_1100.Write("")

if h_tZbW_1100 in templates:

    h_tZbW_1100.SetName("ee"+name+"__TT_tZbW_M1100"+sub)
    h_tZbW_1100.Write("")

if h_tZtZ_1300 in templates:

    h_tZtZ_1300.SetName("ee"+name+"__TT_tZtZ_M1300"+sub)
    h_tZtZ_1300.Write("")

if h_tZtH_1300 in templates:

    h_tZtH_1300.SetName("ee"+name+"__TT_tZtH_M1300"+sub)
    h_tZtH_1300.Write("")

if h_tZbW_1300 in templates:

    h_tZbW_1300.SetName("ee"+name+"__TT_tZbW_M1300"+sub)
    h_tZbW_1300.Write("")

if h_tZtZ_1400 in templates:

    h_tZtZ_1400.SetName("ee"+name+"__TT_tZtZ_M1400"+sub)
    h_tZtZ_1400.Write("")

if h_tZtH_1400 in templates:

    h_tZtH_1400.SetName("ee"+name+"__TT_tZtH_M1400"+sub)
    h_tZtH_1400.Write("")

if h_tZbW_1400 in templates:

    h_tZbW_1400.SetName("ee"+name+"__TT_tZbW_M1400"+sub)
    h_tZbW_1400.Write("")
'''
if h_tZtZ_1500 in templates:

    h_tZtZ_1500.SetName("ee"+name+"__TT_tZtZ_M1500"+sub)
    h_tZtZ_1500.Write("")

if h_tZtH_1500 in templates:

    h_tZtH_1500.SetName("ee"+name+"__TT_tZtH_M1500"+sub)
    h_tZtH_1500.Write("")

if h_tZbW_1500 in templates:

    h_tZbW_1500.SetName("ee"+name+"__TT_tZbW_M1500"+sub)
    h_tZbW_1500.Write("")
'''
if h_tZtZ_1600 in templates:

    h_tZtZ_1600.SetName("ee"+name+"__TT_tZtZ_M1600"+sub)
    h_tZtZ_1600.Write("")

if h_tZtH_1600 in templates:

    h_tZtH_1600.SetName("ee"+name+"__TT_tZtH_M1600"+sub)
    h_tZtH_1600.Write("")

if h_tZbW_1600 in templates:

    h_tZbW_1600.SetName("ee"+name+"__TT_tZbW_M1600"+sub)
    h_tZbW_1600.Write("")

if h_tZtZ_1700 in templates:

    h_tZtZ_1700.SetName("ee"+name+"__TT_tZtZ_M1700"+sub)
    h_tZtZ_1700.Write("")

if h_tZtH_1700 in templates:

    h_tZtH_1700.SetName("ee"+name+"__TT_tZtH_M1700"+sub)
    h_tZtH_1700.Write("")

if h_tZbW_1700 in templates:

    h_tZbW_1700.SetName("ee"+name+"__TT_tZbW_M1700"+sub)
    h_tZbW_1700.Write("")

if h_tZtZ_1800 in templates:

    h_tZtZ_1800.SetName("ee"+name+"__TT_tZtZ_M1800"+sub)
    h_tZtZ_1800.Write("")

if h_tZtH_1800 in templates:

    h_tZtH_1800.SetName("ee"+name+"__TT_tZtH_M1800"+sub)
    h_tZtH_1800.Write("")

if h_tZbW_1800 in templates:

    h_tZbW_1800.SetName("ee"+name+"__TT_tZbW_M1800"+sub)
    h_tZbW_1800.Write("")

'''


f.Close()


#raw_input("hold on")
