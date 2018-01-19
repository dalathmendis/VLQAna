#! /usr/bin/env python

# =====================================================
#  INPUTS		
# =====================================================
#path = '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Histo8ELa/' 
#path = '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Histo5MUb/'
#path = '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Histo3MUb/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/Histo2MUb2g/'
#path= '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/Histo3ELb2gwMC/'
#path = '/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/miniaod/Histo1EL/'

#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/Histo1MUnody/'
#path ='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/Histo2MUsysnody/'
#path ='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/Histo3ELsysnody/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/Ele32_1/'

#DY Inc
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/Histo2ELb2g/'
#path ='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/DYIncstudues/Histo2MUsysnody1/'                                           
path ='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/DYIncstudues/Histo3ELsysnody2/'

#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/followup20.11.17/munodynobug_2/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/followup_doug/Histo1ELwithphoton/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/followup_doug/Histo1ELnophoton/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/Histo2ELcleaning1/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/newera/followuoMU2/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/Histo1MUnody/'
#Histo3ELAndrewremini/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/batch2/Histo1MU/'
#path='/uscms_data/d3/dmendis/80xB2Gana7/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Reminiaod/runonskims/Histo1MU/'  
#path = '/uscms_data/d3/dmendis/80xB2Gana6/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/Histo14MUe/'
#path = '/uscms_data/d3/dmendis/80xB2Gana6/CMSSW_8_0_20/src/Analysis/VLQAna/test/Condor/cutflow/electrons/'
#pathR = '/uscms_data/d3/dmendis/Rachitha2/CMSSW_7_4_16_patch2/src/Analysis/VLQAna/test/CRAB_0n_Skim/Histo/'
#path = '/Users/skhalil/Desktop/Analysis/OSDL/Macro/Analysis/Histo/'
ch = 'Zelel'
channel = 'single'
variable=''
#f_Data_Oct2015 = TFile(path+'DoubleEG-Run2015D-05Oct2015-v1_os2lana_v8_'+ch+'.root')
if ch == 'Zmumu':
    if channel == 'double':
        f_Data1 = TFile(path+'DoubleMu1-Run2016-PromptReco.root') #electron                                                                                            
        f_Data2 = TFile(path+'DoubleMu2-Run2016-PromptReco.root') #electron                                                                                            
        f_Data3 = TFile(path+'DoubleMu3-Run2016-PromptReco.root') #electron    
    else:
        if variable == '':
            f_Data1 = TFile(path+'SingleMu1-Run2016-PromptReco.root') #electron                                               
            f_Data1p1 = TFile(path+'SingleMu1p1-Run2016-PromptReco.root') #electron 
            f_Data2 = TFile(path+'SingleMu2-Run2016-PromptReco.root') #electron                                               
            f_Data3 = TFile(path+'SingleMu3-Run2016-PromptReco.root') #electron        
            f_Data4 = TFile(path+'SingleMu4-Run2016-PromptReco.root') #electron                                               
            f_Data5 = TFile(path+'SingleMu5-Run2016-PromptReco.root') #elec
            f_Data6 = TFile(path+'SingleMu6-Run2016-PromptReco.root') #electron                                              
            f_Data7 = TFile(path+'SingleMu7-Run2016-PromptReco.root') #electron     
            f_Data7p1 = TFile(path+'SingleMu7p1-Run2016-PromptReco.root')
            f_Data7p2 = TFile(path+'SingleMu7p2-Run2016-PromptReco.root')
        else:
            f_Data1 = TFile(path+'SingleMu1-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data1p1 = TFile(path+'SingleMu1p1-Run2016-PromptRecocutflow.root') #electron                                    
            f_Data2 = TFile(path+'SingleMu2-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data3 = TFile(path+'SingleMu3-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data4 = TFile(path+'SingleMu4-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data5 = TFile(path+'SingleMu5-Run2016-PromptRecocutflow.root') #elec                                            
            f_Data6 = TFile(path+'SingleMu6-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data7 = TFile(path+'SingleMu7-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data7p1 = TFile(path+'SingleMu7p1-Run2016-PromptRecocutflow.root')
            f_Data7p2 = TFile(path+'SingleMu7p2-Run2016-PromptRecocutflow.root')

else:
    if channel == 'double':
        f_Data1 = TFile(path+'DoubleEG1-Run2016-PromptReco.root') #electron                                                                                            
        f_Data1p1 = TFile(path+'DoubleEG1p1-Run2016-PromptReco.root') #electron   
        f_Data2 = TFile(path+'DoubleEG2-Run2016-PromptReco.root') #electron                                                                                             
        f_Data3 = TFile(path+'DoubleEG3-Run2016-PromptReco.root') #electron       
        f_Data4 = TFile(path+'DoubleEG4-Run2016-PromptReco.root') 
        f_Data5 = TFile(path+'DoubleEG5-Run2016-PromptReco.root') 
        f_Data6 = TFile(path+'DoubleEG6-Run2016-PromptReco.root') 
        f_Data7 = TFile(path+'DoubleEG7-Run2016-PromptReco.root') 
        f_Data7p1 = TFile(path+'DoubleEG7p1-Run2016-PromptReco.root') 
        f_Data7p2 = TFile(path+'DoubleEG7p2-Run2016-PromptReco.root') 
    else:
        if variable == '':
            f_Data1 = TFile(path+'SingleEG1-Run2016-PromptReco.root') #electron                                             
            f_Data1p1 = TFile(path+'SingleEG1p1-Run2016-PromptReco.root') #electron                                   
            f_Data2 = TFile(path+'SingleEG2-Run2016-PromptReco.root') #electron                                         
            f_Data3 = TFile(path+'SingleEG3-Run2016-PromptReco.root') #electron                                     
            f_Data4 = TFile(path+'SingleEG4-Run2016-PromptReco.root')
            f_Data5 = TFile(path+'SingleEG5-Run2016-PromptReco.root')
            f_Data6 = TFile(path+'SingleEG6-Run2016-PromptReco.root')
            f_Data7 = TFile(path+'SingleEG7-Run2016-PromptReco.root')
            f_Data7p1 = TFile(path+'SingleEG7p1-Run2016-PromptReco.root')
            f_Data7p2 = TFile(path+'SingleEG7p2-Run2016-PromptReco.root')
        else:
            f_Data1 = TFile(path+'SingleEG1-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data1p1 = TFile(path+'SingleEG1p1-Run2016-PromptRecocutflow.root') #electron                                    
            f_Data2 = TFile(path+'SingleEG2-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data3 = TFile(path+'SingleEG3-Run2016-PromptRecocutflow.root') #electron                                        
            f_Data4 = TFile(path+'SingleEG4-Run2016-PromptRecocutflow.root')
            f_Data5 = TFile(path+'SingleEG5-Run2016-PromptRecocutflow.root')
            f_Data6 = TFile(path+'SingleEG6-Run2016-PromptRecocutflow.root')
            f_Data7 = TFile(path+'SingleEG7-Run2016-PromptRecocutflow.root')
            f_Data7p1 = TFile(path+'SingleEG7p1-Run2016-PromptRecocutflow.root')
            f_Data7p2 = TFile(path+'SingleEG7p2-Run2016-PromptRecocutflow.root')



#f_DY100to200 = TFile(path+'DYJetsToLL_M-50_ht-100to200.root')
#f_DY200to400 = TFile(path+'DYJetsToLL_M-50_ht-200to400.root')
#f_DY400to600 = TFile(path+'DYJetsToLL_M-50_ht-400to600.root')
#f_DY600to800 = TFile(path+'DYJetsToLL_M-50_ht-600to800.root')
#f_DY800to1200 = TFile(path+'DYJetsToLL_M-50_ht-800to1200.root')
#f_DY1200to2500 = TFile(path+'DYJetsToLL_M-50_ht-1200to2500.root')
#f_DY2500toInf = TFile(path+'DYJetsToLL_M-50_ht-2500toInf.root')

if variable == '':
   # f_DYmcnlo    = TFile(path+'DYInc.root')    
    f_DYmcnlo    = TFile(path+'DYJetsToLL_M-50_Inc.root')    
    f_DY100to250 = TFile(path+'DYJetsToLL_M-50_pt-100to2501.root')
    f_DY250to400 = TFile(path+'DYJetsToLL_M-50_pt-250to4001.root')
    f_DY400to650 = TFile(path+'DYJetsToLL_M-50_pt-400to650.root')
    f_DY650toInf = TFile(path+'DYJetsToLL_M-50_pt-650toInf.root')
else:
    f_DY100to250 = TFile(path+'DYJetsToLL_M-50_pt-100to2501cutflow.root')
    f_DY250to400 = TFile(path+'DYJetsToLL_M-50_pt-250to4001cutflow.root')
    f_DY400to650 = TFile(path+'DYJetsToLL_M-50_pt-400to650cutflow.root')
    f_DY650toInf = TFile(path+'DYJetsToLL_M-50_pt-650toInfcutflow.root')



#f_DYmcnlo    = TFile(path+'DY_amcatnlo_os2lana_v2_'+ch+'.root')
#f_DYmad     = TFile(path+'DY_madgraph_os2lana_v4_'+ch+'.root')

#f_WJ100to200 = TFile(path+'WJetsToLNu_HT-100To200_os2lana_v1_'+ch+'.root')
#f_WJ200to400 = TFile(path+'WJetsToLNu_HT-200To400_os2lana_v1_'+ch+'.root')
#f_WJ400to600 = TFile(path+'WJetsToLNu_HT-400To600_os2lana_v1_'+ch+'.root')
#f_WJ600to800 = TFile(path+'WJetsToLNu_HT-600To800_os2lana_v1_'+ch+'.root')
#f_WJ800to1200 = TFile(path+'WJetsToLNu_HT-800To1200_os2lana_v1_'+ch+'.root')
#f_WJ1200to2500 = TFile(path+'WJetsToLNu_HT-1200To2500_os2lana_v1_'+ch+'.root')
#f_WJ2500toInf = TFile(path+'WJetsToLNu_HT-2500ToInf_os2lana_v1_'+ch+'.root')

#f_ST_tW_top     = TFile(path+'ST_tW_5f_top_powheg-pythia8_os2lana_v1_'+ch+'.root')
#f_ST_tW_antitop = TFile(path+'ST_tW_5f_antitop_powheg-pythia8_os2lana_v1_'+ch+'.root')
#f_ST_t          = TFile(path+'ST_t_4f_amcatnlo-pythia8_os2lana_v1_'+ch+'.root')
#f_ST_t_ex1      = TFile(path+'ST_t_4f_amcatnlo-pythia8_ext1_os2lana_v1_'+ch+'.root')
#f_ST_s          = TFile(path+'ST_s_4f_amcatnlo-pythia8_os2lana_v1_'+ch+'.root')

#f_ZZTo2L2Nu     = TFile(path+'ZZTo2L2Nu_powheg_pythia8_os2lana_v1_'+ch+'.root')
#f_WZTo2L2Q      = TFile(path+'WZTo2L2Q_amcatnlo_os2lana_v1_'+ch+'.root')
#f_WWTo2L2Nu     = TFile(path+'WWTo2L2Nu_powheg_os2lana_v1_'+ch+'.root')
if variable == '':
    f_ttbar         = TFile(path+'TT_powheg1.root')
else:
    f_ttbar         = TFile(path+'TT_powheg1cutflow.root')
#f_TpTp_tZtZ_700 = TFile(path+'TpTp_tZtZ_M-700_os2lana_v1_'+ch+'.root')
#f_TpTp_tZbW_700 = TFile(path+'TpTp_tZbW_M-700_os2lana_v1_'+ch+'.root')
#f_TpTp_tZtH_700 = TFile(path+'TpTp_tZtH_M-700_os2lana_v1_'+ch+'.root')



#f_TpTp_tZtZ_800 = TFile(path+'TpTp_tZtZ_M-800.root')
#f_TpTp_tZbW_800 = TFile(path+'TpTp_tZbW_M-800.root')
#f_TpTp_tZtH_800 = TFile(path+'TpTp_tZtH_M-800.root')

#f_TpTp_tZtZ_900 = TFile(path+'TpTp_tZtZ_M-900_os2lana_v1_'+ch+'.root')
#f_TpTp_tZbW_900 = TFile(path+'TpTp_tZbW_M-900_os2lana_v1_'+ch+'.root')
#f_TpTp_tZtH_900 = TFile(path+'TpTp_tZtH_M-900_os2lana_v1_'+ch+'.root')

f_TpTp_tZtZ_1000 = TFile(path+'TpTp_tZtZ_M-1000.root')
f_TpTp_tZbW_1000 = TFile(path+'TpTp_tZbW_M-1000.root')
f_TpTp_tZtH_1000 = TFile(path+'TpTp_tZtH_M-1000.root')

#f_TpTp_tZ_1000 = TFile(path+'TpTp_tZ_M-1000.root')

#f_TpTp_tZtZ_1200 = TFile(path+'TpTp_tZtZ_M-1200.root')
#f_TpTp_tZbW_1200 = TFile(path+'TpTp_tZbW_M-1200.root')
#f_TpTp_tZtH_1200 = TFile(path+'TpTp_tZtH_M-1200.root')

#f_TpTp_tZtZ_900 = TFile(path+'TpTp_tZtZ_M-900.root')
#f_TpTp_tZbW_900 = TFile(path+'TpTp_tZbW_M-900.root')
#f_TpTp_tZtH_900 = TFile(path+'TpTp_tZtH_M-900.root')

#f_TpTp_tZtZ_1100 = TFile(path+'TpTp_tZtZ_M-1100.root')
#f_TpTp_tZbW_1100 = TFile(path+'TpTp_tZbW_M-1100.root')
#f_TpTp_tZtH_1100 = TFile(path+'TpTp_tZtH_M-1100.root')

#f_TpTp_tZtZ_1300 = TFile(path+'TpTp_tZtZ_M-1300.root')
#f_TpTp_tZbW_1300 = TFile(path+'TpTp_tZbW_M-1300.root')
#f_TpTp_tZtH_1300 = TFile(path+'TpTp_tZtH_M-1300.root')

#f_TpTp_tZtZ_1400 = TFile(path+'TpTp_tZtZ_M-1400.root')
#f_TpTp_tZbW_1400 = TFile(path+'TpTp_tZbW_M-1400.root')
#f_TpTp_tZtH_1400 = TFile(path+'TpTp_tZtH_M-1400.root')

#f_TpTp_tZtZ_1500 = TFile(path+'TpTp_tZtZ_M-1500.root')
#f_TpTp_tZbW_1500 = TFile(path+'TpTp_tZbW_M-1500.root')
#f_TpTp_tZtH_1500 = TFile(path+'TpTp_tZtH_M-1500.root')

#f_TpTp_tZtZ_1600 = TFile(path+'TpTp_tZtZ_M-1600.root')
#f_TpTp_tZbW_1600 = TFile(path+'TpTp_tZbW_M-1600.root')
#f_TpTp_tZtH_1600 = TFile(path+'TpTp_tZtH_M-1600.root')

#f_TpTp_tZtZ_1700 = TFile(path+'TpTp_tZtZ_M-1700.root')
#f_TpTp_tZbW_1700 = TFile(path+'TpTp_tZbW_M-1700.root')
#f_TpTp_tZtH_1700 = TFile(path+'TpTp_tZtH_M-1700.root')

#f_TpTp_tZtZ_1800 = TFile(path+'TpTp_tZtZ_M-1800.root')
#f_TpTp_tZbW_1800 = TFile(path+'TpTp_tZbW_M-1800.root')
#f_TpTp_tZtH_1800 = TFile(path+'TpTp_tZtH_M-1800.root')

#Bprime
#f_BpBp_bZbZ_700 = TFile(path+'BpBp_bZbZ_M-700_os2lana_v1_'+ch+'.root')
#f_BpTp_bZtW_700 = TFile(path+'BpBp_bZtW_M-700_os2lana_v1_'+ch+'.root')
#f_BpBp_bZbH_700 = TFile(path+'BpBp_bZbH_M-700_os2lana_v1_'+ch+'.root')

#f_BpBp_bZbZ_800 = TFile(path+'BpBp_bZbZ_M-800_os2lana_v4_'+ch+'.root')
#f_BpBp_bZtW_800 = TFile(path+'BpBp_bZtW_M-800_os2lana_v4_'+ch+'.root')
#f_BpBp_bZbH_800 = TFile(path+'BpBp_bZbH_M-800_os2lana_v4_'+ch+'.root')

#f_BpBp_bZbZ_1000 = TFile(path+'BpBp_bZbZ_M-1000_os2lana_v4_'+ch+'.root')                                                                                                                   
#f_BpBp_bZtW_1000 = TFile(path+'BpBp_bZtW_M-1000_os2lana_v4_'+ch+'.root')                                                                                                                    
#f_BpBp_bZbH_1000 = TFile(path+'BpBp_bZbH_M-1000_os2lana_v4_'+ch+'.root')


#f_BpBp_bZbZ_1200 = TFile(path+'BpBp_bZbZ_M-1200_os2lana_v4_'+ch+'.root')
#f_BpBp_bZtW_1200 = TFile(path+'BpBp_bZtW_M-1200_os2lana_v4_'+ch+'.root')
#f_BpBp_bZbH_1200 = TFile(path+'BpBp_bZbH_M-1200_os2lana_v4_'+ch+'.root')


#Diboson backgrou
if variable == '':
    f_WW = TFile(path+'WW.root')
    f_WZ = TFile(path+'WZ.root')
    f_ZZ = TFile(path+'ZZ.root')
else:
    f_WW = TFile(path+'WWcutflow.root')
    f_WZ = TFile(path+'WZcutflow.root')
    f_ZZ = TFile(path+'ZZcutflow.root')

#f_WZto2 = TFile(path+'WZto2.root')
#f_WZto3 = TFile(path+'WZto3.root')


#f_tant = TFile(path+'tant.root')
#f_tt = TFile(path+'tt.root')
#f_twant = TFile(path+'twant.root')
#f_twt = TFile(path+'twt.root')
#f_s = TFile(path+'s.root')
#===== cross sections (pb)==========

Top_xs            = 831.76  *gSF
#Top_xs            = 831.76  *gSF
#DY100to200_xs     = 147.4   *gSF *1.23
#DY200to400_xs     = 40.99   *gSF *1.23
#DY400to600_xs     = 5.678   *gSF *1.23
#DY600to800_xs     = 1.363   *gSF *1.23
#DY800to1200_xs    = 0.6759  *gSF *1.23
#DY1200to2500_xs   = 0.116   *gSF *1.23
#DY2500toInf_xs    = 0.002592*gSF *1.23


DY100to250_xs     = 83.12   *gSF 
DY250to400_xs     = 3.047   *gSF 
DY400to650_xs     = 0.3921   *gSF 
DY650toInf_xs     = 0.03636   *gSF 


#DY100to250_xs     = 86.524   *gSF
#DY250to400_xs     = 3.3247  *gSF
#DY400to650_xs     = 0.4491   *gSF
#DY650toInf_xs     = 0.0422   *gSF


#DY_xs             = 6104.0   *gSF 
DY_xs             = 6025.2  *gSF
#DY_xs             = 1964.0  *gSF
WJ100to200_xs     = 1345.0  *gSF *1.21 
WJ200to400_xs     = 359.7   *gSF *1.21
WJ400to600_xs     = 48.9    *gSF *1.21
WJ600to800_xs     = 12.05   *gSF *1.21
WJ800to1200_xs    = 5.501   *gSF *1.21
WJ1200to2500_xs   = 1.329   *gSF *1.21
WJ2500toInf_xs    = 0.03216 *gSF *1.21
ST_tW_top_xs      = 35.6    *gSF
ST_tW_antitop_xs  = 35.6    *gSF 
ST_t_xs           = 70.69   *gSF
ST_s_xs           = 3.36    *gSF 
ZZTo2L2Nu_xs      = 0.564   *gSF
WZTo2L2Q_xs       = 3.22    *gSF
WWTo2L2Nu_xs      = 12.178  *gSF
TpTp700_xs        = 0.455   *gSF
TpTp800_xs        = 0.196   *gSF  
TpTp900_xs        = 0.0903  *gSF
TpTp1000_xs       = 0.044   *gSF
#TpTp1000_xs       = 1.0   *gSF  

Tp1000_xs       = 0.4875   * 0.25 

'''
TpTp800_xs        = 1
TpTp1000_xs       = 1
TpTp1200_xs       = 1

TpTp900_xs        = 1
TpTp1100_xs        = 1
TpTp1300_xs        = 1
TpTp1400_xs        = 1
TpTp1500_xs        = 1
TpTp1600_xs        = 1
TpTp1700_xs        = 1
TpTp1800_xs        = 1
'''
TpTp1200_xs       = 0.0118  *gSF
TpTp1500_xs       = 0.00200 *gSF
TpTp1400_xs       = 0.003544  *gSF

BpBp700_xs        = 0.455   *gSF
BpBp800_xs        = 0.196   *gSF
BpBp900_xs        = 0.0903  *gSF
BpBp1000_xs       = 0.044   *gSF
BpBp1200_xs       = 0.0118  *gSF
BpBp1500_xs       = 0.00200 *gSF

#WW_xs             = 63.21 *gSF
#WZ_xs             = 22.82 *gSF
#ZZ_xs             = 10.32 *gSF

WW_xs             = 118.7 *gSF                                                                                                                                                                           
WZ_xs             = 46.74 *gSF                                                                                                                                                                          
ZZ_xs             = 16.91 *gSF   



tant_xs          = 26.38  *gSF 
tt_xs            = 44.33  *gSF 
twant_xs         = 35.85  *gSF 
twt_xs           = 35.85  *gSF 
s_xs             = 3.36   *gSF 


#WW_xs             = 12.178  *gSF
#WZto2_xs          = 5.595   *gSF
#WZto3_xs          = 4.42965 *gSF
#ZZto2_xs          = 0.564   *gSF
#ZZto4_xs          = 1.212   *gSF   

#===== Number of generated events ======
'''
print "f1 num ",f_Data1.Get("allEvents/hEventCount_nowt").GetBinContent(1) 
print "f1p1 num ",f_Data1p1.Get("allEvents/hEventCount_nowt").GetBinContent(1)
print "f2 num ",f_Data2.Get("allEvents/hEventCount_nowt").GetBinContent(1)
print "f3 num ",f_Data3.Get("allEvents/hEventCount_nowt").GetBinContent(1)
print "f4 num ",f_Data4.Get("allEvents/hEventCount_nowt").GetBinContent(1)
print "f5 num ",f_Data5.Get("allEvents/hEventCount_nowt").GetBinContent(1)
print "f6 num ",f_Data6.Get("allEvents/hEventCount_nowt").GetBinContent(1)
print "f7 num ",f_Data7.Get("allEvents/hEventCount_nowt").GetBinContent(1)
print "f7p1 num ",f_Data7p1.Get("allEvents/hEventCount_nowt").GetBinContent(1)
print "f7p2 num ",f_Data7p2.Get("allEvents/hEventCount_nowt").GetBinContent(1)
'''
Top_num          =  155235652.
#Top_num          =  77229341. #correct one
#DY100to200_num   =  8434125. #2725655. (numbers updated for 7_6_x)
#DY200to400_num   =  8683719.
#DY400to600_num   =  8317271.
#DY600to800_num   =  8232153.
#DY800to1200_num  =  2650775.
#DY1200to2500_num =  616612.
#DY2500toInf_num   =  376260.
#DYmcnlo_num      = f_DYmcnlo .Get("allEvents/hEventCount_nowt").GetBinContent(1) *2
#print "DYmcnlo_num  = " ,DYmcnlo_num 
#DYmcnlo_num      =  23632312.0 
DYmcnlo_num      =  28968252.0 *2
#DYmcnlo_num      = 15828104.0 * 1.5
#DYmcnlo_num      = 28861410.
#DYmcnlo_num      = 28892340.
#DYmcnlo_num      =  28285317. * 1.0163566 #28825132. (this 0.997.. factor is coming from ratio of igors new file to original file)
#DYmad_num        =  9042031.

#Top_num          =  f_ttbar.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY100to200_num   = f_DY100to200.Get("allEvents/hEventCount_nowt").GetBinContent(1)                                                                                                               
#DY200to400_num   = f_DY200to400.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY400to600_num   = f_DY400to600.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY600to800_num   = f_DY600to800.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY800to1200_num  = f_DY800to1200.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY1200to2500_num = f_DY1200to2500.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY2500toInf_num  = f_DY2500toInf.Get("allEvents/hEventCount_nowt").GetBinContent(1)


#DY100to250_num   = f_DY100to250.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY250to400_num   = f_DY250to400.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY400to650_num   = f_DY400to650.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#DY650toInf_num   = f_DY650toInf.Get("allEvents/hEventCount_nowt").GetBinContent(1)

DY100to250_num   = 81645578.
DY250to400_num   = 20761069.
#DY100to250_num   = 5942627. #correct one 
#DY250to400_num   = 1185123. #correct one 


#DY400to650_num   = 604038.
#DY650toInf_num   = 1197191.
#DY400to650_num  = 1193880. after readin hadd files
DY400to650_num  = 1125455.  #correct one (without 1 file)
DY650toInf_num   = 1112171. #correct one (without 1 file)
#DY400to650_num  = 1193880
#DY650toInf_num   = 1197191
#print "Top_num = ", Top_num 
 

#DY650toInf_num   = 599665

#print "DY100to250_num", DY100to250_num 
#print "DY250to400_num", DY250to400_num 
#print "DY400to650_num", DY400to650_num
#print "DY650toInf_num", DY650toInf_num


WJ100to200_num   =  10152718.
WJ200to400_num   =  5221599.
WJ400to600_num   =  1745914.
WJ600to800_num   =  4041997.
WJ800to1200_num  =  1574633.
WJ1200to2500_num =  255637.
WJ2500toInf_num  =  253036.
ST_tW_top_num    =  995600.
ST_tW_antitop_num=  988500.
ST_t_num         =  19904330.
ST_t_ex1_num     =  29954054. 
ST_s_num         =  984400.  
ZZTo2L2Nu_num    =  8719200.
WZTo2L2Q_num     =  31394787.
WWTo2L2Nu_num    =  1965200

TpTp700_num      =  1597200./9.
#TpTp800_num      =  815600./9.

#TpTp800tZtZ_num      =  89484.
#TpTp800tZtH_num      =  177866.
#TpTp800tZbW_num      =  167369.

#TpTp900_num      =  832800./9.
#TpTp1000_num     =  822800./9.

#TpTp1000tZtZ_num      =92376.
#TpTp1000tZtH_num      =184918.
#TpTp1000tZbW_num      =185329.

#TpTp1200tZtZ_num      =91936.
#TpTp1200tZtH_num      =184813.
#TpTp1200tZbW_num      =184972.

#TpTp800tZtZ_num      = f_TpTp_tZtZ_800.Get("ana/signalEvts").GetBinContent(1)
#print "tZtZ800 num ", TpTp800tZtZ_num
#TpTp800tZtH_num      = f_TpTp_tZtH_800.Get("ana/signalEvts").GetBinContent(1)
#TpTp800tZbW_num      = f_TpTp_tZbW_800.Get("ana/signalEvts").GetBinContent(1)


#print "sum of all M800 ", TpTp800tZtZ_num+TpTp800tZtH_num+TpTp800tZbW_num
#print "tZtZ", TpTp800tZtZ_num
#print "tZtH", TpTp800tZtH_num
#print "tZbW",TpTp800tZbW_num


#TpTp900_num      =  832800./9.                                                                                                                                           
#TpTp1000_num     =  822800./9.                                                                                                                                           

TpTp1000tZtZ_num      =f_TpTp_tZtZ_1000.Get("ana/signalEvts").GetBinContent(1)
TpTp1000tZtH_num      =f_TpTp_tZtH_1000.Get("ana/signalEvts").GetBinContent(1)
TpTp1000tZbW_num      =f_TpTp_tZbW_1000.Get("ana/signalEvts").GetBinContent(1)

#TpTp1000tZ_num        = f_TpTp_tZ_1000 .Get("allEvents/hEventCount_nowt").GetBinContent(1) 

#print "tZ " , TpTp1000tZ_num    
print "tZtZ", TpTp1000tZtZ_num                                                                                                                                                                    
print "tZtH", TpTp1000tZtH_num                                                                                                                                                                   
print "tZbW",  TpTp1000tZbW_num  
#print "tZbW",1000tZbW_num 

#TpTp1200tZtZ_num      =f_TpTp_tZtZ_1200.Get("ana/signalEvts").GetBinContent(1)
#TpTp1200tZtH_num      =f_TpTp_tZtH_1200.Get("ana/signalEvts").GetBinContent(1)
#TpTp1200tZbW_num      =f_TpTp_tZbW_1200.Get("ana/signalEvts").GetBinContent(1)
'''
TpTp900tZtZ_num      = f_TpTp_tZtZ_900.Get("ana/signalEvts").GetBinContent(1) 
TpTp900tZtH_num      = f_TpTp_tZtH_900.Get("ana/signalEvts").GetBinContent(1)
TpTp900tZbW_num      = f_TpTp_tZbW_900.Get("ana/signalEvts").GetBinContent(1)

TpTp1100tZtZ_num      = f_TpTp_tZtZ_1100.Get("ana/signalEvts").GetBinContent(1)
TpTp1100tZtH_num      = f_TpTp_tZtH_1100.Get("ana/signalEvts").GetBinContent(1)
TpTp1100tZbW_num      = f_TpTp_tZbW_1100.Get("ana/signalEvts").GetBinContent(1)

TpTp1300tZtZ_num      = f_TpTp_tZtZ_1300.Get("ana/signalEvts").GetBinContent(1)
TpTp1300tZtH_num      = f_TpTp_tZtH_1300.Get("ana/signalEvts").GetBinContent(1)
TpTp1300tZbW_num      = f_TpTp_tZbW_1300.Get("ana/signalEvts").GetBinContent(1)
'''
#TpTp1400tZtZ_num      = f_TpTp_tZtZ_1400.Get("ana/signalEvts").GetBinContent(1)
#TpTp1400tZtH_num      = f_TpTp_tZtH_1400.Get("ana/signalEvts").GetBinContent(1)
#TpTp1400tZbW_num      = f_TpTp_tZbW_1400.Get("ana/signalEvts").GetBinContent(1)
'''
TpTp1500tZtZ_num      = f_TpTp_tZtZ_1500.Get("ana/signalEvts").GetBinContent(1)
TpTp1500tZtH_num      = f_TpTp_tZtH_1500.Get("ana/signalEvts").GetBinContent(1)
TpTp1500tZbW_num      = f_TpTp_tZbW_1500.Get("ana/signalEvts").GetBinContent(1)

TpTp1600tZtZ_num      = f_TpTp_tZtZ_1600.Get("ana/signalEvts").GetBinContent(1)
TpTp1600tZtH_num      = f_TpTp_tZtH_1600.Get("ana/signalEvts").GetBinContent(1)
TpTp1600tZbW_num      = f_TpTp_tZbW_1600.Get("ana/signalEvts").GetBinContent(1)

TpTp1700tZtZ_num      = f_TpTp_tZtZ_1700.Get("ana/signalEvts").GetBinContent(1)
TpTp1700tZtH_num      = f_TpTp_tZtH_1700.Get("ana/signalEvts").GetBinContent(1)
TpTp1700tZbW_num      = f_TpTp_tZbW_1700.Get("ana/signalEvts").GetBinContent(1)


TpTp1800tZtZ_num      = f_TpTp_tZtZ_1800.Get("ana/signalEvts").GetBinContent(1)
TpTp1800tZtH_num      = f_TpTp_tZtH_1800.Get("ana/signalEvts").GetBinContent(1)
TpTp1800tZbW_num      = f_TpTp_tZbW_1800.Get("ana/signalEvts").GetBinContent(1)
'''

#TpTp1200_num     =  832800./9.
TpTp1500_num     =  812800./9.

BpBp700_num      =  832200./9.
BpBp800_num      =  818200./9.
BpBp900_num      =  832800./9.
BpBp1000_num     =  822800./9.
BpBp1200_num     =  828600./9.
BpBp1500_num     =  812800./9.


#WW_num           = f_WW.Get("allEvents/hEventCount_nowt").GetBinContent(1) 
#WZ_num           = f_WZ.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#ZZ_num           = f_ZZ.Get("allEvents/hEventCount_nowt").GetBinContent(1)

#WW_num = 6987124.
WZ_num = 2995828. #correct one
#ZZ_num = 1988098.

WW_num = 6876112.  #correct one
ZZ_num = 1893274.  #correct one
#print "WZ_num = " ,WZ_num  
#print "ZZ_num ="   , ZZ_num
#print "WW_num =" , WW_num

#tant_num           = f_tant.Get("allEvents/hEventCount_nowt").GetBinContent(1) 
#tt_num           = f_tt.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#twant_num           = f_twant.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#twt_num           = f_twt.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#s_num           = f_s.Get("allEvents/hEventCount_nowt").GetBinContent(1)
#WW_num           = 1979988.
#WZto2_num        = 25704660.
#WZto3_num        = 2000000.
#ZZto2_num        = 8785050.
#ZZto4_num        = 10710290.

# Legend
leg = TLegend(0.76,0.88,0.94,0.50)
#leg = TLegend(0.8,0.88,0.94,0.50)  
#leg = TLegend(0.12,0.88,0.45,0.70)                                                                                                                              
#leg.SetNColumns(2)            

leg.SetBorderSize(0)
leg.SetFillColor(10)
leg.SetLineColor(10)
leg.SetLineWidth(0)


# =====================================================
#  FUNCTIONS		
# =====================================================

def setTitle(hs,xTitle):
    y = hs.GetYaxis()
    x = hs.GetXaxis()
    y.SetTitle("Events / Bin")
    x.SetTitle(xTitle)
    y.SetLabelSize(0.05)
    y.SetTitleSize(0.07)
    y.SetTitleOffset(0.6)
    y.SetTitleFont(42)
    x.SetTitleSize(0.05)
    x.SetTitleFont(42)

def setTitle1(hs,xTitle):
    y = hs.GetYaxis()
    x = hs.GetXaxis()
    y.SetTitle("Events / Bin")
    x.SetTitle(xTitle)
    y.SetLabelSize(0.05)
    y.SetTitleSize(0.07)
    y.SetTitleOffset(0.6)

   # y.SetTitleSize(0.005)
    #y.SetTitleOffset(0.05)    

    y.SetTitleFont(42)
    x.SetTitleSize(0.05)
    x.SetTitleFont(42)
    x.SetLabelSize(0.05)
    x.SetTitleSize(0.07)


def prepareRatio(h_ratio, h_ratio_bkg, scale, xTitle):
    h_ratio.SetTitle("")
    h_ratio.GetYaxis().SetTitle("Data / Bkg")
    h_ratio.GetXaxis().SetTitle(xTitle)   
    h_ratio.SetMarkerStyle(8) 
    h_ratio.SetMaximum(2)#3 was here
    h_ratio.SetMinimum(0)#-1 was here
    h_ratio.GetYaxis().SetLabelSize(0.06*scale)
    h_ratio.GetYaxis().SetTitleOffset(1.00/scale*0.5)
    h_ratio.GetYaxis().SetTitleSize(0.08*scale)
    h_ratio.GetYaxis().SetTitleFont(42)
    h_ratio.GetXaxis().SetLabelSize(0.06*scale)
    h_ratio.GetXaxis().SetTitleOffset(0.45*scale)
    h_ratio.GetXaxis().SetTitleSize(0.09*scale)
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetXaxis().SetNdivisions(510)
    h_ratio.SetTickLength(0.06,"X")
    h_ratio.SetTickLength(0.05,"Y")

    ## The uncertainty band
    h_ratio_bkg.SetMarkerSize(0)
    h_ratio_bkg.SetFillColor(kGray+1)
    h_ratio_bkg.GetYaxis().SetLabelSize(0.6*scale)
    h_ratio_bkg.GetYaxis().SetTitleOffset(1.00/scale*0.6)
    h_ratio_bkg.GetYaxis().SetTitleSize(0.08*scale)
    h_ratio_bkg.GetYaxis().SetTitleFont(42)
    h_ratio_bkg.GetXaxis().SetLabelSize(0.08*scale)
    h_ratio_bkg.GetXaxis().SetTitleOffset(0.45*scale)
    h_ratio_bkg.GetXaxis().SetTitleSize(0.09*scale)
    h_ratio_bkg.GetYaxis().SetNdivisions(505)
    h_ratio_bkg.GetXaxis().SetNdivisions(510)
    h_ratio_bkg.SetTickLength(0.05,"X")
    h_ratio_bkg.SetTickLength(0.05,"y")
    h_ratio_bkg.SetTitle("")    
    h_ratio_bkg.SetMaximum(2)
    h_ratio_bkg.SetMinimum(0)
    
def overUnderFlow(hist):
    xbins = hist.GetNbinsX()
    hist.SetBinContent(xbins, hist.GetBinContent(xbins)+hist.GetBinContent(xbins+1))
    hist.SetBinContent(1, hist.GetBinContent(0)+hist.GetBinContent(1))
    hist.SetBinError(xbins, TMath.Sqrt(TMath.Power(hist.GetBinError(xbins),2)+TMath.Power(hist.GetBinError(xbins+1),2)))
    hist.SetBinError(1, TMath.Sqrt(TMath.Power(hist.GetBinError(0),2)+TMath.Power(hist.GetBinError(1),2)))
    hist.SetBinContent(xbins+1, 0.)
    hist.SetBinContent(0, 0.)
    hist.SetBinError(xbins+1, 0.)
    hist.SetBinError(0, 0.)
    
def setCosmetics(hist, legname, hname, color):
    hist.Rebin(rebinS)
    hist.SetLineColor(color)
    hist.SetName(hname)
    if 'Data' in hname:
        leg.AddEntry(hist, legname, 'pl')
        hist.SetMarkerStyle(8)
    elif 'tZ' in hname or 'bZ' in hname:          
        hist.SetLineWidth(2)
        leg.AddEntry(hist, legname, 'l')
    else:
        hist.SetFillColor(color)
        leg.AddEntry(hist, legname, 'f')

        
def getHisto( label, leg, dir, var, Samples, color, verbose) :
    histos = []
    for iSample in Samples :
        ifile = iSample[0]
        xs = iSample[1]
        nevt = iSample[2]
        lumi = iSample[3]
        readname = dir+'/'+var
      #  hist  = ifile.Get( readname ).Clone()
        if var =="catA_cnt":
            hist  = ifile.Get(dir+'/ST_cntT1Z1H1b1' ).Clone()
            hist1  = ifile.Get(dir+'/ST_cntT1Z1H0b1' ).Clone()
            hist.Add(hist1,1)
            
        elif var =="catB_cnt":
            hist  = ifile.Get(dir+'/ST_cntT1Z1H1b2' ).Clone()
            hist1  = ifile.Get(dir+'/ST_cntT1Z1H0b2' ).Clone()
            hist.Add(hist1,1)

        elif var =="catC_cnt":
            hist  = ifile.Get(dir+'/ST_cntT1Z0H0b1' ).Clone()
            hist1  = ifile.Get(dir+'/ST_cntT0Z1H0b2' ).Clone()
            hist.Add(hist1,1)
  


        elif var =="catD_cnt":
            hist  = ifile.Get(dir+'/ST_cntT1Z0H1b2' ).Clone()
            hist1  = ifile.Get(dir+'/ST_cntT1Z0H0b2' ).Clone()
            hist2  = ifile.Get(dir+'/ST_cntT0Z1H1b2' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var =="catA_sig":
            hist  = ifile.Get(dir+'/ST_sigT1Z1H1b1' ).Clone()
            hist1  = ifile.Get(dir+'/ST_sigT1Z1H0b1' ).Clone()
            hist.Add(hist1,1)

        elif var =="catB_sig":
            hist  = ifile.Get(dir+'/ST_sigT1Z1H1b2' ).Clone()
            hist1  = ifile.Get(dir+'/ST_sigT1Z1H0b2' ).Clone()
            hist.Add(hist1,1)

        elif var =="catC_sig":
            hist  = ifile.Get(dir+'/ST_sigT1Z0H0b1' ).Clone()
            hist1  = ifile.Get(dir+'/ST_sigT0Z1H0b2' ).Clone()
            hist.Add(hist1,1)



        elif var =="catD_sig":
            hist  = ifile.Get(dir+'/ST_sigT1Z0H1b2' ).Clone()
            hist1  = ifile.Get(dir+'/ST_sigT1Z0H0b2' ).Clone()
            hist2  = ifile.Get(dir+'/ST_sigT0Z1H1b2' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)


        elif var == "ST_cnt_CatA":
            hist  = ifile.Get(dir+'/ST_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/ST_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "HT_cnt_CatA":
            hist  = ifile.Get(dir+'/HT_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/HT_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "met_cnt_CatA":
            hist  = ifile.Get(dir+'/met_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/met_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "ptak4jet1_cnt_CatA":
            hist  = ifile.Get(dir+'/ptak4jet1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "ptak4jet2_cnt_CatA":
            hist  = ifile.Get(dir+'/ptak4jet2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "ptak4jet3_cnt_CatA":
            hist  = ifile.Get(dir+'/ptak4jet3_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet3_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "etaak4jet1_cnt_CatA":
            hist  = ifile.Get(dir+'/etaak4jet1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "etaak4jet2_cnt_CatA":
            hist  = ifile.Get(dir+'/etaak4jet2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "etaak4jet3_cnt_CatA":
            hist  = ifile.Get(dir+'/etaak4jet3_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet3_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "phiak4jet1_cnt_CatA":
            hist  = ifile.Get(dir+'/phiak4jet1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "phiak4jet2_cnt_CatA":
            hist  = ifile.Get(dir+'/phiak4jet2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "phiak4jet3_cnt_CatA":
            hist  = ifile.Get(dir+'/phiak4jet3_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet3_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "pt_zelel_cnt_CatA":
            hist  = ifile.Get(dir+'/pt_zelel_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_zelel_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_zelel_cnt_CatA":
            hist  = ifile.Get(dir+'/eta_zelel_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_zelel_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_zelel_cnt_CatA":
            hist  = ifile.Get(dir+'/phi_zelel_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_zelel_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "dr_elel_cnt_CatA":
            hist  = ifile.Get(dir+'/dr_elel_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/dr_elel_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "pt_el1_cnt_CatA":
            hist  = ifile.Get(dir+'/pt_el1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_el1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_el1_cnt_CatA":
            hist  = ifile.Get(dir+'/eta_el1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_el1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_el1_cnt_CatA":
            hist  = ifile.Get(dir+'/phi_el1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_el1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "pt_el2_cnt_CatA":
            hist  = ifile.Get(dir+'/pt_el2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_el2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_el2_cnt_CatA":
            hist  = ifile.Get(dir+'/eta_el2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_el2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_el2_cnt_CatA":
            hist  = ifile.Get(dir+'/phi_el2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_el2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "pt_zmumu_cnt_CatA":
            hist  = ifile.Get(dir+'/pt_zmumu_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_zmumu_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_zmumu_cnt_CatA":
            hist  = ifile.Get(dir+'/eta_zmumu_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_zmumu_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_zmumu_cnt_CatA":
            hist  = ifile.Get(dir+'/phi_zmumu_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_zmumu_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "dr_mumu_cnt_CatA":
            hist  = ifile.Get(dir+'/dr_mumu_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/dr_mumu_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "pt_mu1_cnt_CatA":
            hist  = ifile.Get(dir+'/pt_mu1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_mu1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_mu1_cnt_CatA":
            hist  = ifile.Get(dir+'/eta_mu1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_mu1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_mu1_cnt_CatA":
            hist  = ifile.Get(dir+'/phi_mu1_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_mu1_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "pt_mu2_cnt_CatA":
            hist  = ifile.Get(dir+'/pt_mu2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_mu2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_mu2_cnt_CatA":
            hist  = ifile.Get(dir+'/eta_mu2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_mu2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_mu2_cnt_CatA":
            hist  = ifile.Get(dir+'/phi_mu2_cntT1Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_mu2_cntT1Z1H0' ).Clone()
            hist.Add(hist1,1)



        elif var == "ST_cnt_CatB":
            hist  = ifile.Get(dir+'/ST_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/ST_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "HT_cnt_CatB":
            hist  = ifile.Get(dir+'/HT_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/HT_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "met_cnt_CatB":
            hist  = ifile.Get(dir+'/met_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/met_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "ptak4jet1_cnt_CatB":
            hist  = ifile.Get(dir+'/ptak4jet1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "ptak4jet2_cnt_CatB":
            hist  = ifile.Get(dir+'/ptak4jet2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "ptak4jet3_cnt_CatB":
            hist  = ifile.Get(dir+'/ptak4jet3_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet3_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "etaak4jet1_cnt_CatB":
            hist  = ifile.Get(dir+'/etaak4jet1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "etaak4jet2_cnt_CatB":
            hist  = ifile.Get(dir+'/etaak4jet2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "etaak4jet3_cnt_CatB":
            hist  = ifile.Get(dir+'/etaak4jet3_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet3_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "phiak4jet1_cnt_CatB":
            hist  = ifile.Get(dir+'/phiak4jet1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "phiak4jet2_cnt_CatB":
            hist  = ifile.Get(dir+'/phiak4jet2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)

        elif var == "phiak4jet3_cnt_CatB":
            hist  = ifile.Get(dir+'/phiak4jet3_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet3_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "pt_zelel_cnt_CatB":
            hist  = ifile.Get(dir+'/pt_zelel_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/pt_zelel_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_zelel_cnt_CatB":
            hist  = ifile.Get(dir+'/eta_zelel_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/eta_zelel_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_zelel_cnt_CatB":
            hist  = ifile.Get(dir+'/phi_zelel_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phi_zelel_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "dr_elel_cnt_CatB":
            hist  = ifile.Get(dir+'/dr_elel_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/dr_elel_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "pt_el1_cnt_CatB":
            hist  = ifile.Get(dir+'/pt_el1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/pt_el1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_el1_cnt_CatB":
            hist  = ifile.Get(dir+'/eta_el1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/eta_el1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_el1_cnt_CatB":
            hist  = ifile.Get(dir+'/phi_el1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phi_el1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "pt_el2_cnt_CatB":
            hist  = ifile.Get(dir+'/pt_el2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/pt_el2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_el2_cnt_CatB":
            hist  = ifile.Get(dir+'/eta_el2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/eta_el2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_el2_cnt_CatB":
            hist  = ifile.Get(dir+'/phi_el2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phi_el2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)


        elif var == "pt_zmumu_cnt_CatB":
            hist  = ifile.Get(dir+'/pt_zmumu_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/pt_zmumu_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_zmumu_cnt_CatB":
            hist  = ifile.Get(dir+'/eta_zmumu_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/eta_zmumu_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_zmumu_cnt_CatB":
            hist  = ifile.Get(dir+'/phi_zmumu_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phi_zmumu_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "dr_mumu_cnt_CatB":
            hist  = ifile.Get(dir+'/dr_mumu_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/dr_mumu_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "pt_mu1_cnt_CatB":
            hist  = ifile.Get(dir+'/pt_mu1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/pt_mu1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_mu1_cnt_CatB":
            hist  = ifile.Get(dir+'/eta_mu1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/eta_mu1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_mu1_cnt_CatB":
            hist  = ifile.Get(dir+'/phi_mu1_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phi_mu1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "pt_mu2_cnt_CatB":
            hist  = ifile.Get(dir+'/pt_mu2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/pt_mu2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "eta_mu2_cnt_CatB":
            hist  = ifile.Get(dir+'/eta_mu2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/eta_mu2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
        elif var == "phi_mu2_cnt_CatB":
            hist  = ifile.Get(dir+'/phi_mu2_cntT0Z1H0' ).Clone()
            hist1  = ifile.Get(dir+'/phi_mu2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)




        elif var == "ST_cnt_CatC":
            hist  = ifile.Get(dir+'/ST_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/ST_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/ST_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "HT_cnt_CatC":
            hist  = ifile.Get(dir+'/HT_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/HT_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/HT_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "met_cnt_CatC":
            hist  = ifile.Get(dir+'/met_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/met_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/met_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "ptak4jet1_cnt_CatC":
            hist  = ifile.Get(dir+'/ptak4jet1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/ptak4jet1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "ptak4jet2_cnt_CatC":
            hist  = ifile.Get(dir+'/ptak4jet2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/ptak4jet2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "ptak4jet3_cnt_CatC":
            hist  = ifile.Get(dir+'/ptak4jet3_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/ptak4jet3_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/ptak4jet3_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "etaak4jet1_cnt_CatC":
            hist  = ifile.Get(dir+'/etaak4jet1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/etaak4jet1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "etaak4jet2_cnt_CatC":
            hist  = ifile.Get(dir+'/etaak4jet2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/etaak4jet2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "etaak4jet3_cnt_CatC":
            hist  = ifile.Get(dir+'/etaak4jet3_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/etaak4jet3_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/etaak4jet3_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phiak4jet1_cnt_CatC":
            hist  = ifile.Get(dir+'/phiak4jet1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phiak4jet1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phiak4jet2_cnt_CatC":
            hist  = ifile.Get(dir+'/phiak4jet2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phiak4jet2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phiak4jet3_cnt_CatC":
            hist  = ifile.Get(dir+'/phiak4jet3_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phiak4jet3_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phiak4jet3_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "pt_zelel_cnt_CatC":
            hist  = ifile.Get(dir+'/pt_zelel_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_zelel_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/pt_zelel_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "eta_zelel_cnt_CatC":
            hist  = ifile.Get(dir+'/eta_zelel_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_zelel_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/eta_zelel_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phi_zelel_cnt_CatC":
            hist  = ifile.Get(dir+'/phi_zelel_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_zelel_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phi_zelel_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)
        elif var == "dr_elel_cnt_CatC":
            hist  = ifile.Get(dir+'/dr_elel_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/dr_elel_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/dr_elel_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "pt_el1_cnt_CatC":
            hist  = ifile.Get(dir+'/pt_el1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_el1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/pt_el1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "eta_el1_cnt_CatC":
            hist  = ifile.Get(dir+'/eta_el1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_el1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/eta_el1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phi_el1_cnt_CatC":
            hist  = ifile.Get(dir+'/phi_el1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_el1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phi_el1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "pt_el2_cnt_CatC":
            hist  = ifile.Get(dir+'/pt_el2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_el2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/pt_el2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "eta_el2_cnt_CatC":
            hist  = ifile.Get(dir+'/eta_el2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_el2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/eta_el2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phi_el2_cnt_CatC":
            hist  = ifile.Get(dir+'/phi_el2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_el2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phi_el2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "pt_zmumu_cnt_CatC":
            hist  = ifile.Get(dir+'/pt_zmumu_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_zmumu_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/pt_zmumu_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "eta_zmumu_cnt_CatC":
            hist  = ifile.Get(dir+'/eta_zmumu_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_zmumu_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/eta_zmumu_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phi_zmumu_cnt_CatC":
            hist  = ifile.Get(dir+'/phi_zmumu_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_zmumu_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phi_zmumu_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "dr_mumu_cnt_CatC":
            hist  = ifile.Get(dir+'/dr_mumu_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/dr_mumu_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/dr_mumu_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "pt_mu1_cnt_CatC":
            hist  = ifile.Get(dir+'/pt_mu1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_mu1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/pt_mu1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "eta_mu1_cnt_CatC":
            hist  = ifile.Get(dir+'/eta_mu1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_mu1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/eta_mu1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phi_mu1_cnt_CatC":
            hist  = ifile.Get(dir+'/phi_mu1_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_mu1_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phi_mu1_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "pt_mu2_cnt_CatC":
            hist  = ifile.Get(dir+'/pt_mu2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/pt_mu2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/pt_mu2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "eta_mu2_cnt_CatC":
            hist  = ifile.Get(dir+'/eta_mu2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/eta_mu2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/eta_mu2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)

        elif var == "phi_mu2_cnt_CatC":
            hist  = ifile.Get(dir+'/phi_mu2_cntT0Z1H1' ).Clone()
            hist1  = ifile.Get(dir+'/phi_mu2_cntT1Z0H1' ).Clone()
            hist2  = ifile.Get(dir+'/phi_mu2_cntT1Z0H0' ).Clone()
            hist.Add(hist1,1)
            hist.Add(hist2,1)











        else:
            hist  = ifile.Get( readname ).Clone()
        if verbose:
            print 'file: {0:<20}, histo:{1:<10}, integral before weight:{2:<3.3f}, nEntries:{3:<3.0f}, weight:{4:<2.3f}'.format(
                ifile.GetName(),    
                hist.GetName(),
                hist.Integral(), hist.GetEntries(), xs * lumi /nevt
                )
        hist.Sumw2()    
        hist.Scale( xs * lumi /nevt)
        #if var == "cutflow":
        #if ifile =='f_Data1' or ifile =='f_Data2'or ifile =='f_Data3':
        #if leg =='Data':
         #   if var == "cutflow":
          #      hist.SetBinContent(10,0)
        #if var == "cutflow":
         #   hist.SetBinContent(1,0)


        #if var == "cutflow2" or var == "cutflow4" :
         #   hist.GetXaxis().SetBinLabel(1, "t1Z1H1b1")
         #   hist.GetXaxis().SetBinLabel(2, "t1Z1H1b2")
         #   hist.GetXaxis().SetBinLabel(3, "t1Z1H0b1")
         #   hist.GetXaxis().SetBinLabel(4, "t1Z1H0b2")
         #   hist.GetXaxis().SetBinLabel(5, "t0Z1H1b1")
         #   hist.GetXaxis().SetBinLabel(6, "t0Z1H1b2")
         #   hist.GetXaxis().SetBinLabel(7, "t0Z1H0b1")
         #   hist.GetXaxis().SetBinLabel(8, "t0Z1H0b2")
         #   hist.GetXaxis().SetBinLabel(9, "t1Z0H1b1")
         #   hist.GetXaxis().SetBinLabel(10, "t1Z0H1b2")
         #   hist.GetXaxis().SetBinLabel(11, "t1Z0H0b1")
         #   hist.GetXaxis().SetBinLabel(12, "t1Z0H0b2")
         #   hist.GetXaxis().SetBinLabel(13, "t0Z0H1b1")
         #   hist.GetXaxis().SetBinLabel(14, "t0Z0H1b2")
         #   hist.GetXaxis().SetBinLabel(15, "t0Z0H0b1")
         #   hist.GetXaxis().SetBinLabel(16, "t0Z0H0b2")
#if var == "mass_Zelel_pre":
       # hist.Rebin(4)
       # hist.GetXaxis().SetRangeUser(75., 105.)
        if var == "ht" or var == "met1" or var == "st":# or var == "1b_ht" or var == "1b_st" or var == "nob_ht" or var == "nob_st" :
            hist.Rebin(2) 
        if var == "mass_Zelel" or var == "mass_Zmumu" :
            hist.Rebin(10)
            hist.GetXaxis().SetRangeUser(75., 105.) 
       # if var == "massak4jet1_pre" or var == "massak4jet2_pre" or "massak4jet1_cnt" or var == "massak4jet2_cnt" or "massak4jet1" or var == "massak4jet2"  :
            # hist.Rebin(10)
        if var == "ht_pre" or var == "met1_pre" or var == "st_pre":
            hist.Rebin(2)
        if var == "mass_Zelel_cnt" or var == "mass_Zmumu_cnt" or var == "ht1_cnt" or var == "met1_cnt" or var == "st1_cnt" :
            hist.Rebin(5)
        if  var == "1b_ht" or var == "1b_st" or var == "nob_ht" or  var == "nob_1000_ht" :                                         
            hist.Rebin(4)
        if var == "nob_st" or var == "nob_1000_st" or var =="ptak4jet3_pre" :
            hist.Rebin(2)
        if var == "ST_sig" or  var == "ST_sigT1Z1H1b1"  or  var == "ST_sigT1Z1H1b2" or  var == "catC" or var == "ST_sigT1Z1H0b2_ak8"  or var == "ST_sigT0Z1H0b2_ak8" or var == "ST_sigT1Z0H0b2_ak8":
            hist.Rebin(4)            
        #if var == "ST_cntT1Z1H1b1"  or  var == "ST_cntT1Z1H1b2" or  var == "catC_cnt":
         #   hist.Rebin(2)

        #if var == "catC_cnt":
         #   hist.Rebin(4)

        if "nak4" in var:
            nX = 7
            xBins = array('d',[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 15.5])
            hist = hist.Rebin(nX, hist.GetName()+"_Rebin", xBins)
        
        if var == "nak41_cat" or var== "nak42_cat":
            nX = 7
            xBins = array('d',[0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 15.5])
            hist = hist.Rebin(nX, hist.GetName()+"_Rebin", xBins)        


        #histos.append( hist )

        if "topmas" in var or "topPt" in var or "Hmass" in var or "HPt" in var or "Zmass" in var or "ZPt" in var or "topmass" in var:
            hist.Rebin(4)
        if "Wpt" in var or "Weta" in var or "Wpruned" in var or "Hpt" in var or "Heta" in var or"Hpruned" in var or "Toppt" in var or "Topeta" in var or"Topsoftdrop" in var:
            hist.Rebin(4)
        if "cntT1Z1H1" in var and not "nbjets" in var and not "dr" in var and not "phi" in var:
            hist.Rebin(5)        
        if "cntT1Z1H0" in var and not "nbjets" in var and not "dr" in var and not "phi" in var:
            hist.Rebin(5)
        if "cntT1Z1H1" in var  and "dr" in var:
            hist.Rebin(4)
        if "cntT1Z1H0" in var and "dr" in var:
            hist.Rebin(4)


 
        #if var == "st1_cat" or var == "ht1_cat" or var == "st2_cat" or var == "ht2_cat" or var == "st1_0cat" or var == "ht1_0cat":
         #  hist.Rebin(2)

        if "met1" in var and "_cat" in var or "eta" in var and "_cat" in var or "cvs" in var and "_cat" in var or "mass_Zelel" in var and "_cat" in var or "mass_Zmumu" in var and "_cat" in var:
            hist.Rebin(4)

        if var == 'catA_cnt' or var == 'catB_cnt' or var == 'catC_cnt' or var == 'catD_cnt':
            hist.Rebin(4)        
        if var == 'catA_sig' or var == 'catB_sig' or var == 'catC_sig' or var == 'catD_sig':
            hist.Rebin(4)
            #hist.Rebin(15)   

        if "_ex" in var and not "phi" in var and not "dr_elel_ex" in var and not "_el_ex" in var:
            hist.Rebin(5)

        if "_el_ex" in var and not "missHits" in var and not "conveto" in var:
            hist.Rebin(2)

        if "_ex" in var and "phi" in var:
            hist.Rebin(1)
        if var =="dr_elel_ex_pre":
            hist.Rebin(1)

        if "ptz_" in var:
            hist.GetXaxis().SetRangeUser(600,1000)


        if "mass_Zelel_pre"in var or "mass_Zmumu_pre"in var :
            hist.Rebin(4)

        if "cnt_CatA" in var and not "phi" in var or "cnt_CatB" in var and not "phi" in var  or "cnt_CatC" in var and not "phi" in var:
            hist.Rebin(5)

        if var == 'pt_zelel_cat' or var == 'pt_zmumu_cat':
            hist.Rebin(2)

        if var == 'ht1_cat' or var == 'st1_cat':
            hist.Rebin(10)

        if var == 'st_cat':
             hist.Rebin(4)


        if 'dEtaInSeed' in var or 'Full5x5siee' in var or 'HoE' in var or 'ooEmooP' in var or 'dPhiIn' in var or 'Iso03' in var or 'Eta_' in var or 'scEta_' in var:
            hist.Rebin(5)

#var == 'mass_Zmumu_cat' or var == 'mass_Zelel_cat':

        #if var == 'pt_zelel_pre':
            #hist.Rebin(2)


        #if var == "pt_zelel_pre" or "pt_zmumu_pre":
         #   nX = 33
          #  xBins = array('d',[0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,700,850,1000])
          #  xBins = array('d',[0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360,380,400,420,440,460,480,500,520,540,560,580,600,680,760,1000])
           # hist = hist.Rebin(nX, hist.GetName()+"_Rebin", xBins)
        #print "***********nbins = ",hist.GetNbinsX()

        histos.append( hist )
        #if "Wpruned" in var or "Hpruned" in var or "Topsoftdrop" in var :
         #   hist.Rebin(5)
        #if var == "ht1_cat" or var == "st1_cat" or var == "ht2_cnt" or var == "st2_cnt" :  
         #   hist.Rebin(5)      
    histo = histos[0]
    setCosmetics(histo, leg, label+var, color) 
    for ihisto in range(1, len(histos) ):
        #print 'ihisto =', ihisto, 'integral', histos[ihisto].Integral(), ', entries', histos[ihisto].GetEntries()
        histo.Add( histos[ihisto] )
        #print 'after addition', histo.Integral()
    if verbose:    
        print 'newName: {0:<5}, Entries:{1:5.2f},  newIntegral: {2:5.2f}'.format(label+var, histo.GetEntries(), histo.Integral() )   
    return histo

