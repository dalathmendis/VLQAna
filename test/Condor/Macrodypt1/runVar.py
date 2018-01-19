#!/bin/pythons

import subprocess

dir = 'pre1'
if dir == 'pre1':
    options =[
        #['st'],
        #['pt_el1'],
        #['pt_el2'],
        #['dr_elel_pre'],  
        #['pt_zelel_pre'],
        ['cutflow'],
        ['nak4_pre'],

        ['st_pre'],
        ['ht_pre'],
        ['mass_Zelel_pre'],
        ['pt_zelel_pre'],
        ['mass_Zmumu_pre'],
        ['pt_zmumu_pre'],
        ['dr_mumu_pre'],
        ['dr_elel_pre'],

        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],
        
        ['pt_mu1_pre'],
        ['eta_mu1_pre'],
        ['pt_mu2_pre'],
        ['eta_mu2_pre'],


        ['cutflow'],
        ['st_pre'],
        ['ht_pre'],
        ['mass_Zelel_pre'],
        ['pt_zelel_pre'],
        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],
        
     






        #['dr_mumu_pre'],
        #['Iso04_mu_pre'],
        #['Iso04_mu1_pre'],
        #['Iso04_mu2_pre'],
        #['pt_el1_pre'],
        #['eta_el1_pre'],
       # ['dr_el2minjet_pre'],
      #  ['dr_mu2minjet_pre'],


       # ['pt_mu1_pre'],
       # ['eta_mu1_pre'],


       # ['st_pre'],
       # ['nak4_pre'],
       # ['cutflow'],
        ['nak4_pre'], 

      #  ['dr_el1minjet_pre'],
      #  ['dr_el2minjet_pre'],
      #  ['dr_mu1minjet_pre'],
      #  ['dr_mu2minjet_pre'],

      #  ['dr_el1jet1_pre'],
      #  ['dr_el1jet2_pre'],
      #  ['dr_el1jet3_pre'],
      #  ['dr_el2jet1_pre'],
      #  ['dr_el2jet2_pre'],
      #  ['dr_el2jet3_pre'],

      #  ['dr_mu1jet1_pre'],
      #  ['dr_mu1jet2_pre'],
      #  ['dr_mu1jet3_pre'],
      #  ['dr_mu2jet1_pre'],
      #  ['dr_mu2jet2_pre'],
      #  ['dr_mu2jet3_pre'],


      #  ['Iso04_mu_pre'],

        ['cutflow'],
        ['st_pre'],
        ['nak4_pre'],

        ['mass_Zmumu_pre'],
        ['mass_Zelel_pre'],
        ['pt_zelel_pre'],
        ['pt_zmumu_pre'],
      

        ['dr_elel_ex_pre'],
        
        ['massz_ex_pre'],
        ['ptak4jet1_ex_pre'],
        ['etaak4jet1_ex_pre'],
        ['phiak4jet1_ex_pre'],
        ['massak4jet1_ex_pre'],
        ['energyak4jet1_ex_pre'],
        ['ptak4jet2_ex_pre'],
        ['etaak4jet2_ex_pre'],
        ['phiak4jet2_ex_pre'],
        ['massak4jet2_ex_pre'],
        ['energyak4jet2_ex_pre'],
        ['ptak4jet3_ex_pre'],
        ['etaak4jet3_ex_pre'],
        ['phiak4jet3_ex_pre'],
        ['massak4jet3_ex_pre'],
        ['energyak4jet3_ex_pre'],

        ['ptmu1_ex_pre'],
        ['etamu1_ex_pre'],
        ['phimu1_ex_pre'],
        ['energymu1_ex_pre'],
        ['ptmu2_ex_pre'],
        ['etamu2_ex_pre'],
        ['phimu2_ex_pre'],
        ['energymu2_ex_pre'],

        ['dphi_mu1_jet1_pre'],
        ['dphi_mu1_jet2_pre'],
        ['dphi_mu1_jet3_pre'],

        ['dphi_mu2_jet1_pre'],
        ['dphi_mu2_jet2_pre'],
        ['dphi_mu2_jet3_pre'],


        ['ptel1_ex_pre'],
        ['etael1_ex_pre'],
        ['phiel1_ex_pre'],
        ['energyel1_ex_pre'],
        ['ptel2_ex_pre'],
        ['etael2_ex_pre'],
        ['phiel2_ex_pre'],
        ['energyel2_ex_pre'],

        ['dphi_el1_jet1_pre'],
        ['dphi_el1_jet2_pre'],
        ['dphi_el1_jet3_pre'],

        ['dphi_el2_jet1_pre'],
        ['dphi_el2_jet2_pre'],
        ['dphi_el2_jet3_pre'],

        ['Eta_EB_el_ex'],                                                                                                                                      
        ['Eta_EE_el_ex'],                                                                                                      
        ['scEta_EB_el_ex'],                                                                                                   
        ['scEta_EE_el_ex'],                                                                                                   
        ['Iso03_EB_el_ex'],                                                                                              
        ['Iso03_EE_el_ex'],                                                                                     
        ['scEta_EB_el_ex'],                                                                                            
        ['scEta_EE_el_ex'],                                                                                           
        ['dEtaInSeed_EB_el_ex'],                                                                                       
        ['dEtaInSeed_EE_el_ex'],                                                                                        
        ['dPhiIn_EB_el_ex'],                                                                                          
        ['dPhiIn_EE_el_ex'],                                                                                                
        ['Dz_EB_el_ex'],                                                                                                    
        ['Dz_EE_el_ex'],                                                                                                
        ['Dxy_EB_el_ex'],                                                                                              
        ['Dxy_EE_el_ex'],                                                                                             
        ['Full5x5siee_EB_el_ex'],                                                                         
        ['Full5x5siee_EE_el_ex'],                                                                                      
        ['HoE_EB_el_ex'],                                                                                             
        ['HoE_EE_el_ex'],                                                                                                
        ['ooEmooP_EB_el_ex'],                                                                                                
        ['ooEmooP_EE_el_ex'],                                                                                               
        ['missHits_EB_el_ex'],                                                                                               
        ['missHits_EE_el_ex'],                                                                                            
        ['conveto_EB_El_ex'],                                                                                                                             
        ['conveto_EE_el_ex'] ,          





        ['Wptleading_pre'],
        ['Wetaleading_pre'],
        ['Wprunedleading_pre'],
        ['Wpt2nd_pre'],
        ['Weta2nd_pre'],
        ['Wpruned2nd_pre'],
        ['Wpt_pre'],
        ['Weta_pre'],
        ['Wpruned_pre'],

        ['Hptleading_pre'],
        ['Hetaleading_pre'],
        ['Hprunedleading_pre'],
        ['Hpt2nd_pre'],
        ['Heta2nd_pre'],
        ['Hpruned2nd_pre'],
        ['Hpt_pre'],
        ['Heta_pre'],
        ['Hpruned_pre'],

        ['Topptleading_pre'],
        ['Topetaleading_pre'],
        ['Topsoftdropleading_pre'],
        ['Toppt2nd_pre'],
        ['Topeta2nd_pre'],
        ['Topsoftdrop2nd_pre'],
        ['Toppt_pre'],
        ['Topeta_pre'],
        ['Topsoftdrop_pre'],

        ['ptbjetleading_pre'],
        ['etabjetleading_pre'],
      
]

if dir == 'pre':
    options = [
        ['mass_Zmumu_pre'],
        ['dr_mumu_pre'],
        ['pt_zmumu_pre'],

        #['mass_zelel_pre'],
        ['mass_Zelel_pre'],
        ['dr_elel_pre'],
        ['pt_zelel_pre'],

      #  ['ptz_ex_pre'],
       # ['massz_ex_pre'],


        #['pt_zelel_pre'],
        #['pt_zmumu_pre'],
        #['st_prime'],
        #['st_noprime'],
        #['st_prime1'],
        #['st_noprime1'],
        #['st1_cat'],
        #['st_sig'],
        #['npv'],
        #['cutflow'],

        #['massak4jet1_pre'],
        #['massak4jet2_pre'],

       # ['cutflow'],
       # ['Iso04_mu_pre'],
       # ['Dz_mu_pre'],
       # ['Dxy_mu_pre'],
       # ['IsPFMuon_mu_pre'],
       # ['IsGlobalMuon_mu_pre'],
       # ['GlbTrkNormChi2_mu_pre'],
       # ['NumberValidMuonHits_mu_pre'],
       # ['NumberMatchedStations_mu_pre'],
       # ['NumberValidPixelHits_mu_pre'],
       # ['NumberTrackerLayers_mu_pre'],
       # ['mass_zelel_pre'],
       # ['mass_Zelel_pre'],
       # ['dr_elel_pre'],
       # ['pt_zelel_pre'],
       # ['pt_el1_pre'],
       # ['eta_el1_pre'],
       # ['pt_el2_pre'],
       # ['eta_el2_pre'],
       # ['massCorr1_pre'],
       # ['masssmear1_pre'],
       # ['ptsmear1_pre'],

       # ['Wptleading1_pre'],
       # ['Wetaleading1_pre'],
       # ['Wprunedleading1_pre'],
       # ['Wpt2nd1_pre'],
       # ['Weta2nd1_pre'],
       # ['Wpruned2nd1_pre'],

        ['ht_pre'],
        ['st_pre'],
        ['nak4_pre'],
 

        #['nak4_pre'],
        ['ptak4jet1_pre'],
        ['ptak4jet2_pre'],
        ['ptak4jet3_pre'],

        #['ht_pre'],
        #['st_pre'],
        ['met1_pre'],
     #   ['ptbjetleading_pre'],
      #  ['etabjetleading_pre'],
       # ['massak4jet1_pre'],
       # ['massak4jet2_pre'],       
       # ['cutflow'],
       # ['nob_st'],

        ['npv_pre'],
       # ['nak4_pre'],
       # ['ht_pre'],
       # ['st_pre'],
       # ['met_pre'],
       # ['met1_pre'],
       # ['metPhi_pre'],
        #['ptak4jet1_pre'],
        ['etaak4jet1_pre'],
      #  ['cvsak4jet1_pre'],
        #['ptak4jet2_pre'],
        ['etaak4jet2_pre'],
       # ['cvsak4jet2_pre'],
      #  ['ptak4jet3_pre'],
        ['etaak4jet3_pre'],
      #  ['cvsak4jet3_pre'],
       # ['phi_jet1MET_pre'],
       
       # ['mass_zmumu_pre'],
        #['mass_Zmumu_pre'],
        #['dr_mumu_pre'],
        #['pt_zmumu_pre'],
        ['pt_mu1_pre'],
        ['eta_mu1_pre'],
        ['pt_mu2_pre'],
        ['eta_mu2_pre'],

       # ['mass_zelel_pre'],
       # ['mass_Zelel_pre'],
       # ['dr_elel_pre'],
        #['pt_zelel_pre'],
        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],


        ['ptz_ex_pre'],
        ['massz_ex_pre'],
        ['ptak4jet1_ex_pre'],
        ['etaak4jet1_ex_pre'],
        ['phiak4jet1_ex_pre'],
        ['massak4jet1_ex_pre'],
        ['energyak4jet1_ex_pre'],
        ['ptak4jet2_ex_pre'],
        ['etaak4jet2_ex_pre'],
        ['phiak4jet2_ex_pre'],
        ['massak4jet2_ex_pre'],
        ['energyak4jet2_ex_pre'],       
        ['ptak4jet3_ex_pre'],
        ['etaak4jet3_ex_pre'],
        ['phiak4jet3_ex_pre'],
        ['massak4jet3_ex_pre'],
        ['energyak4jet3_ex_pre'],
      
        ['ptmu1_ex_pre'],
        ['etamu1_ex_pre'],
        ['phimu1_ex_pre'],
        ['energymu1_ex_pre'],
        ['ptmu2_ex_pre'],
        ['etamu2_ex_pre'],
        ['phimu2_ex_pre'],
        ['energymu2_ex_pre'],

        ['dphi_mu1_jet1'],
        ['dphi_mu1_jet2'],
        ['dphi_mu1_jet3'],
        
        ['dphi_mu2_jet1'],
        ['dphi_mu2_jet2'],
        ['dphi_mu2_jet3'],


        ['ptel1_ex_pre'],
        ['etael1_ex_pre'],
        ['phiel1_ex_pre'],
        ['energyel1_ex_pre'],
        ['ptel2_ex_pre'],
        ['etael2_ex_pre'],
        ['phiel2_ex_pre'],
        ['energyel2_ex_pre'],

        ['dphi_el1_jet1'],
        ['dphi_el1_jet2'],
        ['dphi_el1_jet3'],

        ['dphi_el2_jet1'],
        ['dphi_el2_jet2'],
        ['dphi_el2_jet3'],


       # ['mass_zelel_pre'],
       # ['mass_Zelel_pre'],
       # ['dr_elel_pre'],
       # ['pt_zelel_pre'],
       # ['pt_el1_pre'],
       # ['eta_el1_pre'],
       # ['pt_el2_pre'],
       # ['eta_el2_pre'],
       
       # ['Eta_EB_el_pre'],
       # ['Eta_EE_el_pre'],
       # ['scEta_EB_el_pre'],
       # ['scEta_EE_el_pre'],
       # ['Iso03_EB_el_pre'],
       # ['Iso03_EE_el_pre'],
       # ['scEta_EB_el_pre'],
       # ['scEta_EE_el_pre'],
       # ['dEtaInSeed_EB_el_pre'],
       # ['dEtaInSeed_EE_el_pre'],
       # ['dPhiIn_EB_el_pre'],
       # ['dPhiIn_EE_el_pre'],
       # ['Dz_EB_el_pre'],
       # ['Dz_EE_el_pre'],
       # ['Dxy_EB_el_pre'],
       # ['Dxy_EE_el_pre'],
       # ['Full5x5siee_EB_el_pre'],
       # ['Full5x5siee_EE_el_pre'],
       # ['HoE_EB_el_pre'],
       # ['HoE_EE_el_pre'],
       # ['ooEmooP_EB_el_pre'],
       # ['ooEmooP_EE_el_pre'],
       # ['missHits_EB_el_pre'],
       # ['missHits_EE_el_pre'],
       # ['conveto_EB_El_pre'],
       # ['conveto_EE_el_pre'],

       # ['Iso04_mu_pre'],
       # ['Dz_mu_pre'],
       # ['Dxy_mu_pre'],
       # ['IsPFMuon_mu_pre'],
       # ['IsGlobalMuon_mu_pre'],
       # ['GlbTrkNormChi2_mu_pre'],
       # ['NumberValidMuonHits_mu_pre'],
       # ['NumberMatchedStations_mu_pre'],
       # ['NumberValidPixelHits_mu_pre'],
       # ['NumberTrackerLayers_mu_pre'],

        ]


elif dir == 'cnt1':
    options =[
        ['nbjets_cnt'],
        ['ht1_cat'],
        ['st1_cat'],
        ['nak41_cat'],


        ['Wptleading_cat'],
        ['Wetaleading_cat'],
        ['Wprunedleading_cat'],
        ['Wpt2nd_cat'],
        ['Weta2nd_cat'],
        ['Wpruned2nd_cat'],
        ['Wpt_cat'],
        ['Weta_cat'],
        ['Wpruned_cat'],
        ['Wpt2nd_cat'],
        ['Weta2nd_cat'],
        ['Wpruned2nd_cat'],



        ['Hptleading_cat'],
        ['Hetaleading_cat'],
        ['Hprunedleading_cat'],
        ['Hpt2nd_cat'],
        ['Heta2nd_cat'],
        ['Hpruned2nd_cat'],
        ['Hpt_cat'],
        ['Heta_cat'],
        ['Hpruned_cat'],

        ['Topptleading_cat'],
        ['Topetaleading_cat'],
        ['Topsoftdropleading_cat'],
        ['Toppt2nd_cat'],
        ['Topeta2nd_cat'],
        ['Topsoftdrop2nd_cat'],
        ['Toppt_cat'],
        ['Topeta_cat'],
        ['Topsoftdrop_cat'],
        ['nak8_pre'],
        ['nwjet_pre'],
        ['nhjet_pre'],
        ['ntjet_pre'],

        ['nak8_cat'],
        ['nwjet_cat'],
        ['nhjet_cat'],
        ['ntjet_cat'],
        ['ht2_cat'],
        ['st2_cat'],
        ['nak42_cat'],

]


elif dir == 'cnt':
    options = [ 
        #['nbjets_cnt'],     
        #['ht1_cat'],
        #['st1_cat'],
        #['nak41_cat'],
        
        #['ht2_cnt'],
        #['st2_cnt'],
        #['nak42_cat'],
        #['nak8_cat'],
        #['nwjet_cat'],
        #['nhjet_cat'],
        #['ntjet_cat'],
        
        #['nak8_pre'],
        #['nwjet_pre'],
        #['nhjet_pre'],
        #['ntjet_pre'],

        #['nak8_cnt'],
        #['nwjet_cnt'],
        #['nhjet_cnt'],
        #['ntjet_cnt'],
       # ['Wptleading1_cnt'],
       # ['Wetaleading1_cnt'],
       # ['Wprunedleading1_cnt'],
       # ['Wpt2nd1_cnt'],
       # ['Weta2nd1_cnt'],
       # ['Wpruned2nd1_cnt'],


      #  ['Wptleading_cnt'],
      #  ['Wetaleading_cnt'],
      #  ['Wprunedleading_cnt'],
      #  ['Wpt2nd_cnt'],
      #  ['Weta2nd_cnt'],
      #  ['Wpruned2nd_cnt'],

      #  ['Hptleading_cnt'],
      #  ['Hetaleading_cnt'],
      #  ['Hprunedleading_cnt'],
      #  ['Hpt2nd_cnt'],
      #  ['Heta2nd_cnt'],
      #  ['Hpruned2nd_cnt'],

      #  ['Topptleading_cnt'],
      #  ['Topetaleading_cnt'],
      #  ['Topsoftdropleading_cnt'],
      #  ['Toppt2nd_cnt'],
      #  ['Topeta2nd_cnt'],
      #  ['Topsoftdrop2nd_cnt'],


       # ['Wptleading_cat'],
      #  ['Wetaleading_cat'],
      #  ['Wprunedleading_cat'],
      #  ['Wpt2nd_cat'],
      #  ['Weta2nd_cat'],
      #  ['Wpruned2nd_cat'],

      #  ['Hptleading_cat'],
      #  ['Hetaleading_cat'],
      #  ['Hprunedleading_cat'],
      #  ['Hpt2nd_cat'],
      #  ['Heta2nd_cat'],
      #  ['Hpruned2nd_cat'],

      #  ['Topptleading_cat'],
      #  ['Topetaleading_cat'],
      #  ['Topsoftdropleading_cat'],
      #  ['Toppt2nd_cat'],
      #  ['Topeta2nd_cat'],
      #  ['Topsoftdrop2nd_cat'],


        ['nob_1000_st'],
        ['nob_1000_ht'],
        ['nob_1000_pt_zelel'],
        ['nob_1000_pt_zmumu'],


        ['mass_Zmumu_cnt'],
        ['mass_Zelel_cnt'],
        ['1b_st'],
        ['nob_st'],
        ['nob_ht'],
        ['1b_ht'],
        ['massak4jet1_cnt'],
        ['massak4jet2_cnt'],
        #['1b_ht'],
        #['1b_st'],



        ['nak4_cnt'],
        ['ptak4jet1_cnt'],
        ['ptak4jet2_cnt'],
        ['ptak4jet3_cnt'],
        ['ht1_cnt'],
        ['st1_cnt'],
        ['nob_1000_st'],
        ['nob_1000_ht'],
        ['nob_1000_zelel'],
        ['nob_1000_zmumu'],
        ['dr_mumu_cnt'],
        ['dr_elel_cnt'],   
        ['pt_zmumu_cnt'],
        ['pt_zelel_cnt'], 


        ['mass_Zmumu_cnt'],
        ['mass_Zelel_cnt'],
        # ['ht1_cnt'],                                                                                                                                                        
        # ['st1_cnt'],  
        ['met1_cnt'],
        #['nob_ht'],
        #['nob_st'],
       # ['1b_ht'],
       # ['1b_st'],
        #  ['nob_1000_st'],
    #   ['nob_1000_ht'],

       # ['massak4jet1_cnt'],
       # ['massak4jet2_cnt'], 
       # ['nob_ht'],
        ['npv_cnt'],
        #['nak4_cnt'],
       # ['ht_cnt'],
       # ['st_cnt'],
        ['met_cnt'],
        #['met1_cnt'],
        ['metPhi_cnt'],
       # ['ptak4jet1_cnt'],
        ['etaak4jet1_cnt'],
        #['cvsak4jet1_cnt'],
       # ['ptak4jet2_cnt'],
        ['etaak4jet2_cnt'],
        #['cvsak4jet2_cnt'],
        #['ptak4jet3_cnt'],
        ['etaak4jet3_cnt'],
        #['cvsak4jet3_cnt'],
        ['phi_jet1MET_cnt'],
       
       # ['mass_zmumu_cnt'],
       # ['mass_zelel_cnt'],
        #['mass_Zmumu_cnt'],
#        ['dr_mumu_cnt'],
#        ['pt_zmumu_cnt'],
        ['pt_mu1_cnt'],
        ['eta_mu1_cnt'],
        ['pt_mu2_cnt'],
        ['eta_mu2_cnt'],
       # ['ht1_cnt'],
       # ['st1_cnt'],
       # ['nob_ht'],
        ['b_pt_zmumu'],        
        ['b_st'],



       # ['dr_elel_cnt'],
      #  ['pt_zelel_cnt'],
        ['pt_el1_cnt'],
        ['eta_el1_cnt'],
        ['pt_el2_cnt'],
        ['eta_el2_cnt'],            
        ['b_pt_elel'],
        ]
elif dir == 'cat5':
    options = [
        ['pt_el1_cat'],
        ['eta_el1_cat'],
        ['pt_el2_cat'],
        ['eta_el2_cat'],
        ['pt_zelel_cat'],
        ['st1_cat'],
        ['ht1_cat'],

       # ['mass_Zelel_cat'],

      #  ['ST_cnt_CatA'],
      #  ['HT_cnt_CatA'],
      #  ['ST_cnt_CatB'],
      #  ['HT_cnt_CatB'],
      #  ['ST_cnt_CatC'],
      #  ['HT_cnt_CatC'],
        
        ['met1_cat'],
      #  ['ptbjetleading_cat'],
      #  ['etabjetleading_cat'],
      #  ['ptbjetsubleading_cat'],
      #  ['etabjetsubleading_cat'],
       # ['st1_cat'],
       # ['st_cat'],
        ['nak4_cat'],
        ['pt_zelel_cat'],
        ['pt_zmumu_cat'],
        #['st1_0cat'],
        #['st1_cat'],
       # ['ht1_cat'],

        ['mass_Zmumu_cat'],
        ['mass_Zelel_cat'],
        ['dr_elel_cat'],
        ['dr_mumu_cat'],
       # ['pt_zelel_cat'],
        ['pt_el1_cat'],
        ['eta_el1_cat'],
        ['pt_el2_cat'],
        ['eta_el2_cat'],
        ['pt_mu1_cat'],
        ['eta_mu1_cat'],
        ['pt_mu2_cat'],
        ['eta_mu2_cat'],
        ['ptak4jet1_cat'],
        ['ptak4jet2_cat'],
        ['ptak4jet3_cat'],
        

        ['nbjets_cat'],
        ['ptbjet_cat'],
        ['etabjet_cat'],        



['dr_elel_cntT1Z1H1b1'],
        ['dr_elel_cntT1Z1H1b2'],
        ['dr_elel_cntT1Z1H0b1'],
        ['dr_elel_cntT1Z1H0b2'],
        ['dr_elel_cntT1Z1H1'],
        ['dr_elel_cntT1Z1H0'],

        ['ht1_0cat'],
        #['st1_0cat'],
        ['nak4_0cat'],
        ['dr_mumu_0cat'],
        ['dr_elel_0cat'],



        ['dr_mumu_cntT1Z1H1b1'],
        ['dr_mumu_cntT1Z1H1b2'],
        ['dr_mumu_cntT1Z1H0b1'],
        ['dr_mumu_cntT1Z1H0b2'],
        ['dr_mumu_cntT1Z1H1'],
        ['dr_mumu_cntT1Z1H0'],




        ['nak4_cat'],
        ['ht_cat'],
        ['st_cat'],
        ['ht1_cat'],
        #['st1_cat'],
        ['met_cat'],
        ['met1_cat'],
        ['metPhi_cat'],
        ['ptak4jet1_cat'],
        ['etaak4jet1_cat'],
        ['cvsak4jet1_cat'],
        ['massak4jet1_cat'],
        ['ptak4jet2_cat'],
        ['etaak4jet2_cat'],
        ['cvsak4jet2_cat'],
        ['massak4jet2_cat'],
        ['ptak4jet3_cat'],
        ['etaak4jet3_cat'],
        ['cvsak4jet3_cat'],
        ['massak4jet3_cat'],



        ['phi_jet1MET_cat'],
        ['mass_Zelel_cat'],
        ['dr_elel_cat'],
        ['pt_zelel_cat'],
        ['pt_el1_cat'],
        ['eta_el1_cat'],
        ['pt_el2_cat'],
        ['eta_el2_cat'],

        #['mass_Zmumu_cat'],
        ['dr_mumu_cat'],
        ['pt_zmumu_cat'],
        ['pt_mu1_cat'],
        ['eta_mu1_cat'],
        ['pt_mu2_cat'],
        ['eta_mu2_cat'],

        ['nbjets_cat'],
        ['ptbjet_cat'],
        ['etabjet_cat'],
        ['ptbjetleading_cat'],
        ['etabjetleading_cat'],
        ['ptbjetsubleading_cat'],
        ['etabjetsubleading_cat'],

]
elif dir == 'sig1':
    options = [
       # ['1b_ht'],
       # ['1b_st'],
        ['nbjets_sig'],
        ['nbjets_cnt'],
        ['nbjets_0btagcnt'],
        ['nbjets_1btagcnt'],
        ['ht_0btagcnt'],
        #['1b_ht'],
        ['ht_1btagcnt'],
        ['lowmet_ht'],
        ['st_0btagcnt'],
        #['1b_st'],
        ['st_1btagcnt'],
        ['lowmet_st'],
        ]

elif dir == 'sig':
    options = [
        
        ['nak4'],
        ['mass_Zmumu'],
        ['mass_Zelel'],
        ['met1'],
        ['pt_zmumu'],
        ['pt_zelel'],
        ['dr_mumu'],
        #['pt_zmumu'],
        ['pt_mu1'],
        ['eta_mu1'],
        ['pt_mu2'],
        ['eta_mu2'],
        ['dr_elel'],
        #['pt_zelel'],
        ['pt_el1'],
        ['eta_el1'],
        ['pt_el2'],
        ['eta_el2'],
        ['nak8'],
        ['nbjets'],
        ['st'],
        


        ['massak4jet1'],
        ['massak4jet2'],
        ['mass_Zmumu'],
        ['mass_Zelel'],
        ['npv'],
        ['nak4'],
        ['ht'],
        ['st'],
        ['met'],
        ['met1'],
        ['metPhi'],
        ['ptak4jet1'],
        ['etaak4jet1'],
        ['ptak4jet2'],
        ['etaak4jet2'],
        ['phi_jet1MET'],

       # ['mass_zmumu'],
       # ['mass_Zmumu'],
        ['dr_mumu'],
        ['pt_zmumu'],
        ['pt_mu1'],
        ['eta_mu1'],
        ['pt_mu2'],
        ['eta_mu2'],
      

       # ['mass_zelel'],
       # ['mass_Zelel'],
        ['dr_elel'],
        ['pt_zelel'],
        ['pt_el1'],
        ['eta_el1'],
        ['pt_el2'],
        ['eta_el2'],


        ['nbjets'],
        ['ptbjetleading'],
        ['etabjetleading'],
        ['nak8'],
        ['nwjet'],
        ['nhjet'],
        ['ntjet'],
        ['ptak8leading'],
        ['etaak8leading'],
        ['mak8leading'],
        ['prunedmak8leading'],
        ['trimmedmak8leading'],
        ['softdropmak8leading'],
        ['ptak82nd'],
        ['etamak82nd'],
        ['mak82nd'],
        ['purnedmak82nd'],
        ['trimmedmak82nd'],
        ['softdropmak82nd'],
       # ['ptTprime'],
       # ['yTprime'],
       # ['mTprime'],
       # ['ptBprime'],
       # ['yBprime'],
       # ['mBprime'],
       # ['ZJetMasslep'],
       # ['chi2_chi'],
       # ['sqrtChi2'],
       # ['chi_mass'],

        ]
elif dir == 'cat1':
    options = [
        #['ST_sig'],
        #['ST_sigT1Z1H1b1'],
        #['ST_sigT1Z1H1b2'],
        #['catC'],
        ['nak4_0b1'],
        ['nak4_0b2'],
        ['nak4_0b3'],

        ['ptak4jet1_0b1'],
        ['ptak4jet2_0b1'],
        ['ptak4jet3_0b1'],
        ['ptak4jet1_0b2'],
        ['ptak4jet2_0b2'],
        ['ptak4jet3_0b2'],
        ['ptak4jet1_0b3'],
        ['ptak4jet2_0b3'],
        ['ptak4jet3_0b3'],
        



       ['ST_cntT1Z1Hprime1b0'],
        ['met_cntT1Z1Hprime1b0'],
        ['ht_cntT1Z1Hprime1b0'],
        ['st_cntT1Z1Hprime1b0'],
        ['ptak4jet1_cntT1Z1Hprime1b0'],
        ['ptak4jet2_cntT1Z1Hprime1b0'],
        ['ptak4jet3_cntT1Z1Hprime1b0'],
        ['pt_zmumu_cntT1Z1Hprime1b0'],
        ['pt_mu1_cntT1Z1Hprime1b0'],
        ['pt_mu2_cntT1Z1Hprime1b0'],
        ['pt_zelel_cntT1Z1Hprime1b0'],
        ['pt_el1_cntT1Z1Hprime1b0'],
        ['pt_el2_cntT1Z1Hprime1b0'],
        ['met_st1000_e0b'],
        ['st_st1000_e0b'],
        ['ht_st1000_e0b'],
        ['ptak4jet1_st1000_e0b'],
        ['ptak4jet2_st1000_e0b'],
        ['ptak4jet3_st1000_e0b'],
        ['pt_zmumu_st1000_e0b'],
        ['pt_mu1_st1000_e0b'],
        ['pt_mu2_st1000_e0b'],
        ['pt_zelel_st1000_e0b'],
        ['pt_el1_st1000_e0b'],
        ['pt_el2_st1000_e0b'],
        ['met_st1000_e1b'],
        ['st_st1000_e1b'],
        ['ht_st1000_e1b'],
        ['ptak4jet1_st1000_e1b'],
        ['ptak4jet2_st1000_e1b'],
        ['ptak4jet3_st1000_e1b'],
        ['pt_zmumu_st1000_e1b'],
        ['pt_mu1_st1000_e1b'],
        ['pt_mu2_st1000_e1b'],
        ['pt_zelel_st1000_e1b'],
        ['pt_el1_st1000_e1b'],
        ['pt_el2_st1000_e1b'],
        ['ht_st1000_1b'],
        ['st_st1000_1b'],
        ['met_st1000_1b'],
        ['ptak4jet1_st1000_1b'],
        ['ptak4jet2_st1000_1b'],
        ['ptak4jet3_st1000_1b'],
        ['pt_zmumu_st1000_1b'],
        ['pt_mu1_st1000_1b'],
        ['pt_mu2_st1000_1b'],
        ['pt_zelel_st1000_1b'],
        ['pt_el1_st1000_1b'],
        ['pt_el2_st1000_1b'],
        ['ht_st1000_2b'],
        ['st_st1000_2b'],
        ['met_st1000_2b'],
        ['ptak4jet1_st1000_2b'],
        ['ptak4jet2_st1000_2b'],
        ['ptak4jet3_st1000_2b'],
        ['pt_zmumu_st1000_2b'],
        ['pt_mu1_st1000_2b'],
        ['pt_mu2_st1000_2b'],
        ['pt_zelel_st1000_2b'],
        ['pt_el1_st1000_2b'],
        ['pt_el2_st1000_2b'],

        ['ptbjetleading_st1000_e1b'],
        ['ptbjetleading_st1000_1b'],
        ['ptbjetleading_st1000_2b'],
        ['ptbjetsubleading_st1000_e1b'],
        ['ptbjetsubleading_st1000_2b'],
        ['HPrimemass-boosted-cnt'],
        ['HPrimePt-boosted-cnt'],
        ['nHPrimecandidate-boosted-cnt'],
        ['HPrimemassnb-cnt'],
        ['HPrimePtnb-cnt'],
        ['nHPrimecandidatesnb-cnt'],
        ['nHPrimecandidates-tot-cnt'],
        ['nHPrimecandidates1-tot-cnt'],
        ]
elif dir == 'cat':
    options = [
        ['cutflow1'],
        ['cutflow2'],

        ['catA_cnt'],
        ['catB_cnt'],
        ['catC_cnt'],
        ['catD_cnt'],        
        ['ST_cntT1Z1H1b1'],
        ['ST_cntT1Z1H1b2'],
        ['ST_cntT1Z1H0b1'],
        ['ST_cntT1Z1H0b2'],
        ['nbjets_cat'],
        ['ptbjet_cat'],
        ['etabjet_cat'],
        ['ptbjetleading_cat'],
        ['etabjetleading_cat'],
        ['ptbjetsubleading_cat'],
        ['etabjetsubleading_cat'],


        # ['cutflow2'],        
       # ['ptak4jet2_cntT1Z1H1b1'],

       # ['cutflow10'],
       # ['cutflow11'],
       # ['cutflow12'],
       # ['cutflow13'],
       # ['cutflow14'],
       # ['cutflow10'],
       # ['cutflow10'],

       # ['ak4pt1_cattest1'], 
       # ['ak4eta1_cattest1'],
       # ['ak4mass1_cattest1'],

       # ['ak4pt2_cattest1'],
       # ['ak4eta2_cattest1'],
       # ['ak4mass2_cattest1'],

        ['ak4pt1_cattest2'],
        ['ak4eta1_cattest2'],
        ['ak4mass1_cattest2'],

        ['ak4pt2_cattest2'],
        ['ak4eta2_cattest2'],
        ['ak4mass2_cattest2'],

        ['dr_Wb_cnt'],
        ['dphi_Wb_cnt'],
        ['mass_ak4matchedak8'],
        ['pt_ak4matchedak8'],
        ['nak4matchedak8'],
        

#['cutflow4'],
       # ['cutflow6'],
        #['ST_sig'],                                                                                                                                                           
        #['ST_cntT1Z1H1b1_A'],                            
        #['ST_cntT1Z1H1b2_A'],
        #['catC_cnt_A'], 
       # ['ST_cntT1Z1H1b1'],
      #  ['ST_cntT1Z1H1b2'],
      #  ['catC_cnt'],

       # ['met_cntT1Z1H1b1'], 
       # ['metPhi_cntT1Z1H1b1'],
       # ['ptak4jet1_cntT1Z1H1b1'],  
       # ['etaak4jet1_cntT1Z1H1b1'],  
       # ['phiak4jet1_cntT1Z1H1b1'],  
        
       # ['ptak4jet2_cntT1Z1H1b1'],
       # ['etaak4jet2_cntT1Z1H1b1'],
       # ['phiak4jet2_cntT1Z1H1b1'],

       # ['ptak4jet3_cntT1Z1H1b1'],
       # ['etaak4jet3_cntT1Z1H1b1'],
       # ['phiak4jet3_cntT1Z1H1b1'],

       # ['pt_zelel_cntT1Z1H1b1'],
       # ['eta_zelel_cntT1Z1H1b1'],
       # ['phi_zelel_cntT1Z1H1b1'],

       # ['pt_el1_cntT1Z1H1b1'],
       # ['eta_el1_cntT1Z1H1b1'],
       # ['phi_el1_cntT1Z1H1b1'],

       # ['pt_el2_cntT1Z1H1b1'],
       # ['eta_el2_cntT1Z1H1b1'],
       # ['phi_el2_cntT1Z1H1b1'],

       # ['ptbjetleading_cntT1Z1H1b1'],  
       # ['etabjetleading_cntT1Z1H1b1'],  
       # ['phibjetleading_cntT1Z1H1b1'],  

        #['ST_sigT1Z1H1b1'],                                                                                                                                                  
        #['ST_sigT1Z1H1b2'],
        #['catC'],
        #['ST_sigT1Z1H0b2_ak8'],  
        #['ST_sigT0Z1H0b2_ak8'],  
        #['ST_sigT1Z0H0b2_ak8'],  
    
       # ['cutflow4'],
          #  ['Hmass-boosted-cnt'],
          #  ['HPt-boosted-cnt'],
          #  ['nHcandidate-boosted-cnt'],
          #  ['Hmassnb-cnt'],
          #  ['HPtnb-cnt'],
          #  ['nHcandidatesnb-cnt'],
          #  ['nHcandidates-tot-cnt'],
          #  ['nHcandidates1-tot-cnt'],

       # ['ST_cntT1Z1H1b1'],  
       # ['ST_cntT1Z1H1b2'], 
       # ['catC_cnt'],
       # ['met_cntT1Z1H1b1'],                                                                    
        # ['met_cntT1Z1H1b2'], 
       #  ['ptak4jet1_cntT1Z1H1b1'],
        # ['ptak4jet1_cntT1Z1H1b2'],
       #  ['ptak4jet2_cntT1Z1H1b1'],
        # ['ptak4jet2_cntT1Z1H1b2'],
       #  ['ptak4jet3_cntT1Z1H1b1'],
        # ['ptak4jet3_cntT1Z1H1b2'],
       #  ['pt_zelel_cntT1Z1H1b1'],
        # ['pt_zelel_cntT1Z1H1b2'],
        # ['pt_zmumu_cntT1Z1H1b1'],
        # ['pt_zmumu_cntT1Z1H1b2'],
        # ['pt_el1_cntT1Z1H1b1'],
        # ['pt_el1_cntT1Z1H1b2'],
        # ['pt_el2_cntT1Z1H1b1'],
        # ['pt_el2_cntT1Z1H1b2'],
        # ['pt_mu1_cntT1Z1H1b1'],
        # ['pt_mu1_cntT1Z1H1b2'],
        # ['pt_mu2_cntT1Z1H1b1'],
        # ['pt_mu2_cntT1Z1H1b2'],
        # ['ptbjetleading_cntT1Z1H1b1'],
        # ['ptbjetleading_cntT1Z1H1b2'],
        # ['ptbjetsubleading_cntT1Z1H1b1'],
        # ['ptbjetsubleading_cntT1Z1H1b2'],

        # ['catC_cnt'],       
        # ['catC_met'],
        # ['catC_ptak4jet1'],
        # ['catC_ptak4jet2'],
        # ['catC_ptak4jet3'],
        # ['catC_pt_zelel'],
        # ['catC_pt_zmumu'],
        # ['catC_pt_el1'],
        # ['catC_pt_el2'],
        # ['catC_pt_mu1'],
        # ['catC_pt_mu2'],
        # ['catC_ptbjetleading'],
        # ['catC_ptbjetsubleading'],

        ['cutflow1'],
        ['cutflow2'],
       # ['cutflow3'],
       # ['cutflow4'],

        #['ST_sig'],       
        # ['ST_sig1b'],
        # ['ST_sig2b'],
        # ['ST_sigT1Z1'],
        # ['ST_sigT0Z1'],
        # ['ST_sigT1Z0'],
        # ['ST_sigT0Z0'],
        

        #  ['ST_sigT1Z1H1'],
        #  ['ST_sigT1Z1H0'],
        #  ['ST_sigT0Z1H1'],
        #  ['ST_sigT0Z1H0'],
        #  ['ST_sigT1Z0H1'],
        #  ['ST_sigT1Z0H0'],
        #  ['ST_sigT0Z0H1'],
        #  ['ST_sigT0Z0H0'],
        
        # ['ST_sigT1Z1H1b1'],
        # ['ST_sigT1Z1H1b2'],
        # ['ST_sigT1Z1H0b1'],
        # ['ST_sigT1Z1H0b2'],
        # ['ST_sigT0Z1H1b1'],
        # ['ST_sigT0Z1H1b2'],
        # ['ST_sigT0Z1H0b1'],
        # ['ST_sigT0Z1H0b2'],
        # ['ST_sigT1Z0H1b1'],
        # ['ST_sigT1Z0H1b2'],
        # ['ST_sigT1Z0H0b1'],
        # ['ST_sigT1Z0H0b2'],
        # ['ST_sigT0Z0H1b1'],
        # ['ST_sigT0Z0H1b2'],
        # ['ST_sigT0Z0H0b1'],
        # ['ST_sigT0Z0H0b2'],
        
        
        # ['ST_sigT1Z1b1'],
        # ['ST_sigT1Z1b2'],
        # ['ST_sigT0Z1b1'],
        # ['ST_sigT0Z1b2'],
        # ['ST_sigT1Z0b1'],
        # ['ST_sigT1Z0b2'],
        # ['ST_sigT0Z0b1'],
        # ['ST_sigT0Z0b2'],
        #  ['ZHmass-boosted'],                                                                                                                                           
        # ['ZHPt-boosted'],                                                                                                                                             
        # ['nZHcandidate-boosted'],                                                                                                                                     
        # ['ZHmassnb'],                                                                                                                                                 
        # ['ZHPtnb'],                                                                                                                                                   
        # ['nZHcandidatesnb'],                                                                                                                                          
        #['nZHcandidates-tot'],  
        # ['nZHcandidates1-tot'],

      #here  
        ['Hmass-boosted-cnt'],
        ['HPt-boosted-cnt'],
        ['nHcandidate-boosted-cnt'],
        ['Hmassnb-cnt'],
        ['HPtnb-cnt'],
        ['nHcandidatesnb-cnt'],
       # ['nHcandidates-tot-cnt'],                                                                                                                                    
       # ['nHcandidates1-tot-cnt']        
        ['Zmass-boosted-cnt'],
        ['ZPt-boosted-cnt'],
        ['nzcandidate-boosted-cnt'],
        
        ['Zmass-cnt'],
        ['ZPt-cnt'],
        ['nzcandidates-cnt'],
       # ['nzcandidates-tot-cnt'],
       # ['nzcandidates1-tot-cnt'], 
   
        ['topmass-D-cnt'],
        ['topPt-D-cnt'],
        ['ntopcandidate-D-cnt'],
        
        
       # ['Wmass-BC-cnt'],
       # ['nWcandidate-BC-cnt'],
       # ['lightjetmass-BC-cnt'],
       # ['nlightjetcandidate-cnt'],
        
        ['topmas-A-cnt'],
        ['topPt-A-cnt'],
        ['ntopcandidate-A-cnt'],
        
        ['topmass-Bc-cnt'],
        ['topPt-BC-cnt'],
        ['ntopcandidate-BC-cnt'],
        
       # ['ntopcandidate-tot-cnt'],
       # ['ntopcandidate1-tot-cnt'],
        
     
        ]


elif dir == 'cat2':
    options =[
]      

elif dir == 'cat3':
    options =[
        ['dr_elel_cntT1Z1H1b1'],
    ['dr_el1jet1_cntT1Z1H1b1'],
    ['dr_el1jet2_cntT1Z1H1b1'],
    ['dr_el1jet3_cntT1Z1H1b1'],
    ['dr_el2jet1_cntT1Z1H1b1'],
    ['dr_el2jet2_cntT1Z1H1b1'],
    ['dr_el2jet3_cntT1Z1H1b1'],

    ['dr_el1bjet1_cntT1Z1H1b1'],
    ['dr_el2bjet1_cntT1Z1H1b1'],
    ['dr_el1Hb1_cntT1Z1H1b1'],
    ['dr_el2Hb1_cntT1Z1H1b1'],
    ['dr_el1Zb1_cntT1Z1H1b1'],
    ['dr_el2Zb1_cntT1Z1H1b1'],
    ['dr_el1tb1_cntT1Z1H1b1'],
    ['dr_el2tb1_cntT1Z1H1b1'],


        ['dr_elel_cntT1Z1H0b1'],
    ['dr_el1jet1_cntT1Z1H0b1'],
    ['dr_el1jet2_cntT1Z1H0b1'],
    ['dr_el1jet3_cntT1Z1H0b1'],
    ['dr_el2jet1_cntT1Z1H0b1'],
    ['dr_el2jet2_cntT1Z1H0b1'],
    ['dr_el2jet3_cntT1Z1H0b1'],

    ['dr_el1bjet1_cntT1Z1H0b1'],
    ['dr_el2bjet1_cntT1Z1H0b1'],
    ['dr_el1Hb1_cntT1Z1H0b1'],
    ['dr_el2Hb1_cntT1Z1H0b1'],
    ['dr_el1Zb1_cntT1Z1H0b1'],
    ['dr_el2Zb1_cntT1Z1H0b1'],
    ['dr_el1tb1_cntT1Z1H0b1'],
    ['dr_el2tb1_cntT1Z1H0b1'],



    ['nbjets_cntT1Z1H1b1'],
    #['nbjets_cntT1Z1H1b2'],
    ['nbjets_cntT1Z1H0b1'],
    #['nbjets_cntT1Z1H0b2'],


    ['ST_cntT1Z1H0b1'],
    ['ptak4jet1_cntT1Z1H0b1'],

    #['ST_cntT1Z1H1b2'],
    #['catC_cnt'],

    ['ST_cntT1Z1H1b1'],
    ['HT_cntT1Z1H1b1'],   
    ['met_cntT1Z1H1b1'],
    ['metPhi_cntT1Z1H1b1'],
    ['ptak4jet1_cntT1Z1H1b1'],
    ['etaak4jet1_cntT1Z1H1b1'],
    ['phiak4jet1_cntT1Z1H1b1'],

    ['ptak4jet2_cntT1Z1H1b1'],
    ['etaak4jet2_cntT1Z1H1b1'],
    ['phiak4jet2_cntT1Z1H1b1'],
    
    ['ptak4jet3_cntT1Z1H1b1'],
    ['etaak4jet3_cntT1Z1H1b1'],
    ['phiak4jet3_cntT1Z1H1b1'],

    ['pt_zelel_cntT1Z1H1b1'],
    ['eta_zelel_cntT1Z1H1b1'],
    ['phi_zelel_cntT1Z1H1b1'],
    
    ['pt_zmumu_cntT1Z1H1b1'],
    ['eta_zmumu_cntT1Z1H1b1'],
    ['phi_zmumu_cntT1Z1H1b1'], 
    
    ['pt_el1_cntT1Z1H1b1'],
    ['eta_el1_cntT1Z1H1b1'],
    ['phi_el1_cntT1Z1H1b1'],
    

    ['pt_mu1_cntT1Z1H1b1'],
    ['eta_mu1_cntT1Z1H1b1'],
    ['phi_mu1_cntT1Z1H1b1'],

    ['pt_el2_cntT1Z1H1b1'],
    ['eta_el2_cntT1Z1H1b1'],
    ['phi_el2_cntT1Z1H1b1'],

    ['pt_mu2_cntT1Z1H1b1'],
    ['eta_mu2_cntT1Z1H1b1'],
    ['phi_mu2_cntT1Z1H1b1'],

    ['ptbjetleading_cntT1Z1H1b1'],
    ['etabjetleading_cntT1Z1H1b1'],
    ['phibjetleading_cntT1Z1H1b1'],

    ['ptbjetsubleading_cntT1Z1H1b1'],
    ['etabjetsubleading_cntT1Z1H1b1'],
    ['phibjetsubleading_cntT1Z1H1b1'],
    ['nbjets_cntT1Z1H1b1'],

    ['dr_elel_cntT1Z1H1b1'],
    ['dr_el1jet1_cntT1Z1H1b1'],
    ['dr_el1jet2_cntT1Z1H1b1'],
    ['dr_el1jet3_cntT1Z1H1b1'],
    ['dr_el2jet1_cntT1Z1H1b1'],
    ['dr_el2jet2_cntT1Z1H1b1'],
    ['dr_el2jet3_cntT1Z1H1b1'],

    ['dr_el1bjet1_cntT1Z1H1b1'],
    ['dr_el2bjet1_cntT1Z1H1b1'],
    ['dr_el1Hb1_cntT1Z1H1b1'],
    ['dr_el2Hb1_cntT1Z1H1b1'],
    ['dr_el1Zb1_cntT1Z1H1b1'],
    ['dr_el2Zb1_cntT1Z1H1b1'],
    ['dr_el1tb1_cntT1Z1H1b1'],
    ['dr_el2tb1_cntT1Z1H1b1'],
  
    # ['ST_cntT1Z1H1b2'],
   # ['HT_cntT1Z1H1b2'],
   # ['met_cntT1Z1H1b2'],
   # ['metPhi_cntT1Z1H1b2'],
   # ['ptak4jet1_cntT1Z1H1b2'],
    # ['etaak4jet1_cntT1Z1H1b2'],
    # ['phiak4jet1_cntT1Z1H1b2'],

    # ['ptak4jet2_cntT1Z1H1b2'],
    # ['etaak4jet2_cntT1Z1H1b2'],
    # ['phiak4jet2_cntT1Z1H1b2'],

    # ['ptak4jet3_cntT1Z1H1b2'],
    # ['etaak4jet3_cntT1Z1H1b2'],
    # ['phiak4jet3_cntT1Z1H1b2'],

    # ['pt_zelel_cntT1Z1H1b2'],
    # ['eta_zelel_cntT1Z1H1b2'],
    # ['phi_zelel_cntT1Z1H1b2'],

    # ['pt_zmumu_cntT1Z1H1b2'],
    # ['eta_zmumu_cntT1Z1H1b2'],
    # ['phi_zmumu_cntT1Z1H1b2'],

    # ['pt_el1_cntT1Z1H1b2'],
    # ['eta_el1_cntT1Z1H1b2'],
    # ['phi_el1_cntT1Z1H1b2'],


    # ['pt_mu1_cntT1Z1H1b2'],
    # ['eta_mu1_cntT1Z1H1b2'],
    # ['phi_mu1_cntT1Z1H1b2'],

    # ['pt_el2_cntT1Z1H1b2'],
    # ['eta_el2_cntT1Z1H1b2'],
    # ['phi_el2_cntT1Z1H1b2'],

    # ['pt_mu2_cntT1Z1H1b2'],
    # ['eta_mu2_cntT1Z1H1b2'],
    # ['phi_mu2_cntT1Z1H1b2'],

    # ['ptbjetleading_cntT1Z1H1b2'],
    # ['etabjetleading_cntT1Z1H1b2'],
    # ['phibjetleading_cntT1Z1H1b2'],

    # ['ptbjetsubleading_cntT1Z1H1b2'],
    # ['etabjetsubleading_cntT1Z1H1b2'],
    # ['phibjetsubleading_cntT1Z1H1b2'],
    # ['nbjets_cntT1Z1H1b2'],



     
    ['HT_cntT1Z1H0b1'],
    ['met_cntT1Z1H0b1'],
    ['metPhi_cntT1Z1H0b1'],
    ['ptak4jet1_cntT1Z1H0b1'],
    ['etaak4jet1_cntT1Z1H0b1'],
    ['phiak4jet1_cntT1Z1H0b1'],

    ['ptak4jet2_cntT1Z1H0b1'],
    ['etaak4jet2_cntT1Z1H0b1'],
    ['phiak4jet2_cntT1Z1H0b1'],

    ['ptak4jet3_cntT1Z1H0b1'],
    ['etaak4jet3_cntT1Z1H0b1'],
    ['phiak4jet3_cntT1Z1H0b1'],

    ['pt_zelel_cntT1Z1H0b1'],
    ['eta_zelel_cntT1Z1H0b1'],
    ['phi_zelel_cntT1Z1H0b1'],

    ['pt_zmumu_cntT1Z1H0b1'],
    ['eta_zmumu_cntT1Z1H0b1'],
    ['phi_zmumu_cntT1Z1H0b1'],

    ['pt_el1_cntT1Z1H0b1'],
    ['eta_el1_cntT1Z1H0b1'],
    ['phi_el1_cntT1Z1H0b1'],


    ['pt_mu1_cntT1Z1H0b1'],
    ['eta_mu1_cntT1Z1H0b1'],
    ['phi_mu1_cntT1Z1H0b1'],

    ['pt_el2_cntT1Z1H0b1'],
    ['eta_el2_cntT1Z1H0b1'],
    ['phi_el2_cntT1Z1H0b1'],

    ['pt_mu2_cntT1Z1H0b1'],
    ['eta_mu2_cntT1Z1H0b1'],
    ['phi_mu2_cntT1Z1H0b1'],

    ['ptbjetleading_cntT1Z1H0b1'],
    ['etabjetleading_cntT1Z1H0b1'],
    ['phibjetleading_cntT1Z1H0b1'],

    ['ptbjetsubleading_cntT1Z1H0b1'],
    ['etabjetsubleading_cntT1Z1H0b1'],
    ['phibjetsubleading_cntT1Z1H0b1'],
    ['nbjets_cntT1Z1H0b1'],


    ['dr_elel_cntT1Z1H0b1'],
    ['dr_el1jet1_cntT1Z1H0b1'],
    ['dr_el1jet2_cntT1Z1H0b1'],
    ['dr_el1jet3_cntT1Z1H0b1'],
    ['dr_el2jet1_cntT1Z1H0b1'],
    ['dr_el2jet2_cntT1Z1H0b1'],
    ['dr_el2jet3_cntT1Z1H0b1'],

    ['dr_el1bjet1_cntT1Z1H0b1'],
    ['dr_el2bjet1_cntT1Z1H0b1'],
    ['dr_el1Hb1_cntT1Z1H0b1'],
    ['dr_el2Hb1_cntT1Z1H0b1'],
    ['dr_el1Zb1_cntT1Z1H0b1'],
    ['dr_el2Zb1_cntT1Z1H0b1'],
    ['dr_el1tb1_cntT1Z1H0b1'],
    ['dr_el2tb1_cntT1Z1H0b1'],




    ['ST_cntT1Z1H0b2'],
    ['HT_cntT1Z1H0b2'],
    ['met_cntT1Z1H0b2'],
    ['metPhi_cntT1Z1H0b2'],
    ['ptak4jet1_cntT1Z1H0b2'],
    ['etaak4jet1_cntT1Z1H0b2'],
    ['phiak4jet1_cntT1Z1H0b2'],

    ['ptak4jet2_cntT1Z1H0b2'],
    ['etaak4jet2_cntT1Z1H0b2'],
    ['phiak4jet2_cntT1Z1H0b2'],

    ['ptak4jet3_cntT1Z1H0b2'],
    ['etaak4jet3_cntT1Z1H0b2'],
    ['phiak4jet3_cntT1Z1H0b2'],

    ['pt_zelel_cntT1Z1H0b2'],
    ['eta_zelel_cntT1Z1H0b2'],
    ['phi_zelel_cntT1Z1H0b2'],

    ['pt_zmumu_cntT1Z1H0b2'],
    ['eta_zmumu_cntT1Z1H0b2'],
    ['phi_zmumu_cntT1Z1H0b2'],

    ['pt_el1_cntT1Z1H0b2'],
    ['eta_el1_cntT1Z1H0b2'],
    ['phi_el1_cntT1Z1H0b2'],


    ['pt_mu1_cntT1Z1H0b2'],
    ['eta_mu1_cntT1Z1H0b2'],
    ['phi_mu1_cntT1Z1H0b2'],

    ['pt_el2_cntT1Z1H0b2'],
    ['eta_el2_cntT1Z1H0b2'],
    ['phi_el2_cntT1Z1H0b2'],

    ['pt_mu2_cntT1Z1H0b2'],
    ['eta_mu2_cntT1Z1H0b2'],
    ['phi_mu2_cntT1Z1H0b2'],

    ['ptbjetleading_cntT1Z1H0b2'],
    ['etabjetleading_cntT1Z1H0b2'],
    ['phibjetleading_cntT1Z1H0b2'],

    ['ptbjetsubleading_cntT1Z1H0b2'],
    ['etabjetsubleading_cntT1Z1H0b2'],
    ['phibjetsubleading_cntT1Z1H0b2'],
    ['nbjets_cntT1Z1H0b2'],





]


elif dir == 'cat4':
    options =[
        ['ST_cnt_CatA'],
        ['HT_cnt_CatA'],
        ['met_cnt_CatA'],
        ['ptak4jet1_cnt_CatA'],
        ['etaak4jet1_cnt_CatA'],
        ['phiak4jet1_cnt_CatA'],
        ['ptak4jet2_cnt_CatA'],
        ['etaak4jet2_cnt_CatA'],
        ['phiak4jet2_cnt_CatA'],
        ['ptak4jet3_cnt_CatA'],
        ['etaak4jet3_cnt_CatA'],
        ['phiak4jet3_cnt_CatA'],

        ['pt_zelel_cnt_CatA'],
        ['eta_zelel_cnt_CatA'],
        ['phi_zelel_cnt_CatA'],
        ['pt_el1_cnt_CatA'],
        ['eta_el1_cnt_CatA'],
        ['phi_el1_cnt_CatA'],
        
        ['pt_el2_cnt_CatA'],
        ['eta_el2_cnt_CatA'],
        ['phi_el2_cnt_CatA'],

        ['pt_zmumu_cnt_CatA'],
        ['eta_zmumu_cnt_CatA'],
        ['phi_zmumu_cnt_CatA'],
        ['pt_mu1_cnt_CatA'],
        ['eta_mu1_cnt_CatA'],
        ['phi_mu1_cnt_CatA'],

        ['pt_mu2_cnt_CatA'],
        ['eta_mu2_cnt_CatA'],
        ['phi_mu2_cnt_CatA'],
       


        ['ST_cnt_CatB'],
        ['HT_cnt_CatB'],
        ['met_cnt_CatB'],
        ['ptak4jet1_cnt_CatB'],
        ['etaak4jet1_cnt_CatB'],
        ['phiak4jet1_cnt_CatB'],
        ['ptak4jet2_cnt_CatB'],
        ['etaak4jet2_cnt_CatB'],
        ['phiak4jet2_cnt_CatB'],
        ['ptak4jet3_cnt_CatB'],
        ['etaak4jet3_cnt_CatB'],
        ['phiak4jet3_cnt_CatB'],

        ['pt_zelel_cnt_CatB'],
        ['eta_zelel_cnt_CatB'],
        ['phi_zelel_cnt_CatB'],
        ['pt_el1_cnt_CatB'],
        ['eta_el1_cnt_CatB'],
        ['phi_el1_cnt_CatB'],

        ['pt_el2_cnt_CatB'],
        ['eta_el2_cnt_CatB'],
        ['phi_el2_cnt_CatB'],

        ['pt_zmumu_cnt_CatB'],
        ['eta_zmumu_cnt_CatB'],
        ['phi_zmumu_cnt_CatB'],
        ['pt_mu1_cnt_CatB'],
        ['eta_mu1_cnt_CatB'],
        ['phi_mu1_cnt_CatB'],

        ['pt_mu2_cnt_CatB'],
        ['eta_mu2_cnt_CatB'],
        ['phi_mu2_cnt_CatB'],



        ['ST_cnt_CatC'],
        ['HT_cnt_CatC'],
        ['met_cnt_CatC'],
        ['ptak4jet1_cnt_CatC'],
        ['etaak4jet1_cnt_CatC'],
        ['phiak4jet1_cnt_CatC'],
        ['ptak4jet2_cnt_CatC'],
        ['etaak4jet2_cnt_CatC'],
        ['phiak4jet2_cnt_CatC'],
        ['ptak4jet3_cnt_CatC'],
        ['etaak4jet3_cnt_CatC'],
        ['phiak4jet3_cnt_CatC'],

        ['pt_zelel_cnt_CatC'],
        ['eta_zelel_cnt_CatC'],
        ['phi_zelel_cnt_CatC'],
        ['pt_el1_cnt_CatC'],
        ['eta_el1_cnt_CatC'],
        ['phi_el1_cnt_CatC'],

        ['pt_el2_cnt_CatC'],
        ['eta_el2_cnt_CatC'],
        ['phi_el2_cnt_CatC'],

        ['pt_zmumu_cnt_CatC'],
        ['eta_zmumu_cnt_CatC'],
        ['phi_zmumu_cnt_CatC'],
        ['pt_mu1_cnt_CatC'],
        ['eta_mu1_cnt_CatC'],
        ['phi_mu1_cnt_CatC'],

        ['pt_mu2_cnt_CatC'],
        ['eta_mu2_cnt_CatC'],
        ['phi_mu2_cnt_CatC'],



        ['ST_cntT1Z1H1'],
        ['HT_cntT1Z1H1'],
        ['met_cntT1Z1H1'],
        ['metPhi_cntT1Z1H1'],
        ['ptak4jet1_cntT1Z1H1'],
        ['etaak4jet1_cntT1Z1H1'],
        ['phiak4jet1_cntT1Z1H1'],

        ['ptak4jet2_cntT1Z1H1'],
        ['etaak4jet2_cntT1Z1H1'],
        ['phiak4jet2_cntT1Z1H1'],

        ['ptak4jet3_cntT1Z1H1'],
        ['etaak4jet3_cntT1Z1H1'],
        ['phiak4jet3_cntT1Z1H1'],

        ['pt_zelel_cntT1Z1H1'],
        ['eta_zelel_cntT1Z1H1'],
        ['phi_zelel_cntT1Z1H1'],

        ['pt_el1_cntT1Z1H1'],
        ['eta_el1_cntT1Z1H1'],
        ['phi_el1_cntT1Z1H1'],

        ['pt_el2_cntT1Z1H1'],
        ['eta_el2_cntT1Z1H1'],
        ['phi_el2_cntT1Z1H1'],

        ['pt_zmumu_cntT1Z1H1'],
        ['eta_zmumu_cntT1Z1H1'],
        ['phi_zmumu_cntT1Z1H1'],

        ['pt_mu1_cntT1Z1H1'],
        ['eta_mu1_cntT1Z1H1'],
        ['phi_mu1_cntT1Z1H1'],

        ['pt_mu2_cntT1Z1H1'],
        ['eta_mu2_cntT1Z1H1'],
        ['phi_mu2_cntT1Z1H1'],



        ['ST_cntT1Z1H0'],
        ['HT_cntT1Z1H0'],
        ['met_cntT1Z1H0'],
        ['metPhi_cntT1Z1H0'],
        ['ptak4jet1_cntT1Z1H0'],
        ['etaak4jet1_cntT1Z1H0'],
        ['phiak4jet1_cntT1Z1H0'],

        ['ptak4jet2_cntT1Z1H0'],
        ['etaak4jet2_cntT1Z1H0'],
        ['phiak4jet2_cntT1Z1H0'],

        ['ptak4jet3_cntT1Z1H0'],
        ['etaak4jet3_cntT1Z1H0'],
        ['phiak4jet3_cntT1Z1H0'],

        ['pt_zelel_cntT1Z1H0'],
        ['eta_zelel_cntT1Z1H0'],
        ['phi_zelel_cntT1Z1H0'],

        ['pt_el1_cntT1Z1H0'],
        ['eta_el1_cntT1Z1H0'],
        ['phi_el1_cntT1Z1H0'],

        ['pt_el2_cntT1Z1H0'],
        ['eta_el2_cntT1Z1H0'],
        ['phi_el2_cntT1Z1H0'],

        ['pt_zmumu_cntT1Z1H0'],
        ['eta_zmumu_cntT1Z1H0'],
        ['phi_zmumu_cntT1Z1H0'],

        ['pt_mu1_cntT1Z1H0'],
        ['eta_mu1_cntT1Z1H0'],
        ['phi_mu1_cntT1Z1H0'],

        ['pt_mu2_cntT1Z1H0'],
        ['eta_mu2_cntT1Z1H0'],
        ['phi_mu2_cntT1Z1H0'],




        ]

elif dir == 'catsig':
    options = [
       # ['ST_sigT1Z0H0b1'],
       # ['ST_sigT0Z1H0b2'],
        ['cutflow4'],
        ['catA_sig'],
        ['catB_sig'],
        ['catC_sig'],
        ['catD_sig'],
       
        ['ST_sig'],
        ['catA_cnt'],
        ['catB_cnt'],
        ['catC_cnt'],
        ['catD_cnt'],
        ['cutflow1'],
        ['cutflow2'],
       
        ['cutflow3'],
        ['cutflow4'],
        ['catA_sig'],
        ['catB_sig'],
        ['catC_sig'],
        ['catD_sig'],   

        ['Hmass-boosted-sig'],
        ['HPt-boosted-sig'],
        ['nHcandidate-boosted-sig'],
        ['Hmassnb-sig'],
        ['HPtnb-sig'],
        ['nHcandidatesnb-sig'],
        ['nHcandidates-tot-sig'],                                                                                                                                       

        ['nHcandidates1-tot-sig'],




        ['Zmass-boosted-sig'],
        ['ZPt-boosted-sig'],
        ['nzcandidate-boosted-sig'],

        ['Zmass-sig'],
        ['ZPt-sig'],
        ['nzcandidates-sig'],
        ['nzcandidates-tot-sig'],
        ['nzcandidates1-tot-sig'],

        ['topmass-D-sig'],
        ['topPt-D-sig'],
        ['ntopcandidate-D-sig'],


        ['Wmass-BC-sig'],
        ['nWcandidate-BC-sig'],
        ['lightjetmass-BC-sig'],
        ['nlightjetcandidate-sig'],

        ['topmas-A-sig'],
        ['topPt-A-sig'],
        ['ntopcandidate-A-sig'],

        ['topmass-Bc-sig'],
        ['topPt-BC-sig'],
        ['ntopcandidate-BC-sig'],

        ['ntopcandidate-tot-sig'],
        ['ntopcandidate1-tot-sig'],

             

        ]

elif dir == 'gen':
    options = [
        ['mass_zelel_pre'],
        ['mass_Zelel_pre'],
        ['dr_elel_pre'],
        ['pt_zelel_pre'],
        ['pt_el1_pre'],
        ['eta_el1_pre'],
        ['pt_el2_pre'],
        ['eta_el2_pre'],

        ['mass_zelel'],
        ['mass_Zelel'],
        ['dr_elel'],
        ['pt_zelel'],
        ['pt_el1'],
        ['eta_el1'],
        ['pt_el2'],
        ['eta_el2'],


        ['htdr_pre'],
        ['stdr_pre'],
        ['drdr_elel_pre'],
        ['htdr_cnt'],
        ['stdr_cnt'],
        ['drdr_elel_cnt'],
        ['htdr'],
        ['stdr'],
        ['drdr_elel'],
        ['st'],
        ]
command = 'python plot.py --var={0:s}'

for option in options :
    s = command.format(
        option[0]
        )

    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo %s"%s,""]                                                                      , shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( ["echo --------------------------------------------------------------------------",""], shell=True)
    subprocess.call( [s, ""]                                                                               , shell=True)
