import FWCore.ParameterSet.Config as cms

from Analysis.VLQAna.ZCandSelector_cfi import *
from Analysis.VLQAna.ApplyLeptonSFs_cfi import *
from Analysis.VLQAna.ElectronMaker_cfi import *
from Analysis.VLQAna.MuonMaker_cfi import *
from Analysis.VLQAna.PickGenPart_cfi import *
from Analysis.VLQAna.JetMaker_cfi import *
from Analysis.VLQAna.METMaker_cfi import *
from Analysis.VLQAna.JetMakerGen_cfi import *
ana = cms.EDFilter("OS2LAna",
    evtno                      = cms.InputTag("evtcleaner","evtno"),   
    runno                      = cms.InputTag("evtcleaner","runno"),   
    lumisec                    = cms.InputTag("evtcleaner","lumisec"),   
    lhewtids                   = cms.InputTag("evtcleaner", "lhewtids"),
    lhewts                     = cms.InputTag("evtcleaner", "lhewts"),
    hltdecision                = cms.InputTag("evtcleaner","hltdecision"),   
    evttype                    = cms.InputTag("evtcleaner","evttype"),
    evtwtGen                   = cms.InputTag("evtcleaner","evtwtGen"),
    evtwtPV                    = cms.InputTag("evtcleaner","evtwtPV"),
    evtwtPVLow                 = cms.InputTag("evtcleaner","evtwtPVHigh"),
    evtwtPVHigh                = cms.InputTag("evtcleaner","evtwtPVLow"),
    npv                        = cms.InputTag("evtcleaner","npv"),
    skim                       = cms.bool(False),
    filterSignal               = cms.bool(False),
    additionalPlots            = cms.bool(False), 
    signalType                 = cms.string(""), 
    zdecayMode                 = cms.string("zelel"),
    optimizeReco               = cms.bool(False),
    categorize                 = cms.bool(False),
    #lhewtids                   = cms.InputTag("evtcleaner", "lhewtids"),
    #lhewts                     = cms.InputTag("evtcleaner", "lhewts"),
    lheId                      = cms.int32(1001),
    pdfId                      = cms.int32(1001),    

    btagsf_bcUp                = cms.bool(False),
    btagsf_bcDown              = cms.bool(False),
    btagsf_lUp                 = cms.bool(False),
    btagsf_lDown               = cms.bool(False),
    PileupUp                   = cms.bool(False),
    PileupDown                 = cms.bool(False),
    ElTrigUp                   = cms.bool(False),
    ElTrigDown                 = cms.bool(False),

    #vlqMass                    = cms.double(1000.),
    #bosMass                    = cms.double(91.2),
    bosonMass                    = cms.double(125.),
    applyLeptonIDSFs           = cms.bool(False), 
    applyLeptonTrigSFs         = cms.bool(False),
    applyBTagSFs               = cms.bool(False),  
    fnamebtagSF                = cms.string('CSVv2_Moriond17_B_H.csv'),
    fnameSJbtagSF              = cms.string('subjet_CSVv2_Moriond17_B_H.csv'),
    btageffmap                 = cms.string('bpbpbZb1200_bTagEff.root'),
    maketree                   = cms.bool(False),
    applyDYNLOCorr             = cms.bool(False),  
    DYDown                     = cms.bool(False),
    tauShift                   = cms.int32(0),
    mistagW                   = cms.int32(0),
    mistagH                   = cms.int32(0),
    mistagT                   = cms.int32(0),

            pdfID_offset               = cms.int32(1),
                   scale_offset = cms.int32(1),
    scaleVal                   = cms.double(1.),
    pdfVal                     = cms.double(1.),
    File_DYNLOCorr             = cms.string('scalefactors_v4.root'),
    Fun_DYNLOCorr              = cms.string('z_ewkcorr/z_ewkcorr_func'),
    #DoPUReweightingNPV         = cms.bool(False),
    DilepCandParams            = defaultZCandSelectionParameters.clone(
        massMin = cms.double(50),
        ptMaxLeadingLep = cms.double(45),#important, otherwise event weight is affected
    ), 
    ZCandParams                = defaultZCandSelectionParameters.clone(
        massMin = cms.double(75),
        massMax = cms.double(105),
        ptMaxLeadingLep = cms.double(45),
        ptMin = cms.double(0.),
        ), 
    #BoostedZCandParams         = defaultZCandSelectionParameters.clone(
    #    massMin = cms.double(75),
    #    massMax = cms.double(105),
    #    ptMaxLeadingLep = cms.double(45),
    #    ptMin = cms.double(200.),
    #    ), 
    #GenHSelParams              = genPartParams.clone(), 
    STMin                      = cms.double(1000.),
    STMaxControl               = cms.double(700.),
    HTMin                      = cms.double(200.),
    NAK4Min                    = cms.uint32 (3),
    lepIdSFsParams             = defaultWP.clone(
        lepidtype = cms.string("TIGHT"),
        zdecayMode = cms.string("zelel"),
        ),
    lepTrigSFsParams             = defaultWP.clone(
        zdecayMode = cms.string("zelel"),
        eltrigeffmap = cms.string("ElectronTriggerSF.root"), 
        ),
    metselParams               = defaultMETMakerParameters,
    muselParams                = defaultMuonMakerParameters, 
    elselParams                = defaultElectronMakerParameters,
    jetAK4GenselParams         = defaultAK4GenJetParameters, 
    jetAK4selParams            = defaultAK4JetSelectionParameters.clone(
        jetAbsEtaMax              = cms.double(2.4),
        ),
    jetAK4BTaggedselParams     = defaultBTaggedAK4JetSelectionParameters,
    jetAK8selParams            = defaultAK8CHSJetSelectionParameters,
    jetHTaggedselParams        = defaultCHSHJetSelectionParameters.clone(
         jetPtMin            = cms.double(300),
         jetPrunedMassMin    = cms.double(105) ,
         jetPrunedMassMax    = cms.double(135) ,
         jettau2Bytau1Min    = cms.double(0.0) ,
         jettau2Bytau1Max    = cms.double(1.0) ,
         #jettau2Bytau1Max    = cms.double(0.6) ,
         subjetCSVMin        = cms.double(CSVv2L) ,
         #subjetCSVMin        = cms.double(-1000000) , 
         jmrShift            = cms.int32(1),
         #  taggingJER          = cms.bool(False),

    ),

    jetWTaggedselParams        = defaultCHSWJetSelectionParameters.clone(
        jetPtMin            = cms.double(200),
        jetPrunedMassMin    = cms.double(65) ,
        jetPrunedMassMax    = cms.double(105) ,
        jettau2Bytau1Min    = cms.double(0.0) ,
        jettau2Bytau1Max    = cms.double(1.0) ,
        #jettau2Bytau1Max    = cms.double(0.6) , 
        jmrShift            = cms.int32(1),
       # taggingJER          = cms.bool(False),
        ),

#    jetwTaggedselParams        = defaultCHSWJetSelectionParameters.clone(
 #       jetPtMin            = cms.double(200),
  #      jetPrunedMassMin    = cms.double(65) ,
   #     jetPrunedMassMax    = cms.double(105) ,
    #    jettau2Bytau1Min    = cms.double(0.0) ,
     #   jettau2Bytau1Max    = cms.double(1.0) ,
       # jmrShift            = cms.int32(0), 

      # taggingJER          = cms.bool(False),                                                                                                                     
     #   ),





    jetTopTaggedselParams      = defaultCHSTJetSelectionParameters.clone(
      jetPtMin            = cms.double(400),
      jettau3Bytau2Min    = cms.double(0.0) ,
      jettau3Bytau2Max    = cms.double(0.81) ,
      #jettau3Bytau2Max    = cms.double(0.67) ,  
      jetSoftDropMassMin  = cms.double(105) ,
      jetSoftDropMassMax  = cms.double(220) ,
      subjetHighestCSVMin = cms.double(-10000) ,
      jmrShift            = cms.int32(0),  
      ), 
    genParams                  = getSelectParams,
    )
