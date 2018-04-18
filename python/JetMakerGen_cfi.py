import FWCore.ParameterSet.Config as cms

defaultAK4GenJetParameters = cms.PSet(
    jetGenJetPtLabel        = cms.InputTag("jetsAK4CHS", "jetAK4CHSGenJetPt"),
    jetGenJetEtaLabel       = cms.InputTag("jetsAK4CHS", "jetAK4CHSGenJetEta"),
    jetGenJetPhiLabel       = cms.InputTag("jetsAK4CHS", "jetAK4CHSGenJetPhi"),
    jetGenJetELabel         = cms.InputTag("jetsAK4CHS", "jetAK4CHSGenJetE"),
    jetGenJetChargeLabel    = cms.InputTag("jetsAK4CHS", "jetAK4CHSGenJetCharge"),
    )

