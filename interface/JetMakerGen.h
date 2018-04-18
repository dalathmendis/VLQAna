#ifndef ANALYSIS_VLQANA_JETMAKERGEN_H
#define ANALYSIS_VLQANA_JETMAKERGEN_H

#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include "PhysicsTools/SelectorUtils/interface/Selector.h"
#include "Analysis/VLQAna/interface/JetID.h"
#include "AnalysisDataFormats/BoostedObjects/interface/Jet.h"

#include <boost/algorithm/string.hpp>
#include <string>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TF1.h>

using namespace std;
using namespace edm ; 

class JetMakerGen { 
  public:
    JetMakerGen () ; 
    ~JetMakerGen () ; 
    JetMakerGen (edm::ParameterSet const&, edm::ConsumesCollector && iC) ;
    void operator()(edm::Event& evt, vlq::JetCollection& jetOut) ;
    void operator()(vlq::JetCollection jetsIn, vlq::JetCollection& jetsOut) ;

  private:
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetPt        ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetEta       ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetPhi       ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetE         ; 
    edm::EDGetTokenT<std::vector<float>> t_jetGenJetCharge         ;
};

#endif 
