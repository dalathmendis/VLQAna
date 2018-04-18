#include "Analysis/VLQAna/interface/JetMakerGen.h"
//#include "Analysis/VLQAna/interface/ApplyJER.h"
#include "Analysis/VLQAna/interface/Utilities.h"

#include <TRandom.h>

#define DEBUGMORE false 
#define DEBUG false  

using namespace std;
using namespace edm ; 

JetMakerGen::JetMakerGen () {} 

JetMakerGen::~JetMakerGen () {} 

JetMakerGen::JetMakerGen (edm::ParameterSet const& iConfig, edm::ConsumesCollector && iC) :
 

  t_jetGenJetPt       ( iC.consumes<vector<float>>(iConfig.getParameter<edm::InputTag> ("jetGenJetPtLabel")      )),
  t_jetGenJetEta      ( iC.consumes<vector<float>>(iConfig.getParameter<edm::InputTag> ("jetGenJetEtaLabel")     )),
  t_jetGenJetPhi      ( iC.consumes<vector<float>>(iConfig.getParameter<edm::InputTag> ("jetGenJetPhiLabel")     )),
  t_jetGenJetE        ( iC.consumes<vector<float>>(iConfig.getParameter<edm::InputTag> ("jetGenJetELabel")       )),
  t_jetGenJetCharge   ( iC.consumes<vector<float>>(iConfig.getParameter<edm::InputTag> ("jetGenJetChargeLabel")  ))
  {
  }

void JetMakerGen::operator()(edm::Event& evt, vlq::JetCollection& jets) {

  Handle <vector<float>>  h_jetGenJetE        ; evt.getByToken (t_jetGenJetE          , h_jetGenJetE        );
  Handle <vector<float>>  h_jetGenJetEta      ; evt.getByToken (t_jetGenJetEta        , h_jetGenJetEta      );
  Handle <vector<float>>  h_jetGenJetPt       ; evt.getByToken (t_jetGenJetPt         , h_jetGenJetPt       );
  Handle <vector<float>>  h_jetGenJetPhi      ; evt.getByToken (t_jetGenJetPhi        , h_jetGenJetPhi      );

  if ((h_jetGenJetPt.product())->size() < 1) return ; 

 
  for (unsigned ijet = 0; ijet < (h_jetGenJetPt.product())->size(); ++ijet) { 

    double jetGenJetPt = (h_jetGenJetPt.product())->at(ijet) ; 
    if (jetGenJetPt == 0.) continue; 
    if (jetGenJetPt < 30.) continue;
    double jetAbsEta = abs((h_jetGenJetEta.product())->at(ijet)) ;
    if (jetAbsEta  > 2.4) continue ; 

    TLorentzVector genJetP4;
    genJetP4.SetPtEtaPhiE((h_jetGenJetPt.product())->at(ijet), 
			  (h_jetGenJetEta.product())->at(ijet),
			  (h_jetGenJetPhi.product())->at(ijet),
			  (h_jetGenJetE.product())->at(ijet));


    //// Jet to put in the jet collection
    vlq::Jet jet ; 
    jet.setIndex        (ijet)  ;
    jet.setP4     (genJetP4);
    jets.push_back(jet) ; 

#if DEBUGMORE
    cout 
      << " \njet pt finally       = " << jet.getPt() 
      << " \njet eta finally       = " << jet.getEta() 
      << " \njet mass finally     = " << jet.getMass() 
      << endl ; 
#endif 
    // cout<< " &&&&&&&&&&&&&&& end of loop &&&&&&&&&&&&&&&&& " << endl;
  } //// loop over jets 

  std::sort(jets.begin(), jets.end(), Utilities::sortByPt<vlq::Jet>) ; 

}
