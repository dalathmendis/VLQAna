#include "DataFormats/ParticleFlowCandidate/interface/PFCandidatePhotonExtra.h"
#include <ostream>
#include <iomanip>

using namespace reco;

PFCandidatePhotonExtra::PFCandidatePhotonExtra(){}

PFCandidatePhotonExtra::PFCandidatePhotonExtra(const reco::SuperClusterRef& scRef) {
  scRef_ = scRef;
}

void PFCandidatePhotonExtra::addSingleLegConvTrackRef(const reco::TrackRef& trackref){
  assoSingleLegRefTrack_.push_back(trackref);
}

void PFCandidatePhotonExtra::addSingleLegConvMva(float& mvasingleleg){
  assoSingleLegMva_.push_back(mvasingleleg);
}

void PFCandidatePhotonExtra::addConversionRef(const reco::ConversionRef& convref){
  assoConversionsRef_.push_back(convref);
}

void  PFCandidatePhotonExtra::addLCorrClusEnergy(float LCorrE){
  LocalCorr_.push_back(LCorrE);
}
