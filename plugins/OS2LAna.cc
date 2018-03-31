// -*- C++ -*-
//
// Package:    Analysis/VLQAna
// Class:      OS2LAna
// 
/**\class VLQAna OS2LAna.cc Analysis/VLQAna/plugins/OS2LAna.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Devdatta Majumder
//         Created:  Fri, 27 Feb 2015 16:09:10 GMT
// Modified: Sadia Khalil
//           25 Mar 2016 17:11 CDT
//

#include <iostream>
#include <memory>
#include <vector>
#include <fstream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "AnalysisDataFormats/BoostedObjects/interface/GenParticleWithDaughters.h"
#include "AnalysisDataFormats/BoostedObjects/interface/ResolvedVjj.h"

#include "Analysis/VLQAna/interface/Utilities.h"
#include "Analysis/VLQAna/interface/DileptonCandsProducer.h"
#include "Analysis/VLQAna/interface/CandidateFilter.h"
#include "Analysis/VLQAna/interface/MuonMaker.h"
#include "Analysis/VLQAna/interface/ElectronMaker.h"
#include "Analysis/VLQAna/interface/JetMaker.h"
#include "Analysis/VLQAna/interface/JetMakerGen.h"
#include "Analysis/VLQAna/interface/HT.h"
#include "Analysis/VLQAna/interface/ApplyLeptonIDSFs.h"
#include "Analysis/VLQAna/interface/CandidateCleaner.h"
#include "Analysis/VLQAna/interface/METMaker.h"
#include "Analysis/VLQAna/interface/PickGenPart.h"
#include "Analysis/VLQAna/interface/JetID.h"
#include "Analysis/VLQAna/interface/BTagSFUtils.h"
#include "Analysis/VLQAna/interface/ApplyLeptonTrigSFs.h"
#include "Analysis/VLQAna/interface/DYNLOEwkKfact.h"
#include "Analysis/VLQAna/interface/OS2LTree.h"


#include "Analysis/VLQAna/interface/TopCandsProducer.h"
#include "Analysis/VLQAna/interface/ZCandsProducer.h"
#include "Analysis/VLQAna/interface/HCandsProducer.h"
//#include "Analysis/VLQAna/interface/HPrimeCandsProducer.h"
//#include "Analysis/VLQAna/interface/ZHCandsProducer.h"

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

//
// class declaration
//

class OS2LAna : public edm::EDFilter {
  public:
    explicit OS2LAna(const edm::ParameterSet&);
    ~OS2LAna();

  private:
    virtual void beginJob() override;
    virtual bool filter(edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;
    void fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt);
  TLorentzVector getGen(GenParticleCollection, int, int) ;
  TLorentzVector getGen(GenParticleCollection, int, int, int) ;
  TLorentzVector getMatchedJet(TLorentzVector, vlq::JetCollection, double) ;
  double findInvMass(TLorentzVector, TLorentzVector) ;
  double findInvMass(TLorentzVector, TLorentzVector, TLorentzVector) ;
  double findInvMass(TLorentzVector) ;
  double findPt(TLorentzVector) ;
  double findPt(TLorentzVector, TLorentzVector) ;
  double findPt(TLorentzVector, TLorentzVector, TLorentzVector) ;
  bool pairCompare(const std::pair<vlq::JetCollection, double > , const std::pair<vlq::JetCollection, double >);
  bool my_compare (const std::pair<vlq::Jet, double > , const std::pair<vlq::Jet, double >);
 

  // ----------member data ---------------------------
    edm::EDGetTokenT<string>   t_evttype         ;
    edm::EDGetTokenT<double>   t_evtwtGen        ;
    edm::EDGetTokenT<double>   t_evtwtPV         ;
    edm::EDGetTokenT<double>   t_evtwtPVLow      ;
    edm::EDGetTokenT<double>   t_evtwtPVHigh     ;
    edm::EDGetTokenT<unsigned> t_npv             ;
    edm::EDGetTokenT<bool>     t_hltdecision     ;
    edm::EDGetTokenT<int>      t_evtno           ; 
    edm::EDGetTokenT<int>      t_runno           ;
    edm::EDGetTokenT<int>      t_lumisec         ;
  edm::EDGetTokenT<vector<int> > t_lhewtids    ;
  edm::EDGetTokenT<vector<double> > t_lhewts   ;

    edm::ParameterSet DilepCandParams_           ; 
    edm::ParameterSet ZCandParams_               ; 
    //edm::ParameterSet BoostedZCandParams_        ; 
    //edm::ParameterSet GenHSelParams_             ;
    edm::ParameterSet genParams_                 ;
    const unsigned int NAK4Min_                  ;
    const double HTMin_                          ;
    const double STMin_                          ; 
    const double STMaxControl_                   ;
    const bool skim_                             ;
    const bool filterSignal_                     ;
    const bool additionalPlots_                  ;
    const std::string signalType_                ;
    const std::string zdecayMode_                ;
    const bool btagsf_bcUp_                      ;
    const bool btagsf_bcDown_                    ;
    const bool btagsf_lUp_                       ;
    const bool btagsf_lDown_                     ;
    const bool PileupUp_                         ;
    const bool PileupDown_                       ;   
    const bool ElTrigUp_                         ;// Electron Trig unc
    const bool ElTrigDown_                       ;//Electron Trig Unc

    const int lheId_                             ;
    const int  pdfId_                            ;  
  const bool categorize_                       ; 

  const bool applyLeptonIDSFs_                   ;
    const bool applyLeptonTrigSFs_               ;
    const bool applyBTagSFs_                     ;
    const bool applyDYNLOCorr_                   ;
  const bool DYDown_                             ;
  const int  tauShift_                           ;
  const int pdfID_offset_                      ;
  const int scale_offset_ ;

  const double scaleVal_                       ;
  const double pdfVal_ ;  
  const std::string fname_DYNLOCorr_           ; 
    const std::string funname_DYNLOCorr_         ; 
    const bool optimizeReco_                        ; //comment them later                                                                                                                
    const double bosonMass_                      ;//comment them later       

    DYNLOEwkKfact dynloewkkfact                  ;
    ApplyLeptonIDSFs lepIdSFs                    ;
    ApplyLeptonTrigSFs lepTrigSFs                ;
    METMaker metmaker                            ;
    MuonMaker muonmaker                          ; 
    ElectronMaker electronmaker                  ; 
    JetMakerGen jetAK4Genmaker                   ;
    JetMaker jetAK4maker                         ; 
    JetMaker jetAK4BTaggedmaker                  ; 
    JetMaker jetAK8maker                         ; 
    JetMaker jetHTaggedmaker                     ; 
    JetMaker jetWTaggedmaker                     ; 
  // JetMaker jetwTaggedmaker                     ;
    JetMaker jetTopTaggedmaker                   ; 
    edm::Service<TFileService> fs                ; 
    std::map<std::string, TH1D*> h1_             ; 
    std::map<std::string, TH2D*> h2_             ; 
    PickGenPart genpart                          ;
    const std::string fnamebtagSF_               ;
    const std::string fnameSJbtagSF_             ;
    const std::string btageffmap_                ;
    std::unique_ptr<BTagSFUtils> btagsfutils_    ; 
    std::unique_ptr<BTagSFUtils> sjbtagsfutils_  ; 
    const bool maketree_ ;
    TTree* tree_ ; 
    os2l::OS2LAnaTree os2ltree_ ; 
};

using namespace std;

// static data member definitions
void OS2LAna::fillAdditionalPlots( vlq::ElectronCollection goodElectrons,double evtwt){

  for  (unsigned int iele=0; iele<goodElectrons.size(); ++iele){
    float scEta = goodElectrons.at(iele).getscEta();
    if(fabs(scEta) <= 1.479){
      h1_["Eta_EB_el_pre"]-> Fill(goodElectrons.at(iele).getEta(), evtwt);
      h1_["Iso03_EB_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
      h1_["dEtaInSeed_EB_el_pre"]->Fill(goodElectrons.at(iele).getdEtaInSeed(), evtwt);
      h1_["dPhiIn_EB_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
      h1_["Dz_EB_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
      h1_["Dxy_EB_el_pre"]->Fill(goodElectrons.at(iele).getDxy(), evtwt);
      h1_["SCETA_EB_el_pre"]->Fill(goodElectrons.at(iele).getscEta(), evtwt);
      h1_["Full5x5siee_EB_el_pre"]->Fill(goodElectrons.at(iele).getfull5x5siee(), evtwt);
      h1_["HoE_EB_el_pre"]->Fill(goodElectrons.at(iele).getHoE(), evtwt);
      h1_["ooEmooP_EB_el_pre"]->Fill(goodElectrons.at(iele).getooEmooP(), evtwt);
      h1_["missHits_EB_el_pre"]->Fill(goodElectrons.at(iele).getmissHits(), evtwt);
      h1_["conveto_EB_el_pre"]->Fill(goodElectrons.at(iele).gethasMatchedConVeto(), evtwt);


    }
    else if  (fabs(scEta) > 1.479 && fabs(scEta) < 2.4){
      h1_["Eta_EE_el_pre"]->Fill(goodElectrons.at(iele).getEta(), evtwt);
      h1_["Iso03_EE_el_pre"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
      h1_["dEtaInSeed_EE_el_pre"]->Fill(goodElectrons.at(iele).getdEtaInSeed(), evtwt);
      h1_["dPhiIn_EE_el_pre"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
      h1_["Dz_EE_el_pre"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
      h1_["Dxy_EE_el_pre"]->Fill(goodElectrons.at(iele).getDxy(), evtwt);
      h1_["SCETA_EE_el_pre"]->Fill(goodElectrons.at(iele).getscEta(), evtwt);
      h1_["Full5x5siee_EE_el_pre"]->Fill(goodElectrons.at(iele).getfull5x5siee(), evtwt);
      h1_["HoE_EE_el_pre"]->Fill(goodElectrons.at(iele).getHoE(), evtwt);
      h1_["ooEmooP_EE_el_pre"]->Fill(goodElectrons.at(iele).getooEmooP(), evtwt);
      h1_["missHits_EE_el_pre"]->Fill(goodElectrons.at(iele).getmissHits(), evtwt);
      h1_["conveto_EE_el_pre"]->Fill(goodElectrons.at(iele).gethasMatchedConVeto(), evtwt);

    }
  }
}

// constructors and destructor
OS2LAna::OS2LAna(const edm::ParameterSet& iConfig) :
  t_evttype               (consumes<string>  (iConfig.getParameter<edm::InputTag>("evttype"))),
  t_evtwtGen              (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtGen"))),
  t_evtwtPV               (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPV"))),
  t_evtwtPVLow            (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPVLow"))),
  t_evtwtPVHigh           (consumes<double>  (iConfig.getParameter<edm::InputTag>("evtwtPVHigh"))),
  t_npv                   (consumes<unsigned>(iConfig.getParameter<edm::InputTag>("npv"))),
  t_hltdecision           (consumes<bool>    (iConfig.getParameter<edm::InputTag>("hltdecision"))),
  t_evtno                 (consumes<int>     (iConfig.getParameter<edm::InputTag>("evtno"))),
  t_runno                 (consumes<int>     (iConfig.getParameter<edm::InputTag>("runno"))),
  t_lumisec               (consumes<int>     (iConfig.getParameter<edm::InputTag>("lumisec"))),
  t_lhewtids              (consumes<vector<int> >     (iConfig.getParameter<edm::InputTag>("lhewtids"))),
  t_lhewts (consumes<vector<double> > (iConfig.getParameter<edm::InputTag>("lhewts"))),
  

// lheId_                  (iConfig.getParameter<int>               ("lheId")),
  //categorize_             (iConfig.getParameter<bool>              ("categorize")),
  //  optimizeReco_           (iConfig.getParameter<bool>              ("optimizeReco")),
  // bosonMass_              (iConfig.getParameter<double>            ("bosonMass")),

  DilepCandParams_        (iConfig.getParameter<edm::ParameterSet> ("DilepCandParams")),
  ZCandParams_            (iConfig.getParameter<edm::ParameterSet> ("ZCandParams")),
  //BoostedZCandParams_     (iConfig.getParameter<edm::ParameterSet> ("BoostedZCandParams")),
  //GenHSelParams_          (iConfig.getParameter<edm::ParameterSet> ("GenHSelParams")),
  genParams_              (iConfig.getParameter<edm::ParameterSet> ("genParams")),
  NAK4Min_                (iConfig.getParameter<unsigned int>      ("NAK4Min")),
  HTMin_                  (iConfig.getParameter<double>            ("HTMin")),
  STMin_                  (iConfig.getParameter<double>            ("STMin")), 
  STMaxControl_           (iConfig.getParameter<double>            ("STMaxControl")), 
  skim_                   (iConfig.getParameter<bool>              ("skim")), 
  filterSignal_           (iConfig.getParameter<bool>              ("filterSignal")), 
  additionalPlots_        (iConfig.getParameter<bool>              ("additionalPlots")), 
  signalType_             (iConfig.getParameter<std::string>       ("signalType")), 
  zdecayMode_             (iConfig.getParameter<std::string>       ("zdecayMode")),
  btagsf_bcUp_            (iConfig.getParameter<bool>              ("btagsf_bcUp")),
  btagsf_bcDown_          (iConfig.getParameter<bool>              ("btagsf_bcDown")),
  btagsf_lUp_             (iConfig.getParameter<bool>              ("btagsf_lUp")),
  btagsf_lDown_           (iConfig.getParameter<bool>              ("btagsf_lDown")),
  PileupUp_               (iConfig.getParameter<bool>              ("PileupUp")),
  PileupDown_             (iConfig.getParameter<bool>              ("PileupDown")),
  ElTrigUp_               (iConfig.getParameter<bool>              ("ElTrigUp")), // EL trig unc
  ElTrigDown_             (iConfig.getParameter<bool>              ("ElTrigDown")),// EL Trig Unc


  lheId_                  (iConfig.getParameter<int>               ("lheId")),
  pdfId_                  (iConfig.getParameter<int>               ("pdfId")),
  categorize_             (iConfig.getParameter<bool>              ("categorize")),

  applyLeptonIDSFs_       (iConfig.getParameter<bool>              ("applyLeptonIDSFs")), 
  applyLeptonTrigSFs_     (iConfig.getParameter<bool>              ("applyLeptonTrigSFs")),
  applyBTagSFs_           (iConfig.getParameter<bool>              ("applyBTagSFs")),
  applyDYNLOCorr_         (iConfig.getParameter<bool>              ("applyDYNLOCorr")),
  DYDown_                 (iConfig.getParameter<bool>              ("DYDown")),
  tauShift_               (iConfig.getParameter<int>               ("tauShift")),
  pdfID_offset_           (iConfig.getParameter<int>               ("pdfID_offset")),
  scale_offset_           (iConfig.getParameter<int> ("scale_offset")),

  scaleVal_               (iConfig.getParameter<double>            ("scaleVal")),
  pdfVal_                 (iConfig.getParameter<double>            ("pdfVal")),

  fname_DYNLOCorr_        (iConfig.getParameter<std::string>       ("File_DYNLOCorr")),
  funname_DYNLOCorr_      (iConfig.getParameter<std::string>       ("Fun_DYNLOCorr")),

  optimizeReco_           (iConfig.getParameter<bool>              ("optimizeReco")),
  bosonMass_              (iConfig.getParameter<double>            ("bosonMass")),

 
 dynloewkkfact           (DYNLOEwkKfact(fname_DYNLOCorr_,funname_DYNLOCorr_)),
  lepIdSFs                (iConfig.getParameter<edm::ParameterSet> ("lepIdSFsParams")),
  lepTrigSFs              (iConfig.getParameter<edm::ParameterSet> ("lepTrigSFsParams")), 
  metmaker                (iConfig.getParameter<edm::ParameterSet> ("metselParams"),consumesCollector()),
  muonmaker               (iConfig.getParameter<edm::ParameterSet> ("muselParams"),consumesCollector()),
  electronmaker           (iConfig.getParameter<edm::ParameterSet> ("elselParams"),consumesCollector()),

  jetAK4Genmaker             (iConfig.getParameter<edm::ParameterSet> ("jetAK4GenselParams"),consumesCollector()),
  jetAK4maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK4selParams"),consumesCollector()),
  jetAK4BTaggedmaker      (iConfig.getParameter<edm::ParameterSet> ("jetAK4BTaggedselParams"),consumesCollector()),
  jetAK8maker             (iConfig.getParameter<edm::ParameterSet> ("jetAK8selParams"),consumesCollector()),
  jetHTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetHTaggedselParams"),consumesCollector()),
  jetWTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetWTaggedselParams"),consumesCollector()),
  //jetwTaggedmaker         (iConfig.getParameter<edm::ParameterSet> ("jetwTaggedselParams"),consumesCollector()),
  jetTopTaggedmaker       (iConfig.getParameter<edm::ParameterSet> ("jetTopTaggedselParams"),consumesCollector()),   
  genpart                 (genParams_, consumesCollector()),
  fnamebtagSF_            (iConfig.getParameter<std::string>       ("fnamebtagSF")),
  fnameSJbtagSF_          (iConfig.getParameter<std::string>       ("fnameSJbtagSF")),
  btageffmap_             (iConfig.getParameter<std::string>       ("btageffmap")),
  btagsfutils_            (new BTagSFUtils(fnamebtagSF_,BTagEntry::OP_MEDIUM,20., 1000., 20., 1000., 20., 1000.,btageffmap_)),
  sjbtagsfutils_          (new BTagSFUtils(fnameSJbtagSF_,BTagEntry::OP_LOOSE,30., 450., 30., 450., 20., 1000.,btageffmap_)),
  maketree_               (iConfig.getParameter<bool>("maketree"))

{
  produces<vlq::JetCollection>("tjets") ; 
  produces<vlq::JetCollection>("wjets") ; 
  produces<vlq::JetCollection>("hjets") ;
  produces<vlq::JetCollection>("bjets") ; 
  produces<vlq::JetCollection>("jets") ; 
  produces<vlq::CandidateCollection>("zllcands") ; 
  produces<double>("PreWeight");
  produces<double>("btagsf");
  produces<double>("btagsfbcUp");
  produces<double>("btagsfbcDown");
  produces<double>("btagsflUp");
  produces<double>("btagsflDown");
  produces<double>("sjbtagsf");
  produces<double>("sjbtagsfbcUp");
  produces<double>("sjbtagsfbcDown");
  produces<double>("sjbtagsflUp");
  produces<double>("sjbtagsflDown");
  produces<double>("finalWeight");
}


OS2LAna::~OS2LAna() {}

bool OS2LAna::filter(edm::Event& evt, const edm::EventSetup& iSetup) {
  using namespace edm;

  Handle<string>   h_evttype     ; evt.getByToken(t_evttype     ,h_evttype    ) ; 
  Handle<double>   h_evtwtGen    ; evt.getByToken(t_evtwtGen    ,h_evtwtGen   ) ; 
  Handle<double>   h_evtwtPV     ; evt.getByToken(t_evtwtPV     ,h_evtwtPV    ) ; 
  Handle<double>   h_evtwtPVLow  ; evt.getByToken(t_evtwtPVLow  ,h_evtwtPVLow ) ;
  Handle<double>   h_evtwtPVHigh ; evt.getByToken(t_evtwtPVHigh ,h_evtwtPVHigh) ;  
  Handle<unsigned> h_npv         ; evt.getByToken(t_npv         ,h_npv        ) ; 
  Handle<bool>     h_hltdecision ; evt.getByToken(t_hltdecision ,h_hltdecision) ; 
  Handle<int>      h_evtno       ; evt.getByToken(t_evtno       ,h_evtno      ) ; 
  Handle<int>      h_runno       ; evt.getByToken(t_runno       ,h_runno      ) ; 
  Handle<int>      h_lumisec     ; evt.getByToken(t_lumisec     ,h_lumisec    ) ; 
  Handle<vector<int> > h_lhewtids; evt.getByToken(t_lhewtids    ,h_lhewtids   ) ;
  Handle<vector<double> > h_lhewts;evt.getByToken(t_lhewts ,h_lhewts ) ;



  const int npv(*h_npv.product());
  const int evtno(*h_evtno.product()) ;
  const int runno(*h_runno.product()) ;
  const int lumisec(*h_lumisec.product()) ;
  const bool isData(evtno > 0 ? true : false) ; 

  vector<pair<int, double> > lhe_id_wts;

  if (!isData){
    for (unsigned i=0; i<(*h_lhewtids.product()).size(); i++){
      int id = (*h_lhewtids.product()).at(i);
      double wt = (*h_lhewts.product()).at(i);
      lhe_id_wts.push_back(make_pair(id, wt));
    }
  }



  int signalType(-1);
  if ( !isData && filterSignal_ ) {
    if( (skim_ || maketree_ ) && signalType_.empty() ){ 
      if      (*h_evttype.product() == "EvtType_MC_bZbZ") signalType = 1; 
      else if (*h_evttype.product() == "EvtType_MC_bZbH") signalType = 2; 
      else if (*h_evttype.product() == "EvtType_MC_bZtW") signalType = 3; 
      else if (*h_evttype.product() == "EvtType_MC_bHbH") signalType = 4; 
      else if (*h_evttype.product() == "EvtType_MC_bHtW") signalType = 5; 
      else if (*h_evttype.product() == "EvtType_MC_tWtW") signalType = 6; 
      else if (*h_evttype.product() == "EvtType_MC_tZtZ") signalType = 7; 
      else if (*h_evttype.product() == "EvtType_MC_tZtH") signalType = 8; 
      else if (*h_evttype.product() == "EvtType_MC_tZbW") signalType = 9; 
      else if (*h_evttype.product() == "EvtType_MC_tHtH") signalType = 10; 
      else if (*h_evttype.product() == "EvtType_MC_tHbW") signalType = 11; 
      else if (*h_evttype.product() == "EvtType_MC_bWbW") signalType = 12; 
      h1_["signalEvts_all"] -> Fill(signalType);
    }
    else{
      if (*h_evttype.product()!=signalType_) return false ;
      else  h1_["signalEvts"] -> Fill(1) ;
    }
  }

  double evtwt(1.0);
  if (PileupUp_)         evtwt = (*h_evtwtGen.product()) * (*h_evtwtPVHigh.product()) ; 
  else if (PileupDown_)  evtwt = (*h_evtwtGen.product()) * (*h_evtwtPVLow.product()) ;
  else                   evtwt = (*h_evtwtGen.product()) * (*h_evtwtPV.product()) ;
 
  if (!isData){
    for (auto& lhe : lhe_id_wts)
      if (lhe.first == lheId_){
        evtwt *= lhe.second;
      }
  }

        
  /*   if (!isData) {
    for (unsigned i = 0; i < 9; i++){
      //  std::cout << lhe_id_wts.at(i+scale_offset_).first << std::endl;
      //  std::cout << "printed scale" << std::endl;
      h1_[Form("pre_scale%d", i+1)] ->Fill(1, lhe_id_wts.at(i+scale_offset_).second);
    }
    for (unsigned i = 0; i < 100; i++) {
      //  std::cout << lhe_id_wts.at(i+pdfID_offset_).first << std::endl;
      h1_[Form("pre_pdf%d", i+1)] -> Fill(1, lhe_id_wts.at(i+pdfID_offset_).second);
    }
    }*/
  
  

 /*
  double pdfShift(1.);
  if (!isData){
    for (auto& lhe : lhe_id_wts){
      if (lhe.first == lheId_)
        evtwt *= lhe.second;
      if (lhe.first == pdfId_) 
        pdfShift = lhe.second;
    }
  }

  evtwt *= pdfShift;
  h1_["pre_pdf"] -> Fill(1, pdfShift);

  evtwt *= scaleVal_;
  evtwt *= pdfVal_;
  */
  h1_["cutflow"] -> Fill(1, evtwt) ;

  const bool hltdecision(*h_hltdecision.product()) ; 
  if ( hltdecision ) h1_["cutflow"] -> Fill(2, evtwt) ;
  else return false; //// Presel: HLT  

  vlq::MuonCollection goodMuons; 
  muonmaker(evt, goodMuons) ; 

  vlq::ElectronCollection goodElectrons; 
  electronmaker(evt, goodElectrons) ;

  vlq::MetCollection goodMet;
  metmaker(evt, goodMet) ;

  vlq::CandidateCollection tops, BC,D1, Z,ZB1, H,Hb1,HPrime,HbPrime, ZH, ZHb;
  vlq::JetCollection Hb,ZB,D, W,B;
  vlq::CandidateCollection dimuons, dielectrons, dileptons;   
  vlq::CandidateCollection zll; //generic collection

  // dilepton properties: M > 50, lead pt > 45, second pt > 25
  DileptonCandsProducer dileptonsprod(DilepCandParams_) ; 
  dileptonsprod.operator()<vlq::MuonCollection>(dimuons, goodMuons); 
  dileptonsprod.operator()<vlq::ElectronCollection>(dielectrons, goodElectrons) ; 

  //================================================================
  //First pre-selection: 1) 2 OS dileptons from boosted Z, >=3 jets
  //================================================================

  //// Dilepton candidate
  if (zdecayMode_ == "zmumu") {dileptons = dimuons; }
  else if (zdecayMode_ == "zelel") {dileptons = dielectrons;}

  if (dileptons.size() > 0) h1_["cutflow"] -> Fill(3, evtwt) ;
  else return false ; //// Presel: dileptons

  //// Get Dy EWK correction
  if ( applyDYNLOCorr_ ) {
    double EWKNLOkfact( dynloewkkfact(dileptons.at(0).getPt()) ) ; ////GetDYNLOCorr(dileptons.at(0).getPt())) ; 
    evtwt *= EWKNLOkfact ;
  }

  //// Get lepton ID and Iso SF
  if  ( !isData ) {

    if (applyLeptonIDSFs_) {
      if ( zdecayMode_ == "zmumu" ) {
        evtwt *= lepIdSFs.IDSF(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IDSF(goodMuons.at(1).getPt(),goodMuons.at(1).getEta()) ;
        evtwt *= lepIdSFs.IsoSF(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IsoSF(goodMuons.at(1).getPt(), goodMuons.at(1).getEta()) ; 
      }
      else if ( zdecayMode_ == "zelel" ) {
        evtwt *= lepIdSFs.IDSF(goodElectrons.at(0).getPt(),goodElectrons.at(0).getEta()) ;
        evtwt *= lepIdSFs.IDSF(goodElectrons.at(1).getPt(), goodElectrons.at(1).getEta()) ;
      }
    }

    if (applyLeptonTrigSFs_) {
      if ( zdecayMode_ == "zmumu" ) evtwt *= lepTrigSFs(goodMuons.at(0).getPt(),goodMuons.at(0).getEta()) ; 
      // else if ( zdecayMode_ == "zelel" ) evtwt *= lepTrigSFs(goodElectrons.at(0).getPt(),goodElectrons.at(0).getEta()) ; 
      else if ( zdecayMode_ == "zelel" ) {
	 double leg1 = lepTrigSFs(goodElectrons.at(0).getPt(),goodElectrons.at(0).getEta()) ;
	 double leg2 = lepTrigSFs(goodElectrons.at(1).getPt(),goodElectrons.at(1).getEta()) ;
	 double elTrigSF = leg1 + leg2 - leg1*leg2 ;
	 evtwt *= elTrigSF;
      }
    }
    /*
    if ( zdecayMode_ == "zelel" ){
     
      if (ElTrigUp_){   
	//	cout << " ElTrigUp_ = " << ElTrigUp_ << endl;
	if (goodElectrons.at(0).getPt()<= 300){evtwt *= 1.02; }
	else {evtwt *= 1.04; }
	
	if (goodElectrons.at(1).getPt()<= 300){evtwt *= 1.02; }
	else {evtwt *= 1.04; }
      }
      if (ElTrigDown_){
	if (goodElectrons.at(0).getPt()<= 300){evtwt *= 0.98; }
	else {evtwt *= 0.96; }
	
	if (goodElectrons.at(1).getPt()<= 300){evtwt *= 0.98; }
	else {evtwt *= 0.96; }
      }

      }*/


  }

  //// Z mass candidate filter: 75 < M < 105, lead pt > 45, 2nd pt > 25, Z pt > 100
  CandidateFilter zllfilter(ZCandParams_) ; 
  zllfilter(dileptons, zll);

  //// jets
  vlq::JetCollection goodAK4Jets;
  jetAK4maker(evt, goodAK4Jets) ;


  // vlq::JetCollection AK4GenJets;
  // jetAK4Genmaker(evt, AK4GenJets) ;
  // cout << "no of AK4 gen jets ="<< AK4GenJets.size()<<endl;
  // for (unsigned i =0; i < AK4GenJets.size();i++){
  //   cout << i <<"th AK4 gen jet pt is =" << AK4GenJets.at(i).getPt()<< endl;
  //  }

  // cout << "****no of AK4  jets ="<< goodAK4Jets.size()<<endl;
  // for (unsigned i =0; i < goodAK4Jets.size();i++){
  //   cout<<"***** " << i <<"th AK4  jet pt is =" << goodAK4Jets.at(i).getPt()<< endl;
  //  }

  // b-tagging:
  vlq::JetCollection goodBTaggedAK4Jets;
  jetAK4BTaggedmaker(evt, goodBTaggedAK4Jets) ; 

  //// jet cleaning w.r.t dileptons
  CandidateCleaner cleanjets(0.4, -1); //// The second argument is for lepton 2D iso, setting to -1 disables it
  CandidateCleaner cleanak8jets(0.8, -1);
  
 if (zdecayMode_ == "zmumu") {
    cleanjets(goodAK4Jets, goodMuons);
    cleanjets(goodBTaggedAK4Jets, goodMuons);
  }
  else if (zdecayMode_ == "zelel") {
    cleanjets(goodAK4Jets, goodElectrons);
    cleanjets(goodBTaggedAK4Jets, goodElectrons);
  } 

 //// Exactly one Z cand in event
  if(zll.size() == 1) {h1_["cutflow"] -> Fill(4, evtwt) ;}
  else return false ; //// Presel : Z->ll
  /*
  //Eelctron id efficiency 
  if (*h_evttype.product() != "EvtType_Data"){                                                                                                                                                 
                                                                                                                                                                                               
    GenParticleCollection genPartsInfo;                                                                                                                                                        
    genPartsInfo = genpart(evt) ;                                                                                                                                                              
                                                                                                                                                                                               
    //two gen level leptons with opposite sign                                                                                                                                                 
    TLorentzVector e1,e2,Zcan,lepton1,lepton2;                                                                                                                                                 
                                                                                                                                                                                               
                                                                                                                                                                                               
    for (auto& gen1 : genPartsInfo){               
      if  ( zdecayMode_ == "zelel") {                                                                                                                                                          
        if (gen1.getPdgID() == 11 && gen1.getMom0PdgID()==23){ e1 = gen1.getP4();}                                                                                                             
                                                                                                                                                                                               
        else if (gen1.getPdgID() == -11 && gen1.getMom0PdgID()==23){ e2= gen1.getP4();}                                                                                                        
                                                                                                                                                                                               
      }           
      lepton1 = e1;                                                                                                                                                                            
      lepton2 =e2;                                                                                                                                                                             
      if (lepton1.Pt()< 25.0) continue;                                                                                                                                                       
      if (lepton1.Eta()<-2.4 && lepton1.Eta()>2.4) continue;                                                                                                                                   
      if (lepton2.Pt()< 25.0) continue;                                                                                                                                                       
      if (lepton2.Eta()< -2.4 && lepton2.Eta()>2.4) continue;   


      if  ( zdecayMode_ == "zelel") {
	if (goodElectrons.size()>0){
	  //  h1_["pt_leadinlepton"]  -> Fill(goodElectrons.at(0).getPt(), evtwt) ;
	  // h1_["eta_leadinglepton"]  -> Fill(goodElectrons.at(0).getEta(), evtwt) ;
	  // h1_["pt_genlepton1"]  -> Fill(lepton1.Pt(), evtwt) ;
	  // h1_["eta_genlepton1"]  -> Fill(lepton1.Eta(), evtwt) ;
	  // h1_["pt_genlepton2"]  -> Fill(lepton2.Pt(), evtwt) ;
	  //  h1_["eta_genlepton2"]  -> Fill(lepton2.Eta(), evtwt) ;

	  if (fabs(goodElectrons.at(0).getscEta())<= 1.479 && fabs(lepton1.Eta())<= 1.479){
	    if ((lepton1.DeltaR(goodElectrons.at(0).getP4()))<0.4 || (lepton2.DeltaR(goodElectrons.at(0).getP4()))<0.4)  {
	      h1_["ElEB1"] -> Fill(1, evtwt) ;
	      h1_["pt_leadinlepton_drmatchedEB"]  -> Fill(goodElectrons.at(0).getPt(), evtwt) ;
	      h1_["eta_leadinglepton_drmatchedEB"]  -> Fill(goodElectrons.at(0).getEta(), evtwt) ;
	     
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308)  {h1_["ElEB1"] -> Fill(2, evtwt) ;}
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816)  {h1_["ElEB1"] -> Fill(3, evtwt) ;}
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816 && goodElectrons.at(0).getfull5x5siee() < 0.00998 )  {h1_["ElEB1"] -> Fill(4, evtwt) ;}
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816 && goodElectrons.at(0).getfull5x5siee() < 0.00998 && goodElectrons.at(0).getHoE() < 0.0414 )  {h1_["ElEB1"] -> Fill(5, evtwt) ;}
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816 && goodElectrons.at(0).getfull5x5siee() < 0.00998 && goodElectrons.at(0).getHoE() < 0.0414 && goodElectrons.at(0).getooEmooP() <  0.0129  )  {h1_["ElEB1"] -> Fill(6, evtwt) ;}
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816 && goodElectrons.at(0).getfull5x5siee() < 0.00998 && goodElectrons.at(0).getHoE() < 0.0414 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0588  )  {h1_["ElEB1"] -> Fill(7, evtwt) ;}
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816 && goodElectrons.at(0).getfull5x5siee() < 0.00998 && goodElectrons.at(0).getHoE() < 0.0414 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0588 && !(goodElectrons.at(0).gethasMatchedConVeto())  )  {h1_["ElEB1"] -> Fill(8, evtwt) ;}
	    
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816 && goodElectrons.at(0).getfull5x5siee() < 0.00998 && goodElectrons.at(0).getHoE() < 0.0414 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0588 && !(goodElectrons.at(0).gethasMatchedConVeto()) && goodElectrons.at(0).getmissHits() <= 1 )  {h1_["ElEB1"] -> Fill(9, evtwt) ;}

	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816 && goodElectrons.at(0).getfull5x5siee() < 0.00998 && goodElectrons.at(0).getHoE() < 0.0414 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0588 && !(goodElectrons.at(0).gethasMatchedConVeto()) && goodElectrons.at(0).getmissHits() <= 1 && fabs(goodElectrons.at(0).getDxy())< 0.05 )  {h1_["ElEB1"] -> Fill(10, evtwt) ;}
	    
	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00308 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0816 && goodElectrons.at(0).getfull5x5siee() < 0.00998 && goodElectrons.at(0).getHoE() < 0.0414 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0588 && !(goodElectrons.at(0).gethasMatchedConVeto()) && goodElectrons.at(0).getmissHits() <= 1 && fabs(goodElectrons.at(0).getDxy())< 0.05 && fabs(goodElectrons.at(0).getDz())<0.10 )  {
		h1_["ElEB1"] -> Fill(11, evtwt) ;
		h1_["pt_leadinlepton_idmatchedEB"]  -> Fill(goodElectrons.at(0).getPt(), evtwt) ;
	       	h1_["eta_leadinglepton_idmatchedEB"]  -> Fill(goodElectrons.at(0).getEta(), evtwt) ;
		//	h1_["dr_elel_idmatchedEB"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
	      }
	    }
	  }

	  else if (fabs(goodElectrons.at(0).getscEta())> 1.479 && fabs(goodElectrons.at(0).getscEta())< 2.5 && fabs(lepton1.Eta())>1.479 && fabs(lepton1.Eta())<2.5){
	    if ((lepton1.DeltaR(goodElectrons.at(0).getP4()))<0.4 || (lepton2.DeltaR(goodElectrons.at(0).getP4()))<0.4){
	      h1_["ElEE1"] -> Fill(1, evtwt) ;
	      // h1_["ElEE2"] -> Fill(1, evtwt) ;
	       h1_["pt_leadinlepton_drmatchedEE"]  -> Fill(goodElectrons.at(0).getPt(), evtwt) ;
	       h1_["eta_leadinglepton_drmatchedEE"]  -> Fill(goodElectrons.at(0).getEta(), evtwt) ;
	      // h1_["dr_elel_drmatchedEE"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );

	      if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605)  {h1_["ElEE1"] -> Fill(2, evtwt) ;}
              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394)  {h1_["ElEE1"] -> Fill(3, evtwt) ;}
              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394 && goodElectrons.at(0).getfull5x5siee() < 0.0292 )  {h1_["ElEE1"] -> Fill(4, evtwt) ;}
              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394 && goodElectrons.at(0).getfull5x5siee() < 0.0292 && goodElectrons.at(0).getHoE() < 0.0641 )  {h1_["ElEE1"] -> Fill(5, evtwt) ;}
              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394 && goodElectrons.at(0).getfull5x5siee() < 0.0292 && goodElectrons.at(0).getHoE() < 0.0641 && goodElectrons.at(0).getooEmooP() <  0.0129  )  {h1_["ElEE1"] -> Fill(6, evtwt) ;}
              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394 && goodElectrons.at(0).getfull5x5siee() < 0.0292 && goodElectrons.at(0).getHoE() < 0.0641 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0571  )  {h1_["ElEE1"] -> Fill(7, evtwt) ;}
              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394 && goodElectrons.at(0).getfull5x5siee() < 0.0292 && goodElectrons.at(0).getHoE() < 0.0641 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0571 && !(goodElectrons.at(0).gethasMatchedConVeto())  )  {h1_["ElEE1"] -> Fill(8, evtwt) ;}

              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394 && goodElectrons.at(0).getfull5x5siee() < 0.0292 && goodElectrons.at(0).getHoE() < 0.0641 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0571 && !(goodElectrons.at(0).gethasMatchedConVeto()) && goodElectrons.at(0).getmissHits() <= 1 )  {h1_["ElEE1"] -> Fill(9, evtwt) ;}

              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394 && goodElectrons.at(0).getfull5x5siee() < 0.0292 && goodElectrons.at(0).getHoE() < 0.0641 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0571 && !(goodElectrons.at(0).gethasMatchedConVeto()) && goodElectrons.at(0).getmissHits() <= 1 && fabs(goodElectrons.at(0).getDxy())< 0.10 )  {h1_["ElEE1"] -> Fill(10, evtwt) ;}

              if (fabs(goodElectrons.at(0).getdEtaInSeed()) < 0.00605 && fabs(goodElectrons.at(0).getdPhiIn())< 0.0394 && goodElectrons.at(0).getfull5x5siee() < 0.0292 && goodElectrons.at(0).getHoE() < 0.0641 && goodElectrons.at(0).getooEmooP() <  0.0129 && goodElectrons.at(0).getIso03() < 0.0571 && !(goodElectrons.at(0).gethasMatchedConVeto()) && goodElectrons.at(0).getmissHits() <= 1 && fabs(goodElectrons.at(0).getDxy())< 0.10 && fabs(goodElectrons.at(0).getDz())<0.20 )  {
		h1_["ElEE1"] -> Fill(11, evtwt) ;
		h1_["pt_leadinlepton_idmatchedEE"]  -> Fill(goodElectrons.at(0).getPt(), evtwt) ;
		h1_["eta_leadinglepton_idmatchedEE"]  -> Fill(goodElectrons.at(0).getEta(), evtwt) ;
		//		h1_["dr_elel_idmatchedEE"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );

	      }
	    }                           
	  }



	}
      }
    }
  }//end of electron id efficiency
*/

  //// HT selection
  HT htak4(goodAK4Jets) ; 
  if ( htak4.getHT() > HTMin_ ){
    h1_["cutflow"] -> Fill(5, evtwt) ;
    // h1_["nak4HT_pre"] -> Fill(goodAK4Jets.size(),evtwt) ;
    // if ( goodBTaggedAK4Jets.size() == 0){
    // if ( (goodAK4Jets.size()>0 && goodAK4Jets.at(0).getPt() > 100)|| (goodAK4Jets.size()>1 && goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50)){
    //	h1_["nak4HT_0b3"] -> Fill(goodAK4Jets.size(), evtwt) ;
    //  }
    // }
  }
  else return false ; //// Presel: HT cut

  //// At least 3 AK4 jets in event
  if (goodAK4Jets.size()  >= NAK4Min_ ) {h1_["cutflow"] -> Fill(6, evtwt) ;} 
  else return false; //// Presel N(AK4) 
  /*  
  //gen AK4 studies 

  if (*h_evttype.product() != "EvtType_Data"){                                                                                                                                               
    GenParticleCollection genPartsInfo;       
    genPartsInfo = genpart(evt) ;                                                                                                                                                            
    //two gen level leptons with opposite sign                                                                                                                                               
    TLorentzVector e1,e2,lepton1,lepton2,Z1,Z2;                                                                                                                                               
    vlq::CandidateCollection ZCan1, ZCan2, ZCan ;
    ZCan1.clear();
    ZCan2.clear();
    for (auto& gen : genPartsInfo){
      if (gen.getPdgID() == 23){
	if (gen.getP4().Pt() <100) continue;
	for (auto& gen1 : genPartsInfo){                                                                                                                                                        
	  if  ( zdecayMode_ == "zelel") {  
	    if (gen1.getPdgID() == 11 && gen1.getMom0PdgID()==23){ e1 = gen1.getP4();}    
	    else if (gen1.getPdgID() == -11 && gen1.getMom0PdgID()==23){ e2= gen1.getP4();}                                         
	  }

	  else if ( zdecayMode_ == "zmumu") {
            if (gen1.getPdgID() == 13 && gen1.getMom0PdgID()==23){ e1 = gen1.getP4();}
            else if (gen1.getPdgID() == -13 && gen1.getMom0PdgID()==23){ e2= gen1.getP4();}
	  }
	  if (e1.Pt()< 25.0) continue;
	  if (e1.Eta()< -2.4 && e1.Eta()>2.4) continue;
	  if (e2.Pt()< 25.0) continue;
          if (e2.Eta()< -2.4 && e2.Eta()>2.4) continue;
	  lepton1 = e1;                                                                                                            
	  lepton2 =e2;                                                                                                           
	  ZCan1.push_back(lepton1);
	  ZCan2.push_back(lepton2);
	}
	if (ZCan1.size()>0 && ZCan2.size()>0)  ZCan.push_back(gen.getP4());
							      
      }
    }
    
    
    // cout <<"can 1 size ="<< ZCan1.size()<<endl;
    // cout <<"can 2 size ="<< ZCan2.size()<<endl;
    // cout <<"Z Can size =" << ZCan.size()<<endl;
    if (ZCan.size()>0){

      for (unsigned i =0; i < ZCan.size();i++){
	//cout << " ZCan.at(i).getPt() = " << ZCan.at(i).getPt()<<endl;
	h1_["ZPtGen_pre"] -> Fill( ZCan.at(i).getPt()) ;
	
      }



      vlq::JetCollection AK4GenJets , matchedAK4GenJets;                                                                                                                                                         
      jetAK4Genmaker(evt, AK4GenJets); 
	
      double htgen=0;
      for (unsigned i =0; i < AK4GenJets.size();i++){                                                                                                                                      
	//cout << i <<"th AK4 gen jet pt is =" << AK4GenJets.at(i).getPt()<< endl;                                                                                                      
	htgen += AK4GenJets.at(i).getPt();
      }      
      //    cout << " ht gen = " << htgen <<endl;
      if ( htgen > 200){
	if (AK4GenJets.size() >= 3){
	  //if ( AK4GenJets.size() >=3 && AK4GenJets.at(0).getPt()> 100 && AK4GenJets.at(1).getPt()> 50 && AK4GenJets.at(2).getPt()> 30){

	  for (unsigned i =0; i < AK4GenJets.size();i++){
	    for (unsigned j =0; j < goodAK4Jets.size();j++){
	      if (AK4GenJets.at(i).getP4().DeltaR(goodAK4Jets.at(j).getP4()) < 0.4){ matchedAK4GenJets.push_back(AK4GenJets.at(i));}
	      // cout << " AK4GenJets.at(i).Pt , goodAK4Jets.at(j).Pt , dR = " << AK4GenJets.at(i).getPt()<<","<<goodAK4Jets.at(j).getPt()<<","<< AK4GenJets.at(i).getP4().DeltaR(goodAK4Jets.at(j).getP4())<< endl;
	    }
	  }
	  if ( matchedAK4GenJets.size()>=3 ){ 
	// cout <<"can 1 size ="<< ZCan1.size()<<endl;
	  // cout <<"can 2 size ="<< ZCan2.size()<<endl;
	  // cout << " ht gen = " << htgen <<endl;
	  // cout << " ak 4 passed = " <<  matchedAK4GenJets.size()<<endl;
	  h1_["nak4Gen_pre"] -> Fill( matchedAK4GenJets.size()) ;
	  h1_["htGen_pre"] -> Fill( htgen) ;
	  // cout << "&&& " << endl;
	  // cout << "matched AK4 gen size = " << matchedAK4GenJets.size()<<endl;
	  ////Gen  Ak4 jet plots               
	  for(unsigned j=0; j<matchedAK4GenJets.size(); ++j){          
	    //  cout<<j << "th matched jet pt = "<< matchedAK4GenJets.at(j).getPt() << endl;
	  }
	  //	  cout << " End event %%%%%%%%%%%%%%%%%%%%%% "<< endl;
     
	  for(unsigned j=0; j<3; ++j){
	    // if ( matchedAK4GenJets.size()> j){
	    h1_[Form("ptak4jetGen%d_pre", j+1)]  -> Fill( matchedAK4GenJets.at(j).getPt()) ;
	    h1_[Form("etaak4jetGen%d_pre", j+1)] -> Fill( matchedAK4GenJets.at(j).getEta()) ;
	    // h1_[Form("cvsak4jet%d_pre", j+1)] -> Fill( matchedAK4GenJets.at(j).getCSV()) ;
	    // h1_[Form("massak4jetGen%d_pre", j+1)] -> Fill( matchedAK4GenJets.at(j).getMass()) ;
	    }
	  }


	}
	
      }
      
    }
  }

  */

  
  vlq::JetCollection goodAK8Jets ;
  jetAK8maker(evt, goodAK8Jets); 
  // cleanjets(goodAK8Jets, goodMuons); 
  // cleanjets(goodAK8Jets, goodElectrons); 
  cleanak8jets(goodAK8Jets, goodMuons); 
  cleanak8jets(goodAK8Jets, goodElectrons); 

  double presel_wt(evtwt);
  double btagsf(1) ;
  double btagsf_bcUp(1) ; 
  double btagsf_bcDown(1) ; 
  double btagsf_lUp(1) ; 
  double btagsf_lDown(1) ; 
  if ( applyBTagSFs_ ) {
    std::vector<double>csvs;
    std::vector<double>pts;
    std::vector<double>etas;
    std::vector<int>   flhads;

    for (vlq::Jet jet : goodAK4Jets) {
      csvs.push_back(jet.getCSV()) ; 
      pts.push_back(jet.getPt()) ; 
      etas.push_back(jet.getEta()) ; 
      flhads.push_back(jet.getHadronFlavour()) ; 
    }

    btagsfutils_->getBTagSFs (csvs, pts, etas, flhads, jetAK4BTaggedmaker.idxjetCSVDiscMin_, btagsf, btagsf_bcUp, btagsf_bcDown, btagsf_lUp, btagsf_lDown) ;

    //// bTag SF along with sys. unc. options
    if (btagsf_bcUp_)
      evtwt *= btagsf_bcUp;
    else if (btagsf_bcDown_)
      evtwt *= btagsf_bcDown;
    else  if (btagsf_lUp_)
      evtwt *= btagsf_lUp;
    else if (btagsf_lDown_)
      evtwt *= btagsf_lDown;
    else
      evtwt *= btagsf;
  }

  double ST(htak4.getHT() + zll.at(0).getPt() + goodMet.at(0).getFullPt()); 
 
  if (goodAK4Jets.at(0).getPt() > 100 ) {
    h1_["cutflow"] -> Fill(7, presel_wt) ; 
    if (goodAK4Jets.at(1).getPt() > 50){
      h1_["cutflow"] -> Fill(8, presel_wt) ; 
      if ( goodBTaggedAK4Jets.size() > 0 ) { 
        h1_["cutflow"] -> Fill(9, evtwt) ;
        if ( ST > STMin_ ) h1_["cutflow"] -> Fill(10, evtwt) ;  
      }
    }
  } //// Completing the cutflow

  vlq::JetCollection  goodWTaggedJets;
  jetWTaggedmaker(evt, goodWTaggedJets);
  //  cleanjets(goodWTaggedJets, goodMuons); 
  // cleanjets(goodWTaggedJets, goodElectrons); 
  cleanak8jets(goodWTaggedJets, goodMuons); 
  cleanak8jets(goodWTaggedJets, goodElectrons); 


  // vlq::JetCollection  goodwTaggedJets;
  // jetwTaggedmaker(evt, goodwTaggedJets);
  // cleanjets(goodwTaggedJets, goodMuons);
  // cleanjets(goodwTaggedJets, goodElectrons);


  vlq::JetCollection goodHTaggedJets; 
  jetHTaggedmaker(evt, goodHTaggedJets);
  //  cleanjets(goodHTaggedJets, goodMuons); 
  //cleanjets(goodHTaggedJets, goodElectrons); 
  cleanak8jets(goodHTaggedJets, goodMuons); 
  cleanak8jets(goodHTaggedJets, goodElectrons); 

  //// http://cms.cern.ch/iCMS/jsp/openfile.jsp?tp=draft&files=AN2016_245_v3.pdf
  //// SF of 0.93+/-0.09 required for top tag WP with mistag rate 1% (no subjet b tag): AN2016-245v3
  vlq::JetCollection goodTopTaggedJets;
  jetTopTaggedmaker(evt, goodTopTaggedJets);
  //cleanjets(goodTopTaggedJets, goodMuons); 
  // cleanjets(goodTopTaggedJets, goodElectrons);
  cleanak8jets(goodTopTaggedJets, goodMuons); 
  cleanak8jets(goodTopTaggedJets, goodElectrons); 

  //cout<< "event weight befor taushift, totags =" << evtwt <<","<< goodTopTaggedJets.size() <<endl; 
  // cout<< "tauShift_ = " << tauShift_ <<endl;
  if (!isData){ 
    for (unsigned i=0; i<goodTopTaggedJets.size(); i++) {
      // evtwt *= ( 1.06 + (tauShift_ * .09));
      evtwt *= ( 1.00 + (tauShift_ * .09)); 
    }
  }
  //cout<< "event weight after taushift =" << evtwt<<endl;
  double sjbtagsf(1) ;
  double sjbtagsf_bcUp(1) ; 
  double sjbtagsf_bcDown(1) ; 
  double sjbtagsf_lUp(1) ; 
  double sjbtagsf_lDown(1) ; 
  if ( applyBTagSFs_ ) {
    std::vector<double>csvs;
    std::vector<double>pts;
    std::vector<double>etas;
    std::vector<int>   flhads;

    for (vlq::Jet jet : goodHTaggedJets) {
      csvs.push_back(jet.getCSVSubjet0()) ; 
      csvs.push_back(jet.getCSVSubjet1()) ; 
      pts.push_back(jet.getPtSubjet0()) ; 
      pts.push_back(jet.getPtSubjet1()) ; 
      etas.push_back(jet.getEtaSubjet0()) ; 
      etas.push_back(jet.getEtaSubjet1()) ; 
      flhads.push_back(jet.getHadronFlavourSubjet0()) ; 
      flhads.push_back(jet.getHadronFlavourSubjet1()) ; 
    }

    sjbtagsfutils_->getBTagSFs (csvs, pts, etas, flhads, jetHTaggedmaker.idxsjCSVMin_, sjbtagsf, sjbtagsf_bcUp, sjbtagsf_bcDown, sjbtagsf_lUp, sjbtagsf_lDown) ;

    //// bTag SF along with sys. unc. options
    if (btagsf_bcUp_)
      evtwt *= sjbtagsf_bcUp;
    else if (btagsf_bcDown_)
      evtwt *= sjbtagsf_bcDown;
    else  if (btagsf_lUp_)
      evtwt *= sjbtagsf_lUp;
    else if (btagsf_lDown_)
      evtwt *= sjbtagsf_lDown;
    else
      evtwt *= sjbtagsf;

  }
  // double w1=1.0,w2=1.0,w3=1.0,w4=1.0,w5=1.0;

    if ( skim_ ) {

    h1_["ht_preSel_bfFit"] -> Fill(htak4.getHT(), presel_wt);
    //fill jet flavor pt and eta for b-tagging efficiencies 
    if (goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){// not sure if we should also cut on ST
      for (vlq::Jet jet : goodAK4Jets) {
        if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
        else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
        else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_all"] -> Fill(jet.getPt(), jet.getEta()) ; 
      }
      if ( goodBTaggedAK4Jets.size() > 0 ){
        for (vlq::Jet jet : goodBTaggedAK4Jets) {
          if ( abs(jet.getHadronFlavour()) == 5) h2_["pt_eta_b_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
          else if ( abs(jet.getHadronFlavour()) == 4) h2_["pt_eta_c_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
          else if ( abs(jet.getHadronFlavour()) == 0) h2_["pt_eta_l_btagged"] -> Fill(jet.getPt(), jet.getEta()) ; 
        }
      }
    }
  
    std::unique_ptr<double> ptr_fevtwt ( new double(evtwt) ) ;
    std::unique_ptr<double> ptr_prewt ( new double(presel_wt) ) ;
    std::unique_ptr<double> ptr_btagsf          ( new double(btagsf         ) ) ;
    std::unique_ptr<double> ptr_btagsf_bcUp     ( new double(btagsf_bcUp    ) ) ;
    std::unique_ptr<double> ptr_btagsf_bcDown   ( new double(btagsf_bcDown  ) ) ;
    std::unique_ptr<double> ptr_btagsf_lUp      ( new double(btagsf_lUp     ) ) ;
    std::unique_ptr<double> ptr_btagsf_lDown    ( new double(btagsf_lDown   ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf        ( new double(sjbtagsf       ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_bcUp   ( new double(sjbtagsf_bcUp  ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_bcDown ( new double(sjbtagsf_bcDown) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_lUp    ( new double(sjbtagsf_lUp   ) ) ;
    std::unique_ptr<double> ptr_sjbtagsf_lDown  ( new double(sjbtagsf_lDown ) ) ;
    evt.put(std::move(ptr_fevtwt), "finalWeight");
    evt.put(std::move(ptr_prewt), "PreWeight");
    evt.put(std::move(ptr_btagsf         ), "btagsf"         );
    evt.put(std::move(ptr_btagsf_bcUp    ), "btagsfbcUp"    );
    evt.put(std::move(ptr_btagsf_bcDown  ), "btagsfbcDown"  );
    evt.put(std::move(ptr_btagsf_lUp     ), "btagsflUp"     );
    evt.put(std::move(ptr_btagsf_lDown   ), "btagsflDown"   );
    evt.put(std::move(ptr_sjbtagsf       ), "sjbtagsf"       );
    evt.put(std::move(ptr_sjbtagsf_bcUp  ), "sjbtagsfbcUp"  );
    evt.put(std::move(ptr_sjbtagsf_bcDown), "sjbtagsfbcDown");
    evt.put(std::move(ptr_sjbtagsf_lUp   ), "sjbtagsflUp"   );
    evt.put(std::move(ptr_sjbtagsf_lDown ), "sjbtagsflDown" );

    if(goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50 && goodBTaggedAK4Jets.size() > 0 && ST > STMin_){
      std::unique_ptr<vlq::JetCollection> ptr_tjets( new vlq::JetCollection(goodTopTaggedJets) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_wjets( new vlq::JetCollection(goodWTaggedJets) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_hjets( new vlq::JetCollection(goodHTaggedJets) ) ;
      std::unique_ptr<vlq::JetCollection> ptr_bjets( new vlq::JetCollection(goodBTaggedAK4Jets ) ) ; 
      std::unique_ptr<vlq::JetCollection> ptr_jets ( new vlq::JetCollection(goodAK4Jets ) ) ; 
      std::unique_ptr<vlq::CandidateCollection> ptr_zllcands ( new vlq::CandidateCollection(zll) ) ; 

      evt.put(std::move(ptr_tjets), "tjets") ; 
      evt.put(std::move(ptr_wjets), "wjets") ; 
      evt.put(std::move(ptr_hjets), "hjets") ; 
      evt.put(std::move(ptr_bjets), "bjets") ; 
      evt.put(std::move(ptr_jets), "jets")  ; 
      evt.put(std::move(ptr_zllcands), "zllcands")  ;
    }   

    } //// if skim 
  else if ( !maketree_ ) { 
    
    if ( applyDYNLOCorr_ ) {
      
	//if ( zdecayMode_ == "zmumu" ){
     //   if      (goodAK4Jets.size()==3){w1=1.0105 ;presel_wt *= w1; evtwt *= w1;}
       // else if (goodAK4Jets.size()==4){w2=1.0604 ;presel_wt *= w2; evtwt *= w2;}
       // else if (goodAK4Jets.size()==5){w3=1.1882 ;presel_wt *= w3; evtwt *= w3;}
       // else if (goodAK4Jets.size()==6){w4=1.3290 ;presel_wt *= w4; evtwt *= w4;}
       // else if (goodAK4Jets.size()>=7){w5=1.6049 ;presel_wt *= w5; evtwt *= w5;}
//	}
//	else if( zdecayMode_ == "zelel" ){
  //      if      (goodAK4Jets.size()==3){w1=0.9730 ; presel_wt *= w1; evtwt *= w1;}
    //    else if (goodAK4Jets.size()==4){w2=1.0248 ; presel_wt *= w2; evtwt *= w2;}
     //   else if (goodAK4Jets.size()==5){w3=1.1125 ; presel_wt *= w3; evtwt *= w3;}
     //   else if (goodAK4Jets.size()==6){w4=1.224 ; presel_wt *= w4; evtwt *= w4;}
      //  else if (goodAK4Jets.size()>=7){ w5=1.364 ; presel_wt *= w5; evtwt *= w5;}
//	}
      
      //if (!DYDown_){
	if (DYDown_){  
	//cout << " Pass " << endl;   
	/*
	if ( zdecayMode_ == "zmumu"){
	  if (goodAK4Jets.size() == 3) { evtwt *= 1.0105; presel_wt *= 1.0105;}
	  else if (goodAK4Jets.size() == 4) { evtwt*= 1.0604; presel_wt *= 1.0604;}
	  else if (goodAK4Jets.size() == 5) { evtwt*= 1.1882; presel_wt *= 1.1882;}
	  else if (goodAK4Jets.size() == 6) { evtwt*= 1.3290; presel_wt *=1.3290;}
	  else if (goodAK4Jets.size() >= 7) { evtwt*= 1.6049; presel_wt *= 1.6049;}
	}
	else if ( zdecayMode_ == "zelel"){
	  if (goodAK4Jets.size() == 3) { evtwt *= 0.9730; presel_wt   *= 0.9730;}
	  else if (goodAK4Jets.size() == 4) { evtwt *= 1.0248; presel_wt  *= 1.0248;}
	  else if (goodAK4Jets.size() == 5) { evtwt *= 1.1125; presel_wt  *= 1.1125;}
	  else if (goodAK4Jets.size() == 6) { evtwt *= 1.224; presel_wt  *= 1.224;}
	  else if (goodAK4Jets.size() >= 7) { evtwt *= 1.364; presel_wt  *= 1.364;}
	}
	
	if ( zdecayMode_ == "zelel"){
          if (goodAK4Jets.size() == 3) { evtwt *= 0.9522; presel_wt *= 0.9522;}
          else if (goodAK4Jets.size() == 4) { evtwt*= 1.0723; presel_wt *= 1.0723;}
          else if (goodAK4Jets.size() == 5) { evtwt*= 1.0928; presel_wt *= 1.0928;}
          else if (goodAK4Jets.size() == 6) { evtwt*= 1.2423; presel_wt *=1.2423;}
          else if (goodAK4Jets.size() >= 7) { evtwt*= 1.4094; presel_wt *= 1.4095;}
        }
        else if ( zdecayMode_ == "zmumu"){
          if (goodAK4Jets.size() == 3) { evtwt *= 1.0201; presel_wt   *= 1.0201;}
          else if (goodAK4Jets.size() == 4) { evtwt *= 1.1148; presel_wt  *= 1.1148;}
          else if (goodAK4Jets.size() == 5) { evtwt *= 1.2352; presel_wt  *= 1.2352;}
          else if (goodAK4Jets.size() == 6) { evtwt *= 1.3741; presel_wt  *= 1.3741;}
          else if (goodAK4Jets.size() >= 7) { evtwt *= 1.6653; presel_wt  *= 1.6653;}
        }
	*/

	if ( zdecayMode_ == "zelel"){
          if (goodAK4Jets.size() == 3) { evtwt *= 0.87; presel_wt *= 0.87;}
          else if (goodAK4Jets.size() == 4) { evtwt*= 1.02; presel_wt *= 1.02;}
          else if (goodAK4Jets.size() == 5) { evtwt*= 1.05; presel_wt *= 1.05;}
          else if (goodAK4Jets.size() == 6) { evtwt*= 1.24; presel_wt *=1.24;}
          else if (goodAK4Jets.size() >= 7) { evtwt*= 1.60; presel_wt *= 1.60;}
        }
        else if ( zdecayMode_ == "zmumu"){
          if (goodAK4Jets.size() == 3) { evtwt *= 0.94; presel_wt   *= 0.94;}
          else if (goodAK4Jets.size() == 4) { evtwt *= 1.05; presel_wt  *= 1.05;}
          else if (goodAK4Jets.size() == 5) { evtwt *= 1.18; presel_wt  *= 1.18;}
          else if (goodAK4Jets.size() == 6) { evtwt *= 1.39; presel_wt  *= 1.39;}
          else if (goodAK4Jets.size() >= 7) { evtwt *= 1.77; presel_wt  *= 1.77;}
        }





	   }
    }
    
    std::string lep("");
    if(zdecayMode_ == "zmumu") {lep = "mu";}
    else if ( zdecayMode_ == "zelel") {lep = "el";}
    else edm::LogError("OS2LAna::filter") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;




    for (auto izll : zll) {
      h1_["mass_z"+lep+lep+"_pre"] -> Fill(izll.getMass(), presel_wt) ;
      h1_["mass_Z"+lep+lep+"_pre"] -> Fill(izll.getMass(), presel_wt) ;
      h1_["pt_z"+lep+lep+"_pre"] -> Fill(izll.getPt(), presel_wt) ; 
    }
      // if (izll.getPt()>700){
    
    //if (zdecayMode_ == "zelel" ) { 
    // if (goodElectrons.at(0).getP4().DeltaR(goodElectrons.at(1).getP4()) < 0.4){
	
    //	for (auto izll : zll) {
    //	  h1_["ptz_ex_pre"] -> Fill(izll.getPt(), presel_wt);
    //	  h1_["ptz1_ex_pre"] -> Fill(izll.getPt(), presel_wt);
    //	}
    //	h1_["ptak4jet1_ex_pre"]  -> Fill(goodAK4Jets.at(0).getPt(), presel_wt) ;
    //	h1_["etaak4jet1_ex_pre"] -> Fill(goodAK4Jets.at(0).getEta(), presel_wt) ;
    //	h1_["massak4jet1_ex_pre"] -> Fill(goodAK4Jets.at(0).getMass(), presel_wt) ;
    //	h1_["phiak4jet1_ex_pre"] -> Fill(goodAK4Jets.at(0).getPhi(), presel_wt) ;
    //	h1_["energyak4jet1_ex_pre"] -> Fill(goodAK4Jets.at(0).getP4().Energy(), presel_wt) ;
	

    //	h1_["ptak4jet2_ex_pre"]  -> Fill(goodAK4Jets.at(1).getPt(), presel_wt) ;
    //	h1_["etaak4jet2_ex_pre"] -> Fill(goodAK4Jets.at(1).getEta(), presel_wt) ;
    //	h1_["massak4jet2_ex_pre"] -> Fill(goodAK4Jets.at(1).getMass(), presel_wt) ;
    //	h1_["phiak4jet2_ex_pre"] -> Fill(goodAK4Jets.at(1).getPhi(), presel_wt) ;
    //	h1_["energyak4jet2_ex_pre"] -> Fill(goodAK4Jets.at(1).getP4().Energy(), presel_wt) ;
	
    //	h1_["ptak4jet3_ex_pre"]  -> Fill(goodAK4Jets.at(2).getPt(), presel_wt) ;
    //	h1_["etaak4jet3_ex_pre"] -> Fill(goodAK4Jets.at(2).getEta(), presel_wt) ;
    //	h1_["massak4jet3_ex_pre"] -> Fill(goodAK4Jets.at(2).getMass(), presel_wt) ;
    //	h1_["phiak4jet3_ex_pre"] -> Fill(goodAK4Jets.at(2).getPhi(), presel_wt) ;
    //	h1_["energyak4jet3_ex_pre"] -> Fill(goodAK4Jets.at(2).getP4().Energy(), presel_wt) ;
	

      /*
      if ( zdecayMode_ == "zmumu" ){
	h1_["ptmu1_ex_pre"]  -> Fill(goodMuons.at(0).getPt(), presel_wt) ;
	h1_["etamu1_ex_pre"]  -> Fill(goodMuons.at(0).getEta(), presel_wt) ;
	h1_["phimu1_ex_pre"]  -> Fill(goodMuons.at(0).getPhi(), presel_wt) ;
	h1_["energymu1_ex_pre"]  -> Fill(goodMuons.at(0).getP4().Energy(), presel_wt) ;
	
	h1_["ptmu2_ex_pre"]  -> Fill(goodMuons.at(1).getPt(), presel_wt) ;
	h1_["etamu2_ex_pre"]  -> Fill(goodMuons.at(1).getEta(), presel_wt) ;
	h1_["phimu2_ex_pre"]  -> Fill(goodMuons.at(1).getPhi(), presel_wt) ;
	h1_["energymu2_ex_pre"]  -> Fill(goodMuons.at(1).getP4().Energy(), presel_wt) ;

	h1_["dphi_mu1_jet1_pre"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMuons.at(0).getP4()), presel_wt);
	h1_["dphi_mu1_jet2_pre"] -> Fill( (goodAK4Jets.at(1).getP4()).DeltaPhi(goodMuons.at(0).getP4()), presel_wt);
	h1_["dphi_mu1_jet3_pre"] -> Fill( (goodAK4Jets.at(2).getP4()).DeltaPhi(goodMuons.at(0).getP4()), presel_wt);
	
	h1_["dphi_mu2_jet1_pre"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMuons.at(1).getP4()), presel_wt);
	h1_["dphi_mu2_jet2_pre"] -> Fill( (goodAK4Jets.at(1).getP4()).DeltaPhi(goodMuons.at(1).getP4()), presel_wt);
	h1_["dphi_mu2_jet3_pre"] -> Fill( (goodAK4Jets.at(2).getP4()).DeltaPhi(goodMuons.at(1).getP4()), presel_wt);



      }
      else if (zdecayMode_ == "zelel" ) {
      */
    //	h1_["ptel1_ex_pre"]  -> Fill(goodElectrons.at(0).getPt(), presel_wt) ;
    //	h1_["etael1_ex_pre"]  -> Fill(goodElectrons.at(0).getEta(), presel_wt) ;
    //	h1_["phiel1_ex_pre"]  -> Fill(goodElectrons.at(0).getPhi(), presel_wt) ;
    //	h1_["energyel1_ex_pre"]  -> Fill(goodElectrons.at(0).getP4().Energy(), presel_wt) ;
	
    //	h1_["ptel2_ex_pre"]  -> Fill(goodElectrons.at(1).getPt(), presel_wt) ;
    //	h1_["etael2_ex_pre"]  -> Fill(goodElectrons.at(1).getEta(), presel_wt) ;
    //	h1_["phiel2_ex_pre"]  -> Fill(goodElectrons.at(1).getPhi(), presel_wt) ;
    //	h1_["energyel2_ex_pre"]  -> Fill(goodElectrons.at(1).getP4().Energy(), presel_wt) ;
	
    //	h1_["dphi_el1_jet1_pre"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodElectrons.at(0).getP4()), presel_wt);
    //	h1_["dphi_el1_jet2_pre"] -> Fill( (goodAK4Jets.at(1).getP4()).DeltaPhi(goodElectrons.at(0).getP4()), presel_wt);
    //	h1_["dphi_el1_jet3_pre"] -> Fill( (goodAK4Jets.at(2).getP4()).DeltaPhi(goodElectrons.at(0).getP4()), presel_wt);
	
    //	h1_["dphi_el2_jet1_pre"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodElectrons.at(1).getP4()), presel_wt);
    //	h1_["dphi_el2_jet2_pre"] -> Fill( (goodAK4Jets.at(1).getP4()).DeltaPhi(goodElectrons.at(1).getP4()), presel_wt);
    //	h1_["dphi_el2_jet3_pre"] -> Fill( (goodAK4Jets.at(2).getP4()).DeltaPhi(goodElectrons.at(1).getP4()), presel_wt);


	/*	
	for  (unsigned int iele=0; iele<goodElectrons.size(); ++iele){
	  ofstream myfile;
	  // if(-2.0< goodElectrons.at(iele).getEta() && goodElectrons.at(iele).getEta()<-1.6){
	  myfile.open ("EleInfo.txt", std::ios::app);
	  myfile << iele+1<< " Leading Electron information \n";
	  for (auto izll : zll) {
	    myfile << "Zpt         = "<< izll.getPt()<<"\n";
	  }
	  myfile <<"Dr(e1,e2)    = "<< goodElectrons.at(0).getP4().DeltaR(goodElectrons.at(1).getP4())<<"\n";
	  myfile << "Pt          = "<< goodElectrons.at(iele).getPt()<<"\n";
	  myfile << "Eta         = "<<goodElectrons.at(iele).getEta()<<"\n";
	  myfile << "Phi         = "<<goodElectrons.at(iele).getPhi()<<"\n";
	  myfile << "Energy      = "<<goodElectrons.at(iele).getP4().Energy()<<"\n";
	  

	  
	  myfile << "Iso03       = "<< goodElectrons.at(iele).getIso03()<<"\n";
	  myfile << "dEtaInSeed  = "<<  goodElectrons.at(iele).getdEtaInSeed()<<"\n";
	  myfile << "dPhiIn      = "<<goodElectrons.at(iele).getdPhiIn()<<"\n";
	  myfile << "Dz          = "<<goodElectrons.at(iele).getDz()<<"\n";
	  myfile << "Dxy         = "<<goodElectrons.at(iele).getDxy()<<"\n";
	  myfile << "scEta       = "<<goodElectrons.at(iele).getscEta()<<"\n";
	  myfile << "full5x5siee = "<<goodElectrons.at(iele).getfull5x5siee()<<"\n";
	  myfile << "H/E         = "<<goodElectrons.at(iele).getHoE()<<"\n";
	  myfile << "1/E - 1/P   = "<<goodElectrons.at(iele).getooEmooP()<<"\n";
	  myfile << "missHits    = "<<goodElectrons.at(iele).getmissHits()<<"\n";
	  myfile << "hasMatchedConVeto = "<<goodElectrons.at(iele).gethasMatchedConVeto()<<"\n";
		  
	  myfile << "***********************\n";
	  myfile.close();
	}
	  
	*/
       
    //	for  (unsigned int iele=0; iele<goodElectrons.size(); ++iele){
    //	  float scEta = goodElectrons.at(iele).getscEta();
    //	  if(fabs(scEta) <= 1.479){
    //	    h1_["Eta_EB_el_ex"]-> Fill(goodElectrons.at(iele).getEta(), evtwt);
    //	    h1_["Iso03_EB_el_ex"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
    //	    h1_["dEtaInSeed_EB_el_ex"]->Fill(goodElectrons.at(iele).getdEtaInSeed(), evtwt);
    //	    h1_["dPhiIn_EB_el_ex"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
    //	    h1_["Dz_EB_el_ex"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
    //	    h1_["Dxy_EB_el_ex"]->Fill(goodElectrons.at(iele).getDxy(), evtwt);
    //	    h1_["SCETA_EB_el_ex"]->Fill(goodElectrons.at(iele).getscEta(), evtwt);
    //	    h1_["Full5x5siee_EB_el_ex"]->Fill(goodElectrons.at(iele).getfull5x5siee(), evtwt);
    //	    h1_["HoE_EB_el_ex"]->Fill(goodElectrons.at(iele).getHoE(), evtwt);
    //	    h1_["ooEmooP_EB_el_ex"]->Fill(goodElectrons.at(iele).getooEmooP(), evtwt);
    //	    h1_["missHits_EB_el_ex"]->Fill(goodElectrons.at(iele).getmissHits(), evtwt);
    //	    h1_["conveto_EB_el_ex"]->Fill(goodElectrons.at(iele).gethasMatchedConVeto(), evtwt);
	      
	      
    //	  }
  //	  else if  (fabs(scEta) > 1.479 && fabs(scEta) < 2.4){
  //	    h1_["Eta_EE_el_ex"]->Fill(goodElectrons.at(iele).getEta(), evtwt);
    //	    h1_["Iso03_EE_el_ex"]->Fill(goodElectrons.at(iele).getIso03(), evtwt);
    //	    h1_["dEtaInSeed_EE_el_ex"]->Fill(goodElectrons.at(iele).getdEtaInSeed(), evtwt);
    //	    h1_["dPhiIn_EE_el_ex"]->Fill(goodElectrons.at(iele).getdPhiIn(), evtwt);
    //	    h1_["Dz_EE_el_ex"]->Fill(goodElectrons.at(iele).getDz(), evtwt);
    //	    h1_["Dxy_EE_el_ex"]->Fill(goodElectrons.at(iele).getDxy(), evtwt);
    //	    h1_["SCETA_EE_el_ex"]->Fill(goodElectrons.at(iele).getscEta(), evtwt);
    //	    h1_["Full5x5siee_EE_el_ex"]->Fill(goodElectrons.at(iele).getfull5x5siee(), evtwt);
    //	    h1_["HoE_EE_el_ex"]->Fill(goodElectrons.at(iele).getHoE(), evtwt);
    //	    h1_["ooEmooP_EE_el_ex"]->Fill(goodElectrons.at(iele).getooEmooP(), evtwt);
    //	    h1_["missHits_EE_el_ex"]->Fill(goodElectrons.at(iele).getmissHits(), evtwt);
    //	    h1_["conveto_EE_el_ex"]->Fill(goodElectrons.at(iele).gethasMatchedConVeto(), evtwt);
    //	    
    //	  }
    //	}






      


    //}
    //  else {
	
    //	for (auto izll : zll) {
    //	  h1_["ptz2_ex_pre"] -> Fill(izll.getPt(), presel_wt);
    //	}
	
	//else if (zdecayMode_ == "zelel" ) {
	
    //	h1_["dr_elel_ex_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), presel_wt);
	  //	}



    // }
    // }
    h1_["nak4_pre"] -> Fill(goodAK4Jets.size(), presel_wt) ;
    h1_["ht_pre"] -> Fill(htak4.getHT(), presel_wt);
    h1_["st_pre"] -> Fill(ST, presel_wt) ;   

    //// Ak4 jet plots
    for(int j=0; j<3; ++j){
      h1_[Form("ptak4jet%d_pre", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), presel_wt) ; 
      h1_[Form("etaak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), presel_wt) ;
      h1_[Form("cvsak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), presel_wt) ;
      h1_[Form("massak4jet%d_pre", j+1)] -> Fill(goodAK4Jets.at(j).getMass(), presel_wt) ;

    }
    h1_["phi_jet1MET_pre"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), presel_wt);

    //// npv
    h1_["npv_noweight_pre"] -> Fill(npv, *h_evtwtGen.product()); 
    h1_["npv_pre"] -> Fill(npv, presel_wt);

    //// met
    h1_["met_pre"] -> Fill(goodMet.at(0).getFullPt(), presel_wt);
    h1_["met1_pre"] -> Fill(goodMet.at(0).getFullPt(),  presel_wt);
    h1_["metPhi_pre"] -> Fill(goodMet.at(0).getFullPhi(), presel_wt);

    //// Lepton specfic properties
    if ( zdecayMode_ == "zmumu" ){       
      for(int l=0; l<2; ++l){
        h1_["pt_"+lep+Form("%d_pre", l+1)]  -> Fill(goodMuons.at(l).getPt(), presel_wt) ; 
        h1_["eta_"+lep+Form("%d_pre", l+1)]  -> Fill(goodMuons.at(l).getEta(), presel_wt) ; 
      } 
      h1_["Iso04_mu1_pre"]->Fill(goodMuons.at(0).getIso04(), presel_wt);
      h1_["Iso04_mu2_pre"]->Fill(goodMuons.at(1).getIso04(), presel_wt);
      for  (unsigned int imu=0; imu<goodMuons.size(); ++imu){
        h1_["Iso04_mu_pre"]->Fill(goodMuons.at(imu).getIso04(), presel_wt);
        h1_["Dz_mu_pre"]->Fill(goodMuons.at(imu).getDz(), presel_wt);
        h1_["Dxy_mu_pre"]->Fill(goodMuons.at(imu).getDxy(), presel_wt);
        h1_["IsGlobalMuon_mu_pre"]->Fill(goodMuons.at(imu).getIsGlobalMuon(), presel_wt);
        h1_["IsPFMuon_mu_pre"]->Fill(goodMuons.at(imu).getIsPFMuon(), presel_wt);
        h1_["GlbTrkNormChi2_mu_pre"]->Fill(goodMuons.at(imu).getGlbTrkNormChi2(), presel_wt);
        h1_["NumberValidMuonHits_mu_pre"]->Fill(goodMuons.at(imu).getNumberValidMuonHits(), presel_wt);
        h1_["NumberMatchedStations_mu_pre"]->Fill(goodMuons.at(imu).getNumberMatchedStations(), presel_wt);
        h1_["NumberValidPixelHits_mu_pre"]->Fill(goodMuons.at(imu).getNumberValidPixelHits(), presel_wt);
        h1_["NumberTrackerLayers_mu_pre"]->Fill(goodMuons.at(imu).getNumberTrackerLayers(), presel_wt);
      }

      h1_["dr_mumu_pre"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), presel_wt);
    }
    else if (zdecayMode_ == "zelel" ) {
      for(int l=0; l<2; ++l){
        h1_["pt_"+lep+Form("%d_pre", l+1)]   -> Fill(goodElectrons.at(l).getPt(), presel_wt) ; 
        h1_["eta_"+lep+Form("%d_pre", l+1)]  -> Fill(goodElectrons.at(l).getEta(), presel_wt) ; 
      } 

      h1_["dr_elel_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), presel_wt);
      if(additionalPlots_) fillAdditionalPlots(goodElectrons, presel_wt);
    }

    //// Z pt
    for (auto izll : zll) h1_["pt_z"+lep+lep+"_pre"] -> Fill(izll.getPt(), presel_wt) ;
    if(zdecayMode_ == "zmumu") {
      h1_["dr_mu1jet1_pre"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodAK4Jets.at(0).getP4()),  presel_wt );
      h1_["dr_mu1jet2_pre"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodAK4Jets.at(1).getP4()),  presel_wt );
      h1_["dr_mu1jet3_pre"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodAK4Jets.at(2).getP4()),  presel_wt );

      h1_["dr_mu2jet1_pre"]-> Fill( (goodMuons.at(1).getP4()).DeltaR(goodAK4Jets.at(0).getP4()),  presel_wt );
      h1_["dr_mu2jet2_pre"]-> Fill( (goodMuons.at(1).getP4()).DeltaR(goodAK4Jets.at(1).getP4()),  presel_wt );
      h1_["dr_mu2jet3_pre"]-> Fill( (goodMuons.at(1).getP4()).DeltaR(goodAK4Jets.at(2).getP4()),  presel_wt );      


      // double dr=1000.0;
      vector<std::pair<vlq::Jet , double> > pfEn , pfEn1;

     for (unsigned i=0; i<goodAK4Jets.size();i++){
       std::pair<vlq::Jet , double> newPair = make_pair(goodAK4Jets.at(i), (goodMuons.at(0).getP4()).DeltaR(goodAK4Jets.at(i).getP4()));
   
       //       cout << "dr, pt,eta = " << (goodMuons.at(0).getP4()).DeltaR(goodAK4Jets.at(i).getP4())<<","<< goodAK4Jets.at(i).getPt()<<","<< goodAK4Jets.at(i).getEta() <<endl;

       pfEn.push_back(newPair);
     }

     //std::sort(pfEn.begin(), pfEn.end(),my_compare );
     //pfEn.sort(my_compare ); 
     double drmin =1000.0;
     vlq::Jet minjet;
     for(auto& cand : pfEn ){
       if (cand.second <drmin) { drmin = cand.second; minjet=cand.first; };
       }
     // cout << "can dr,pt,eta =" << drmin<<","<< minjet.getPt()<<","<< minjet.getEta() <<endl;

     //    cout << "******* W.r.t 2n lepton "<<endl;   
   

     for (unsigned i=0; i<goodAK4Jets.size();i++){
       std::pair<vlq::Jet , double> newPair1 = make_pair(goodAK4Jets.at(i), (goodMuons.at(1).getP4()).DeltaR(goodAK4Jets.at(i).getP4()));

       //  cout << "dr, pt,eta = " << (goodMuons.at(1).getP4()).DeltaR(goodAK4Jets.at(i).getP4())<<","<< goodAK4Jets.at(i).getPt()<<","<< goodAK4Jets.at(i).getEta() <<endl;                     

       pfEn1.push_back(newPair1);
     }
 
     double drmin1 =1000.0;
     vlq::Jet minjet1;
     for(auto& cand : pfEn1 ){
       if (cand.second <drmin1) { drmin1 = cand.second; minjet1=cand.first; };
     }
     //  cout << "can dr,pt,eta =" << drmin1<<","<< minjet1.getPt()<<","<< minjet1.getEta() <<endl;                                                                                                  

 
     //  cout << "%%%%%%%% end event %%%%%%%% " << endl;
     h1_["dr_mu1minjet_pre"]-> Fill( drmin,  presel_wt );
     h1_["dr_mu2minjet_pre"]-> Fill( drmin1,  presel_wt );


    }
    else if (zdecayMode_ == "zelel" ) {
      // h1_["dr_elel_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()),  presel_wt );
      h1_["dr_el1jet1_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(0).getP4()),  presel_wt );
      h1_["dr_el1jet2_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(1).getP4()),  presel_wt );
      h1_["dr_el1jet3_pre"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(2).getP4()),  presel_wt );

      h1_["dr_el2jet1_pre"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(0).getP4()),  presel_wt );
      h1_["dr_el2jet2_pre"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(1).getP4()),  presel_wt );
      h1_["dr_el2jet3_pre"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(2).getP4()),  presel_wt );


      vector<std::pair<vlq::Jet , double> > pfEn , pfEn1;

      for (unsigned i=0; i<goodAK4Jets.size();i++){
	std::pair<vlq::Jet , double> newPair = make_pair(goodAK4Jets.at(i), (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(i).getP4()));

	//	cout << "dr, pt,eta = " << (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(i).getP4())<<","<< goodAK4Jets.at(i).getPt()<<","<< goodAK4Jets.at(i).getEta() <<endl;

	pfEn.push_back(newPair);
      }

      //std::sort(pfEn.begin(), pfEn.end(),my_compare );                                                                                                                                                                                                           
      //pfEn.sort(my_compare );                                                                                                                                                                                                                                    
      double drmin =1000.0;
      vlq::Jet minjet;
      for(auto& cand : pfEn ){
	if (cand.second <drmin) { drmin = cand.second; minjet=cand.first; };
      }
      // cout << "can dr,pt,eta =" << drmin<<","<< minjet.getPt()<<","<< minjet.getEta() <<endl;

      // cout << "******* W.r.t 2n lepton "<<endl;

      for (unsigned i=0; i<goodAK4Jets.size();i++){
	std::pair<vlq::Jet , double> newPair1 = make_pair(goodAK4Jets.at(i), (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(i).getP4()));

	//	cout << "dr, pt,eta = " << (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(i).getP4())<<","<< goodAK4Jets.at(i).getPt()<<","<< goodAK4Jets.at(i).getEta() <<endl;

	pfEn1.push_back(newPair1);
      }

      double drmin1 =1000.0;
      vlq::Jet minjet1;
      for(auto& cand : pfEn1 ){
	if (cand.second <drmin1) { drmin1 = cand.second; minjet1=cand.first; };
      }
      //  cout << "can dr,pt,eta =" << drmin1<<","<< minjet1.getPt()<<","<< minjet1.getEta() <<endl;


      //  cout << "%%%%%%%% end event %%%%%%%% " << endl;
      h1_["dr_el1minjet_pre"]-> Fill( drmin,  presel_wt );
      h1_["dr_el2minjet_pre"]-> Fill( drmin1,  presel_wt );

      //      cout << "%%%%%%%% end event %%%%%%%% " << endl;




    }


    if (goodWTaggedJets.size() > 0) {
      h1_["Wptleading_pre"] -> Fill((goodWTaggedJets.at(0)).getPt(), presel_wt) ;
      h1_["Wetaleading_pre"] -> Fill((goodWTaggedJets.at(0)).getEta(), presel_wt) ;
      h1_["Wprunedleading_pre"] -> Fill((goodWTaggedJets.at(0)).getPrunedMass(), presel_wt) ;
    }
    if (goodWTaggedJets.size() > 1) {
      h1_["Wpt2nd_pre"] -> Fill((goodWTaggedJets.at(1)).getPt(), presel_wt) ;
      h1_["Weta2nd_pre"] -> Fill((goodWTaggedJets.at(1)).getEta(), presel_wt) ;
      h1_["Wpruned2nd_pre"] -> Fill((goodWTaggedJets.at(1)).getPrunedMass(), presel_wt) ;
    }

    //   for  (unsigned int i=0; i<goodAK8Jets.size(); i++){
    //  cout << " AK8d jet size = " << goodAK8Jets.size() <<endl;
    //  cout << " i th index = " << goodAK8Jets.at(i).getIndex()<<endl;
    
    //  }
    // cout << " ** end AK8 **"<< endl;

    for  (unsigned int i=0; i<goodWTaggedJets.size(); i++){
      // cout << " W tagged jet size = " << goodWTaggedJets.size() <<endl;
      // cout << " i th index = " << goodWTaggedJets.at(i).getIndex()<<endl;
      h1_["Wpt_pre"] -> Fill((goodWTaggedJets.at(i)).getPt(), presel_wt) ;
      h1_["Weta_pre"] -> Fill((goodWTaggedJets.at(i)).getEta(), presel_wt) ;
      h1_["Wpruned_pre"] -> Fill((goodWTaggedJets.at(i)).getPrunedMass(), presel_wt) ;
    }
    // cout << " ** end W  **"<< endl;



    /*   
    if (goodwTaggedJets.size() > 0) {
      h1_["Wptleading1_pre"] -> Fill((goodwTaggedJets.at(0)).getPt(), presel_wt) ;
      h1_["Wetaleading1_pre"] -> Fill((goodwTaggedJets.at(0)).getEta(), presel_wt) ;
      h1_["Wprunedleading1_pre"] -> Fill((goodwTaggedJets.at(0)).getPrunedMass(), presel_wt) ;
    }
    if (goodwTaggedJets.size() > 1) {
      h1_["Wpt2nd1_pre"] -> Fill((goodwTaggedJets.at(1)).getPt(), presel_wt) ;
      h1_["Weta2nd1_pre"] -> Fill((goodwTaggedJets.at(1)).getEta(), presel_wt) ;
      h1_["Wpruned2nd1_pre"] -> Fill((goodwTaggedJets.at(1)).getPrunedMass(), presel_wt) ;
    }
   
    for(unsigned i=0; i < goodWTaggedJets.size(); i++){
      h1_["ptsmear1_pre"] -> Fill((goodWTaggedJets.at(i)).getptsmear() , presel_wt) ;
      h1_["massCorr1_pre"] -> Fill((goodWTaggedJets.at(i)).getmassCorr () , presel_wt) ;
      h1_["masssmear1_pre"] -> Fill((goodWTaggedJets.at(i)).getmasssmear() , presel_wt) ;
    }
   
    for(unsigned i=0; i < goodwTaggedJets.size(); i++){
      h1_["ptsmear2_pre"] -> Fill((goodwTaggedJets.at(i)).getptsmear() , presel_wt) ;
      h1_["massCorr2_pre"] -> Fill((goodwTaggedJets.at(i)).getmassCorr () , presel_wt) ;
      h1_["masssmear2_pre"] -> Fill((goodwTaggedJets.at(i)).getmasssmear() , presel_wt) ;
      h1_["masssmear21_pre"] -> Fill((goodwTaggedJets.at(i)).getmasssmear() , presel_wt) ;

    }
    */


    if (goodHTaggedJets.size() > 0) {
      h1_["Hptleading_pre"] -> Fill((goodHTaggedJets.at(0)).getPt(), presel_wt) ;
      h1_["Hetaleading_pre"] -> Fill((goodHTaggedJets.at(0)).getEta(), presel_wt) ;
      h1_["Hprunedleading_pre"] -> Fill((goodHTaggedJets.at(0)).getPrunedMass(), presel_wt) ;
    }
    if (goodHTaggedJets.size() > 1) {
      h1_["Hpt2nd_pre"] -> Fill((goodHTaggedJets.at(1)).getPt(), presel_wt) ;
      h1_["Heta2nd_pre"] -> Fill((goodHTaggedJets.at(1)).getEta(), presel_wt) ;
      h1_["Hpruned2nd_pre"] -> Fill((goodHTaggedJets.at(1)).getPrunedMass(), presel_wt) ;
    }

    for  (unsigned int i=0; i<goodHTaggedJets.size(); i++){
      h1_["Hpt_pre"] -> Fill((goodHTaggedJets.at(i)).getPt(), presel_wt) ;
      h1_["Heta_pre"] -> Fill((goodHTaggedJets.at(i)).getEta(), presel_wt) ;
      h1_["Hpruned_pre"] -> Fill((goodHTaggedJets.at(i)).getPrunedMass(), presel_wt) ;
      // cout << " H tagged jet size = " << goodHTaggedJets.size() <<endl;
      //  cout << " i th index = " << goodHTaggedJets.at(i).getIndex()<<endl;
      
    }
    // cout << " ** end H **"<< endl;


    if (goodTopTaggedJets.size() > 0) {
      h1_["Topptleading_pre"] -> Fill((goodTopTaggedJets.at(0)).getPt(), presel_wt) ;
      h1_["Topetaleading_pre"] -> Fill((goodTopTaggedJets.at(0)).getEta(), presel_wt) ;
      h1_["Topsoftdropleading_pre"] -> Fill((goodTopTaggedJets.at(0)).getSoftDropMass(), presel_wt) ;
    }
    if (goodTopTaggedJets.size() > 1) {
      h1_["Toppt2nd_pre"] -> Fill((goodTopTaggedJets.at(1)).getPt(), presel_wt) ;
      h1_["Topeta2nd_pre"] -> Fill((goodTopTaggedJets.at(1)).getEta(), presel_wt) ;
      h1_["Topsoftdrop2nd_pre"] -> Fill((goodTopTaggedJets.at(1)).getSoftDropMass(), presel_wt) ;
    }

    for  (unsigned int i=0; i<goodTopTaggedJets.size(); i++){
      h1_["Toppt_pre"] -> Fill((goodTopTaggedJets.at(i)).getPt(), presel_wt) ;
      h1_["Topeta_pre"] -> Fill((goodTopTaggedJets.at(i)).getEta(), presel_wt) ;
      h1_["Topsoftdrop_pre"] -> Fill((goodTopTaggedJets.at(i)).getSoftDropMass(), presel_wt) ;
      //  cout << " Top tagged jet size = " << goodTopTaggedJets.size() <<endl;
      // cout << " i th index = " << goodTopTaggedJets.at(i).getIndex()<<endl;

   }

    // cout << " ** end Top **"<< endl;




    h1_["nak8_pre"] -> Fill(goodAK8Jets.size(), evtwt) ;
    h1_["nwjet_pre"] -> Fill(goodWTaggedJets.size(), evtwt) ;
    h1_["nhjet_pre"] -> Fill(goodHTaggedJets.size(), evtwt) ;
    h1_["ntjet_pre"] -> Fill(goodTopTaggedJets.size(), evtwt) ;





    //b jets pt plots                                                                                             
    if (goodAK4Jets.at(0).getPt() > 100){
      if (goodAK4Jets.at(1).getPt() > 50){
        if ( goodBTaggedAK4Jets.size() > 0 ) {

          h1_["ptbjetleading_pre"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
          h1_["etabjetleading_pre"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;

          //  if (goodBTaggedAK4Jets.size() >= 1){                                                                
          // h1_["ptbjetsubleading_pre"] -> Fill(goodBTaggedAK4Jets.at(1).getPt(), evtwt) ;                       
          //  h1_["etabjetsubleading_pre"] -> Fill(goodBTaggedAK4Jets.at(1).getEta(), evtwt) ;                    
          // }                                                                                                    
        }
      }
    }


    //========================================================
    // Preselection done, proceeding with control selections
    //========================================================

    if ( goodBTaggedAK4Jets.size() == 0) {
      for (auto izll : zll) {
	//	h1_["nob_pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ;
      }

      h1_["nob_ht"] ->Fill(htak4.getHT(), evtwt);
      h1_["nob_st"] ->Fill(ST, evtwt);

      if (ST > 1000 && goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){
	h1_["nob_1000_ht"] ->Fill(htak4.getHT(), evtwt);
	h1_["nob_1000_st"] ->Fill(ST, evtwt);


	for (unsigned int i=0; i< goodAK8Jets.size();i++){
	  h1_["ptak8_nob1000"] -> Fill((goodAK8Jets.at(i)).getPt(), evtwt) ;
	}

	for (unsigned int i=0; i< goodWTaggedJets.size();i++){
          h1_["ptW_nob1000"] -> Fill((goodWTaggedJets.at(i)).getPt(), evtwt) ;
        }
	for (unsigned int i=0; i< goodHTaggedJets.size();i++){
          h1_["ptH_nob1000"] -> Fill((goodHTaggedJets.at(i)).getPt(), evtwt) ;
        }
	for (unsigned int i=0; i< goodTopTaggedJets.size();i++){
          h1_["ptT_nob1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
        }


	for (auto izll : zll) {
	  if (zdecayMode_ == "zelel" ) h1_["nob_1000_pt_zelel"] -> Fill(izll.getPt(), evtwt) ;
	  else if (zdecayMode_ == "zmumu" ) h1_["nob_1000_pt_zmumu"] -> Fill(izll.getPt(), evtwt) ;
	}

      }
      if (goodMet.at(0).getFullPt()<60){
	//	h1_["nbjets_met_0btagcnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
	//	h1_["ht_met_0btagcnt"] ->Fill(htak4.getHT(), evtwt);
	//	h1_["st_met_0btagcnt"] ->Fill(ST, evtwt);
      }
    }

    if( goodBTaggedAK4Jets.size() == 0){
      h1_["nak4_0b1"] -> Fill(goodAK4Jets.size(), evtwt) ;
      for(int j=0; j<3; ++j){
	h1_[Form("ptak4jet%d_0b1", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
      }

      if(goodAK4Jets.at(0).getPt() > 100){// not sure if we should also cut on ST                                                        
	h1_["nak4_0b2"] -> Fill(goodAK4Jets.size(), evtwt) ;
	for(int j=0; j<3; ++j){
	  h1_[Form("ptak4jet%d_0b2", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
	}

	if (goodAK4Jets.at(1).getPt() > 50){
	  h1_["nak4_0b3"] -> Fill(goodAK4Jets.size(), evtwt) ;
	  for(int j=0; j<3; ++j){
	    h1_[Form("ptak4jet%d_0b3", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
	  }
	  if ( ST <= 1000) {h1_["nak4_0b3_l1000"] -> Fill(goodAK4Jets.size(), evtwt) ;}
	  else if ( ST > 1000){ h1_["nak4_0b3_g1000"] -> Fill(goodAK4Jets.size(), evtwt) ;}
	  
	}
      }
    }
      
    else if (goodBTaggedAK4Jets.size() > 0){
      for (auto izll : zll) {
	//	h1_["b_pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ;
	//	h1_["b_st"] ->Fill(ST, evtwt);
      }
      h1_["1b_ht"] ->Fill(htak4.getHT(), evtwt);
      h1_["1b_st"] ->Fill(ST, evtwt);
      if (ST > 1000 && goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){
        h1_["1b_1000_ht"] ->Fill(htak4.getHT(), evtwt);
        h1_["1b_1000_st"] ->Fill(ST, evtwt);

      }

      if (goodMet.at(0).getFullPt()<60){
	//	h1_["nbjets_met_1btagcnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
	//	h1_["ht_met_1btagcnt"] ->Fill(htak4.getHT(), evtwt);
	//	h1_["st_met_1btagcnt"] ->Fill(ST, evtwt);
      }
    }
    if (goodMet.at(0).getFullPt()<60){
      // h1_["nbjets_met_cnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
      // h1_["lowmet_ht"] ->Fill(htak4.getHT(), evtwt);
      //  h1_["lowmet_st"] ->Fill(ST, evtwt);
    }





    //fill control plots
    if ( goodBTaggedAK4Jets.size() > 0 && ST < STMaxControl_  ) {
      for (auto izll : zll) {
        h1_["mass_z"+lep+lep+"_cnt"] -> Fill(izll.getMass(), evtwt) ;  
	h1_["mass_Z"+lep+lep+"_cnt"] -> Fill(izll.getMass(), evtwt) ;
        h1_["pt_z"+lep+lep+"_cnt"] -> Fill(izll.getPt(), evtwt) ; 
      }
      h1_["nak4_cnt"] -> Fill(goodAK4Jets.size(), evtwt) ;
      h1_["ht_cnt"] -> Fill(htak4.getHT(), evtwt) ;
      h1_["st_cnt"] -> Fill(ST, evtwt) ;   
      h1_["ht1_cnt"] -> Fill(htak4.getHT(), evtwt) ;
      h1_["st1_cnt"] -> Fill(ST, evtwt) ;
      h1_["npv_noweight_cnt"] -> Fill(npv, *h_evtwtGen.product()); 
      h1_["npv_cnt"] -> Fill(npv, evtwt);
      h1_["met_cnt"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["met1_cnt"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["metPhi_cnt"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

      //lepton specfic properties
      if ( zdecayMode_ == "zmumu" ){       
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_mumu_cnt"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
      }
      else if (zdecayMode_ == "zelel" ) {
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d_cnt", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d_cnt", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_elel_cnt"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
      }

      //ak4 jet plots
      for(int j=0; j<3; ++j){
        h1_[Form("ptak4jet%d_cnt", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
        h1_[Form("etaak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
        h1_[Form("cvsak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
	h1_[Form("massak4jet%d_cnt", j+1)] -> Fill(goodAK4Jets.at(j).getMass(), evtwt) ;
      }
      
      //cout<< " nbjets_cnt = " << goodBTaggedAK4Jets.size() <<endl;
      h1_["phi_jet1MET_cnt"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);
      h1_["nbjets_cnt"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;


      if (goodWTaggedJets.size() > 0) {
        h1_["Wptleading_cnt"] -> Fill((goodWTaggedJets.at(0)).getPt(), evtwt) ;
        h1_["Wetaleading_cnt"] -> Fill((goodWTaggedJets.at(0)).getEta(), evtwt) ;
        h1_["Wprunedleading_cnt"] -> Fill((goodWTaggedJets.at(0)).getPrunedMass(), evtwt) ;
      }
      if (goodWTaggedJets.size() > 1) {
        h1_["Wpt2nd_cnt"] -> Fill((goodWTaggedJets.at(1)).getPt(), evtwt) ;
        h1_["Weta2nd_cnt"] -> Fill((goodWTaggedJets.at(1)).getEta(), evtwt) ;
        h1_["Wpruned2nd_cnt"] -> Fill((goodWTaggedJets.at(1)).getPrunedMass(), evtwt) ;
      }

      for  (unsigned int i=0; i<goodWTaggedJets.size(); i++){
	h1_["Wpt_cnt"] -> Fill((goodWTaggedJets.at(i)).getPt(), evtwt) ;
        h1_["Weta_cnt"] -> Fill((goodWTaggedJets.at(i)).getEta(), evtwt) ;
        h1_["Wpruned_cnt"] -> Fill((goodWTaggedJets.at(i)).getPrunedMass(), evtwt) ;
      }


      /*   
      if (goodwTaggedJets.size() > 0) {
	h1_["Wptleading1_cnt"] -> Fill((goodwTaggedJets.at(0)).getPt(), evtwt) ;
	h1_["Wetaleading1_cnt"] -> Fill((goodwTaggedJets.at(0)).getEta(), evtwt) ;
	h1_["Wprunedleading1_cnt"] -> Fill((goodwTaggedJets.at(0)).getPrunedMass(), evtwt) ;
      }
      if (goodwTaggedJets.size() > 1) {
	h1_["Wpt2nd1_cnt"] -> Fill((goodwTaggedJets.at(1)).getPt(), evtwt) ;
	h1_["Weta2nd1_cnt"] -> Fill((goodwTaggedJets.at(1)).getEta(), evtwt) ;
	h1_["Wpruned2nd1_cnt"] -> Fill((goodwTaggedJets.at(1)).getPrunedMass(), evtwt) ;
      }

      */
      if (goodHTaggedJets.size() > 0) {
	h1_["Hptleading_cnt"] -> Fill((goodHTaggedJets.at(0)).getPt(), evtwt) ;
	h1_["Hetaleading_cnt"] -> Fill((goodHTaggedJets.at(0)).getEta(), evtwt) ;
	h1_["Hprunedleading_cnt"] -> Fill((goodHTaggedJets.at(0)).getPrunedMass(), evtwt) ;
      }
      if (goodHTaggedJets.size() > 1) {
	h1_["Hpt2nd_cnt"] -> Fill((goodHTaggedJets.at(1)).getPt(), evtwt) ;
	h1_["Heta2nd_cnt"] -> Fill((goodHTaggedJets.at(1)).getEta(), evtwt) ;
	h1_["Hpruned2nd_cnt"] -> Fill((goodHTaggedJets.at(1)).getPrunedMass(), evtwt) ;
      }

      for  (unsigned int i=0; i<goodHTaggedJets.size(); i++){
	h1_["Hpt_cnt"] -> Fill((goodHTaggedJets.at(i)).getPt(), evtwt) ;
        h1_["Heta_cnt"] -> Fill((goodHTaggedJets.at(i)).getEta(), evtwt) ;
        h1_["Hpruned_cnt"] -> Fill((goodHTaggedJets.at(i)).getPrunedMass(), evtwt) ;
      }


      if (goodTopTaggedJets.size() > 0) {
	h1_["Topptleading_cnt"] -> Fill((goodTopTaggedJets.at(0)).getPt(), evtwt) ;
	h1_["Topetaleading_cnt"] -> Fill((goodTopTaggedJets.at(0)).getEta(), evtwt) ;
	h1_["Topsoftdropleading_cnt"] -> Fill((goodTopTaggedJets.at(0)).getSoftDropMass(), evtwt) ;
      }
      if (goodTopTaggedJets.size() > 1) {
	h1_["Toppt2nd_cnt"] -> Fill((goodTopTaggedJets.at(1)).getPt(), evtwt) ;
	h1_["Topeta2nd_cnt"] -> Fill((goodTopTaggedJets.at(1)).getEta(), evtwt) ;
	h1_["Topsoftdrop2nd_cnt"] -> Fill((goodTopTaggedJets.at(1)).getSoftDropMass(), evtwt) ;
      }
      for  (unsigned int i=0; i<goodTopTaggedJets.size(); i++){
	h1_["Toppt_cnt"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
        h1_["Topeta_cnt"] -> Fill((goodTopTaggedJets.at(i)).getEta(), evtwt) ;
        h1_["Topsoftdrop_cnt"] -> Fill((goodTopTaggedJets.at(i)).getSoftDropMass(), evtwt) ;
      }



      h1_["nak8_cnt"] -> Fill(goodAK8Jets.size(), evtwt) ;
      h1_["nwjet_cnt"] -> Fill(goodWTaggedJets.size(), evtwt) ;
      h1_["nhjet_cnt"] -> Fill(goodHTaggedJets.size(), evtwt) ;
      h1_["ntjet_cnt"] -> Fill(goodTopTaggedJets.size(), evtwt) ;



    } //// Control region 

    /*
    //fill categorization control plots
    if (ST < 1000){
      if ( goodBTaggedAK4Jets.size() == 0) {
	h1_["ht_st1000_e0b"] ->Fill(htak4.getHT(), evtwt);
	h1_["st_st1000_e0b"] ->Fill(ST, evtwt);

	for (auto izll : zll) {
	  h1_["pt_z"+lep+lep+"_st1000_e0b"] -> Fill(izll.getPt(), evtwt) ;
	}
	if ( zdecayMode_ == "zmumu" ){
	  for(int l=0; l<2; ++l){
	    h1_["pt_"+lep+Form("%d_st1000_e0b", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
	    h1_["eta_"+lep+Form("%d_st1000_e0b", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
	  }
	}
	else if (zdecayMode_ == "zelel" ) {
	  for(int l=0; l<2; ++l){
	    h1_["pt_"+lep+Form("%d_st1000_e0b", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
	    h1_["eta_"+lep+Form("%d_st1000_e0b", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
	  }
	}
	h1_["met_st1000_e0b"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	
	for(int j=0; j<3; ++j){
	  h1_[Form("ptak4jet%d_st1000_e0b", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
	}
      }
      else if( goodBTaggedAK4Jets.size() > 0) {
	h1_["ht_st1000_e1b"] ->Fill(htak4.getHT(), evtwt);
        h1_["st_st1000_e1b"] ->Fill(ST, evtwt);

	for (auto izll : zll) {
          h1_["pt_z"+lep+lep+"_st1000_e1b"] -> Fill(izll.getPt(), evtwt) ;
        }
        if ( zdecayMode_ == "zmumu" ){
          for(int l=0; l<2; ++l){
            h1_["pt_"+lep+Form("%d_st1000_e1b", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
            h1_["eta_"+lep+Form("%d_st1000_e1b", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
          }
        }
        else if (zdecayMode_ == "zelel" ) {
          for(int l=0; l<2; ++l){
            h1_["pt_"+lep+Form("%d_st1000_e1b", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
            h1_["eta_"+lep+Form("%d_st1000_e1b", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
          }
        }
        h1_["met_st1000_e1b"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	
        for(int j=0; j<3; ++j){
          h1_[Form("ptak4jet%d_st1000_e1b", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;

        }
        h1_["ptbjetleading_st1000_e1b"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
	if( goodBTaggedAK4Jets.size() > 1) {
	  h1_["ptbjetsubleading_st1000_e1b"] -> Fill(goodBTaggedAK4Jets.at(1).getPt(), evtwt) ;
	}
	if ( goodBTaggedAK4Jets.size() == 1){
	  h1_["ht_st1000_1b"] ->Fill(htak4.getHT(), evtwt);
	  h1_["st_st1000_1b"] ->Fill(ST, evtwt);

	  for (auto izll : zll) {
	    h1_["pt_z"+lep+lep+"_st1000_1b"] -> Fill(izll.getPt(), evtwt) ;
	  }
	  if ( zdecayMode_ == "zmumu" ){
	    for(int l=0; l<2; ++l){
	      h1_["pt_"+lep+Form("%d_st1000_1b", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
	      h1_["eta_"+lep+Form("%d_st1000_1b", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
	    }
	  }
	  else if (zdecayMode_ == "zelel" ) {
	    for(int l=0; l<2; ++l){
	      h1_["pt_"+lep+Form("%d_st1000_1b", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
	      h1_["eta_"+lep+Form("%d_st1000_1b", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
	    }
	  }
	  h1_["met_st1000_1b"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	  
	  for(int j=0; j<3; ++j){
	    h1_[Form("ptak4jet%d_st1000_1b", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;

	  }
	  h1_["ptbjetleading_st1000_1b"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >=2 ){
	  h1_["ht_st1000_2b"] ->Fill(htak4.getHT(), evtwt);
          h1_["st_st1000_2b"] ->Fill(ST, evtwt);

	  for (auto izll : zll) {
	    h1_["pt_z"+lep+lep+"_st1000_2b"] -> Fill(izll.getPt(), evtwt) ;
	  }
	  if ( zdecayMode_ == "zmumu" ){
	    for(int l=0; l<2; ++l){
	      h1_["pt_"+lep+Form("%d_st1000_2b", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
	      h1_["eta_"+lep+Form("%d_st1000_2b", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
	    }
	  }
	  else if (zdecayMode_ == "zelel" ) {
	    for(int l=0; l<2; ++l){
	      h1_["pt_"+lep+Form("%d_st1000_2b", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
	      h1_["eta_"+lep+Form("%d_st1000_2b", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
	    }
	  }
	  h1_["met_st1000_2b"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	  
	  for(int j=0; j<3; ++j){
	    h1_[Form("ptak4jet%d_st1000_2b", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
	  
	  }
	  h1_["ptbjetleading_st1000_2b"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
	  h1_["ptbjetsubleading_st1000_2b"] -> Fill(goodBTaggedAK4Jets.at(1).getPt(), evtwt) ;
	}

      }

    }

    */
 
    vlq::JetCollection goodAK4Jetscleaned,ak4matchedak8, ak4nonmatched1, ak4nonmatched2, ak4nonmatched3;
    /*
    for (unsigned i=0; i<goodAK4Jets.size(); i++) {
      goodAK4Jetscleaned.push_back(goodAK4Jets[i]);
    }
    // cout << " %%%%%%% start %%%%%%%%% " << endl;
    // for (unsigned i=0; i<goodAK4Jetscleaned.size(); i++) {
    //  cout << " before cleaning = " << goodAK4Jetscleaned.at(i).getMass() <<","<< goodAK4Jetscleaned.at(i).getPt()<< ","<<goodAK4Jetscleaned.at(i).getEta()<<","<<goodAK4Jetscleaned.at(i).getPhi() << endl;
    // }

    for (unsigned i=0; i<goodHTaggedJets.size(); i++) {
      for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
	//	cout << " dr   = " << goodHTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4())<<endl;
	if (goodHTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){
	  //  cout <<" matched jet = " << goodAK4Jetscleaned.at(j).getMass() <<","<< goodAK4Jetscleaned.at(j).getPt()<< ","<<goodAK4Jetscleaned.at(j).getEta()<<","<<goodAK4Jetscleaned.at(j).getPhi() << endl;

	  ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
	  goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
	}
	else{
	  ak4nonmatched1.push_back(goodAK4Jetscleaned.at(j));
	}
      }
    }

    //   cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% cleaned %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" <<endl;
    //    for (unsigned i=0; i<goodAK4Jetscleaned.size(); i++) {
    //	cout << " before cleaning = " << goodAK4Jetscleaned.at(i).getMass() <<","<< goodAK4Jetscleaned.at(i).getPt()<< ","<<goodAK4Jetscleaned.at(i).getEta()<<","<<goodAK4Jetscleaned.at(i).getPhi() << endl;
    //  }
      // cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& END &&&&&&&&&&&&&&&&&&&&&&&&&&&" << endl;

    for (unsigned i=0; i<goodWTaggedJets.size(); i++) {
      for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
	if (goodWTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){
	  ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
	  goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
	} 
	else{
          ak4nonmatched2.push_back(goodAK4Jetscleaned.at(j));
        }
      }
    }

    for (unsigned i=0; i<goodTopTaggedJets.size(); i++) {
      for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
	if (goodTopTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){
	  ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
	  goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
	}
	else{
	  ak4nonmatched3.push_back(goodAK4Jetscleaned.at(j));
	}
      }
    }
    */

    if (categorize_){
      if ( goodBTaggedAK4Jets.size() == 0  && goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50  && ST < 1000) {
	h1_["nak4_0cat"] -> Fill(goodAK4Jets.size(), evtwt) ;
        h1_["ht1_0cat"] -> Fill(htak4.getHT(), evtwt) ;
        h1_["st1_0cat"] -> Fill(ST, evtwt) ;
	if ( zdecayMode_ == "zmumu" ){
	  h1_["dr_mumu_0cat"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
	}
	else if (zdecayMode_ == "zelel" ) {
	  h1_["dr_elel_0cat"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
        }


      }
      
      if ( goodBTaggedAK4Jets.size() >=1  && goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50  && ST < 1000) {

	//if ( goodBTaggedAK4Jets.size() ==0  && goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50  && ST > 1000) {
	
	for (auto izll : zll) {
	  h1_["mass_z"+lep+lep+"_cat"] -> Fill(izll.getMass(), evtwt) ;
	  h1_["mass_Z"+lep+lep+"_cat"] -> Fill(izll.getMass(), evtwt) ;
	  h1_["pt_z"+lep+lep+"_cat"] -> Fill(izll.getPt(), evtwt) ;
	}
	h1_["nak4_cat"] -> Fill(goodAK4Jets.size(), evtwt) ;
	h1_["ht_cat"] -> Fill(htak4.getHT(), evtwt) ;
	h1_["st_cat"] -> Fill(ST, evtwt) ;
	h1_["ht1_cat"] -> Fill(htak4.getHT(), evtwt) ;
	h1_["st1_cat"] -> Fill(ST, evtwt) ;
	h1_["npv_noweight_cat"] -> Fill(npv, *h_evtwtGen.product());
	h1_["npv_cat"] -> Fill(npv, evtwt);
	h1_["met_cat"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	h1_["met1_cat"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	h1_["metPhi_cat"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

	//lepton specfic properties                                                                                                                                                       
	if ( zdecayMode_ == "zmumu" ){
	  for(int l=0; l<2; ++l){
	    h1_["pt_"+lep+Form("%d_cat", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
	    h1_["eta_"+lep+Form("%d_cat", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
	  }
	  h1_["dr_mumu_cat"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
	}
	else if (zdecayMode_ == "zelel" ) {
	  for(int l=0; l<2; ++l){
	    h1_["pt_"+lep+Form("%d_cat", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
	    h1_["eta_"+lep+Form("%d_cat", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
	  }
	  h1_["dr_elel_cat"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
	}

	//ak4 jet plots                                                                                                                                                                   
	for(int j=0; j<3; ++j){
	  h1_[Form("ptak4jet%d_cat", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
	  h1_[Form("etaak4jet%d_cat", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
	  h1_[Form("cvsak4jet%d_cat", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
	  h1_[Form("massak4jet%d_cat", j+1)] -> Fill(goodAK4Jets.at(j).getMass(), evtwt) ;
	}
	h1_["phi_jet1MET_cat"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);
	
	//	h1_["nbjets_cat"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;













	//	if (goodWTaggedJets.size() > 0) {
	//  h1_["Wptleading_cat"] -> Fill((goodWTaggedJets.at(0)).getPt(), evtwt) ;
	//  h1_["Wetaleading_cat"] -> Fill((goodWTaggedJets.at(0)).getEta(), evtwt) ;
	//  h1_["Wprunedleading_cat"] -> Fill((goodWTaggedJets.at(0)).getPrunedMass(), evtwt) ;
	//}
	//if (goodWTaggedJets.size() > 1) {
	//  h1_["Wpt2nd_cat"] -> Fill((goodWTaggedJets.at(1)).getPt(), evtwt) ;
	//  h1_["Weta2nd_cat"] -> Fill((goodWTaggedJets.at(1)).getEta(), evtwt) ;
	//  h1_["Wpruned2nd_cat"] -> Fill((goodWTaggedJets.at(1)).getPrunedMass(), evtwt) ;
	//	}

	//	for  (unsigned int i=0; i<goodWTaggedJets.size(); i++){
	//  h1_["Wpt_cat"] -> Fill((goodWTaggedJets.at(i)).getPt(), evtwt) ;
	// h1_["Weta_cat"] -> Fill((goodWTaggedJets.at(i)).getEta(), evtwt) ;
        //  h1_["Wpruned_cat"] -> Fill((goodWTaggedJets.at(i)).getPrunedMass(), evtwt) ;
	// }

	/*
	if (goodwTaggedJets.size() > 0) {
          h1_["Wptleading1_cat"] -> Fill((goodwTaggedJets.at(0)).getPt(), evtwt) ;
          h1_["Wetaleading1_cat"] -> Fill((goodwTaggedJets.at(0)).getEta(), evtwt) ;
          h1_["Wprunedleading1_cat"] -> Fill((goodwTaggedJets.at(0)).getPrunedMass(), evtwt) ;
        }
        if (goodwTaggedJets.size() > 1) {
          h1_["Wpt2nd1_cat"] -> Fill((goodwTaggedJets.at(1)).getPt(), evtwt) ;
          h1_["Weta2nd1_cat"] -> Fill((goodwTaggedJets.at(1)).getEta(), evtwt) ;
          h1_["Wpruned2nd1_cat"] -> Fill((goodwTaggedJets.at(1)).getPrunedMass(), evtwt) ;
        }
	*/
	


	//	if (goodHTaggedJets.size() > 0) {
	//  h1_["Hptleading_cat"] -> Fill((goodHTaggedJets.at(0)).getPt(), evtwt) ;
	//  h1_["Hetaleading_cat"] -> Fill((goodHTaggedJets.at(0)).getEta(), evtwt) ;
	//  h1_["Hprunedleading_cat"] -> Fill((goodHTaggedJets.at(0)).getPrunedMass(), evtwt) ;
	//	}
	//	if (goodHTaggedJets.size() > 1) {
	//  h1_["Hpt2nd_cat"] -> Fill((goodHTaggedJets.at(1)).getPt(), evtwt) ;
	//  h1_["Heta2nd_cat"] -> Fill((goodHTaggedJets.at(1)).getEta(), evtwt) ;
	//  h1_["Hpruned2nd_cat"] -> Fill((goodHTaggedJets.at(1)).getPrunedMass(), evtwt) ;
	//	}

	//	for  (unsigned int i=0; i<goodHTaggedJets.size(); i++){
	//  h1_["Hpt_cat"] -> Fill((goodHTaggedJets.at(i)).getPt(), evtwt) ;
        //  h1_["Heta_cat"] -> Fill((goodHTaggedJets.at(i)).getEta(), evtwt) ;
        //  h1_["Hpruned_cat"] -> Fill((goodHTaggedJets.at(i)).getPrunedMass(), evtwt) ;
	//	}

	//	if (goodTopTaggedJets.size() > 0) {
	//  h1_["Topptleading_cat"] -> Fill((goodTopTaggedJets.at(0)).getPt(), evtwt) ;
	//  h1_["Topetaleading_cat"] -> Fill((goodTopTaggedJets.at(0)).getEta(), evtwt) ;
	//  h1_["Topsoftdropleading_cat"] -> Fill((goodTopTaggedJets.at(0)).getSoftDropMass(), evtwt) ;
	//	}
	//if (goodTopTaggedJets.size() > 1) {
	//  h1_["Toppt2nd_cat"] -> Fill((goodTopTaggedJets.at(1)).getPt(), evtwt) ;
	//  h1_["Topeta2nd_cat"] -> Fill((goodTopTaggedJets.at(1)).getEta(), evtwt) ;
	//  h1_["Topsoftdrop2nd_cat"] -> Fill((goodTopTaggedJets.at(1)).getSoftDropMass(), evtwt) ;
	//	}

	//for  (unsigned int i=0; i<goodTopTaggedJets.size(); i++){
	//  h1_["Toppt_cat"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
        //  h1_["Topeta_cat"] -> Fill((goodTopTaggedJets.at(i)).getEta(), evtwt) ;
        //  h1_["Topsoftdrop_cat"] -> Fill((goodTopTaggedJets.at(i)).getSoftDropMass(), evtwt) ;
	// }


	h1_["nak8_cat"] -> Fill(goodAK8Jets.size(), evtwt) ;
	h1_["nwjet_cat"] -> Fill(goodWTaggedJets.size(), evtwt) ;
	h1_["nhjet_cat"] -> Fill(goodHTaggedJets.size(), evtwt) ;
	h1_["ntjet_cat"] -> Fill(goodTopTaggedJets.size(), evtwt) ;

	//// Ak4 jet plots                                                                                                                                                    

	//	if (goodAK4Jets.size()>0){
	//  h1_["ak4pt1_cattest1"]  -> Fill(goodAK4Jets.at(0).getPt(), evtwt) ;
        //  h1_["ak4eta1_cattest1"] -> Fill(goodAK4Jets.at(0).getEta(), evtwt) ;
        //  h1_["ak4mass1_cattest1"] -> Fill(goodAK4Jets.at(0).getMass(), evtwt) ;
	//	}
	//	if (goodAK4Jets.size()>1){
        //  h1_["ak4pt2_cattest1"]  -> Fill(goodAK4Jets.at(1).getPt(), evtwt) ;
        //  h1_["ak4eta2_cattest1"] -> Fill(goodAK4Jets.at(1).getEta(), evtwt) ;
        //  h1_["ak4mass2_cattest1"] -> Fill(goodAK4Jets.at(1).getMass(), evtwt) ;
	//	}
	//	h1_["ht1_cat"] -> Fill(htak4.getHT(), evtwt) ;
	//h1_["st1_cat"] -> Fill(ST, evtwt) ;
	h1_["nak41_cat"] -> Fill(goodAK4Jets.size(), evtwt) ;


	for (unsigned i=0; i<goodAK4Jets.size(); i++) {
	  goodAK4Jetscleaned.push_back(goodAK4Jets[i]);
	}
	for (unsigned i=0; i<goodHTaggedJets.size(); i++) {
	  for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
	    if (goodHTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){	
	      ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
	      goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
	    }
	    else{
	      ak4nonmatched1.push_back(goodAK4Jetscleaned.at(j));
	    }
	  }
	}
	for (unsigned i=0; i<goodWTaggedJets.size(); i++) {
	  for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
	    if (goodWTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){
	      ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
	      goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
	    }
	    else{
	      ak4nonmatched2.push_back(goodAK4Jetscleaned.at(j));
	    }
	  }
	}

	for (unsigned i=0; i<goodTopTaggedJets.size(); i++) {
	  for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
	    if (goodTopTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){
	      ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
	      goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
	    }
	    else{
	      ak4nonmatched3.push_back(goodAK4Jetscleaned.at(j));
	    }
	  }
	}

	h1_["nbjets_cat"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;

	if ( goodBTaggedAK4Jets.size() > 0 ) {
          h1_["ptbjetleading_cat"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
          h1_["etabjetleading_cat"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;
	}
	
	if (goodBTaggedAK4Jets.size() > 1){                                                                                                                                       
	  h1_["ptbjetsubleading_cat"] -> Fill(goodBTaggedAK4Jets.at(1).getPt(), evtwt) ;                                                                                              
          h1_["etabjetsubleading_cat"] -> Fill(goodBTaggedAK4Jets.at(1).getEta(), evtwt) ;                                                                                           
	}                                                                                                                                                                           
        
	for (unsigned i=0; i<goodBTaggedAK4Jets.size(); i++) {
	  h1_["ptbjet_cat"] -> Fill(goodBTaggedAK4Jets.at(i).getPt(), evtwt) ;
          h1_["etabjet_cat"] -> Fill(goodBTaggedAK4Jets.at(i).getEta(), evtwt) ;
        }





	if (goodAK4Jetscleaned.size()>0){
          h1_["ak4pt1_cattest2"]  -> Fill(goodAK4Jetscleaned.at(0).getPt(), evtwt) ;
          h1_["ak4eta1_cattest2"] -> Fill(goodAK4Jetscleaned.at(0).getEta(), evtwt) ;
          h1_["ak4mass1_cattest2"] -> Fill(goodAK4Jetscleaned.at(0).getMass(), evtwt) ;
	}
        if (goodAK4Jetscleaned.size()>1){
          h1_["ak4pt2_cattest2"]  -> Fill(goodAK4Jetscleaned.at(1).getPt(), evtwt) ;
          h1_["ak4eta2_cattest2"] -> Fill(goodAK4Jetscleaned.at(1).getEta(), evtwt) ;
          h1_["ak4mass2_cattest2"] -> Fill(goodAK4Jetscleaned.at(1).getMass(), evtwt) ;
        }

	h1_["ht2_cat"] -> Fill(htak4.getHT(), evtwt) ;
        h1_["st2_cat"] -> Fill(ST, evtwt) ;
	h1_["nak42_cat"] -> Fill(goodAK4Jetscleaned.size(), evtwt) ;




      
	// HCandidate cut(boosted)                                                                             
	for (unsigned i=0; i<goodHTaggedJets.size(); i++) {
	  Hb.push_back(goodHTaggedJets.at(i)); 
	}
   
	//HCandidates (nonboosted)
	HCandsProducer h;
	if (goodAK4Jetscleaned.size() >= 2){
	  h.operator()(goodAK4Jetscleaned.size(), 2, goodAK4Jetscleaned,H)  ;
	}
	//Zcandidate cut (boosted)                                     
	//      for (unsigned i=0; i<goodWTaggedJets.size(); i++) {
	for (unsigned i=0; i<goodWTaggedJets.size(); i++) { 	
	  ZB.push_back(goodWTaggedJets.at(i));   
	}

	//Z candidate cut (non boosted)                                                                                                            
	ZCandsProducer z;
	if (goodAK4Jetscleaned.size() >= 2){
	  z.operator()(goodAK4Jetscleaned.size(), 2, goodAK4Jetscleaned,Z) ;
	}   

	//Top Candidates (Category D)                                                                                                                                                       
	for (unsigned i=0; i<goodTopTaggedJets.size(); i++) {
	  D.push_back(goodTopTaggedJets.at(i));
	}
	TopCandsProducer top,w;
           
	// Category BC                                                                                                                                       
	if (goodWTaggedJets.size()>0 && goodAK4Jetscleaned.size()>0){               
	  for (unsigned i=0; i<goodWTaggedJets.size(); i++) {
	    for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
	      TLorentzVector bc1;
	      bc1= goodWTaggedJets.at(i).getP4()+goodAK4Jetscleaned.at(j).getP4();
	      if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
		vlq::Candidate bc2(bc1);
		h1_["dr_Wb_cnt"]-> Fill( (goodWTaggedJets.at(i).getP4()).DeltaR(goodAK4Jetscleaned.at(j).getP4()), evtwt);
		h1_["dphi_Wb_cnt"]-> Fill( (goodWTaggedJets.at(i).getP4()).DeltaPhi(goodAK4Jetscleaned.at(j).getP4()), evtwt);
		W.push_back(goodWTaggedJets.at(i));
		B.push_back(goodAK4Jetscleaned.at(j));
		BC.push_back(bc2);
	      }
         
	    }
	  }
	}
        
	// category A                                                                                                                                                                       
	if (goodAK4Jetscleaned.size() >= 3){
	  top.operator()(goodAK4Jetscleaned.size(), 3, goodAK4Jetscleaned,tops) ;
	}


	//plots
	//ak4 jets matched to ak8
	for (unsigned i=0; i<ak4matchedak8.size(); i++) {
	  h1_["mass_ak4matchedak8"] -> Fill(ak4matchedak8.at(i).getMass(), evtwt) ;
	  h1_["pt_ak4matchedak8"] -> Fill(ak4matchedak8.at(i).getPt(), evtwt) ;
	}
	h1_["nak4matchedak8"] -> Fill(ak4matchedak8.size(), evtwt) ;


	// HCandidate cut(boosted)
	for (unsigned i=0; i<Hb.size(); i++) {
	  h1_["H_mass_b_cnt"] -> Fill(Hb.at(i).getPrunedMass(), evtwt) ;
	  h1_["H_Pt_b_cnt"] -> Fill(Hb.at(i).getPt(), evtwt) ;
	}
	h1_["nHcandidatejets_b_cnt"] -> Fill(Hb.size(), evtwt) ;
      
	// if(Hb.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}
      
	//Hcandidate cut(nonboosted)
	for (unsigned i=0; i<H.size(); i++) {
	  h1_["H_mass_nb_cnt"] -> Fill(H.at(i).getMass(), evtwt) ;
	  h1_["H_Pt_nb_cnt"] -> Fill(H.at(i).getPt(), evtwt) ;
	}
      
	h1_["nHcandidatejets_nb_cnt"] -> Fill(H.size(), evtwt) ;

	//if(H.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
      
	double  nHcandidates=0.0;
	double  nHcandidates1=0.0;
	if(Hb.size()>0 ||H.size()>0){
	  //h1_["cutflow"] -> Fill(16, evtwt) ;
	  nHcandidates = Hb.size()+H.size();
	  //nHcandidates = H.size();
	  //nHcandidates = Hb.size(); 
	  h1_["nHcandidatejets_cnt"] -> Fill(nHcandidates, evtwt) ;
      
	}

	nHcandidates1 = Hb.size()+H.size();
	// nHcandidates1 = H.size();

	h1_["nHcandidatejets1_cnt"] -> Fill(nHcandidates1, evtwt) ;

	//Zcandidate cut (boosted)
	for (unsigned i=0; i<ZB.size(); i++) {     
	  h1_["Z_mass_a_cnt"] -> Fill(ZB.at(i).getPrunedMass(), evtwt) ;
	  h1_["Z_Pt_a_cnt"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
	}                                                                                                                                                                                       
	h1_["nzcandidatejets_a_cnt"] -> Fill(ZB.size(), evtwt) ;                                                                                                                          
   
	//if(ZB.size()>0)  { h1_["cutflow"] -> Fill(10, evtwt) ;}     
      
	//Z candidate cut (non boosted)                                                                                                                                                         
	for (unsigned i=0; i<Z.size(); i++) {
	  h1_["Z_mass_b_cnt"] -> Fill(Z.at(i).getMass(), evtwt) ;
	  h1_["Z_Pt_b_cnt"] -> Fill(Z.at(i).getPt(), evtwt) ;
	}

	h1_["nzcandidatejets_b_cnt"] -> Fill(Z.size(), evtwt) ;
      
	// if(Z.size()>0)  { h1_["cutflow"] -> Fill(11, evtwt) ;}
 
	double nzcandidates=0.0;
	double nzcandidates1=0.0;
	if ( ZB.size() || Z.size()){
	  //h1_["cutflow"] -> Fill(12, evtwt);
	  nzcandidates = ZB.size()+ Z.size();
	  //nzcandidates = Z.size();
	  // nzcandidates = ZB.size();
	  h1_["nzcandidatejets_tot_cnt"] -> Fill(nzcandidates, evtwt) ;
	
	}
	nzcandidates1 = ZB.size()+ Z.size();
	// nzcandidates1 = Z.size();
	h1_["nzcandidatejets1_tot_cnt"] -> Fill(nzcandidates1, evtwt) ;

      
	// Category D    
	for (unsigned i=0; i<D.size(); i++) {
	  h1_["top_mass_d_cnt"] -> Fill(D.at(i).getSoftDropMass(), evtwt) ;
	  h1_["top_Pt_d_cnt"] -> Fill(D.at(i).getPt(), evtwt) ;
	}
	h1_["ntopcandidatejets_d_cnt"] -> Fill(D.size(), evtwt) ;
      
	// if(D.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}


	// Category BC
	// for (unsigned i=0; i<W.size(); i++) {
	//	h1_["W_mass_bc_cnt"] -> Fill(W.at(i).getMass(), evtwt) ;
	// }
	//  h1_["nWcandidatejets_bc_cnt"] -> Fill(W.size(), evtwt) ;
      
	//for (unsigned i=0; i<B.size(); i++) {
	//	h1_["lightjet_mass_bc_cnt"] -> Fill(B.at(i).getMass(), evtwt) ; 
	// }
	//h1_["nlightjetcandidatejets_bc_cnt"] -> Fill(B.size(), evtwt) ;
      
	for (unsigned i=0; i<W.size(); i++) {
	  h1_["W_mass_bc_cnt"] -> Fill(W.at(i).getPrunedMass(), evtwt) ;
	}
	h1_["nWcandidatejets_bc_cnt"] -> Fill(W.size(), evtwt) ;

	for (unsigned i=0; i<B.size(); i++) {                                                                                                                             
	  h1_["lightjet_mass_bc_cnt"] -> Fill(B.at(i).getMass(), evtwt) ;                                                                                                
	}                                                                                                                                                                                
	h1_["nlightjetcandidatejets_bc_cnt"] -> Fill(B.size(), evtwt) ;    

	for (unsigned i=0; i<BC.size(); i++) {
	  h1_["top_mass_bc_cnt"] -> Fill(BC.at(i).getMass(), evtwt) ;
	  h1_["top_Pt_bc_cnt"] -> Fill(BC.at(i).getPt(), evtwt) ;
	}
	h1_["ntopcandidatejets_bc_cnt"] -> Fill(BC.size(), evtwt) ;
      
	// if(BC.size()>0)  { h1_["cutflow"] -> Fill(14, evtwt) ;}


	for (unsigned i=0; i<tops.size(); i++) {
	  h1_["top_mass_a_cnt"] -> Fill(tops.at(i).getMass(), evtwt) ;
	  h1_["top_Pt_a_cnt"] -> Fill(tops.at(i).getPt(), evtwt) ;
	}
	h1_["ntopcandidatejets_a_cnt"] -> Fill(tops.size(), evtwt) ;
      
	// if(tops.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
      
	double  ntopcandidates=0.0;
	double  ntopcandidates1=0.0;
	if(D.size() || BC.size() ||tops.size()){
	  //   h1_["cutflow"] -> Fill(16, evtwt) ;
	  ntopcandidates = D.size()+BC.size()+tops.size();
	  //ntopcandidates = tops.size();    
	  //ntopcandidates = D.size()+BC.size();
	  h1_["ntopcandidatejets_cnt"] -> Fill(ntopcandidates, evtwt) ;  
	
	} 
	// ntopcandidates1 = D.size()+BC.size()+tops.size();
	ntopcandidates1 = D.size()+BC.size()+tops.size(); 
	// ntopcandidates1 = tops.size();

	h1_["ntopcandidatejets1_cnt"] -> Fill(ntopcandidates1, evtwt) ;
      
	//cout << " pass 2 " <<endl;
	//Z and top corelations and ST tempelates

	//h1_["st_cnt"] -> Fill(ST,evtwt);
	h1_["cutflow1"] -> Fill(1, evtwt) ;
	if (goodBTaggedAK4Jets.size() == 1){
	  h1_["cutflow1"] -> Fill(2, evtwt) ;   
	  //h1_["st_cnt1b"] -> Fill(ST,evtwt);
	}
      
	if (goodBTaggedAK4Jets.size() >=2){
	  //h1_["st_cnt2b"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(3, evtwt) ;
	}
	//   cout<<" eventweight after multiplying  **** background ######################################## 22222222222222 =  " << evtwt<<endl;

	//n,Z,H,B
	//(1)
	if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
	  //h1_["st_cntT1Z1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(4, evtwt) ;
	  if (nHcandidates >= 1.0){
	    // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	    h1_["cutflow1"] -> Fill(8, evtwt) ;

	    h1_["st_cntT1Z1H1"] -> Fill(ST, evtwt) ;
	    h1_["ht_cntT1Z1H1"] -> Fill(htak4.getHT(), evtwt) ;




	    //  for (auto izll : zll) {
	    //  h1_["pt_z"+lep+lep+"_cntT1Z1H1"] -> Fill(izll.getPt(), evtwt) ;
	    //  h1_["eta_z"+lep+lep+"_cntT1Z1H1"] -> Fill(izll.getEta(), evtwt) ;
	    //  h1_["phi_z"+lep+lep+"_cntT1Z1H1"] -> Fill(izll.getPhi(), evtwt) ;

	    //  }
	    // if ( zdecayMode_ == "zmumu" ){
	    //  h1_["dr_mumu_cntT1Z1H1"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
	    //  for(int l=0; l<2; ++l){
	    //	h1_["pt_"+lep+Form("%d_cntT1Z1H1", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
	    //	h1_["eta_"+lep+Form("%d_cntT1Z1H1", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
	    //	h1_["phi_"+lep+Form("%d_cntT1Z1H1", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
	    //  }
	    // }
	    // else if (zdecayMode_ == "zelel" ) {
	    //  h1_["dr_elel_cntT1Z1H1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
	    //  for(int l=0; l<2; ++l){
	    //		h1_["pt_"+lep+Form("%d_cntT1Z1H1", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
	    //	h1_["eta_"+lep+Form("%d_cntT1Z1H1", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
	    //	h1_["phi_"+lep+Form("%d_cntT1Z1H1", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
	    //  }
	    // }
	    // h1_["met_cntT1Z1H1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	    // h1_["metPhi_cntT1Z1H1"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
	    //  h1_["metEta_cntT1Z1H1"] -> Fill(goodMet.at(0).getFullEta(), evtwt); //do not uncomment this line                                                                                              

	    //  for(int j=0; j<3; ++j){
	    //  h1_[Form("ptak4jet%d_cntT1Z1H1", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
	    //  h1_[Form("etaak4jet%d_cntT1Z1H1", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
	    //  h1_[Form("phiak4jet%d_cntT1Z1H1", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;

	    //   }







	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow2"] -> Fill(1, evtwt) ;
	      h1_["st_cntT1Z1H1b1"] -> Fill(ST, evtwt) ;
	      h1_["ht_cntT1Z1H1b1"] -> Fill(htak4.getHT(), evtwt) ;
	      h1_["nbjets_cntT1Z1H1b1"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;

	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT1Z1H1b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT1Z1H1b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



	      /*
	      if (D.size() >=1.0){
		if (ZB.size() >=1.0){
		  if (Hb.size() >=1.0){
		    h1_["cutflow10"] -> Fill(1, evtwt);
		    h1_["st_cntD1ZB1Hb1b1"] -> Fill(ST, evtwt) ;
		    //h1_["ht_cntD1ZB1Hb1b1"] -> Fill(htak4.getHT(), evtwt) ;
		  }
		  else if(H.size() >=1.0){
                    h1_["cutflow10"] -> Fill(2, evtwt);
                    h1_["st_cntD1ZB1H1b1"] -> Fill(ST, evtwt) ;
		    // h1_["ht_cntD1ZB1H1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
		}
		else if (Z.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow10"] -> Fill(3, evtwt);
                    h1_["st_cntD1Z1Hb1b1"] -> Fill(ST, evtwt) ;
		    //   h1_["ht_cntD1Z1Hb1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow10"] -> Fill(4, evtwt);
                    h1_["st_cntD1Z1H1b1"] -> Fill(ST, evtwt) ;
		    //   h1_["ht_cntD1Z1H1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
		}
	      }

	      else if (BC.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow10"] -> Fill(5, evtwt);
                    h1_["st_cntBC1ZB1Hb1b1"] -> Fill(ST, evtwt) ;
		    //   h1_["ht_cntBC1ZB1Hb1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow10"] -> Fill(6, evtwt);
                    h1_["st_cntBC1ZB1H1b1"] -> Fill(ST, evtwt) ;
		    //  h1_["ht_cntBC1ZB1H1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow10"] -> Fill(7, evtwt);
                    h1_["st_cntBC1Z1Hb1b1"] -> Fill(ST, evtwt) ;
		    //  h1_["ht_cntBC1Z1Hb1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow10"] -> Fill(8, evtwt);
                    h1_["st_cntBC1Z1H1b1"] -> Fill(ST, evtwt) ;
                    h1_["ht_cntBC1Z1H1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                }
              }

	      else if (tops.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow10"] -> Fill(9, evtwt);
                    h1_["st_cntt1ZB1Hb1b1"] -> Fill(ST, evtwt) ;
		    //  h1_["ht_cntt1ZB1Hb1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow10"] -> Fill(10, evtwt);
                    h1_["st_cntt1ZB1H1b1"] -> Fill(ST, evtwt) ;
		    //   h1_["ht_cntt1ZB1H1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow10"] -> Fill(11, evtwt);
                    h1_["st_cntt1Z1Hb1b1"] -> Fill(ST, evtwt) ;
		    // h1_["ht_cntt1Z1Hb1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow10"] -> Fill(12, evtwt);
                    h1_["st_cntt1Z1H1b1"] -> Fill(ST, evtwt) ;
		    // h1_["ht_cntt1Z1H1b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                }
              }
	      */


	         for (auto izll : zll) {
	      	h1_["pt_z"+lep+lep+"_cntT1Z1H1b1"] -> Fill(izll.getPt(), evtwt) ;
	      	h1_["eta_z"+lep+lep+"_cntT1Z1H1b1"] -> Fill(izll.getEta(), evtwt) ;
	      	h1_["phi_z"+lep+lep+"_cntT1Z1H1b1"] -> Fill(izll.getPhi(), evtwt) ;
		
	            }
	      if ( zdecayMode_ == "zmumu" ){
	      	h1_["dr_mumu_cntT1Z1H1b1"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
	      //	for(int l=0; l<2; ++l){
	      //	  h1_["pt_"+lep+Form("%d_cntT1Z1H1b1", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
	      //	  h1_["eta_"+lep+Form("%d_cntT1Z1H1b1", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
	      //	  h1_["phi_"+lep+Form("%d_cntT1Z1H1b1", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
	      //	}
	      }
	      else if (zdecayMode_ == "zelel" ) {
	      	h1_["dr_elel_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
	      	h1_["dr_el1jet1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(0).getP4()), evtwt );
	      	h1_["dr_el1jet2_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(1).getP4()), evtwt );
	      	h1_["dr_el1jet3_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(2).getP4()), evtwt );

	      	h1_["dr_el2jet1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(0).getP4()), evtwt );
	      	h1_["dr_el2jet2_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(1).getP4()), evtwt );
	      	h1_["dr_el2jet3_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(2).getP4()), evtwt );

		h1_["dr_el1bjet1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodBTaggedAK4Jets.at(0).getP4()), evtwt );
	      	h1_["dr_el2bjet1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodBTaggedAK4Jets.at(0).getP4()), evtwt );


	      	if (Hb.size() >0){
	      	  h1_["dr_el1Hb1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(Hb.at(0).getP4()), evtwt );
	      	  h1_["dr_el2Hb1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(Hb.at(0).getP4()), evtwt );
	      	}

	      	if (ZB.size() >0){
		  h1_["dr_el1Zb1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(ZB.at(0).getP4()), evtwt );
	      	  h1_["dr_el2Zb1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(ZB.at(0).getP4()), evtwt );

		}

		if (D.size() >0){
		  h1_["dr_el1tb1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(D.at(0).getP4()), evtwt );
	      	  h1_["dr_el2tb1_cntT1Z1H1b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(D.at(0).getP4()), evtwt );
		}

	      

		for(int l=0; l<2; ++l){
	      	  h1_["pt_"+lep+Form("%d_cntT1Z1H1b1", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
	      	  h1_["eta_"+lep+Form("%d_cntT1Z1H1b1", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
	      	  h1_["phi_"+lep+Form("%d_cntT1Z1H1b1", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
		}
	      }
	      // h1_["met_cntT1Z1H1b1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	      // h1_["metPhi_cntT1Z1H1b1"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
	      //  h1_["metEta_cntT1Z1H1b1"] -> Fill(goodMet.at(0).getFullEta(), evtwt);
	      
	      //  for(int j=0; j<3; ++j){
	      //	h1_[Form("ptak4jet%d_cntT1Z1H1b1", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
	      //h1_[Form("etaak4jet%d_cntT1Z1H1b1", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
	      //h1_[Form("phiak4jet%d_cntT1Z1H1b1", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;	 
		
	      //      }
	      // h1_["ptbjetleading_cntT1Z1H1b1"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
	      // h1_["etabjetleading_cntT1Z1H1b1"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;
	      // h1_["phibjetleading_cntT1Z1H1b1"] -> Fill(goodBTaggedAK4Jets.at(0).getPhi(), evtwt) ;
	      
	      // if(goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){
	      //		h1_["st_cntT1Z1H1b1_A"] -> Fill(ST, evtwt) ;
	      //		h1_["ht_cntT1Z1H1b1_A"] -> Fill(htak4.getHT(), evtwt) ;
	      // }
	      
	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow2"] -> Fill(2, evtwt) ;
	      h1_["st_cntT1Z1H1b2"] -> Fill(ST, evtwt) ;
	      h1_["ht_cntT1Z1H1b2"] -> Fill(htak4.getHT(), evtwt) ;
	      h1_["nbjets_cntT1Z1H1b2"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;

	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT1Z1H1b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT1Z1H1b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



	      /*
	      if (D.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow11"] -> Fill(1, evtwt);
                    h1_["st_cntD1ZB1Hb1b2"] -> Fill(ST, evtwt) ;
                    //h1_["ht_cntD1ZB1Hb1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                           
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow11"] -> Fill(2, evtwt);
                    h1_["st_cntD1ZB1H1b2"] -> Fill(ST, evtwt) ;
                    // h1_["ht_cntD1ZB1H1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                           
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow11"] -> Fill(3, evtwt);
                    h1_["st_cntD1Z1Hb1b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntD1Z1Hb1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                         
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow11"] -> Fill(4, evtwt);
                    h1_["st_cntD1Z1H1b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntD1Z1H1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                          
                  }
                }
              }
	      
              else if (BC.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow11"] -> Fill(5, evtwt);
                    h1_["st_cntBC1ZB1Hb1b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntBC1ZB1Hb1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                       
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow11"] -> Fill(6, evtwt);
                    h1_["st_cntBC1ZB1H1b2"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntBC1ZB1H1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                         
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() >=1.0){
                    h1_["cutflow11"] -> Fill(7, evtwt);
                    h1_["st_cntBC1Z1Hb1b2"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntBC1Z1Hb1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                         
                  }
                  else if(H.size() >=1.0){
                    h1_["cutflow11"] -> Fill(8, evtwt);
                    h1_["st_cntBC1Z1H1b2"] -> Fill(ST, evtwt) ;
                    h1_["ht_cntBC1Z1H1b2"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                }
              }
		else if (tops.size() >=1.0){
		  if (ZB.size() >=1.0){
		    if (Hb.size() >=1.0){
		      h1_["cutflow11"] -> Fill(9, evtwt);
		      h1_["st_cntt1ZB1Hb1b2"] -> Fill(ST, evtwt) ;
		      //  h1_["ht_cntt1ZB1Hb1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                         
		    }
		    else if(H.size() >=1.0){
		      h1_["cutflow11"] -> Fill(10, evtwt);
		      h1_["st_cntt1ZB1H1b2"] -> Fill(ST, evtwt) ;
		      //   h1_["ht_cntt1ZB1H1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                         
		    }
		  }
		  else if (Z.size() >=1.0){
		    if (Hb.size() >=1.0){
		      h1_["cutflow11"] -> Fill(11, evtwt);
		      h1_["st_cntt1Z1Hb1b2"] -> Fill(ST, evtwt) ;
		      // h1_["ht_cntt1Z1Hb1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                           
		    }
		    else if(H.size() >=1.0){
		      h1_["cutflow11"] -> Fill(12, evtwt);
		      h1_["st_cntt1Z1H1b2"] -> Fill(ST, evtwt) ;
		      // h1_["ht_cntt1Z1H1b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                            
		    }
		  }
		}
	    
	      
	      */
	      


	      /*  
	      for (auto izll : zll) {
		h1_["pt_z"+lep+lep+"_cntT1Z1H1b2"] -> Fill(izll.getPt(), evtwt) ;
		h1_["eta_z"+lep+lep+"_cntT1Z1H1b2"] -> Fill(izll.getEta(), evtwt) ;
		h1_["phi_z"+lep+lep+"_cntT1Z1H1b2"] -> Fill(izll.getPhi(), evtwt) ;

	      }
	      if ( zdecayMode_ == "zmumu" ){
		h1_["dr_mumu_cntT1Z1H1b2"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
		for(int l=0; l<2; ++l){
		  h1_["pt_"+lep+Form("%d_cntT1Z1H1b2", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
		  h1_["eta_"+lep+Form("%d_cntT1Z1H1b2", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
		  h1_["phi_"+lep+Form("%d_cntT1Z1H1b2", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
		}
	      }
	      else if (zdecayMode_ == "zelel" ) {
		for(int l=0; l<2; ++l){
		  h1_["dr_elel_cntT1Z1H1b2"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
		  h1_["pt_"+lep+Form("%d_cntT1Z1H1b2", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
		  h1_["eta_"+lep+Form("%d_cntT1Z1H1b2", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
		  h1_["phi_"+lep+Form("%d_cntT1Z1H1b2", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
		}
	      }
	      h1_["met_cntT1Z1H1b2"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	      h1_["metPhi_cntT1Z1H1b2"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
	      for(int j=0; j<3; ++j){
		h1_[Form("ptak4jet%d_cntT1Z1H1b2", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
		h1_[Form("etaak4jet%d_cntT1Z1H1b2", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
                h1_[Form("phiak4jet%d_cntT1Z1H1b2", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;		

	      }
		h1_["ptbjetleading_cntT1Z1H1b2"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
		h1_["etabjetleading_cntT1Z1H1b2"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;
		h1_["phibjetleading_cntT1Z1H1b2"] -> Fill(goodBTaggedAK4Jets.at(0).getPhi(), evtwt) ;
		h1_["ptbjetsubleading_cntT1Z1H1b2"] -> Fill(goodBTaggedAK4Jets.at(1).getPt(), evtwt) ;
		h1_["etabjetsubleading_cntT1Z1H1b2"] -> Fill(goodBTaggedAK4Jets.at(1).getEta(), evtwt) ;
		h1_["phibjetsubleading_cntT1Z1H1b2"] -> Fill(goodBTaggedAK4Jets.at(1).getPhi(), evtwt) ;	



	//	if(goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){
		// h1_["st_cntT1Z1H1b2_A"] -> Fill(ST, evtwt) ;
		// h1_["ht_cntT1Z1H1b2_A"] -> Fill(htak4.getHT(), evtwt) ;
		//	}

		*/
	      
	    }
	  }

	  else if (nHcandidates == 0.0){
	    //h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	    h1_["cutflow1"] -> Fill(9, evtwt) ;

	    h1_["st_cntT1Z1H0"] -> Fill(ST, evtwt) ;
            h1_["ht_cntT1Z1H0"] -> Fill(htak4.getHT(), evtwt) ;


	    /*	    
            for (auto izll : zll) {
              h1_["pt_z"+lep+lep+"_cntT1Z1H0"] -> Fill(izll.getPt(), evtwt) ;
              h1_["eta_z"+lep+lep+"_cntT1Z1H0"] -> Fill(izll.getEta(), evtwt) ;
              h1_["phi_z"+lep+lep+"_cntT1Z1H0"] -> Fill(izll.getPhi(), evtwt) ;

            }
            if ( zdecayMode_ == "zmumu" ){
	      h1_["dr_mumu_cntT1Z1H0"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT1Z1H0", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT1Z1H0", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT1Z1H0", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
              }
            }
            else if (zdecayMode_ == "zelel" ) {
	      h1_["dr_elel_cntT1Z1H0"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT1Z1H0", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT1Z1H0", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT1Z1H0", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
              }
            }
            h1_["met_cntT1Z1H0"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
            h1_["metPhi_cntT1Z1H0"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
            //  h1_["metEta_cntT1Z1H0"] -> Fill(goodMet.at(0).getFullEta(), evtwt);                                                                                                      

            for(int j=0; j<3; ++j){
              h1_[Form("ptak4jet%d_cntT1Z1H0", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
              h1_[Form("etaak4jet%d_cntT1Z1H0", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
              h1_[Form("phiak4jet%d_cntT1Z1H0", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;

            }



	    */



	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow2"] -> Fill(3, evtwt) ;
	      h1_["st_cntT1Z1H0b1"] -> Fill(ST, evtwt) ;
	      h1_["ht_cntT1Z1H0b1"] -> Fill(htak4.getHT(), evtwt) ;
	      h1_["nbjets_cntT1Z1H0b1"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
	    
	      /* if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT1Z1H0b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT1Z1H0b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/




	      /*
	      if (D.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow12"] -> Fill(1, evtwt);
                    h1_["st_cntD1ZB1Hb0b1"] -> Fill(ST, evtwt) ;
                    //h1_["ht_cntD1ZB1Hb0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                  
                  }
                  else if(H.size() ==0.0){
                    h1_["cutflow12"] -> Fill(2, evtwt);
                    h1_["st_cntD1ZB1H0b1"] -> Fill(ST, evtwt) ;
                    // h1_["ht_cntD1ZB1H0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                  
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow12"] -> Fill(3, evtwt);
                    h1_["st_cntD1Z1Hb0b1"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntD1Z1Hb0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow12"] -> Fill(4, evtwt);
                    h1_["st_cntD1Z1H0b1"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntD1Z1H0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                 
                  }
                }
              }

              else if (BC.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow12"] -> Fill(5, evtwt);
                    h1_["st_cntBC1ZB1Hb0b1"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntBC1ZB1Hb0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                              
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow12"] -> Fill(6, evtwt);
                    h1_["st_cntBC1ZB1H0b1"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntBC1ZB1H0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow12"] -> Fill(7, evtwt);
                    h1_["st_cntBC1Z1Hb0b1"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntBC1Z1Hb0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow12"] -> Fill(8, evtwt);
                    h1_["st_cntBC1Z1H0b1"] -> Fill(ST, evtwt) ;
                    h1_["ht_cntBC1Z1H0b1"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                }
              }

	      else if (tops.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow12"] -> Fill(9, evtwt);
                    h1_["st_cntt1ZB1Hb0b1"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntt1ZB1Hb0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow12"] -> Fill(10, evtwt);
                    h1_["st_cntt1ZB1H0b1"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntt1ZB1H0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow12"] -> Fill(11, evtwt);
                    h1_["st_cntt1Z1Hb0b1"] -> Fill(ST, evtwt) ;
                    // h1_["ht_cntt1Z1Hb0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                  
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow12"] -> Fill(12, evtwt);
                    h1_["st_cntt1Z1H0b1"] -> Fill(ST, evtwt) ;
                    // h1_["ht_cntt1Z1H0b1"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                   
                  }
                }
              }



	      */

	      
	      for (auto izll : zll) {
		h1_["pt_z"+lep+lep+"_cntT1Z1H0b1"] -> Fill(izll.getPt(), evtwt) ;
		h1_["eta_z"+lep+lep+"_cntT1Z1H0b1"] -> Fill(izll.getEta(), evtwt) ;
                h1_["phi_z"+lep+lep+"_cntT1Z1H0b1"] -> Fill(izll.getPhi(), evtwt) ;
	      }
	      
	      if ( zdecayMode_ == "zmumu" ){
		h1_["dr_mumu_cntT1Z1H0b1"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
		for(int l=0; l<2; ++l){
		  h1_["pt_"+lep+Form("%d_cntT1Z1H0b1", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
		  h1_["eta_"+lep+Form("%d_cntT1Z1H0b1", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
		  h1_["phi_"+lep+Form("%d_cntT1Z1H0b1", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
		}
	      }
	      else if (zdecayMode_ == "zelel" ) {
		h1_["dr_elel_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );

		h1_["dr_el1jet1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(0).getP4()), evtwt );
                h1_["dr_el1jet2_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(1).getP4()), evtwt );
                h1_["dr_el1jet3_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodAK4Jets.at(2).getP4()), evtwt );

                h1_["dr_el2jet1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(0).getP4()), evtwt );
		h1_["dr_el2jet2_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(1).getP4()), evtwt );
		h1_["dr_el2jet3_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodAK4Jets.at(2).getP4()), evtwt );

                h1_["dr_el1bjet1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodBTaggedAK4Jets.at(0).getP4()), evtwt );
                h1_["dr_el2bjet1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(goodBTaggedAK4Jets.at(0).getP4()), evtwt );


                if (Hb.size() >0){
                  h1_["dr_el1Hb1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(Hb.at(0).getP4()), evtwt );
                  h1_["dr_el2Hb1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(Hb.at(0).getP4()), evtwt );
                }

                if (ZB.size() >0){
                  h1_["dr_el1Zb1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(ZB.at(0).getP4()), evtwt );
                  h1_["dr_el2Zb1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(ZB.at(0).getP4()), evtwt );

                }

                if (D.size() >0){
                  h1_["dr_el1tb1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(D.at(0).getP4()), evtwt );
                  h1_["dr_el2tb1_cntT1Z1H0b1"]-> Fill( (goodElectrons.at(1).getP4()).DeltaR(D.at(0).getP4()), evtwt );
                }


	      
		
		for(int l=0; l<2; ++l){
		  h1_["pt_"+lep+Form("%d_cntT1Z1H0b1", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
		  h1_["eta_"+lep+Form("%d_cntT1Z1H0b1", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
		  h1_["phi_"+lep+Form("%d_cntT1Z1H0b1", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
		}
	      }
	      /*
	      h1_["met_cntT1Z1H0b1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	      h1_["metPhi_cntT1Z1H0b1"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
	      for(int j=0; j<3; ++j){
		h1_[Form("ptak4jet%d_cntT1Z1H0b1", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
		h1_[Form("etaak4jet%d_cntT1Z1H0b1", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
                h1_[Form("phiak4jet%d_cntT1Z1H0b1", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;
	      }
	      h1_["ptbjetleading_cntT1Z1H0b1"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
	      h1_["etabjetleading_cntT1Z1H0b1"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;
	      h1_["phibjetleading_cntT1Z1H0b1"] -> Fill(goodBTaggedAK4Jets.at(0).getPhi(), evtwt) ;
	      if(goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){
		h1_["st_cntT1Z1H0b1_A"] -> Fill(ST, evtwt) ;
		h1_["ht_cntT1Z1H0b1_A"] -> Fill(htak4.getHT(), evtwt) ;
	      }
	      */
	     
		}
	    
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow2"] -> Fill(4, evtwt) ;
	      h1_["st_cntT1Z1H0b2"] -> Fill(ST, evtwt) ;
	      h1_["ht_cntT1Z1H0b2"] -> Fill(htak4.getHT(), evtwt) ;
	      h1_["nbjets_cntT1Z1H0b2"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
	     
	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT1Z1H0b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT1Z1H0b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/


	      /*
	      if (D.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow13"] -> Fill(1, evtwt);
                    h1_["st_cntD1ZB1Hb0b2"] -> Fill(ST, evtwt) ;
                    //h1_["ht_cntD1ZB1Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                  
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow13"] -> Fill(2, evtwt);
                    h1_["st_cntD1ZB1H0b2"] -> Fill(ST, evtwt) ;
                    // h1_["ht_cntD1ZB1H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                  
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow13"] -> Fill(3, evtwt);
                    h1_["st_cntD1Z1Hb0b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntD1Z1Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow13"] -> Fill(4, evtwt);
                    h1_["st_cntD1Z1H0b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntD1Z1H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                 
                  }
                }
              }

              else if (BC.size() >=1.0){
                if (ZB.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow13"] -> Fill(5, evtwt);
                    h1_["st_cntBC1ZB1Hb0b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntBC1ZB1Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                              
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow13"] -> Fill(6, evtwt);
                    h1_["st_cntBC1ZB1H0b2"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntBC1ZB1H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
                  }
                }
                else if (Z.size() >=1.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow13"] -> Fill(7, evtwt);
                    h1_["st_cntBC1Z1Hb0b2"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntBC1Z1Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow13"] -> Fill(8, evtwt);
                    h1_["st_cntBC1Z1H0b2"] -> Fill(ST, evtwt) ;
                    h1_["ht_cntBC1Z1H0b2"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                }
              }

	      else if (tops.size() >=1.0){
		if (ZB.size() >=1.0){
		  if (Hb.size() == 0.0){
		    h1_["cutflow13"] -> Fill(9, evtwt);
		    h1_["st_cntt1ZB1Hb0b2"] -> Fill(ST, evtwt) ;
		    //  h1_["ht_cntt1ZB1Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                              
		  }
		  else if(H.size() == 0.0){
		    h1_["cutflow13"] -> Fill(10, evtwt);
		    h1_["st_cntt1ZB1H0b2"] -> Fill(ST, evtwt) ;
		    //   h1_["ht_cntt1ZB1H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                              
		  }
		}
		else if (Z.size() >=1.0){
		  if (Hb.size() == 0.0){
		    h1_["cutflow13"] -> Fill(11, evtwt);
		    h1_["st_cntt1Z1Hb0b2"] -> Fill(ST, evtwt) ;
		    // h1_["ht_cntt1Z1Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                
		  }
		  else if(H.size() == 0.0){
		    h1_["cutflow13"] -> Fill(12, evtwt);
		    h1_["st_cntt1Z1H0b2"] -> Fill(ST, evtwt) ;
		    // h1_["ht_cntt1Z1H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                 
		  }
		}
	      }
	      */




	      /*
	      
	      for (auto izll : zll) {
		h1_["pt_z"+lep+lep+"_cntT1Z1H0b2"] -> Fill(izll.getPt(), evtwt) ;
		h1_["eta_z"+lep+lep+"_cntT1Z1H0b2"] -> Fill(izll.getEta(), evtwt) ;
                h1_["phi_z"+lep+lep+"_cntT1Z1H0b2"] -> Fill(izll.getPhi(), evtwt) ;

	      }
	      if ( zdecayMode_ == "zmumu" ){
		h1_["dr_mumu_cntT1Z1H0b2"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
		for(int l=0; l<2; ++l){
		  h1_["pt_"+lep+Form("%d_cntT1Z1H0b2", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
		  h1_["eta_"+lep+Form("%d_cntT1Z1H0b2", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
		  h1_["phi_"+lep+Form("%d_cntT1Z1H0b2", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
		}
	      }
	      else if (zdecayMode_ == "zelel" ) {
		h1_["dr_elel_cntT1Z1H0b2"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
		for(int l=0; l<2; ++l){
		  h1_["pt_"+lep+Form("%d_cntT1Z1H0b2", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
		  h1_["eta_"+lep+Form("%d_cntT1Z1H0b2", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
		  h1_["phi_"+lep+Form("%d_cntT1Z1H0b2", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
		}
	      }
	      h1_["met_cntT1Z1H0b2"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
	      h1_["metPhi_cntT1Z1H0b2"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
	      for(int j=0; j<3; ++j){
		h1_[Form("ptak4jet%d_cntT1Z1H0b2", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
		h1_[Form("etaak4jet%d_cntT1Z1H0b2", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
                h1_[Form("phiak4jet%d_cntT1Z1H0b2", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;
	      }
	      h1_["ptbjetleading_cntT1Z1H0b2"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
	      h1_["etabjetleading_cntT1Z1H0b2"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;
	      h1_["phibjetleading_cntT1Z1H0b2"] -> Fill(goodBTaggedAK4Jets.at(0).getPhi(), evtwt) ;
	      h1_["ptbjetsubleading_cntT1Z1H0b2"] -> Fill(goodBTaggedAK4Jets.at(1).getPt(), evtwt) ;
	      h1_["etabjetsubleading_cntT1Z1H0b2"] -> Fill(goodBTaggedAK4Jets.at(1).getEta(), evtwt) ;
	      h1_["phibjetsubleading_cntT1Z1H0b2"] -> Fill(goodBTaggedAK4Jets.at(1).getPhi(), evtwt) ;
	      if(goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){
		h1_["st_cntT1Z1H0b2_A"] -> Fill(ST, evtwt) ;
		h1_["ht_cntT1Z1H0b2_A"] -> Fill(htak4.getHT(), evtwt) ;
	      }
	      */
	    }

	  }
	}

	/*
	//HPrime Candidates
	// cout << "ntopcandidates =" <<ntopcandidates <<endl;
	//  cout <<"nzcandidates = "<<nzcandidates <<endl;
	//  cout<<"nHPrimecandidates =" << nHPrimecandidates <<endl;
	//  cout<<"goodBTaggedAK4Jets.size() = " << goodBTaggedAK4Jets.size() <<endl;
	if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
	if (nHPrimecandidates >= 1.0){
	// h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);                                                                                                                                 
	// h1_["cutflow1"] -> Fill(8, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 0 ){
	//	    cout<<"st T1Z1H1b0 " << ST <<endl;
	//  h1_["cutflow2"] -> Fill(1, evtwt) ;
	h1_["st_cntT1Z1Hprime1b0"] -> Fill(ST, evtwt) ;

	for (auto izll : zll) {
	h1_["pt_z"+lep+lep+"_cntT1Z1Hprime1b0"] -> Fill(izll.getPt(), evtwt) ;
	}
	if ( zdecayMode_ == "zmumu" ){
	for(int l=0; l<2; ++l){
	h1_["pt_"+lep+Form("%d_cntT1Z1Hprime1b0", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
	h1_["eta_"+lep+Form("%d_cntT1Z1Hprime1b0", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
	}
	}
	else if (zdecayMode_ == "zelel" ) {
	for(int l=0; l<2; ++l){
	h1_["pt_"+lep+Form("%d_cntT1Z1Hprime1b0", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
	h1_["eta_"+lep+Form("%d_cntT1Z1Hprime1b0", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
	}
	}
	h1_["met_cntT1Z1Hprime1b0"] -> Fill(goodMet.at(0).getFullPt(), evtwt);

	for(int j=0; j<3; ++j){
	h1_[Form("ptak4jet%d_cntT1Z1Hprime1b0", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;

	}
     
	}

	}
  
	}
	*/
	//(2)                                                                                                                                                                                        
	if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
	  //h1_["st_cntT0Z1"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(5, evtwt) ;
	  if (nHcandidates >= 1.0){
	    // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	    h1_["cutflow1"] -> Fill(10, evtwt) ;
	   

	    h1_["st_cntT0Z1H1"] -> Fill(ST, evtwt) ;
            h1_["ht_cntT0Z1H1"] -> Fill(htak4.getHT(), evtwt) ;
	    /*
            for (auto izll : zll) {
              h1_["pt_z"+lep+lep+"_cntT0Z1H1"] -> Fill(izll.getPt(), evtwt) ;
              h1_["eta_z"+lep+lep+"_cntT0Z1H1"] -> Fill(izll.getEta(), evtwt) ;
              h1_["phi_z"+lep+lep+"_cntT0Z1H1"] -> Fill(izll.getPhi(), evtwt) ;

            }
            if ( zdecayMode_ == "zmumu" ){
              h1_["dr_mumu_cntT0Z1H1"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z1H1", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z1H1", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT0Z1H1", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
              }
            }
            else if (zdecayMode_ == "zelel" ) {
              h1_["dr_elel_cntT0Z1H1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z1H1", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z1H1", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT0Z1H1", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
              }
            }
            h1_["met_cntT0Z1H1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
            h1_["metPhi_cntT0Z1H1"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
            //  h1_["metEta_cntT0Z1H1"] -> Fill(goodMet.at(0).getFullEta(), evtwt);                                                                                                                                           

            for(int j=0; j<3; ++j){
              h1_[Form("ptak4jet%d_cntT0Z1H1", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
              h1_[Form("etaak4jet%d_cntT0Z1H1", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
              h1_[Form("phiak4jet%d_cntT0Z1H1", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;

            }
	    */




	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow2"] -> Fill(5, evtwt) ;
	      h1_["st_cntT0Z1H1b1"] -> Fill(ST, evtwt) ;
	    
	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT0Z1H1b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT0Z1H1b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow2"] -> Fill(6, evtwt) ;
	      h1_["st_cntT0Z1H1b2"] -> Fill(ST, evtwt) ;
	      h1_["ht_cntT0Z1H1b2"] -> Fill(htak4.getHT(), evtwt) ;
	
	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT0Z1H1b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT0Z1H1b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



      /*
		for (auto izll : zll) {
		h1_["pt_z"+lep+lep+"_cntT0Z1H1b2"] -> Fill(izll.getPt(), evtwt) ;
		}
		if ( zdecayMode_ == "zmumu" ){
		for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z1H1b2", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z1H1b2", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
		}
		}
		else if (zdecayMode_ == "zelel" ) {
		for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z1H1b2", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z1H1b2", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
		}
		}
		h1_["met_cntT0Z1H1b2"] -> Fill(goodMet.at(0).getFullPt(), evtwt);

		for(int j=0; j<3; ++j){
		h1_[Form("ptak4jet%d_cntT0Z1H1b2", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;

		}
		h1_["ptbjetleading_cntT0Z1H1b2"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
		h1_["ptbjetsubleading_cntT0Z1H1b2"] -> Fill(goodBTaggedAK4Jets.at(1).getPt(), evtwt) ;
		if(goodAK4Jets.at(0).getPt() > 100 && goodAK4Jets.at(1).getPt() > 50){
		h1_["st_cntT0Z1H1b2_A"] -> Fill(ST, evtwt) ;
		h1_["ht_cntT0Z1H1b2_A"] -> Fill(htak4.getHT(), evtwt) ;
		}

	      */
	    }
	  }
	  else if (nHcandidates == 0.0){


	    // h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	    h1_["cutflow1"] -> Fill(11, evtwt) ;

	    h1_["st_cntT0Z1H0"] -> Fill(ST, evtwt) ;
            h1_["ht_cntT0Z1H0"] -> Fill(htak4.getHT(), evtwt) ;
	    /*
            for (auto izll : zll) {
              h1_["pt_z"+lep+lep+"_cntT0Z1H0"] -> Fill(izll.getPt(), evtwt) ;
              h1_["eta_z"+lep+lep+"_cntT0Z1H0"] -> Fill(izll.getEta(), evtwt) ;
              h1_["phi_z"+lep+lep+"_cntT0Z1H0"] -> Fill(izll.getPhi(), evtwt) ;

            }
            if ( zdecayMode_ == "zmumu" ){
              h1_["dr_mumu_cntT0Z1H0"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z1H0", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z1H0", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT0Z1H0", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
              }
            }
            else if (zdecayMode_ == "zelel" ) {
              h1_["dr_elel_cntT0Z1H0"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z1H0", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z1H0", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT0Z1H0", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
              }
            }
            h1_["met_cntT0Z1H0"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
            h1_["metPhi_cntT0Z1H0"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
            //  h1_["metEta_cntT0Z1H0"] -> Fill(goodMet.at(0).getFullEta(), evtwt);                                                                                                                                          \
                                                                                                                                                                                                                              

            for(int j=0; j<3; ++j){
              h1_[Form("ptak4jet%d_cntT0Z1H0", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
              h1_[Form("etaak4jet%d_cntT0Z1H0", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
              h1_[Form("phiak4jet%d_cntT0Z1H0", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;

            }
	    */

	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow2"] -> Fill(7, evtwt) ;
	      h1_["st_cntT0Z1H0b1"] -> Fill(ST, evtwt) ;

	      /*   if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT0Z1H0b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT0Z1H0b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/




	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow2"] -> Fill(8, evtwt) ;
	      h1_["st_cntT0Z1H0b2"] -> Fill(ST, evtwt) ;

	      /* if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT0Z1H0b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT0Z1H0b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



	    }
	  }
	}

	//(3)                                                                                                                                                                                      
	if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
	  //h1_["st_cntT1Z0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(6, evtwt) ;
	  if (nHcandidates >= 1.0){
	    // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	    h1_["cutflow1"] -> Fill(12, evtwt) ;

	    h1_["st_cntT1Z0H1"] -> Fill(ST, evtwt) ;
            h1_["ht_cntT1Z0H1"] -> Fill(htak4.getHT(), evtwt) ;
	    /*
            for (auto izll : zll) {
              h1_["pt_z"+lep+lep+"_cntT1Z0H1"] -> Fill(izll.getPt(), evtwt) ;
              h1_["eta_z"+lep+lep+"_cntT1Z0H1"] -> Fill(izll.getEta(), evtwt) ;
              h1_["phi_z"+lep+lep+"_cntT1Z0H1"] -> Fill(izll.getPhi(), evtwt) ;

            }
            if ( zdecayMode_ == "zmumu" ){
              h1_["dr_mumu_cntT1Z0H1"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT1Z0H1", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT1Z0H1", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT1Z0H1", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
              }
            }
            else if (zdecayMode_ == "zelel" ) {
              h1_["dr_elel_cntT1Z0H1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT1Z0H1", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT1Z0H1", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT1Z0H1", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
              }
            }
            h1_["met_cntT1Z0H1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
            h1_["metPhi_cntT1Z0H1"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
            //  h1_["metEta_cntT1Z0H1"] -> Fill(goodMet.at(0).getFullEta(), evtwt);                                                                                                                                          \
                                                                                                                                                                                                                              

            for(int j=0; j<3; ++j){
              h1_[Form("ptak4jet%d_cntT1Z0H1", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
              h1_[Form("etaak4jet%d_cntT1Z0H1", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
              h1_[Form("phiak4jet%d_cntT1Z0H1", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;

            }
	    */
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow2"] -> Fill(9, evtwt) ;
	      h1_["st_cntT1Z0H1b1"] -> Fill(ST, evtwt) ;

	      /* if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT1Z0H1b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT1Z0H1b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow2"] -> Fill(10, evtwt) ;
	      h1_["st_cntT1Z0H1b2"] -> Fill(ST, evtwt) ;
	 
	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT1Z0H1b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT1Z0H1b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/




	    }
	  }
	  else if (nHcandidates == 0.0){
	    // h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	    h1_["cutflow1"] -> Fill(13, evtwt) ;

	    h1_["st_cntT1Z0H0"] -> Fill(ST, evtwt) ;
            h1_["ht_cntT1Z0H0"] -> Fill(htak4.getHT(), evtwt) ;
	    /*
            for (auto izll : zll) {
              h1_["pt_z"+lep+lep+"_cntT1Z0H0"] -> Fill(izll.getPt(), evtwt) ;
              h1_["eta_z"+lep+lep+"_cntT1Z0H0"] -> Fill(izll.getEta(), evtwt) ;
              h1_["phi_z"+lep+lep+"_cntT1Z0H0"] -> Fill(izll.getPhi(), evtwt) ;

            }
            if ( zdecayMode_ == "zmumu" ){
              h1_["dr_mumu_cntT1Z0H0"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT1Z0H0", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT1Z0H0", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT1Z0H0", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
              }
            }
            else if (zdecayMode_ == "zelel" ) {
              h1_["dr_elel_cntT1Z0H0"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT1Z0H0", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT1Z0H0", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT1Z0H0", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
              }
            }
            h1_["met_cntT1Z0H0"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
            h1_["metPhi_cntT1Z0H0"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
            //  h1_["metEta_cntT1Z0H0"] -> Fill(goodMet.at(0).getFullEta(), evtwt);                                                                                                                                        
                                                                                                                                                                                                                           
            for(int j=0; j<3; ++j){
              h1_[Form("ptak4jet%d_cntT1Z0H0", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
              h1_[Form("etaak4jet%d_cntT1Z0H0", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
              h1_[Form("phiak4jet%d_cntT1Z0H0", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;

            }
	    */
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow2"] -> Fill(11, evtwt) ;
	      h1_["st_cntT1Z0H0b1"] -> Fill(ST, evtwt) ;

	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT1Z0H0b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT1Z0H0b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow2"] -> Fill(12, evtwt) ;
	      h1_["st_cntT1Z0H0b2"] -> Fill(ST, evtwt) ;
	  
	      /* if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT1Z0H0b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT1Z0H0b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/




	      /*

	      if (D.size() >=1.0){
                if (ZB.size() == 0.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow14"] -> Fill(1, evtwt);
                    h1_["st_cntD1ZB0Hb0b2"] -> Fill(ST, evtwt) ;
                    //h1_["ht_cntD1ZB0Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                       
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow14"] -> Fill(2, evtwt);
                    h1_["st_cntD1ZB0H0b2"] -> Fill(ST, evtwt) ;
                    // h1_["ht_cntD1ZB0H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                       
                  }
                }
                else if (Z.size() == 0.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow14"] -> Fill(3, evtwt);
                    h1_["st_cntD1Z0Hb0b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntD1Z0Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                     
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow14"] -> Fill(4, evtwt);
                    h1_["st_cntD1Z0H0b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntD1Z0H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                      
                  }
                }
              }

              else if (BC.size() >=1.0){
                if (ZB.size() == 0.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow14"] -> Fill(5, evtwt);
                    h1_["st_cntBC1ZB0Hb0b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntBC1ZB0Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                   
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow14"] -> Fill(6, evtwt);
                    h1_["st_cntBC1ZB0H0b2"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntBC1ZB0H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                     
                  }
                }
                else if (Z.size() == 0.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow14"] -> Fill(7, evtwt);
                    h1_["st_cntBC1Z0Hb0b2"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntBC1Z0Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                     
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow14"] -> Fill(8, evtwt);
                    h1_["st_cntBC1Z0H0b2"] -> Fill(ST, evtwt) ;
                    h1_["ht_cntBC1Z0H0b2"] -> Fill(htak4.getHT(), evtwt) ;
                  }
                }
              }
	      else if (tops.size() >=1.0){
                if (ZB.size() == 0.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow14"] -> Fill(9, evtwt);
                    h1_["st_cntt1ZB0Hb0b2"] -> Fill(ST, evtwt) ;
                    //  h1_["ht_cntt1ZB0Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                     
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow14"] -> Fill(10, evtwt);
                    h1_["st_cntt1ZB0H0b2"] -> Fill(ST, evtwt) ;
                    //   h1_["ht_cntt1ZB0H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                     
                  }
                }
                else if (Z.size() == 0.0){
                  if (Hb.size() == 0.0){
                    h1_["cutflow14"] -> Fill(11, evtwt);
                    h1_["st_cntt1Z0Hb0b2"] -> Fill(ST, evtwt) ;
                    // h1_["ht_cntt1Z0Hb0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                       
                  }
                  else if(H.size() == 0.0){
                    h1_["cutflow14"] -> Fill(12, evtwt);
                    h1_["st_cntt1Z0H0b2"] -> Fill(ST, evtwt) ;
                    // h1_["ht_cntt1Z0H0b2"] -> Fill(htak4.getHT(), evtwt) ;                                                                                                                                                                                                                        
                  }
                }
              }

	      */




	    }
	  }
	}

	//(4)                                                                                                                                                                                        
	if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
	  //h1_`["st_cntT0Z0"] -> Fill(ST,evtwt);
	  h1_["cutflow1"] -> Fill(7, evtwt) ;
	  if (nHcandidates >= 1.0){
	    // h1_["st_cntT1Z1H1"] -> Fill(ST,evtwt);
	    h1_["cutflow1"] -> Fill(14, evtwt) ;

	    h1_["st_cntT0Z0H1"] -> Fill(ST, evtwt) ;
            h1_["ht_cntT0Z0H1"] -> Fill(htak4.getHT(), evtwt) ;

            for (auto izll : zll) {
              h1_["pt_z"+lep+lep+"_cntT0Z0H1"] -> Fill(izll.getPt(), evtwt) ;
              h1_["eta_z"+lep+lep+"_cntT0Z0H1"] -> Fill(izll.getEta(), evtwt) ;
              h1_["phi_z"+lep+lep+"_cntT0Z0H1"] -> Fill(izll.getPhi(), evtwt) ;

            }
            if ( zdecayMode_ == "zmumu" ){
              h1_["dr_mumu_cntT0Z0H1"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z0H1", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z0H1", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT0Z0H1", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
              }
            }
            else if (zdecayMode_ == "zelel" ) {
              h1_["dr_elel_cntT0Z0H1"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z0H1", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z0H1", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT0Z0H1", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
              }
            }
            h1_["met_cntT0Z0H1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
            h1_["metPhi_cntT0Z0H1"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
            //  h1_["metEta_cntT0Z0H1"] -> Fill(goodMet.at(0).getFullEta(), evtwt);                                                                                                                                          \
                                                                                                                                                                                                                              

            for(int j=0; j<3; ++j){
              h1_[Form("ptak4jet%d_cntT0Z0H1", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
              h1_[Form("etaak4jet%d_cntT0Z0H1", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
              h1_[Form("phiak4jet%d_cntT0Z0H1", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;

            }

	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow2"] -> Fill(13, evtwt) ;
	      h1_["st_cntT0Z0H1b1"] -> Fill(ST, evtwt) ;
	    
	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT0Z0H1b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT0Z0H1b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow2"] -> Fill(14, evtwt) ;
	      h1_["st_cntT0Z0H1b2"] -> Fill(ST, evtwt) ;

	      /*   if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT0Z0H1b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT0Z0H1b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/

	    }
	  }
	  else if (nHcandidates == 0.0){
	    //h1_["st_cntT1Z1H0"] -> Fill(ST,evtwt);
	    h1_["cutflow1"] -> Fill(15, evtwt) ;

	    h1_["st_cntT0Z0H0"] -> Fill(ST, evtwt) ;
            h1_["ht_cntT0Z0H0"] -> Fill(htak4.getHT(), evtwt) ;

            for (auto izll : zll) {
              h1_["pt_z"+lep+lep+"_cntT0Z0H0"] -> Fill(izll.getPt(), evtwt) ;
              h1_["eta_z"+lep+lep+"_cntT0Z0H0"] -> Fill(izll.getEta(), evtwt) ;
              h1_["phi_z"+lep+lep+"_cntT0Z0H0"] -> Fill(izll.getPhi(), evtwt) ;

            }
            if ( zdecayMode_ == "zmumu" ){
              h1_["dr_mumu_cntT0Z0H0"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z0H0", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z0H0", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT0Z0H0", l+1)]  -> Fill(goodMuons.at(l).getPhi(), evtwt) ;
              }
            }
            else if (zdecayMode_ == "zelel" ) {
              h1_["dr_elel_cntT0Z0H0"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
              for(int l=0; l<2; ++l){
                h1_["pt_"+lep+Form("%d_cntT0Z0H0", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ;
                h1_["eta_"+lep+Form("%d_cntT0Z0H0", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ;
                h1_["phi_"+lep+Form("%d_cntT0Z0H0", l+1)]  -> Fill(goodElectrons.at(l).getPhi(), evtwt) ;
              }
            }
            h1_["met_cntT0Z0H0"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
            h1_["metPhi_cntT0Z0H0"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);
            //  h1_["metEta_cntT0Z0H0"] -> Fill(goodMet.at(0).getFullEta(), evtwt);                                                                                                                                                                                                                                                                                                                                                                        

            for(int j=0; j<3; ++j){
              h1_[Form("ptak4jet%d_cntT0Z0H0", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ;
              h1_[Form("etaak4jet%d_cntT0Z0H0", j+1)]  -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
              h1_[Form("phiak4jet%d_cntT0Z0H0", j+1)]  -> Fill(goodAK4Jets.at(j).getPhi(), evtwt) ;

            }

	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow2"] -> Fill(15, evtwt) ;
	      h1_["st_cntT0Z0H0b1"] -> Fill(ST, evtwt) ;
	   
	      /* if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT0Z0H0b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT0Z0H0b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/



	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow2"] -> Fill(16, evtwt) ;
	      h1_["st_cntT0Z0H0b2"] -> Fill(ST, evtwt) ;
	   
	      /*  if (!isData) {
                for (unsigned i = 0; i < 9; i++) {
                  h1_[Form("st_cntT0Z0H0b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
                for (unsigned i = 0; i < 100; i++) {
                  h1_[Form("st_cntT0Z0H0b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
                }
		}*/




	    }
	  }
	}
      }
    }
    
    Hb.clear();
    H.clear();
    ZB.clear();
    Z.clear();
    D.clear();
    BC.clear();
    tops.clear();
    HbPrime.clear();
    HPrime.clear();
    W.clear();
    B.clear();
 
    goodAK4Jetscleaned.clear();
    ak4matchedak8.clear();
    ak4nonmatched1.clear(); 
    ak4nonmatched2.clear();
    ak4nonmatched3.clear();





    //if ( goodAK4Jets.at(0).getPt() > 100
    //   && goodAK4Jets.at(1).getPt() > 50
    //   && goodBTaggedAK4Jets.size() > 0){

    //  h1_["st_noprime1"] -> Fill(ST, evtwt) ;
    //   h1_["st_prime1"] -> Fill(ST, evtwt) ;

    //  if (ST >2500.00) ST = 2498.00;
    //  if (ST >2500.00) ST = 2498.00;
    //  h1_["st_noprime"] -> Fill(ST, evtwt) ;
      //  h1_["st_noprime1"] -> Fill(ST, evtwt) ;
    //  h1_["st_prime"] -> Fill(ST, evtwt) ;
      // h1_["st_prime1"] -> Fill(ST, evtwt) ;
    // }



    if ( goodAK4Jets.at(0).getPt() > 100 
	 && goodAK4Jets.at(1).getPt() > 50
	 && goodBTaggedAK4Jets.size() > 0
	 && ST > STMin_   
	 ) { //// Signal region 

      // h1_["post_pdf"] -> Fill(1, pdfShift);
      //// fill all the plots in signal region
          
      /* if (!isData) {
        for (unsigned i = 0; i < 9; i++) {
          h1_[Form("st_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
        }
        for (unsigned i = 0; i < 100; i++) {
          h1_[Form("st_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
        }
	}*/
            
      //   goodTopTaggedJets.at(i).getIndex()<<endl;
      // goodWTaggedJets.size()
      // mistag rates
      for (unsigned int i=0; i< goodAK8Jets.size();i++){
	h1_["ptak8_st1000"] -> Fill((goodAK8Jets.at(i)).getPt(), evtwt) ;
      }
      
      for (unsigned int i=0; i< goodWTaggedJets.size();i++){
	h1_["ptW_st1000"] -> Fill((goodWTaggedJets.at(i)).getPt(), evtwt) ;
      }
      for (unsigned int i=0; i< goodHTaggedJets.size();i++){
	h1_["ptH_st1000"] -> Fill((goodHTaggedJets.at(i)).getPt(), evtwt) ;
      }
      for (unsigned int i=0; i< goodTopTaggedJets.size();i++){
	h1_["ptT_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
      }
      for (unsigned int i=0; i< goodTopTaggedJets.size();i++){
	for (unsigned int j=0; j< goodHTaggedJets.size();j++){
	  if (goodTopTaggedJets.at(i).getIndex()== goodHTaggedJets.at(j).getIndex()){
	    double evtwtTop = evtwt;
	    // double evtwt2 = evtwt;

	    h1_["ptTHnowwt_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;

	    evtwtTop *= 1.06;
	    h1_["ptTHwt_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwtTop) ;
	  }
	}
      }


      //signal mistags                                                                                                                                                     
      if (!isData  && filterSignal_ ){                                                                                                                          
  	GenParticleCollection genPartsInfo;                                                                                                                                 
	genPartsInfo = genpart(evt) ;                                                                                                                                       	//  if (gen1.getPdgID() == -6 && gen1.getMom0PdgID()== -8001){ Topgen1 = gen1.getP4();}
	//  else if (gen1.getPdgID() == 6 && gen1.getMom0PdgID()== 8001){ Topgen2 = gen1.getP4();}  
	  //  else if (gen1.getPdgID() == -24 && gen1.getMom0PdgID()== -8001 || gen1.getPdgID() == 24 && gen1.getMom0PdgID()== 8001 || gen1.getPdgID() == -23 && gen1.getMom0PdgID()== -8001 || gen1.getPdgID() == 23 && gen1.getMom0PdgID()== 8001){Zgen = gen1.getP4();}                                                
	  //  else if (gen1.getPdgID() == -25 && gen1.getMom0PdgID()== -8001 || gen1.getPdgID() == 25 && gen1.getMom0PdgID()== 8001){ Hgen = gen1.getP4();}
                                                                                                                                                                            
	  for (unsigned int i=0; i< goodTopTaggedJets.size();i++){
	    TLorentzVector Topgen1,Topgen2,Hgen1, Hgen2;;
	    //bool matchedtop = false;
	    // bool matchedH = false;
	    if (goodHTaggedJets.size()>0){
	      for (unsigned int j=0; j< goodHTaggedJets.size();j++){
		bool matchedtop = false;
		bool matchedH = false;
		if (goodTopTaggedJets.at(i).getIndex()== goodHTaggedJets.at(j).getIndex()){
		  for (auto& gen1 : genPartsInfo){
		    if ((gen1.getPdgID() == -6 && gen1.getMom0PdgID()== -8000001) || (gen1.getPdgID() == 6 && gen1.getMom0PdgID()== -8000001)){ Topgen1 = gen1.getP4();}
		    else if ((gen1.getPdgID() == 6 && gen1.getMom0PdgID()== 8000001) || (gen1.getPdgID() == -6 && gen1.getMom0PdgID()== 8000001) ){ Topgen2 = gen1.getP4();}
		    else if ((gen1.getPdgID() == -25 && gen1.getMom0PdgID()== -8000001)|| (gen1.getPdgID() == 25 && gen1.getMom0PdgID()== -8000001)){ Hgen1 = gen1.getP4();}
		    else if ((gen1.getPdgID() == 25 && gen1.getMom0PdgID()== 8000001)|| (gen1.getPdgID() == -25 && gen1.getMom0PdgID()== 8000001)){ Hgen2 = gen1.getP4();}
		    

		    if (goodTopTaggedJets.at(i).getP4().DeltaR(Topgen1) < 0.8 || goodTopTaggedJets.at(i).getP4().DeltaR(Topgen2)< 0.8){
		      matchedtop = true;
		      // h1_["ptTmatched_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
		      //break;
		    }
		    else if (goodHTaggedJets.at(j).getP4().DeltaR(Hgen1) < 0.8 || goodHTaggedJets.at(j).getP4().DeltaR(Hgen2)< 0.8){
		      matchedH = true;
		      //break;
		     
		    }
		   
		    cout<< " matchedtop ="<< matchedtop <<endl;
		    cout<< " matchedH ="<< matchedH <<endl;
		    cout <<" end ***** " <<endl;
		    if (matchedtop && !matchedH){h1_["ptTmatched_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ; break;}
		    else if (!matchedtop && matchedH){break;}
		  }
		}

		else{
		  for (auto& gen1 : genPartsInfo){
		    if ((gen1.getPdgID() == -6 && gen1.getMom0PdgID()== -8000001) || (gen1.getPdgID() == 6 && gen1.getMom0PdgID()== -8000001)){ Topgen1 = gen1.getP4();}
		    else if ((gen1.getPdgID() == 6 && gen1.getMom0PdgID()== 8000001) || (gen1.getPdgID() == -6 && gen1.getMom0PdgID()== 8000001) ){ Topgen2 = gen1.getP4();}
		    if (goodTopTaggedJets.at(i).getP4().DeltaR(Topgen1) < 0.8 || goodTopTaggedJets.at(i).getP4().DeltaR(Topgen2)< 0.8){
		      h1_["ptTmatched_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
		      matchedtop = true;
		      break;
		    }
		  }
		  if (matchedtop == true){h1_["ptTmatched_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;}
		  else { h1_["ptTnonmatched_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;}
		}
	      }
	    }
	    else{
	      bool matchedtop = false;
		for (auto& gen1 : genPartsInfo){
		  if ((gen1.getPdgID() == -6 && gen1.getMom0PdgID()== -8000001) || (gen1.getPdgID() == 6 && gen1.getMom0PdgID()== -8000001)){ Topgen1 = gen1.getP4();}
		  else if ((gen1.getPdgID() == 6 && gen1.getMom0PdgID()== 8000001) || (gen1.getPdgID() == -6 && gen1.getMom0PdgID()== 8000001) ){ Topgen2 = gen1.getP4();}
		  if (goodTopTaggedJets.at(i).getP4().DeltaR(Topgen1) < 0.8 || goodTopTaggedJets.at(i).getP4().DeltaR(Topgen2)< 0.8){
		    matchedtop = true;
		    break;
		  }
		}
		if (matchedtop == true){h1_["ptTmatched_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;}
		else { h1_["ptTnonmatched_st1000"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;}
	    }//end toptagged
	  }
      
	  for (unsigned int i=0; i< goodHTaggedJets.size();i++){
            TLorentzVector Hgen1, Hgen2;
            bool matchedH = false;
            for (auto& gen1 : genPartsInfo){
              if ((gen1.getPdgID() == -25 && gen1.getMom0PdgID()== -8000001)|| (gen1.getPdgID() == 25 && gen1.getMom0PdgID()== -8000001)){ Hgen1 = gen1.getP4();}
              else if ((gen1.getPdgID() == 25 && gen1.getMom0PdgID()== 8000001)|| (gen1.getPdgID() == -25 && gen1.getMom0PdgID()== 8000001)){ Hgen2 = gen1.getP4();}
              if (goodHTaggedJets.at(i).getP4().DeltaR(Hgen1) < 0.8 || goodHTaggedJets.at(i).getP4().DeltaR(Hgen2)< 0.8){
                matchedH = true;
                break;
	      } 
            }
            if (matchedH == true){h1_["ptHmatched_st1000"] -> Fill((goodHTaggedJets.at(i)).getPt(), evtwt) ;}
            else { h1_["ptHnonmatched_st1000"] -> Fill((goodHTaggedJets.at(i)).getPt(), evtwt) ;}
          }//end Htagged     

	  
	  for (unsigned int i=0; i< goodWTaggedJets.size();i++){
            TLorentzVector Wgen1, Wgen2, Wgen3,Wgen4, Wgen5,Wgen6;
            bool matchedW = false;
            for (auto& gen1 : genPartsInfo){
              if ((gen1.getPdgID() == -24 && gen1.getMom0PdgID()== -8000001) || (gen1.getPdgID() == 24 && gen1.getMom0PdgID()== -8000001)){ Wgen1 = gen1.getP4();}
              else if ((gen1.getPdgID() == 24 && gen1.getMom0PdgID()== 8000001) || (gen1.getPdgID() == -24 && gen1.getMom0PdgID()== 8000001)){ Wgen2 = gen1.getP4();}
	      else if ((gen1.getPdgID() == 23 && gen1.getMom0PdgID()== 8000001) || (gen1.getPdgID() == -23 && gen1.getMom0PdgID()== 8000001)){ Wgen3 = gen1.getP4();}
	      else if ((gen1.getPdgID() == -23 && gen1.getMom0PdgID()== -8000001) || (gen1.getPdgID() == 23 && gen1.getMom0PdgID()== -8000001)){ Wgen4 = gen1.getP4();}
              else if ((gen1.getPdgID() == -24 && gen1.getMom0PdgID()== 6) || (gen1.getPdgID() == 24 && gen1.getMom0PdgID()== 6)){ Wgen5 = gen1.getP4();}
              else if ((gen1.getPdgID() == 24 && gen1.getMom0PdgID()== -6) || (gen1.getPdgID() == -24 && gen1.getMom0PdgID()== -6)){ Wgen6 = gen1.getP4();}
            
	      if (goodWTaggedJets.at(i).getP4().DeltaR(Wgen1) < 0.8 || goodWTaggedJets.at(i).getP4().DeltaR(Wgen2)< 0.8 || goodWTaggedJets.at(i).getP4().DeltaR(Wgen3) < 0.8 || goodWTaggedJets.at(i).getP4().DeltaR(Wgen4)< 0.8 || goodWTaggedJets.at(i).getP4().DeltaR(Wgen5) < 0.8 || goodWTaggedJets.at(i).getP4().DeltaR(Wgen6)< 0.8 ){
                matchedW = true;
                break;
              }
            }
            if (matchedW == true){h1_["ptWmatched_st1000"] -> Fill((goodWTaggedJets.at(i)).getPt(), evtwt) ;}
            else { h1_["ptWnonmatched_st1000"] -> Fill((goodWTaggedJets.at(i)).getPt(), evtwt) ;}
          }//end Wtagged        	          





                                                                                                        
      }

    
    
      for (auto izll : zll) {
        h1_["mass_z"+lep+lep] -> Fill(izll.getMass(), evtwt) ;  
	h1_["mass_Z"+lep+lep] -> Fill(izll.getMass(), evtwt) ;
        h1_["pt_z"+lep+lep] -> Fill(izll.getPt(), evtwt) ; 
      }

      h1_["nak4"] -> Fill(goodAK4Jets.size(), evtwt) ;
      h1_["ht"] -> Fill(htak4.getHT(), evtwt) ;
      h1_["st"] -> Fill(ST, evtwt) ;   
      h1_["npv_noweight"] -> Fill(npv, *h_evtwtGen.product()); 
      h1_["npv"] -> Fill(npv, evtwt);
      h1_["met"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["met1"] -> Fill(goodMet.at(0).getFullPt(), evtwt);
      h1_["metPhi"] -> Fill(goodMet.at(0).getFullPhi(), evtwt);

      //// lepton specfic properties
      if ( zdecayMode_ == "zmumu" ){       
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d", l+1)]  -> Fill(goodMuons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d", l+1)]  -> Fill(goodMuons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_mumu"]-> Fill( (goodMuons.at(0).getP4()).DeltaR(goodMuons.at(1).getP4()), evtwt );
      }
      else if (zdecayMode_ == "zelel" ) {
        for(int l=0; l<2; ++l){
          h1_["pt_"+lep+Form("%d", l+1)]   -> Fill(goodElectrons.at(l).getPt(), evtwt) ; 
          h1_["eta_"+lep+Form("%d", l+1)]  -> Fill(goodElectrons.at(l).getEta(), evtwt) ; 
        } 
        h1_["dr_elel"]-> Fill( (goodElectrons.at(0).getP4()).DeltaR(goodElectrons.at(1).getP4()), evtwt );
      }
      //// ak4 jet plots
      for(int j=0; j<3; ++j){
        h1_[Form("ptak4jet%d", j+1)]  -> Fill(goodAK4Jets.at(j).getPt(), evtwt) ; 
        h1_[Form("etaak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getEta(), evtwt) ;
        h1_[Form("cvsak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getCSV(), evtwt) ;
	h1_[Form("massak4jet%d", j+1)] -> Fill(goodAK4Jets.at(j).getMass(), evtwt) ;
      }
      h1_["phi_jet1MET"] -> Fill( (goodAK4Jets.at(0).getP4()).DeltaPhi(goodMet.at(0).getP4()), evtwt);

      //// fill the b-tagging plots and efficiency maps
      h1_["nbjets"] -> Fill(goodBTaggedAK4Jets.size(), evtwt) ;
      h1_["ptbjetleading"] -> Fill(goodBTaggedAK4Jets.at(0).getPt(), evtwt) ;
      h1_["etabjetleading"] -> Fill(goodBTaggedAK4Jets.at(0).getEta(), evtwt) ;

      //// fill the additional plots
      h1_["nak8"] -> Fill(goodAK8Jets.size(), evtwt) ;
      h1_["nwjet"] -> Fill(goodWTaggedJets.size(), evtwt) ; 
      h1_["nhjet"] -> Fill(goodHTaggedJets.size(), evtwt) ; 
      h1_["ntjet"] -> Fill(goodTopTaggedJets.size(), evtwt) ; 

      if (goodAK8Jets.size() > 0) {
        h1_["ptak8leading"] -> Fill((goodAK8Jets.at(0)).getPt(), evtwt) ; 
        h1_["etaak8leading"] -> Fill((goodAK8Jets.at(0)).getEta(), evtwt) ;
        //h1_["mak8leading"] -> Fill((goodAK8Jets.at(0)).getMass(), evtwt) ; 
        //h1_["trimmedmak8leading"] -> Fill((goodAK8Jets.at(0)).getTrimmedMass(), evtwt) ;
        //h1_["prunedmak8leading"] -> Fill((goodAK8Jets.at(0)).getPrunedMass(), evtwt) ;
        h1_["softdropmak8leading"] -> Fill((goodAK8Jets.at(0)).getSoftDropMass(), evtwt) ;
      }
      if (goodAK8Jets.size() > 1) {
        h1_["ptak82nd"] -> Fill((goodAK8Jets.at(1)).getPt(), evtwt) ; 
        h1_["etaak82nd"] -> Fill((goodAK8Jets.at(1)).getEta(), evtwt) ;
        //h1_["mak82nd"] -> Fill((goodAK8Jets.at(1)).getMass(), evtwt) ; 
        //h1_["trimmedmak82nd"] -> Fill((goodAK8Jets.at(1)).getTrimmedMass(), evtwt) ;
        //h1_["prunedmak82nd"] -> Fill((goodAK8Jets.at(1)).getPrunedMass(), evtwt) ;
        h1_["softdropmak82nd"] -> Fill((goodAK8Jets.at(1)).getSoftDropMass(), evtwt) ;
      }


      //Do mass reconstruction TPrime
      if (categorize_){

	for (unsigned i=0; i<goodAK4Jets.size(); i++) {
          goodAK4Jetscleaned.push_back(goodAK4Jets[i]);
        }
        for (unsigned i=0; i<goodHTaggedJets.size(); i++) {
          for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
            if (goodHTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){
              ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
              goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
            }
            else{
              ak4nonmatched1.push_back(goodAK4Jetscleaned.at(j));
            }
          }
        }
        for (unsigned i=0; i<goodWTaggedJets.size(); i++) {
          for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
            if (goodWTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){
              ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
              goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
            }
            else{
              ak4nonmatched2.push_back(goodAK4Jetscleaned.at(j));
            }
          }
        }

        for (unsigned i=0; i<goodTopTaggedJets.size(); i++) {
          for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
            if (goodTopTaggedJets.at(i).getP4().DeltaR(goodAK4Jetscleaned.at(j).getP4()) < 0.8){
              ak4matchedak8.push_back(goodAK4Jetscleaned.at(j));
              goodAK4Jetscleaned.erase(goodAK4Jetscleaned.begin()+j);
            }
            else{
              ak4nonmatched3.push_back(goodAK4Jetscleaned.at(j));
            }
          }
        }



	// HCandidate cut(boosted)                                                                             
	// for (unsigned i=0; i<goodAK4Jets.size(); i++) { 
	// if (goodAK4Jets.at(i).getMass()>=80 && goodAK4Jets.at(i).getMass()<= 160 && goodAK4Jets.at(i).getPt()> 450 && goodAK4Jets.at(i).getCSV()>0.800){
	for (unsigned i=0; i<goodHTaggedJets.size(); i++) {
	  // if (goodHTaggedJets.at(i).getPt()> 450){
	  //TLorentzVector h;
	  //h= goodAK4Jets.at(i).getP4();
	  //	vlq::Candidate h1(h);
	  //Hb.push_back(h1);
	  Hb.push_back(goodHTaggedJets.at(i)); 
	  // }
	}
	//HCandidates (nonboosted)
	HCandsProducer h;
	if (goodAK4Jetscleaned.size() >= 2){
	  h.operator()(goodAK4Jetscleaned.size(), 2, goodAK4Jetscleaned,H)  ;
	}
	//Zcandidate cut (boosted)                                     
	//for (unsigned i=0; i<goodAK4Jets.size(); i++) {                                                  
	//if (goodAK4Jets.at(i).getMass()>= 70 && goodAK4Jets.at(i).getMass()<= 120 && goodAK4Jets.at(i).getPt()> 300){  
	for (unsigned i=0; i<goodWTaggedJets.size(); i++) {
	  //  if (goodWTaggedJets.at(i).getPt()> 300){
	  //	TLorentzVector zb;    
	  //zb= goodAK4Jets.at(i).getP4();
	  //vlq::Candidate zb1(zb);       
	  //	ZB.push_back(zb1);
	  ZB.push_back(goodWTaggedJets.at(i));
	  // }                                                                                                                                       
	}

	//Z candidate cut (non boosted)                                                                                                            
	ZCandsProducer z;
	if (goodAK4Jetscleaned.size() >= 2){
	  z.operator()(goodAK4Jetscleaned.size(), 2, goodAK4Jetscleaned,Z) ;
	}
	//Top Candidates (Category D)    
	// for (unsigned i=0; i<goodAK4Jets.size(); i++) {
	//if (goodAK4Jets.at(i).getMass()>= 140 && goodAK4Jets.at(i).getMass()<= 200 && goodAK4Jets.at(i).getPt()> 600){
	for (unsigned i=0; i<goodTopTaggedJets.size(); i++) {
	  //if (goodTopTaggedJets.at(i).getPt()> 600){
	  //TLorentzVector d;
	  //	d= goodAK4Jets.at(i).getP4();
	  //	vlq::Candidate d1(d);
	  //	D.push_back(d1);
	  D.push_back(goodTopTaggedJets.at(i)); 
	  //      }
	}

	TopCandsProducer top,w;
	// Category BC
	// w.operator()(goodAK4Jets, W,B) ;
	// for (unsigned i=0; i<W.size(); i++) {
	// for (unsigned j=0; j<B.size(); j++) {
	// cout << " goodWTaggedJets.size() =" << goodWTaggedJets.size() <<endl;
	// cout<< "  goodAK4Jetscleaned.size() =" <<  goodAK4Jetscleaned.size()<< endl;
	// if (goodWTaggedJets.size()>0 && goodAK4Jetscleaned.size()>0){
	for (unsigned i=0; i<goodWTaggedJets.size(); i++) {                                                                                  
	  for (unsigned j=0; j<goodAK4Jetscleaned.size(); j++) {
	    TLorentzVector bc1;
	    bc1= goodWTaggedJets.at(i).getP4()+goodAK4Jetscleaned.at(j).getP4();
	    if (bc1.Mag() >= 120 && bc1.Mag() <= 240 && bc1.Pt() >= 150){
	      vlq::Candidate bc2(bc1);
	      h1_["dr_Wb_sig"]-> Fill( (goodWTaggedJets.at(i).getP4()).DeltaR(goodAK4Jetscleaned.at(j).getP4()), evtwt);
	      h1_["dphi_Wb_sig"]-> Fill( (goodWTaggedJets.at(i).getP4()).DeltaPhi(goodAK4Jetscleaned.at(j).getP4()), evtwt);
	      W.push_back(goodWTaggedJets.at(i));
	      B.push_back(goodAK4Jetscleaned.at(j));
	      BC.push_back(bc2); 
	    }  
	  }
	}
      																       
		  
	// category A
	if (goodAK4Jetscleaned.size() >= 3){
	  top.operator()(goodAK4Jetscleaned.size(), 3, goodAK4Jetscleaned,tops) ;
	}
	for (unsigned i=0; i<Hb.size(); i++) {
	  h1_["H_mass_b_sig"] -> Fill(Hb.at(i).getPrunedMass(), evtwt) ;
	  h1_["H_Pt_b_sig"] -> Fill(Hb.at(i).getPt(), evtwt) ;
	}
	h1_["nHcandidatejets_b_sig"] -> Fill(Hb.size(), evtwt) ;
    
	// if(Hb.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}

	//Hcandidate cut(nonboosted)
	for (unsigned i=0; i<H.size(); i++) {
	  h1_["H_mass_nb_sig"] -> Fill(H.at(i).getMass(), evtwt) ;
	  h1_["H_Pt_nb_sig"] -> Fill(H.at(i).getPt(), evtwt) ;
	}
	h1_["nHcandidatejets_nb_sig"] -> Fill(H.size(), evtwt) ;

	//if(H.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
    
	double  nHcandidates=0.0;
	double  nHcandidates1=0.0;
	if(Hb.size()>0 ||H.size()>0){
	  //h1_["cutflow"] -> Fill(16, evtwt) ;
	  nHcandidates = Hb.size()+H.size();
	  // nHcandidates = H.size();
	  h1_["nHcandidatejets_sig"] -> Fill(nHcandidates, evtwt) ;
	}

	nHcandidates1 = Hb.size()+H.size();
	// nHcandidates1 = H.size();
	h1_["nHcandidatejets1_sig"] -> Fill(nHcandidates1, evtwt) ;
    
	//Zcandidate cut (boosted)
	for (unsigned i=0; i<ZB.size(); i++) {     
	  h1_["Z_mass_a_sig"] -> Fill(ZB.at(i).getPrunedMass(), evtwt) ;
	  h1_["Z_Pt_a_sig"] -> Fill(ZB.at(i).getPt(), evtwt) ;  
	}
	h1_["nzcandidatejets_a_sig"] -> Fill(ZB.size(), evtwt) ;
    
	//if(ZB.size()>0)  { h1_["cutflow"] -> Fill(10, evtwt) ;}     
    
	//Z candidate cut (non boosted)
	for (unsigned i=0; i<Z.size(); i++) {
	  h1_["Z_mass_b_sig"] -> Fill(Z.at(i).getMass(), evtwt) ;
	  h1_["Z_Pt_b_sig"] -> Fill(Z.at(i).getPt(), evtwt) ;
	}
	h1_["nzcandidatejets_b_sig"] -> Fill(Z.size(), evtwt) ;

	// if(Z.size()>0)  { h1_["cutflow"] -> Fill(11, evtwt) ;}
    
	double nzcandidates=0.0;
	double nzcandidates1=0.0;
	if ( ZB.size() || Z.size()){
	  // h1_["cutflow"] -> Fill(12, evtwt);
	  nzcandidates = ZB.size()+ Z.size();
	  // nzcandidates = Z.size();  

	  h1_["nzcandidatejets_tot_sig"] -> Fill(nzcandidates, evtwt) ;
	}
	nzcandidates1 = ZB.size()+ Z.size();
	// nzcandidates1 = Z.size();
	h1_["nzcandidatejets1_tot_sig"] -> Fill(nzcandidates1, evtwt) ;

    
	// Category D    
	for (unsigned i=0; i<D.size(); i++) {
	  h1_["top_mass_d_sig"] -> Fill(D.at(i).getSoftDropMass(), evtwt) ;
	  h1_["top_Pt_d_sig"] -> Fill(D.at(i).getPt(), evtwt) ;
	}
	h1_["ntopcandidatejets_d_sig"] -> Fill(D.size(), evtwt) ;
    
	// if(D.size()>0)  { h1_["cutflow"] -> Fill(13, evtwt) ;}


	// Category BC
	// for (unsigned i=0; i<W.size(); i++) {
	//  h1_["W_mass_bc_sig"] -> Fill(W.at(i).getMass(), evtwt) ;
	// }
	// h1_["nWcandidatejets_bc_sig"] -> Fill(W.size(), evtwt) ;
    
	// for (unsigned i=0; i<B.size(); i++) {
	//  h1_["lightjet_mass_bc_sig"] -> Fill(B.at(i).getMass(), evtwt) ; 
	// }
	// h1_["nlightjetcandidatejets_bc_sig"] -> Fill(B.size(), evtwt) ;
	//if (W.size()>0){
	for (unsigned i=0; i<W.size(); i++) {
	  h1_["W_mass_bc_sig"] -> Fill(W.at(i).getPrunedMass(), evtwt) ;
	}
	// }
	h1_["nWcandidatejets_bc_sig"] -> Fill(W.size(), evtwt) ;
	// if (B.size()>0){
	for (unsigned i=0; i<B.size(); i++) {                                                                                                
	  h1_["lightjet_mass_bc_sig"] -> Fill(B.at(i).getMass(), evtwt) ;
	}
	//  }
	h1_["nlightjetcandidatejets_bc_sig"] -> Fill(B.size(), evtwt) ;



	for (unsigned i=0; i<BC.size(); i++) {
	  h1_["top_mass_bc_sig"] -> Fill(BC.at(i).getMass(), evtwt) ;
	  h1_["top_Pt_bc_sig"] -> Fill(BC.at(i).getPt(), evtwt) ;
	}
	h1_["ntopcandidatejets_bc_sig"] -> Fill(BC.size(), evtwt) ;

	// if(BC.size()>0)  { h1_["cutflow"] -> Fill(14, evtwt) ;}
      
	for (unsigned i=0; i<tops.size(); i++) {
	  h1_["top_mass_a_sig"] -> Fill(tops.at(i).getMass(), evtwt) ;
	  h1_["top_Pt_a_sig"] -> Fill(tops.at(i).getPt(), evtwt) ;
	}
	h1_["ntopcandidatejets_a_sig"] -> Fill(tops.size(), evtwt) ;
    
	// if(tops.size()>0)  { h1_["cutflow"] -> Fill(15, evtwt) ;}
    
	double  ntopcandidates=0.0;
	double  ntopcandidates1=0.0;
	if(D.size() || BC.size() ||tops.size()){
	  //   h1_["cutflow"] -> Fill(16, evtwt) ;
	  ntopcandidates = D.size()+BC.size()+tops.size();
	  //ntopcandidates = tops.size();  
	  h1_["ntopcandidatejets_sig"] -> Fill(ntopcandidates, evtwt) ;  
	} 
	ntopcandidates1 = D.size()+BC.size()+tops.size();
	// ntopcandidates1 = tops.size();
	h1_["ntopcandidatejets1_sig"] -> Fill(ntopcandidates1, evtwt) ;

	// cout << " st befor rebin " << ST <<endl;
	//Z and top corelations and ST tempelates
	if (ST >2500.00) ST = 2498.00;    
	//cout << " st after rebin ****" << ST <<endl;
	h1_["st_sig"] -> Fill(ST,evtwt);
	h1_["cutflow3"] -> Fill(1, evtwt) ;
	if (goodBTaggedAK4Jets.size() == 1){
	  h1_["cutflow3"] -> Fill(2, evtwt) ;   
	  h1_["st_sig1b"] -> Fill(ST,evtwt);
	}
	if (goodBTaggedAK4Jets.size() >=2){
	  h1_["st_sig2b"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(3, evtwt) ;
	}

	//n,Z,H,B
	//(1)
	if (ntopcandidates >=1.0    &&   nzcandidates>=1.0){
	  h1_["st_sigT1Z1"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(4, evtwt) ;
	  if (nHcandidates >= 1.0){
	    h1_["st_sigT1Z1H1"] -> Fill(ST,evtwt);
	    h1_["cutflow3"] -> Fill(8, evtwt) ;
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow4"] -> Fill(1, evtwt) ;
	      h1_["st_sigT1Z1H1b1"] -> Fill(ST, evtwt) ;
	      
	      for (unsigned int i=0; i< goodTopTaggedJets.size();i++){
		h1_["ptT_st1000_t1Z1H1b1"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
		for (unsigned int j=0; j< goodHTaggedJets.size();j++){
		  if (goodTopTaggedJets.at(i).getIndex()== goodHTaggedJets.at(j).getIndex()){
		    double evtwtTop = evtwt;
		    h1_["ptTHnowwt_st1000_t1Z1H1b1"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
		    evtwtTop *= 1.06;
		    h1_["ptTHwt_st1000_t1Z1H1b1"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwtTop) ;
		  }
		}
	      }




	           /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT1Z1H1b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT1Z1H1b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      

	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow4"] -> Fill(2, evtwt) ;
	      h1_["st_sigT1Z1H1b2"] -> Fill(ST, evtwt) ;

	      for (unsigned int i=0; i< goodTopTaggedJets.size();i++){
		h1_["ptT_st1000_t1Z1H1b2"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
                for (unsigned int j=0; j< goodHTaggedJets.size();j++){
                  if (goodTopTaggedJets.at(i).getIndex()== goodHTaggedJets.at(j).getIndex()){
                    double evtwtTop = evtwt;
                    h1_["ptTHnowwt_st1000_t1Z1H1b2"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
                    evtwtTop *= 1.06;
                    h1_["ptTHwt_st1000_t1Z1H1b2"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwtTop) ;
                  }
                }
	      } 



	      
	        /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT1Z1H1b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT1Z1H1b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      

	    }
	  }
	  else if (nHcandidates == 0.0){
	    h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	    h1_["cutflow3"] -> Fill(9, evtwt) ;
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow4"] -> Fill(3, evtwt) ;
	      h1_["st_sigT1Z1H0b1"] -> Fill(ST, evtwt) ;
	      
	       /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT1Z1H0b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT1Z1H0b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		} */
	      
	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow4"] -> Fill(4, evtwt) ;
	      h1_["st_sigT1Z1H0b2"] -> Fill(ST, evtwt) ;
	      
	        /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT1Z1H0b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT1Z1H0b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	  }  
	}

	//(2)
	if (ntopcandidates ==0.0    &&   nzcandidates>=1.0){
	  h1_["st_sigT0Z1"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(5, evtwt) ;
	  if (nHcandidates >= 1.0){
	    h1_["st_sigT0Z1H1"] -> Fill(ST,evtwt);
	    h1_["cutflow3"] -> Fill(10, evtwt) ;
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow4"] -> Fill(5, evtwt) ;
	      h1_["st_sigT0Z1H1b1"] -> Fill(ST, evtwt) ;
	      
	         /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT0Z1H1b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT0Z1H1b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow4"] -> Fill(6, evtwt) ;
	      h1_["st_sigT0Z1H1b2"] -> Fill(ST, evtwt) ;
	      
	         /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT0Z1H1b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT0Z1H1b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	  }
	  else if (nHcandidates == 0.0){
	    h1_["st_sigT0Z1H0"] -> Fill(ST,evtwt);
	    h1_["cutflow3"] -> Fill(11, evtwt) ;
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow4"] -> Fill(7, evtwt) ;
	      h1_["st_sigT0Z1H0b1"] -> Fill(ST, evtwt) ;
	      
	        /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT0Z1H0b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT0Z1H0b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow4"] -> Fill(8, evtwt) ;
	      h1_["st_sigT0Z1H0b2"] -> Fill(ST, evtwt) ;
	      
	        /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT0Z1H0b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT0Z1H0b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	  }
	}

	//(3)
	if (ntopcandidates >=1.0    &&   nzcandidates==0.0){
	  h1_["st_sigT1Z0"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(6, evtwt) ;
	  if (nHcandidates >= 1.0){
	    h1_["st_sigT1Z0H1"] -> Fill(ST,evtwt);
	    h1_["cutflow3"] -> Fill(12, evtwt) ;
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow4"] -> Fill(9, evtwt) ;
	      h1_["st_sigT1Z0H1b1"] -> Fill(ST, evtwt) ;

	      for (unsigned int i=0; i< goodTopTaggedJets.size();i++){
		h1_["ptT_st1000_t1Z0H1b2"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
                for (unsigned int j=0; j< goodHTaggedJets.size();j++){
                  if (goodTopTaggedJets.at(i).getIndex()== goodHTaggedJets.at(j).getIndex()){
                    double evtwtTop = evtwt;
                    h1_["ptTHnowwt_st1000_t1Z0H1b2"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwt) ;
                    evtwtTop *= 1.06;
                    h1_["ptTHwt_st1000_t1Z0H1b2"] -> Fill((goodTopTaggedJets.at(i)).getPt(), evtwtTop) ;
                  }
                }
	      } 





	         /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT1Z0H1b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT1Z0H1b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow4"] -> Fill(10, evtwt) ;
	      h1_["st_sigT1Z0H1b2"] -> Fill(ST, evtwt) ;
	      
	        /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT1Z0H1b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT1Z0H1b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	  }
	  else if (nHcandidates == 0.0){
	    h1_["st_sigT1Z0H0"] -> Fill(ST,evtwt);
	    h1_["cutflow3"] -> Fill(13, evtwt) ;
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow4"] -> Fill(11, evtwt) ;
	      h1_["st_sigT1Z0H0b1"] -> Fill(ST, evtwt) ;
	      
	       /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT1Z0H0b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT1Z0H0b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      

	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow4"] -> Fill(12, evtwt) ;
	      h1_["st_sigT1Z0H0b2"] -> Fill(ST, evtwt) ;
	      
	        /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT1Z0H0b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT1Z0H0b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      

	    }
	  }
	}

	//(4)
	if (ntopcandidates == 0.0    &&   nzcandidates ==0.0){
	  h1_["st_sigT0Z0"] -> Fill(ST,evtwt);
	  h1_["cutflow3"] -> Fill(7, evtwt) ;
	  if (nHcandidates >= 1.0){
	    h1_["st_sigT0Z0H1"] -> Fill(ST,evtwt);
	    h1_["cutflow3"] -> Fill(14, evtwt) ;
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow4"] -> Fill(13, evtwt) ;
	      h1_["st_sigT0Z0H1b1"] -> Fill(ST, evtwt) ;
	      
	       /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT0Z0H1b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT0Z0H1b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow4"] -> Fill(14, evtwt) ;
	      h1_["st_sigT0Z0H1b2"] -> Fill(ST, evtwt) ;
	      
	        /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT0Z0H1b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT0Z0H1b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	  }
	  else if (nHcandidates == 0.0){
	    h1_["st_sigT0Z0H0"] -> Fill(ST,evtwt);
	    h1_["cutflow3"] -> Fill(15, evtwt) ;
	    if( goodBTaggedAK4Jets.size() == 1 ){
	      h1_["cutflow4"] -> Fill(15, evtwt) ;
	      h1_["st_sigT0Z0H0b1"] -> Fill(ST, evtwt) ;
	      
	        /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT0Z0H0b1_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT0Z0H0b1_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	    else if( goodBTaggedAK4Jets.size() >= 2 ){
	      h1_["cutflow4"] -> Fill(16, evtwt) ;
	      h1_["st_sigT0Z0H0b2"] -> Fill(ST, evtwt) ;
	      
	       /* if (!isData) {
		for (unsigned i = 0; i < 9; i++) {
		  h1_[Form("st_sigT0Z0H0b2_scale%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+scale_offset_).second );
		}
		for (unsigned i = 0; i < 100; i++) {
		  h1_[Form("st_sigT0Z0H0b2_pdf%d", i+1)] -> Fill(ST, evtwt*lhe_id_wts.at(i+pdfID_offset_).second );
		}
		}*/
	      
	    }
	  }
	}
  
	/*
	//TPrime categorization using AK8 jets

	if (ST >2500.00) ST = 2498.00;    
	//n,Z,H,B
	//(1)
	if (goodTopTaggedJets.size() >=1.0    &&   goodWTaggedJets.size()>=1.0){
	// h1_["st_sigT1Z1"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(4, evtwt) ;
	if (goodHTaggedJets.size() >= 1.0){
	//	h1_["st_sigT1Z1H1_ak8"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(8, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow6"] -> Fill(1, evtwt) ;
	h1_["st_sigT1Z1H1b1_ak8"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow6"] -> Fill(2, evtwt) ;
	h1_["st_sigT1Z1H1b2_ak8"] -> Fill(ST, evtwt) ;
	}
	}
	else if (goodHTaggedJets.size() == 0.0){
	//	h1_["st_sigT1Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(9, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow6"] -> Fill(3, evtwt) ;
	h1_["st_sigT1Z1H0b1_ak8"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow6"] -> Fill(4, evtwt) ;
	h1_["st_sigT1Z1H0b2_ak8"] -> Fill(ST, evtwt) ;
	}
	}  
	}

	//(2)
	if (goodTopTaggedJets.size() ==0.0    &&   goodWTaggedJets.size() >=1.0){
	// h1_["st_sigT0Z1"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(5, evtwt) ;
	if (goodHTaggedJets.size() >= 1.0){
	//h1_["st_sigT0Z1H1"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(10, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow6"] -> Fill(5, evtwt) ;
	h1_["st_sigT0Z1H1b1_ak8"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow6"] -> Fill(6, evtwt) ;
	h1_["st_sigT0Z1H1b2_ak8"] -> Fill(ST, evtwt) ;
	}
	}
	else if (goodHTaggedJets.size() == 0.0){
	//h1_["st_sigT0Z1H0"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(11, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow6"] -> Fill(7, evtwt) ;
	h1_["st_sigT0Z1H0b1_ak8"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow6"] -> Fill(8, evtwt) ;
	h1_["st_sigT0Z1H0b2_ak8"] -> Fill(ST, evtwt) ;
	}
	}
	}

	//(3)
	if (goodTopTaggedJets.size() >=1.0    &&   goodWTaggedJets.size()==0.0){
	// h1_["st_sigT1Z0"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(6, evtwt) ;
	if (goodHTaggedJets.size() >= 1.0){
	//	h1_["st_sigT1Z0H1"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(12, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow6"] -> Fill(9, evtwt) ;
	h1_["st_sigT1Z0H1b1_ak8"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow6"] -> Fill(10, evtwt) ;
	h1_["st_sigT1Z0H1b2_ak8"] -> Fill(ST, evtwt) ;
	}
	}
	else if (goodHTaggedJets.size() == 0.0){
	//	h1_["st_sigT1Z0H0"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(13, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow6"] -> Fill(11, evtwt) ;
	h1_["st_sigT1Z0H0b1_ak8"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow6"] -> Fill(12, evtwt) ;
	h1_["st_sigT1Z0H0b2_ak8"] -> Fill(ST, evtwt) ;
	}
	}
	}

	//(4)
	if (goodTopTaggedJets.size() == 0.0    &&   goodWTaggedJets.size() ==0.0){
	// h1_["st_sigT0Z0"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(7, evtwt) ;
	if (goodHTaggedJets.size() >= 1.0){
	//h1_["st_sigT0Z0H1"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(14, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow6"] -> Fill(13, evtwt) ;
	h1_["st_sigT0Z0H1b1_ak8"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow6"] -> Fill(14, evtwt) ;
	h1_["st_sigT0Z0H1b2_ak8"] -> Fill(ST, evtwt) ;
	}
	}
	else if (goodHTaggedJets.size() == 0.0){
	//	h1_["st_sigT0Z0H0"] -> Fill(ST,evtwt);
	h1_["cutflow5"] -> Fill(15, evtwt) ;
	if( goodBTaggedAK4Jets.size() == 1 ){
	h1_["cutflow6"] -> Fill(15, evtwt) ;
	h1_["st_sigT0Z0H0b1_ak8"] -> Fill(ST, evtwt) ;
	}
	else if( goodBTaggedAK4Jets.size() >= 2 ){
	h1_["cutflow6"] -> Fill(16, evtwt) ;
	h1_["st_sigT0Z0H0b2_ak8"] -> Fill(ST, evtwt) ;
	}
	}
	}


	*/












      }
  


      //Do mass reconstruction TPrime
  
      //  MassReco reco;
      //  CandidateFilter boostedzllfilter(BoostedZCandParams_) ;
      /*
	TLorentzVector lep1, lep2;
	if (zdecayMode_ == "zelel"){
	lep1 = goodElectrons.at(0).getP4();
	lep2 = goodElectrons.at(1).getP4();
	}
	else if (zdecayMode_ == "zmumu" ){
	lep1 = goodMuons.at(0).getP4();
	lep2 = goodMuons.at(1).getP4();
	}
	//  cout << " mass of lepton 1 = " << lep1.M()<<endl;
	// cout << " mass of lepton 2 = " << lep2.M()<<endl;
	// cout << " mass of Z boson " << (lep1+lep2).M()<<endl; 
  
	if (optimizeReco_ && *h_evttype.product() != "EvtType_Data"){
	vlq::GenParticleCollection genPartsInfo;
	genPartsInfo = genpart(evt) ;
       
	// Declare TLorentzVectors to fill with genParticles                                                                                                                                      
	TLorentzVector tGen, tbarGen, q1, q2;// Z1, Z2;                                                                                    
	TLorentzVector qJet, qbarJet, tJet, tbarJet , WJet , WbarJet, bJet;
	TLorentzVector W1Jet,Wbar1Jet,W2Jet,Wbar2Jet,b1Jet,b2Jet;  
	TLorentzVector had_tjet, lep_tjet, had_tGen, lep_tGen;

	if (bosonMass_ == 91.2){
	q1 = getGen(genPartsInfo, 1, 5, 23);
	q2 = getGen(genPartsInfo, -5, -1, 23);
	}
	else{
	q1 = getGen(genPartsInfo, 1, 5, 25);
	q2 = getGen(genPartsInfo, -5, -1, 25);
	}
	qJet = getMatchedJet(q1, goodAK4Jets, 0.3);
	qbarJet = getMatchedJet(q2, goodAK4Jets, 0.3);
 
    
	//top reconstruction from W and b
	TLorentzVector q3, q4,q5,q6, b1,b2;
	for (auto& gen1 : genPartsInfo){     
	if (gen1.getPdgID() == 24 && gen1.getMom0PdgID()==6){
	for (auto& gen : genPartsInfo){	
	if (gen.getPdgID() >=1 && gen.getPdgID() <= 4 && gen.getMom0PdgID() == 24){q3 = gen.getP4();}
	else  if (gen.getPdgID() >=-4 && gen.getPdgID() <= -1 && gen.getMom0PdgID() == 24){q4 = gen.getP4();}
	else if (gen.getPdgID() == 5 && gen.getMom0PdgID() == 6){ b1 = gen.getP4();}
	W1Jet = getMatchedJet(q3, goodAK4Jets, 0.3);
	Wbar1Jet = getMatchedJet(q4, goodAK4Jets, 0.3);
	b1Jet = getMatchedJet(b1, goodAK4Jets, 0.3);
	  
	if (b1Jet.M() == 0 && b1Jet.Pt()==0 && W1Jet.M() == 0 && W1Jet.Pt()==0 && Wbar1Jet.M() == 0 && Wbar1Jet.Pt()==0) {continue;}	
	// else if (W1Jet.M() == Wbar1Jet.M()) {continue;}
	tGen = q3+q4+b1;
	if (W1Jet.M() != Wbar1Jet.M()){tJet =W1Jet+Wbar1Jet+ b1Jet;}
	if (tGen.M() == 0 && tGen.Pt()==0 && tJet.M() == 0 && tJet.Pt()==0) {continue;}
	}
	}
    
        
	else if (gen1.getPdgID() == -24 && gen1.getMom0PdgID() == -6){
	for (auto& gen : genPartsInfo){
	if (gen.getPdgID() >=1 && gen.getPdgID() <= 4 && gen.getMom0PdgID() == -24){ q5 = gen.getP4();}
	else if (gen.getPdgID() >=-4 && gen.getPdgID() <= -1 && gen.getMom0PdgID() ==- 24){ q6 = gen.getP4();}
	else if (gen.getPdgID() == -5 && gen.getMom0PdgID() == -6){ b2 = gen.getP4();}
	W2Jet = getMatchedJet(q5, goodAK4Jets, 0.3);
	Wbar2Jet = getMatchedJet(q6, goodAK4Jets, 0.3);
	b2Jet = getMatchedJet(b2, goodAK4Jets, 0.3);
	  
	if (b2Jet.M() == 0 && b2Jet.Pt()==0 && W2Jet.M() == 0 && W2Jet.Pt()==0 && Wbar2Jet.M() == 0 && Wbar2Jet.Pt()==0) {continue;}
	// else if (W2Jet.M() == Wbar2Jet.M()) {continue;}
	tbarGen = q5+q6+b2;
	if (W2Jet.M() != Wbar2Jet.M()){tbarJet =W2Jet+Wbar2Jet+ b2Jet;}
	  
	if (tbarGen.M() == 0 && tbarGen.Pt()==0 && tbarJet.M() == 0 && tbarJet.Pt()==0 && tGen.M() == 0 && tGen.Pt()==0 && tJet.M() == 0 && tJet.Pt()==0) {continue;}
	   
	}
	}
	}
 
    
	//Higgs and Z candidates combines , only difference from Z candidate is that we have expand the mass window
	//Higgs candidates
	//case (1)non boosted

	if (q1.M() >0 && q2.M()> 0  && q1 != q2 &&  q1.M() != q2.M() && q1.Pt() > 0 && q2.Pt()>0  && (q1+q2).M()> 0 && qJet.M()> 0 && qJet.Pt()>0 && qbarJet.M()> 0 && qbarJet.Pt()>0 && qJet.M() != qbarJet.M()){ 
	if (bosonMass_ == 91.2){

	double hadgenZ = findInvMass(q1, q2);
	double hadZJet = findInvMass(qJet, qbarJet);
	double hadZJetPt = findPt(qJet, qbarJet);
      
	h1_["hadgenZ_nonb"]->Fill(hadgenZ, evtwt);
	h1_["hadZJetMass_nonb"]->Fill(hadZJet,evtwt);
	h1_["hadZJetPt_nonb"]->Fill(hadZJetPt,evtwt);
	}
	else {
	double hadgenH = findInvMass(q1, q2);
	double hadHJet = findInvMass(qJet, qbarJet);
	double hadHJetPt =findPt(qJet, qbarJet);
	h1_["hadgenH_nonb"]->Fill(hadgenH, evtwt);
	h1_["hadHJetMass_nonb"]->Fill(hadHJet,evtwt);
	h1_["hadHJetPt_nonb"]->Fill(hadHJetPt,evtwt);  
	}    
	}
	//case (2)boosted
	else if (q1.M() >0 && q2.M()> 0  && q1.Pt() > 0 && q2.Pt()>0  && (q1+q2).M()> 0 && qJet.M()> 0 && qJet.Pt()>0 && qbarJet.M()> 0 && qbarJet.Pt()>0 && qJet.M() == qbarJet.M()){ 
	if (bosonMass_ == 91.2){
	double hadgenZ = findInvMass(q1, q2);
	double hadZJet = findInvMass(qJet);
	double hadZJetPt = findPt(qJet);
	h1_["hadgenZ_b"]->Fill(hadgenZ, evtwt);
	h1_["hadZJetMass_b"]->Fill(hadZJet,evtwt);
	h1_["hadZJetPt_b"]->Fill(hadZJetPt,evtwt);
	}
	else {
	double hadgenH = findInvMass(q1, q2);
	double hadHJet = findInvMass(qJet);
	double hadHJetPt =findPt(qJet);
	h1_["hadgenH_b"]->Fill(hadgenH, evtwt);
	h1_["hadHJetMass_b"]->Fill(hadHJet,evtwt);
	h1_["hadHJetPt_b"]->Fill(hadHJetPt,evtwt);  
	}
	}
      
	//top candidates
	//Resolved Case ( 3 distinct reconstructed AK4 jets matched to 3 gene level prodicts)   
	//  if (q3.M() >0 && q4.M()> 0 && b1.M()>0 && q5.M() >0 && q6.M()> 0 && b2.M()>0 && q3 != q4 && q5 != q6 && q3.M()!= q4.M() && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()>0 &&  q5.M() != q6.M() && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0 && Wbar2Jet.Pt()>0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0 && Wbar1Jet.Pt()>0){
	if (q5.M() >0 && q6.M()> 0 && b2.M()>0  && q5 != q6 &&  q5.M() != q6.M() && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0 && W2Jet.M() != Wbar2Jet.M()  && W2Jet.M() != b2Jet.M() &&  Wbar2Jet.M() != b2Jet.M()){   
      	  
	had_tGen = tbarGen;
	had_tjet = tbarJet;    
      
	//  double hadgenT = findInvMass(q1, q2, had_tGen);
	double hadgent = findInvMass(had_tGen);
	// double hadTJet = findInvMass(qJet, qbarJet, had_tjet);
	double hadtJet = findInvMass(had_tjet);
	double hadtJetPt = findPt(had_tjet);
  
	if (bosonMass_ == 91.2){
	double hadgenZ = findInvMass(q1, q2);
	double hadZJet = findInvMass(qJet, qbarJet);
	double hadZJetPt = findPt(qJet, qbarJet);
	//	double hadTJetPt = findPt(qJet, qbarJet, had_tjet);
  
	h1_["hadgenZ"]->Fill(hadgenZ, evtwt);
	h1_["hadZJetMass"]->Fill(hadZJet,evtwt);
	h1_["hadZJetPt"]->Fill(hadZJetPt,evtwt);
	}
	else {
	double hadgenH = findInvMass(q1, q2);
	double hadHJet = findInvMass(qJet, qbarJet);
	double hadHJetPt =findPt(qJet, qbarJet);
	//	double hadTJetPt = findPt(qJet, qbarJet, had_tjet);
	h1_["hadgenH"]->Fill(hadgenH, evtwt);
	h1_["hadHJetMass"]->Fill(hadHJet,evtwt);
	h1_["hadHJetPt"]->Fill(hadHJetPt,evtwt);  
	}

	h1_["hadgentMass"] -> Fill(hadgent, evtwt);      
	h1_["hadtJetMass"] ->Fill(hadtJet, evtwt);
	h1_["hadtJetPt"] ->Fill(hadtJetPt, evtwt);
	h2_["genhadtJetMasshad"] ->Fill(hadgent,hadtJet, evtwt);
	}
    
	else if (q3.M() >0 && q4.M()> 0 && b1.M()>0  && q3 != q4 &&  q3.M() != q4.M() && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()> 0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0 && W1Jet.M() != Wbar1Jet.M() && W1Jet.M() != b1Jet.M() &&  Wbar1Jet.M() != b1Jet.M()){
      
	lep_tGen = tGen;
	lep_tjet = tJet;
      
	// double lepgenT = findInvMass(lep1, lep2, lep_tGen);
	double lepgent = findInvMass(lep_tGen);
	// double lepTJet = findInvMass(lep1, lep2, lep_tjet);
	double leptJet = findInvMass(lep_tjet);
	double leptJetPt = findPt(lep_tjet);
	if (bosonMass_ == 91.2){
	double lepZ = findInvMass(lep1, lep2);
	double lepZPt = findPt(lep1, lep2);
	//	double lepTJetPt = findPt(lep1, lep2, lep_tjet);

	h1_["lepZ"]->Fill(lepZ, evtwt);
	h1_["lepZPt"]->Fill(lepZPt,evtwt);
	}
	else {
	double lepH = findInvMass(lep1, lep2);
	double lepHPt =findPt(lep1, lep2);
	//	double lepTJetPt = findPt(lep1, lep2, lep_tjet);
	h1_["lepH"]->Fill(lepH, evtwt);
	h1_["lepHPt"]->Fill(lepHPt, evtwt);
	}
   
	h1_["lepgentMass"] -> Fill(lepgent, evtwt);
	h1_["leptJetMass"] ->Fill(leptJet, evtwt);
	h1_["leptJetPt"] ->Fill(leptJetPt, evtwt);      
	h2_["genhadtJetMasslep"] ->Fill(lepgent,leptJet, evtwt);
	}
    
	// W jets ( Category B)
	else if (q5.M() >0 && q6.M()> 0 && b2.M()>0  && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0  && W2Jet.M() == Wbar2Jet.M()  &&  W2Jet.M() != b2Jet.M() &&  Wbar2Jet.M() != b2Jet.M()){ 
      
	tbarJet =W2Jet+ b2Jet;
	had_tGen = tbarGen;
	had_tjet = tbarJet;    
      
	//  double hadgenT = findInvMass(q1, q2, had_tGen);
	double hadgent = findInvMass(had_tGen);
	// double hadTJet = findInvMass(qJet, qbarJet, had_tjet);
	double hadtJet = findInvMass(had_tjet);
	double hadWJetmass = W2Jet.M();
	double hadWJetpt = W2Jet.Pt();
	double hadtJetPt = findPt(had_tjet);
	 
	if (bosonMass_ == 91.2){
	    
	double hadgenZ = findInvMass(q1, q2);
	double hadZJet = findInvMass(qJet, qbarJet);
	double hadZJetPt = findPt(qJet, qbarJet);
	// double hadTJetPt = findPt(qJet, qbarJet, had_tjet);
	h1_["hadgenZ_B"]->Fill(hadgenZ, evtwt);
	h1_["hadZJetMass_B"]->Fill(hadZJet,evtwt);
	h1_["hadZJetPt_B"]->Fill(hadZJetPt,evtwt);
	}
	else {
	double hadgenH = findInvMass(q1, q2);
	double hadHJet = findInvMass(qJet, qbarJet);
	double hadHJetPt =findPt(qJet, qbarJet);
	//  double hadTJetPt = findPt(qJet, qbarJet, had_tjet);
	h1_["hadgenH_B"]->Fill(hadgenH, evtwt);
	h1_["hadHJetMass_B"]->Fill(hadHJet,evtwt);
	h1_["hadHJetPt_B"]->Fill(hadHJetPt,evtwt);  
	}

	// h1_["hadgenTMass_B"] -> Fill(hadgenT, evtwt);
	h1_["hadgentMass_B"] -> Fill(hadgent, evtwt);   
	h1_["hadtJetMass_B"] ->Fill(hadtJet, evtwt);
	h1_["hadWJetMass_B"] ->Fill(hadWJetmass, evtwt);
	h1_["hadWJetPt_B"] ->Fill(hadWJetpt, evtwt);
	h1_["hadtJetPt_B"] ->Fill(hadtJetPt, evtwt);
	h2_["genhadtJetMasshad_B"] ->Fill(hadgent,hadtJet, evtwt);
	h2_["hadWjetmassWjetpt_B"] ->Fill(W2Jet.M(),W2Jet.Pt(), evtwt);
	}
    
	else if (q3.M() >0 && q4.M()> 0 && b1.M()>0  && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()> 0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0 && W1Jet.M() == Wbar1Jet.M()  &&  W1Jet.M() != b1Jet.M() &&  Wbar1Jet.M() != b1Jet.M()){
      
	tbarJet =W1Jet+ b1Jet;     
	lep_tGen = tGen;
	lep_tjet = tJet;
      
	// double lepgenT = findInvMass(lep1, lep2, lep_tGen);
	double lepgent = findInvMass(lep_tGen);
	//  double lepTJet = findInvMass(lep1, lep2, lep_tjet);
	double leptJet = findInvMass(lep_tjet);
	double lepWJetmass = W1Jet.M();
	double lepWJetpt = W1Jet.Pt();
	double leptJetPt = findPt(lep_tjet);
       
	if (bosonMass_ == 91.2){
	
	double lepZ = findInvMass(lep1, lep2);
	double lepZPt = findPt(lep1, lep2);
	//	double lepTJetPt = findPt(lep1, lep2, lep_tjet);
	
	h1_["lepZ_B"]->Fill(lepZ, evtwt);
	h1_["lepZPt_B"]->Fill(lepZPt,evtwt);
	}
	else {	
	double lepH = findInvMass(lep1, lep2);
	double lepHPt =findPt(lep1, lep2);
	//	double lepTJetPt = findPt(lep1, lep2, lep_tjet);
	h1_["lepH_B"]->Fill(lepH, evtwt);
	h1_["lepHPt_B"]->Fill(lepHPt, evtwt);
	}
	h1_["lepgentMass_B"] -> Fill(lepgent, evtwt);
	h1_["lepWJetMass_B"] ->Fill(lepWJetmass, evtwt);
	h1_["lepWJetPt_B"] ->Fill(lepWJetpt, evtwt);
	h1_["leptJetMass_B"] ->Fill(leptJet, evtwt);
	h1_["leptJetPt_B"] ->Fill(leptJetPt, evtwt);
	h2_["genhadtJetMasslep_B"] ->Fill(lepgent,leptJet, evtwt);
	h2_["lepWjetmassWjetpt_B"] ->Fill(W1Jet.M(),W1Jet.Pt(), evtwt);
	}

	// Category C ( one of q jet from W overlaps with b jet, no W jet category here)

	else if (q5.M() >0 && q6.M()> 0 && b2.M()>0  && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0  &&  W2Jet.M() != Wbar2Jet.M() && ( W2Jet.M() == b2Jet.M() ||  Wbar2Jet.M() == b2Jet.M())){

	tbarJet =W2Jet+ Wbar2Jet; 
	had_tGen = tbarGen;
	had_tjet = tbarJet;

	// double hadgenT = findInvMass(q1, q2, had_tGen);
	double hadgent = findInvMass(had_tGen);
	// double hadTJet = findInvMass(qJet, qbarJet, had_tjet);
	double hadtJet = findInvMass(had_tjet);
	double hadtJetPt = findPt(had_tjet);
	double hadmergedjetmass = b2Jet.M();
	double hadmergedjetpt = b2Jet.Pt();
      
	if (bosonMass_ == 91.2){
	double hadgenZ = findInvMass(q1, q2);
	double hadZJet = findInvMass(qJet, qbarJet);
	double hadZJetPt = findPt(qJet, qbarJet);
	//	double hadTJetPt = findPt(qJet, qbarJet, had_tjet);
	h1_["hadgenZ_C"]->Fill(hadgenZ, evtwt);
	h1_["hadZJetMass_C"]->Fill(hadZJet,evtwt);
	h1_["hadZJetPt_C"]->Fill(hadZJetPt,evtwt);
	}
	else {
	double hadgenH = findInvMass(q1, q2);
	double hadHJet = findInvMass(qJet, qbarJet);
	double hadHJetPt =findPt(qJet, qbarJet);
	//	double hadTJetPt = findPt(qJet, qbarJet, had_tjet);
	h1_["hadgenH_C"]->Fill(hadgenH, evtwt);
	h1_["hadHJetMass_C"]->Fill(hadHJet,evtwt);
	h1_["hadHJetPt_C"]->Fill(hadHJetPt,evtwt);
	}
	h1_["hadgentMass_C"] -> Fill(hadgent, evtwt);	
	h1_["hadtJetMass_C"] ->Fill(hadtJet, evtwt);
	h1_["hadtJetPt_C"] ->Fill(hadtJetPt, evtwt);
	h1_["hadmergedJetMass_C"] ->Fill(hadmergedjetmass, evtwt);
	h1_["hadmergedJetPt_C"] ->Fill(hadmergedjetpt , evtwt);
	h2_["genhadtJetMasshad_C"] ->Fill(hadgent,hadtJet, evtwt);
	h2_["hadWjetmassWjetpt_C"] ->Fill(W2Jet.M(),W2Jet.Pt(), evtwt);
	}
 
	else if (q3.M() >0 && q4.M()> 0 && b1.M()>0  && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()> 0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0  &&  W1Jet.M()!= Wbar1Jet.M() && ( W1Jet.M() == b1Jet.M() ||  Wbar1Jet.M() == b1Jet.M())){

	tJet =W1Jet+ Wbar1Jet; 
	lep_tGen = tGen;
	lep_tjet = tJet;
	// double lepgenT = findInvMass(lep1, lep2, lep_tGen);
	double lepgent = findInvMass(lep_tGen);
	// double lepTJet = findInvMass(lep1, lep2, lep_tjet);
	double leptJet = findInvMass(lep_tjet);
	double leptJetPt = findPt(lep_tjet);
	double lepmergedjetmass = b1Jet.M();
	double lepmergedjetpt = b1Jet.Pt();

	if (bosonMass_ == 91.2){
        double lepZ = findInvMass(lep1, lep2);
        double lepZPt = findPt(lep1, lep2);
	// double lepTJetPt = findPt(lep1, lep2, lep_tjet);
        h1_["lepZ_C"]->Fill(lepZ, evtwt);
        h1_["lepZPt_C"]->Fill(lepZPt,evtwt);
	}
	else {
        double lepH = findInvMass(lep1, lep2);
        double lepHPt =findPt(lep1, lep2);
	// double lepTJetPt = findPt(lep1, lep2, lep_tjet);
        h1_["lepH_C"]->Fill(lepH, evtwt);
        h1_["lepHPt_C"]->Fill(lepHPt, evtwt);
	}
	h1_["lepgentMass_C"] -> Fill(lepgent, evtwt);
	h1_["leptJetMass_C"] ->Fill(leptJet, evtwt);
	h1_["leptJetPt_C"] ->Fill(leptJetPt, evtwt);
	h1_["lepmergedJetMass_C"] ->Fill(lepmergedjetmass, evtwt);
	h1_["lepmergedJetPt_C"] ->Fill(lepmergedjetpt , evtwt);
	h2_["genhadtJetMasslep_C"] ->Fill(lepgent,leptJet, evtwt);
	h2_["lepWjetmassWjetpt_C"] ->Fill(W1Jet.M(),W1Jet.Pt(), evtwt);
	}
 
	// Category D ( All three jet masses are equal , i.e. fully merged case)
	else if (q5.M() >0 && q6.M()> 0 && b2.M()>0  && q5.Pt() > 0 && q6.Pt()>0 && b2.Pt()>0 && (q5+q6).M()> 0 && b2Jet.M() > 0 && b2Jet.Pt()>0 && W2Jet.M()> 0 && W2Jet.Pt()>0 && Wbar2Jet.M()> 0  &&  W2Jet.M() == Wbar2Jet.M() &&  W2Jet.M() == b2Jet.M() &&  Wbar2Jet.M() == b2Jet.M()){
  
	tbarJet =W2Jet;
	had_tGen = tbarGen;
	had_tjet = tbarJet;

	// double hadgenT = findInvMass(q1, q2, had_tGen);
	double hadgent = findInvMass(had_tGen);
	//  double hadTJet = findInvMass(qJet, qbarJet, had_tjet);
	double hadtJet = findInvMass(had_tjet);
	double hadtJetPt = findPt(had_tjet);

	if (bosonMass_ == 91.2){
	
	double hadgenZ = findInvMass(q1, q2);
	double hadZJet = findInvMass(qJet, qbarJet);
	double hadZJetPt = findPt(qJet, qbarJet);
	//	double hadTJetPt = findPt(qJet, qbarJet, had_tjet);

	h1_["hadgenZ_D"]->Fill(hadgenZ, evtwt);
	h1_["hadZJetMass_D"]->Fill(hadZJet,evtwt);
	h1_["hadZJetPt_D"]->Fill(hadZJetPt,evtwt);
	}
	else {
	double hadgenH = findInvMass(q1, q2);
	double hadHJet = findInvMass(qJet, qbarJet);
	double hadHJetPt =findPt(qJet, qbarJet);
	//	double hadTJetPt = findPt(qJet, qbarJet, had_tjet);
	h1_["hadgenH_D"]->Fill(hadgenH, evtwt);
	h1_["hadHJetMass_D"]->Fill(hadHJet,evtwt);
	h1_["hadHJetPt_D"]->Fill(hadHJetPt,evtwt);
	}
	// h1_["hadgenTMass_D"] -> Fill(hadgenT, evtwt);
	h1_["hadgentMass_D"] -> Fill(hadgent, evtwt);
	h1_["hadtJetMass_D"] ->Fill(hadtJet, evtwt);
	h1_["hadtJetPt_D"] ->Fill(hadtJetPt, evtwt);
	h2_["genhadtJetMasshad_D"] ->Fill(hadgent,hadtJet, evtwt);
	}
    
	else if (q3.M() >0 && q4.M()> 0 && b1.M()>0  && q3.Pt() > 0 && q4.Pt()>0 && b1.Pt()>0 && (q3+q4).M()> 0 && b1Jet.M() > 0 && b1Jet.Pt()>0 && W1Jet.M()> 0 && W1Jet.Pt()>0 && Wbar1Jet.M()> 0  &&  W1Jet.M()== Wbar1Jet.M() &&  W1Jet.M() == b1Jet.M() &&  Wbar1Jet.M() == b1Jet.M()){
         
	tJet =W1Jet;
	lep_tGen = tGen;
	lep_tjet = tJet;
    
	//double lepgenT = findInvMass(lep1, lep2, lep_tGen);
	double lepgent = findInvMass(lep_tGen);
	// double lepTJet = findInvMass(lep1, lep2, lep_tjet);
	double leptJet = findInvMass(lep_tjet);
	double leptJetPt = findPt(lep_tjet);

	if (bosonMass_ == 91.2){

        double lepZ = findInvMass(lep1, lep2);
        double lepZPt = findPt(lep1, lep2);
	// double lepTJetPt = findPt(lep1, lep2, lep_tjet);

        h1_["lepZ_D"]->Fill(lepZ, evtwt);
        h1_["lepZPt_D"]->Fill(lepZPt,evtwt);
	}   
	else {
        double lepH = findInvMass(lep1, lep2);
        double lepHPt =findPt(lep1, lep2);
	// double lepTJetPt = findPt(lep1, lep2, lep_tjet);
        h1_["lepH_D"]->Fill(lepH, evtwt);
        h1_["lepHPt_D"]->Fill(lepHPt, evtwt);
	}
	h1_["lepgentMass_D"] -> Fill(lepgent, evtwt);
	h1_["leptJetMass_D"] ->Fill(leptJet, evtwt);
	h1_["leptJetPt_D"] ->Fill(leptJetPt, evtwt);
	h2_["genhadtJetMasslep_D"] ->Fill(lepgent,leptJet, evtwt);
	}
	else if (q1.M() >0 && q2.M()> 0 && q1.Pt() > 0 && q2.Pt()>0 && (q1+q2).M()> 0  && qJet.M()==qbarJet.M()){
	if (bosonMass_ == 91.2){
        double hadgenZ = findInvMass(q1, q2);
        double hadZJet = findInvMass(qJet);
        double hadZJetPt = findPt(qJet);
	// double hadTJetPt = findPt(qJet, had_tjet);
     
        h1_["hadgenZ_E"]->Fill(hadgenZ, evtwt);
        h1_["hadZJetMass_E"]->Fill(hadZJet,evtwt);
        h1_["hadZJetPt_E"]->Fill(hadZJetPt,evtwt);
	}
	}
	}

      */






    } //// Signal region 
    else return false;

  } //// if not skim  and not maketree

  if ( maketree_ ) {

    os2ltree_.clearTreeVectors();

    os2ltree_.t_signalType = signalType ; 

    os2ltree_.ta_npv = npv ; 
    os2ltree_.t_evtInfoEventNumber =  evtno ; 
    os2ltree_.t_evtInfoRunNumber = runno ; 
    os2ltree_.t_evtInfoLumiBlock = lumisec ; 

    for (vlq::Electron e : goodElectrons) {
      os2ltree_.t_elPt        .push_back(e.getPt()) ; 
      os2ltree_.t_elPhi       .push_back(e.getPhi());
      os2ltree_.t_elEta       .push_back(e.getEta());
      os2ltree_.t_elE         .push_back(e.getE());
      os2ltree_.t_elCharge    .push_back(e.getCharge());
      os2ltree_.t_elIso03     .push_back(e.getIso03());
      os2ltree_.t_elecIsLoose .push_back(e.getisLoose()); 
      os2ltree_.t_elecIsMedium.push_back(e.getisMedium()); 
      os2ltree_.t_elecIsTight .push_back(e.getisTight()); 
      os2ltree_.t_elecIsVeto  .push_back(e.getisVeto());
    }

    for (vlq::Muon m : goodMuons) {
      os2ltree_.t_muPt         .push_back(m.getPt());
      os2ltree_.t_muPhi        .push_back(m.getPhi());
      os2ltree_.t_muEta        .push_back(m.getEta());
      os2ltree_.t_muE          .push_back(m.getE());
      os2ltree_.t_muCharge     .push_back(m.getCharge());
      os2ltree_.t_muIso04      .push_back(m.getIso04());
      os2ltree_.t_muonIsTight  .push_back(m.getIsTightMuon()); 
      os2ltree_.t_muonIsLoose  .push_back(m.getIsLooseMuon()); 
      os2ltree_.t_muonIsGlobal .push_back(m.getIsGlobalMuon()); 
      os2ltree_.t_muonIsPFMuon .push_back(m.getIsPFMuon()); 
      os2ltree_.t_muonIsTracker.push_back(m.getIsTrackerMuon()); 
    }

    os2ltree_.t_ZllPt   .push_back(zll.at(0).getPt())   ;
    os2ltree_.t_ZllEta  .push_back(zll.at(0).getEta())   ;
    os2ltree_.t_ZllPhi  .push_back(zll.at(0).getPhi())   ;
    os2ltree_.t_ZllE    .push_back(zll.at(0).getEnergy())   ;
    os2ltree_.t_ZllMass .push_back(zll.at(0).getMass())   ;

    for (vlq::Jet jet : goodAK4Jets) {
      os2ltree_.t_jetAK4Pt           .push_back(jet.getPt());
      os2ltree_.t_jetAK4Phi          .push_back(jet.getPhi());
      os2ltree_.t_jetAK4Eta          .push_back(jet.getEta());
      os2ltree_.t_jetAK4E            .push_back((jet.getP4()).E());
      os2ltree_.t_jetAK4CSV          .push_back(jet.getCSV());
      os2ltree_.t_jetAK4Mass         .push_back(jet.getMass());
      os2ltree_.t_jetAK4HadronFlavour.push_back(jet.getHadronFlavour());
      os2ltree_.t_jetAK4PartonFlavour.push_back(jet.getPartonFlavour());
    }

    for (vlq::Jet bjet : goodBTaggedAK4Jets) {
      os2ltree_.t_jetAK4BPt            .push_back(bjet.getPt());
      os2ltree_.t_jetAK4BPhi           .push_back(bjet.getPhi());
      os2ltree_.t_jetAK4BEta           .push_back(bjet.getEta());
      os2ltree_.t_jetAK4BE             .push_back((bjet.getP4()).E());
      os2ltree_.t_jetAK4BCSV           .push_back(bjet.getCSV());
      os2ltree_.t_jetAK4BMass          .push_back(bjet.getMass());
      os2ltree_.t_jetAK4BHadronFlavour .push_back(bjet.getHadronFlavour());
      os2ltree_.t_jetAK4BPartonFlavour .push_back(bjet.getPartonFlavour());
    }

    for (vlq::Jet ak8 : goodAK8Jets) {
      os2ltree_.t_jetAK8Pt            .push_back(ak8.getPt());
      os2ltree_.t_jetAK8Phi           .push_back(ak8.getPhi());
      os2ltree_.t_jetAK8Eta           .push_back(ak8.getEta());
      os2ltree_.t_jetAK8E             .push_back((ak8.getP4()).E());
      os2ltree_.t_jetAK8CSV           .push_back(ak8.getCSV());
      os2ltree_.t_jetAK8Mass          .push_back(ak8.getMass());
      os2ltree_.t_jetAK8HadronFlavour .push_back(ak8.getHadronFlavour());
      os2ltree_.t_jetAK8PartonFlavour .push_back(ak8.getPartonFlavour());
      os2ltree_.t_jetAK8_tau1         .push_back(ak8.getTau1());
      os2ltree_.t_jetAK8_tau2         .push_back(ak8.getTau2());
      os2ltree_.t_jetAK8_tau3         .push_back(ak8.getTau3());
      os2ltree_.t_jetAK8_MassPruned   .push_back(ak8.getPrunedMass());
      os2ltree_.t_jetAK8_SoftDropMass .push_back(ak8.getSoftDropMass());
      os2ltree_.t_jetAK8_NSubJets     .push_back(ak8.getNSubjets());
    }

    for (vlq::Jet wjet : goodWTaggedJets) {
      os2ltree_.t_jetWJetPt            .push_back(wjet.getPt());
      os2ltree_.t_jetWJetPhi           .push_back(wjet.getPhi());
      os2ltree_.t_jetWJetEta           .push_back(wjet.getEta());
      os2ltree_.t_jetWJetE             .push_back((wjet.getP4()).E());
      os2ltree_.t_jetWJetCSV           .push_back(wjet.getCSV());
      os2ltree_.t_jetWJetMass          .push_back(wjet.getMass());
      os2ltree_.t_jetWJetHadronFlavour .push_back(wjet.getHadronFlavour());
      os2ltree_.t_jetWJetPartonFlavour .push_back(wjet.getPartonFlavour());
      os2ltree_.t_jetWJet_tau1         .push_back(wjet.getTau1());
      os2ltree_.t_jetWJet_tau2         .push_back(wjet.getTau2());
      os2ltree_.t_jetWJet_tau3         .push_back(wjet.getTau3());
      os2ltree_.t_jetWJet_MassPruned   .push_back(wjet.getPrunedMass());
      os2ltree_.t_jetWJet_SoftDropMass .push_back(wjet.getSoftDropMass());
      os2ltree_.t_jetWJet_NSubJets     .push_back(wjet.getNSubjets());
    }

    for (vlq::Jet hjet : goodHTaggedJets) {
      os2ltree_.t_jetHJetPt            .push_back(hjet.getPt());
      os2ltree_.t_jetHJetPhi           .push_back(hjet.getPhi());
      os2ltree_.t_jetHJetEta           .push_back(hjet.getEta());
      os2ltree_.t_jetHJetE             .push_back((hjet.getP4()).E());
      os2ltree_.t_jetHJetCSV           .push_back(hjet.getCSV());
      os2ltree_.t_jetHJetMass          .push_back(hjet.getMass());
      os2ltree_.t_jetHJetHadronFlavour .push_back(hjet.getHadronFlavour());
      os2ltree_.t_jetHJetPartonFlavour .push_back(hjet.getPartonFlavour());
      os2ltree_.t_jetHJet_tau1         .push_back(hjet.getTau1());
      os2ltree_.t_jetHJet_tau2         .push_back(hjet.getTau2());
      os2ltree_.t_jetHJet_tau3         .push_back(hjet.getTau3());
      os2ltree_.t_jetHJet_MassPruned   .push_back(hjet.getPrunedMass());
      os2ltree_.t_jetHJet_SoftDropMass .push_back(hjet.getSoftDropMass());
      os2ltree_.t_jetHJet_NSubJets     .push_back(hjet.getNSubjets());
    }

    for (vlq::Jet tjet : goodTopTaggedJets) {
      os2ltree_.t_jetTopJetPt            .push_back(tjet.getPt());
      os2ltree_.t_jetTopJetPhi           .push_back(tjet.getPhi());
      os2ltree_.t_jetTopJetEta           .push_back(tjet.getEta());
      os2ltree_.t_jetTopJetE             .push_back((tjet.getP4()).E());
      os2ltree_.t_jetTopJetCSV           .push_back(tjet.getCSV());
      os2ltree_.t_jetTopJetMass          .push_back(tjet.getMass());
      os2ltree_.t_jetTopJetHadronFlavour .push_back(tjet.getHadronFlavour());
      os2ltree_.t_jetTopJetPartonFlavour .push_back(tjet.getPartonFlavour());
      os2ltree_.t_jetTopJet_tau1         .push_back(tjet.getTau1());
      os2ltree_.t_jetTopJet_tau2         .push_back(tjet.getTau2());
      os2ltree_.t_jetTopJet_tau3         .push_back(tjet.getTau3());
      os2ltree_.t_jetTopJet_MassPruned   .push_back(tjet.getPrunedMass());
      os2ltree_.t_jetTopJet_SoftDropMass .push_back(tjet.getSoftDropMass());
      os2ltree_.t_jetTopJet_NSubJets     .push_back(tjet.getNSubjets());
    }

    os2ltree_.t_HT = htak4.getHT();
    os2ltree_.t_ST = ST;

    os2ltree_.t_presel_wt = presel_wt;
    os2ltree_.t_evtwt = evtwt;
    os2ltree_.t_btagsf = btagsf;
    os2ltree_.t_btagsf_bcUp = btagsf_bcUp;
    os2ltree_.t_btagsf_bcDown = btagsf_bcDown;
    os2ltree_.t_btagsf_lUp = btagsf_lUp;
    os2ltree_.t_btagsf_lDown = btagsf_lDown;
    os2ltree_.t_sjbtagsf = sjbtagsf;
    os2ltree_.t_sjbtagsf_bcUp = sjbtagsf_bcUp;
    os2ltree_.t_sjbtagsf_bcDown = sjbtagsf_bcDown;
    os2ltree_.t_sjbtagsf_lUp = sjbtagsf_lUp;
    os2ltree_.t_sjbtagsf_lDown = sjbtagsf_lDown;

    if ( !isData && filterSignal_ ) {
      vlq::GenParticleCollection vlqGen = genpart(evt) ;
      for(vlq::GenParticle p : vlqGen) { 
        os2ltree_.t_genPartPt        .push_back(p.getP4().Pt());
        os2ltree_.t_genPartPhi       .push_back(p.getP4().Phi());
        os2ltree_.t_genPartEta       .push_back(p.getP4().Eta());
        os2ltree_.t_genPartE         .push_back(p.getP4().E());
        os2ltree_.t_genPartID        .push_back(p.getPdgID());
        os2ltree_.t_genPartStatus    .push_back(p.getStatus());
        os2ltree_.t_genPartMom1ID    .push_back(p.getMom0PdgID());
        os2ltree_.t_genPartMom2ID    .push_back(p.getMom1PdgID());
        os2ltree_.t_genPartDau1ID    .push_back(p.getDau0PdgID());
        os2ltree_.t_genPartDau2ID    .push_back(p.getDau1PdgID());
        os2ltree_.t_genPartMom1Status.push_back(p.getMom0Status());
        os2ltree_.t_genPartMom2Status.push_back(p.getMom1Status());
        os2ltree_.t_genPartDau1Status.push_back(p.getDau0Status());
        os2ltree_.t_genPartDau2Status.push_back(p.getDau1Status());
      }
    }

    for( vlq::Met m : goodMet ) { 
      os2ltree_.t_metPt .push_back(m.getPt());
      os2ltree_.t_metPhi.push_back(m.getPhi());
      os2ltree_.t_metEta.push_back(m.getEta());
      os2ltree_.t_metE  .push_back(m.getP4().E());
    }

    tree_->Fill(); 
  } //// maketree

  return true ; 
}

// ------------ method called once each job just before starting event loop  ------------
void OS2LAna::beginJob() {
 std::string lep(""); 
  if(zdecayMode_ == "zmumu") {lep = "mu";}
  else if ( zdecayMode_ == "zelel") {lep = "el";}
  else edm::LogError("OS2LAna::beginJob") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;

  if (filterSignal_){
    if(skim_ || maketree_){
      const int nCh = 12;
      const char *channel[nCh] = {"bZbZ", "bZbH", "bZtW", "bHbH", "bHtW", "tWtW",
        "tZtZ", "tZtH", "tZbW", "tHtH", "tHbW", "bWbW"};
      h1_["signalEvts_all"] = fs->make<TH1D>("signalEvts_all", "All signal events", 12, 0.5, 12.5) ;
      for (int i=1;i<=nCh;i++) h1_["signalEvts_all"]->GetXaxis()->SetBinLabel(i,channel[i-1]);
    }
    else{
      h1_["signalEvts"] = fs->make<TH1D>("signalEvts", "signal events", 2, 0.5, 2.5) ;
    }
  }

  h1_["cutflow"] = fs->make<TH1D>("cutflow", "cut flow", 10, 0.5, 10.5) ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(1 , "All") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(2  , "Trig.") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(3  , "l^{+}l^{-}") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(4 , "75 #lt M(l^{+}l^{-}) #lt 105") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(5 , "H_{T} #geq 200") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(6 , "N(AK4) #geq 3") ;
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(7 , "leading jet pt > 100") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(8 , "2nd jet pt > 50") ;  
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(9 , "N(b jet) #geq 1") ; 
  h1_["cutflow"] -> GetXaxis() -> SetBinLabel(10, "S_{T} #geq 1000") ; 

  if(skim_){
    // b-tagging efficiency maps:
    h2_["pt_eta_b_all"] = fs->make<TH2D>("pt_eta_b_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_all"] = fs->make<TH2D>("pt_eta_c_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_all"] = fs->make<TH2D>("pt_eta_l_all", "b flavoured jets;p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 

    h2_["pt_eta_b_btagged"] = fs->make<TH2D>("pt_eta_b_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_c_btagged"] = fs->make<TH2D>("pt_eta_c_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ; 
    h2_["pt_eta_l_btagged"] = fs->make<TH2D>("pt_eta_l_btagged", "b flavoured jets (b-tagged);p_{T} [GeV];#eta;", 50, 0., 1000., 80, -4, 4) ;

    h1_["ht_preSel_bfFit"] = fs->make<TH1D>("ht_preSel_bfFit", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ; 
  }


  h1_["pre_pdf"] = fs->make<TH1D>("pre_pdf", "Preselection Event Yield", 1, 0.5, 1.5);

  h1_["post_pdf"] = fs->make<TH1D>("post_pdf", "Signal Event Yield", 1, 0.5, 1.5);


  if (!skim_ and !maketree_) {
    TFileDirectory pre = fs->mkdir ("pre");
    TFileDirectory sig = fs->mkdir ("sig");
    TFileDirectory cnt = fs->mkdir ("cnt");
    TFileDirectory cat = fs->mkdir ("cat");
    TFileDirectory cat1 = fs->mkdir ("cat1");

    TFileDirectory *bookDir[3]; bookDir[0] = &pre; bookDir[1] = &cnt; bookDir[2] = &sig;  bookDir[3] = &cat; bookDir[4] = &cat1;
    std::vector<string> suffix = {"_pre", "_cnt",""};
    /*    
    for (unsigned i = 0; i < 9; i++) {
      string preName_scale = Form("pre_scale%d", i+1);
      string STName_scale = Form("st_scale%d", i+1);
      string STName_scale1 = Form("st_sigT1Z1H1b1_scale%d", i+1);
      string STName_scale2 = Form("st_sigT1Z1H1b2_scale%d", i+1);
      string STName_scale3 = Form("st_sigT1Z1H0b1_scale%d", i+1);
      string STName_scale4 = Form("st_sigT1Z1H0b2_scale%d", i+1);
      string STName_scale5 = Form("st_sigT0Z1H1b1_scale%d", i+1);
      string STName_scale6 = Form("st_sigT0Z1H1b2_scale%d", i+1);
      string STName_scale7 = Form("st_sigT0Z1H0b1_scale%d", i+1);
      string STName_scale8 = Form("st_sigT0Z1H0b2_scale%d", i+1);
      string STName_scale9 = Form("st_sigT1Z0H1b1_scale%d", i+1);
      string STName_scale10 = Form("st_sigT1Z0H1b2_scale%d", i+1);
      string STName_scale11 = Form("st_sigT1Z0H0b1_scale%d", i+1);
      string STName_scale12 = Form("st_sigT1Z0H0b2_scale%d", i+1);
      string STName_scale13 = Form("st_sigT0Z0H1b1_scale%d", i+1);
      string STName_scale14 = Form("st_sigT0Z0H1b2_scale%d", i+1);
      string STName_scale15 = Form("st_sigT0Z0H0b1_scale%d", i+1);
      string STName_scale16 = Form("st_sigT0Z0H0b2_scale%d", i+1);

      //  string STName_scale1a = Form("st_cntT1Z1H1b1_scale%d", i+1);
      //  string STName_scale2a = Form("st_cntT1Z1H1b2_scale%d", i+1);
      // string STName_scale3a = Form("st_cntT1Z1H0b1_scale%d", i+1);
      // string STName_scale4a = Form("st_cntT1Z1H0b2_scale%d", i+1);
      // string STName_scale5a = Form("st_cntT0Z1H1b1_scale%d", i+1);
      // string STName_scale6a = Form("st_cntT0Z1H1b2_scale%d", i+1);
      // string STName_scale7a = Form("st_cntT0Z1H0b1_scale%d", i+1);
      // string STName_scale8a = Form("st_cntT0Z1H0b2_scale%d", i+1);
      // string STName_scale9a = Form("st_cntT1Z0H1b1_scale%d", i+1);
      // string STName_scale10a = Form("st_cntT1Z0H1b2_scale%d", i+1);
      // string STName_scale11a = Form("st_cntT1Z0H0b1_scale%d", i+1);
      //string STName_scale12a = Form("st_cntT1Z0H0b2_scale%d", i+1);
      //string STName_scale13a = Form("st_cntT0Z0H1b1_scale%d", i+1);
      //string STName_scale14a = Form("st_cntT0Z0H1b2_scale%d", i+1);
      //string STName_scale15a = Form("st_cntT0Z0H0b1_scale%d", i+1);
      //string STName_scale16a = Form("st_cntT0Z0H0b2_scale%d", i+1);






      h1_[preName_scale.c_str()] = cat1.make<TH1D>(preName_scale.c_str(), "preScale", 100, 0., 4000.);
      h1_[STName_scale.c_str()] = cat1.make<TH1D>(STName_scale.c_str(), "scaleST", 100, 0., 4000.);
      h1_[STName_scale1.c_str()] = cat1.make<TH1D>(STName_scale1.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale2.c_str()] = cat1.make<TH1D>(STName_scale2.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale3.c_str()] = cat1.make<TH1D>(STName_scale3.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale4.c_str()] = cat1.make<TH1D>(STName_scale4.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale5.c_str()] = cat1.make<TH1D>(STName_scale5.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale6.c_str()] = cat1.make<TH1D>(STName_scale6.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale7.c_str()] = cat1.make<TH1D>(STName_scale7.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale8.c_str()] = cat1.make<TH1D>(STName_scale8.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale9.c_str()] = cat1.make<TH1D>(STName_scale9.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale10.c_str()] = cat1.make<TH1D>(STName_scale10.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale11.c_str()] = cat1.make<TH1D>(STName_scale11.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale12.c_str()] = cat1.make<TH1D>(STName_scale12.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale13.c_str()] = cat1.make<TH1D>(STName_scale13.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale14.c_str()] = cat1.make<TH1D>(STName_scale14.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale15.c_str()] = cat1.make<TH1D>(STName_scale15.c_str(), "scaleST",50, 1000.,2500.);
      h1_[STName_scale16.c_str()] = cat1.make<TH1D>(STName_scale16.c_str(), "scaleST",50, 1000.,2500.);


      //   h1_[STName_scale1a.c_str()] = cat1.make<TH1D>(STName_scale1a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale2a.c_str()] = cat1.make<TH1D>(STName_scale2a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale3a.c_str()] = cat1.make<TH1D>(STName_scale3a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale4a.c_str()] = cat1.make<TH1D>(STName_scale4a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale5a.c_str()] = cat1.make<TH1D>(STName_scale5a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale6a.c_str()] = cat1.make<TH1D>(STName_scale6a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale7a.c_str()] = cat1.make<TH1D>(STName_scale7a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale8a.c_str()] = cat1.make<TH1D>(STName_scale8a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale9a.c_str()] = cat1.make<TH1D>(STName_scale9a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale10a.c_str()] = cat1.make<TH1D>(STName_scale10a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale11a.c_str()] = cat1.make<TH1D>(STName_scale11a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale12a.c_str()] = cat1.make<TH1D>(STName_scale12a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale13a.c_str()] = cat1.make<TH1D>(STName_scale13a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale14a.c_str()] = cat1.make<TH1D>(STName_scale14a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale15a.c_str()] = cat1.make<TH1D>(STName_scale15a.c_str(), "scaleST",100,0.,1000.);
      // h1_[STName_scale16a.c_str()] = cat1.make<TH1D>(STName_scale16a.c_str(), "scaleST",100,0.,1000.);




   }

    for (unsigned i = 0; i < 101; i++) {
      string preName_pdf = Form("pre_pdf%d", i+1);
      string STName_pdf = Form("st_pdf%d", i+1);

      string STName_pdf1 = Form("st_sigT1Z1H1b1_pdf%d", i+1);
      string STName_pdf2 = Form("st_sigT1Z1H1b2_pdf%d", i+1);
      string STName_pdf3 = Form("st_sigT1Z1H0b1_pdf%d", i+1);
      string STName_pdf4 = Form("st_sigT1Z1H0b2_pdf%d", i+1);
      string STName_pdf5 = Form("st_sigT0Z1H1b1_pdf%d", i+1);
      string STName_pdf6 = Form("st_sigT0Z1H1b2_pdf%d", i+1);
      string STName_pdf7 = Form("st_sigT0Z1H0b1_pdf%d", i+1);
      string STName_pdf8 = Form("st_sigT0Z1H0b2_pdf%d", i+1);
      string STName_pdf9 = Form("st_sigT1Z0H1b1_pdf%d", i+1);
      string STName_pdf10 = Form("st_sigT1Z0H1b2_pdf%d", i+1);
      string STName_pdf11 = Form("st_sigT1Z0H0b1_pdf%d", i+1);
      string STName_pdf12 = Form("st_sigT1Z0H0b2_pdf%d", i+1);
      string STName_pdf13 = Form("st_sigT0Z0H1b1_pdf%d", i+1);
      string STName_pdf14 = Form("st_sigT0Z0H1b2_pdf%d", i+1);
      string STName_pdf15 = Form("st_sigT0Z0H0b1_pdf%d", i+1);
      string STName_pdf16 = Form("st_sigT0Z0H0b2_pdf%d", i+1);

      // string STName_pdf1a = Form("st_cntT1Z1H1b1_pdf%d", i+1);
      // string STName_pdf2a = Form("st_cntT1Z1H1b2_pdf%d", i+1);
      // string STName_pdf3a = Form("st_cntT1Z1H0b1_pdf%d", i+1);
      // string STName_pdf4a = Form("st_cntT1Z1H0b2_pdf%d", i+1);
      // string STName_pdf5a = Form("st_cntT0Z1H1b1_pdf%d", i+1);
      // string STName_pdf6a = Form("st_cntT0Z1H1b2_pdf%d", i+1);
      // string STName_pdf7a = Form("st_cntT0Z1H0b1_pdf%d", i+1);
      // string STName_pdf8a = Form("st_cntT0Z1H0b2_pdf%d", i+1);
      // string STName_pdf9a = Form("st_cntT1Z0H1b1_pdf%d", i+1);
      // string STName_pdf10a = Form("st_cntT1Z0H1b2_pdf%d", i+1);
      // string STName_pdf11a = Form("st_cntT1Z0H0b1_pdf%d", i+1);
      // string STName_pdf12a = Form("st_cntT1Z0H0b2_pdf%d", i+1);
      // string STName_pdf13a = Form("st_cntT0Z0H1b1_pdf%d", i+1);
      // string STName_pdf14a = Form("st_cntT0Z0H1b2_pdf%d", i+1);
      // string STName_pdf15a = Form("st_cntT0Z0H0b1_pdf%d", i+1);
      // string STName_pdf16a = Form("st_cntT0Z0H0b2_pdf%d", i+1);




      h1_[preName_pdf.c_str()] = cat1.make<TH1D>(preName_pdf.c_str(), "prePDF", 100, 0., 4000.);
      h1_[STName_pdf.c_str()] = cat1.make<TH1D>(STName_pdf.c_str(), "pdfST", 100, 0., 4000.);

      h1_[STName_pdf1.c_str()] = cat1.make<TH1D>(STName_pdf1.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf2.c_str()] = cat1.make<TH1D>(STName_pdf2.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf3.c_str()] = cat1.make<TH1D>(STName_pdf3.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf4.c_str()] = cat1.make<TH1D>(STName_pdf4.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf5.c_str()] = cat1.make<TH1D>(STName_pdf5.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf6.c_str()] = cat1.make<TH1D>(STName_pdf6.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf7.c_str()] = cat1.make<TH1D>(STName_pdf7.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf8.c_str()] = cat1.make<TH1D>(STName_pdf8.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf9.c_str()] = cat1.make<TH1D>(STName_pdf9.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf10.c_str()] = cat1.make<TH1D>(STName_pdf10.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf11.c_str()] = cat1.make<TH1D>(STName_pdf11.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf12.c_str()] = cat1.make<TH1D>(STName_pdf12.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf13.c_str()] = cat1.make<TH1D>(STName_pdf13.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf14.c_str()] = cat1.make<TH1D>(STName_pdf14.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf15.c_str()] = cat1.make<TH1D>(STName_pdf15.c_str(), "pdfST",50, 1000.,2500.);
      h1_[STName_pdf16.c_str()] = cat1.make<TH1D>(STName_pdf16.c_str(), "pdfST",50, 1000.,2500.);

      // h1_[STName_pdf1a.c_str()] = cat1.make<TH1D>(STName_pdf1a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf2a.c_str()] = cat1.make<TH1D>(STName_pdf2a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf3a.c_str()] = cat1.make<TH1D>(STName_pdf3a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf4a.c_str()] = cat1.make<TH1D>(STName_pdf4a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf5a.c_str()] = cat1.make<TH1D>(STName_pdf5a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf6a.c_str()] = cat1.make<TH1D>(STName_pdf6a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf7a.c_str()] = cat1.make<TH1D>(STName_pdf7a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf8a.c_str()] = cat1.make<TH1D>(STName_pdf8a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf9a.c_str()] = cat1.make<TH1D>(STName_pdf9a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf10a.c_str()] = cat1.make<TH1D>(STName_pdf10a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf11a.c_str()] = cat1.make<TH1D>(STName_pdf11a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf12a.c_str()] = cat1.make<TH1D>(STName_pdf12a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf13a.c_str()] = cat1.make<TH1D>(STName_pdf13a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf14a.c_str()] = cat1.make<TH1D>(STName_pdf14a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf15a.c_str()] = cat1.make<TH1D>(STName_pdf15a.c_str(), "pdfST",100,0.,1000.);
      // h1_[STName_pdf16a.c_str()] = cat1.make<TH1D>(STName_pdf16a.c_str(), "pdfST",100,0.,1000.);





    }
    
    */
    


    for (int i=0; i<3; i++){
      h1_[("npv_noweight"+suffix[i]).c_str()] = bookDir[i]->make<TH1D>( ("npv_noweight"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("npv"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("npv"+suffix[i]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ; 
      h1_[("nak4"+suffix[i]).c_str()] =  bookDir[i]->make<TH1D>( ("nak4"+suffix[i]).c_str(), ";N(AK4 jets);;" , 21, -0.5, 20.5) ;
      h1_[("nak4Gen"+suffix[i]).c_str()] =  bookDir[i]->make<TH1D>( ("nak4Gen"+suffix[i]).c_str(), ";N(matched AK4 Gen jets);;" , 21, -0.5, 20.5) ; //Gen
      h1_[("htGen"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("htGen"+suffix[i]).c_str(), ";H_{T} (matched AK4 Gen jets) [GeV]", 100, 0., 4000.) ; //Gen
      h1_[("ht"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("ht"+suffix[i]).c_str(), ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
      h1_[("st"+suffix[i]).c_str()]   =  bookDir[i]->make<TH1D>( ("st"+suffix[i]).c_str() ,";S_{T} [GeV]", 100, 0., 4000.) ;
      h1_[("met"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("met"+suffix[i]).c_str(), "MET [GeV]", 100, 0., 1000.);
      h1_[("met1"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("met1"+suffix[i]).c_str(), ";MET [GeV]", 100, 0., 200.);   
      h1_[("metPhi"+suffix[i]).c_str()]  =  bookDir[i]->make<TH1D>( ("metPhi"+suffix[i]).c_str(), "#Phi(MET)", 20, -5., 5.);

      //jets
      for(int j=1; j<4; ++j){
        string jetPtName = Form("ptak4jet%d", j)+suffix[i]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
        h1_[jetPtName.c_str()] = bookDir[i]->make<TH1D>(jetPtName.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
        string jetEtaName = Form("etaak4jet%d", j)+suffix[i]; string jetEtaTitle  = Form(";#eta(%d leading AK4 jet) ;;",j);
        h1_[jetEtaName.c_str()] = bookDir[i]->make<TH1D>(jetEtaName.c_str(), jetEtaTitle.c_str(), 80 ,-4. ,4.) ;
        string jetCVSName = Form("cvsak4jet%d", j)+suffix[i]; string jetCVSTitle  = Form(";CVS(%d leading AK4 jet) ;;",j); 
        h1_[jetCVSName.c_str()] = bookDir[i]->make<TH1D>(jetCVSName.c_str(), jetCVSTitle.c_str(), 50 ,0. ,1.) ;
	string jetMassName = Form("massak4jet%d", j)+suffix[i]; string jetMassTitle  = Form(";Mass(%d leading AK4 jet) ;;",j);
        h1_[jetMassName.c_str()] = bookDir[i]->make<TH1D>(jetMassName.c_str(), jetMassTitle.c_str(), 100 ,0. ,1200.) ;

	//Gen
	string GenjetPtName = Form("ptak4jetGen%d", j)+suffix[i]; string GenjetPtTitle  = Form(";p_{T}(%d leading Gen AK4 jet) [GeV];;",j);
        h1_[GenjetPtName.c_str()] = bookDir[i]->make<TH1D>(GenjetPtName.c_str(), GenjetPtTitle.c_str(), 50, 0., 1000.) ;
        string GenjetEtaName = Form("etaak4jetGen%d", j)+suffix[i]; string GenjetEtaTitle  = Form(";#eta(%d leading Gen AK4 jet) ;;",j);
        h1_[GenjetEtaName.c_str()] = bookDir[i]->make<TH1D>(GenjetEtaName.c_str(), GenjetEtaTitle.c_str(), 80 ,-4. ,4.) ;


      }
      string jet1METPhiName = "phi_jet1MET"+suffix[i];
      h1_[jet1METPhiName.c_str()] = bookDir[i]->make<TH1D>(jet1METPhiName.c_str(), ";#Phi(leading jet, MET)", 20, -5., 5.) ;

      //leptons
      //  std::string lep("");
      // if(zdecayMode_ == "zmumu") {lep = "mu";}
      // else if ( zdecayMode_ == "zelel") {lep = "el";}
      // else edm::LogError("OS2LAna::beginJob") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;
      string mass_Z = "mass_z"+lep+lep+suffix[i];
      h1_[mass_Z.c_str()] = bookDir[i]->make<TH1D>(mass_Z.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 100, 20., 220.) ;
      string dr_ll = "dr_"+lep+lep+suffix[i];  
      string mass_Z1 = "mass_Z"+lep+lep+suffix[i];
      h1_[mass_Z1.c_str()] = bookDir[i]->make<TH1D>(mass_Z1.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 200, 60., 110.) ;

      h1_[dr_ll.c_str()] = bookDir[i]->make<TH1D>(dr_ll.c_str(), ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
      string pt_Z = "pt_z"+lep+lep+suffix[i];
      h1_[pt_Z.c_str()] = bookDir[i]->make<TH1D>(pt_Z.c_str(), ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ; 
      for(int l=1; l<3; ++l){
        string lepPtName = "pt_"+lep+Form("%d",l)+suffix[i]; string lepPtTitle = Form(";p_{T}(%d leading lepton) [GeV];;",l);
        h1_[lepPtName.c_str()] = bookDir[i]->make<TH1D>(lepPtName.c_str(), lepPtTitle.c_str(), 50, 0., 500.) ;
        string lepEtaName = "eta_"+lep+Form("%d",l)+suffix[i]; string lepEtaTitle  = Form(";#eta(%d leading lepton) ;;",l);
        h1_[lepEtaName.c_str()] = bookDir[i]->make<TH1D>(lepEtaName.c_str(), lepEtaTitle.c_str(), 80, -4., 4.) ;
      }
    }

    std::vector<string> suffix4 = {"_cat"};
    h1_[("npv_noweight"+suffix4[0]).c_str()] = cat.make<TH1D>( ("npv_noweight"+suffix4[0]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ;
    h1_[("npv"+suffix4[0]).c_str()]  =  cat.make<TH1D>( ("npv"+suffix4[0]).c_str(), ";N(PV);;", 51, -0.5, 50.5) ;
    h1_[("nak4"+suffix4[0]).c_str()] =  cat.make<TH1D>( ("nak4"+suffix4[0]).c_str(), ";N(AK4 jets);;" , 21, -0.5, 20.5) ;
    h1_[("ht"+suffix4[0]).c_str()]   =  cat.make<TH1D>( ("ht"+suffix4[0]).c_str(), ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_[("st"+suffix4[0]).c_str()]   =  cat.make<TH1D>( ("st"+suffix4[0]).c_str() ,";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_[("met"+suffix4[0]).c_str()]  =  cat.make<TH1D>( ("met"+suffix4[0]).c_str(), "MET [GeV]", 100, 0., 1000.);
    h1_[("met1"+suffix4[0]).c_str()]  =  cat.make<TH1D>( ("met1"+suffix4[0]).c_str(), ";MET [GeV]", 100, 0., 200.);
    h1_[("metPhi"+suffix4[0]).c_str()]  =  cat.make<TH1D>( ("metPhi"+suffix4[0]).c_str(), "#Phi(MET)", 20, -5., 5.);

    //jets                                                                                                                                                                                                                                                        
    for(int j=1; j<4; ++j){
      string jetPtName = Form("ptak4jet%d", j)+suffix4[0]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
      h1_[jetPtName.c_str()] = cat.make<TH1D>(jetPtName.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
      string jetEtaName = Form("etaak4jet%d", j)+suffix4[0]; string jetEtaTitle  = Form(";#eta(%d leading AK4 jet) ;;",j);
      h1_[jetEtaName.c_str()] = cat.make<TH1D>(jetEtaName.c_str(), jetEtaTitle.c_str(), 80 ,-4. ,4.) ;
      string jetCVSName = Form("cvsak4jet%d", j)+suffix4[0]; string jetCVSTitle  = Form(";CVS(%d leading AK4 jet) ;;",j);
      h1_[jetCVSName.c_str()] = cat.make<TH1D>(jetCVSName.c_str(), jetCVSTitle.c_str(), 50 ,0. ,1.) ;
      string jetMassName = Form("massak4jet%d", j)+suffix4[0]; string jetMassTitle  = Form(";Mass(%d leading AK4 jet) ;;",j);
      h1_[jetMassName.c_str()] = cat.make<TH1D>(jetMassName.c_str(), jetMassTitle.c_str(), 100 ,0. ,1200.) ;

    }
    string jet1METPhiName = "phi_jet1MET"+suffix4[0];
    h1_[jet1METPhiName.c_str()] = cat.make<TH1D>(jet1METPhiName.c_str(), ";#Phi(leading jet, MET)", 20, -5., 5.) ;

    //leptons                                                                                                                                                                                                                                                     
    //  std::string lep("");                                                                                                                                                                                                                                      
    // if(zdecayMode_ == "zmumu") {lep = "mu";}                                                                                                                                                                                                                   
    // else if ( zdecayMode_ == "zelel") {lep = "el";}                                                                                                                                                                                                            
    // else edm::LogError("OS2LAna::beginJob") << " >>>> WrongleptonType: " << lep << " Check lep name !!!" ;                                                                                                                                                     
    string mass_Z = "mass_z"+lep+lep+suffix4[0];
    h1_[mass_Z.c_str()] = cat.make<TH1D>(mass_Z.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 100, 20., 220.) ;
    string dr_ll = "dr_"+lep+lep+suffix4[0];
    string mass_Z1 = "mass_Z"+lep+lep+suffix4[0];
    h1_[mass_Z1.c_str()] = cat.make<TH1D>(mass_Z1.c_str(), ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 200, 60., 110.) ;

    h1_[dr_ll.c_str()] = cat.make<TH1D>(dr_ll.c_str(), ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
    string pt_Z = "pt_z"+lep+lep+suffix4[0];
    h1_[pt_Z.c_str()] = cat.make<TH1D>(pt_Z.c_str(), ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
    for(int l=1; l<3; ++l){
      string lepPtName = "pt_"+lep+Form("%d",l)+suffix4[0]; string lepPtTitle = Form(";p_{T}(%d leading lepton) [GeV];;",l);
      h1_[lepPtName.c_str()] = cat.make<TH1D>(lepPtName.c_str(), lepPtTitle.c_str(), 50, 0., 500.) ;
      string lepEtaName = "eta_"+lep+Form("%d",l)+suffix4[0]; string lepEtaTitle  = Form(";#eta(%d leading lepton) ;;",l);
      h1_[lepEtaName.c_str()] = cat.make<TH1D>(lepEtaName.c_str(), lepEtaTitle.c_str(), 80, -4., 4.) ;
    }


    // Preselection Z pt excess plots

    h1_["massz_ex_pre"] = pre.make<TH1D>("massz_ex_pre", ";M(Z#rightarrow l^{+}l^{-}) [GeV]", 200, 60., 110.) ;
    h1_["ptz_ex_pre"] = pre.make<TH1D>("ptz_ex_pre", ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 100, 500., 1000.) ;
    h1_["ptz1_ex_pre"] = pre.make<TH1D>("ptz1_ex_pre", ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
    h1_["ptz2_ex_pre"] = pre.make<TH1D>("ptz2_ex_pre", ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
    h1_["dr_elel_ex_pre"] = pre.make<TH1D>("dr_elel_ex_pre", ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
    h1_["ZPtGen_pre"] = pre.make<TH1D>("ZPtGen_pre", ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;



    h1_["ptak4jet1_ex_pre"] = pre.make<TH1D>("ptak4jet1_ex_pre", ";p_{T}(1 st AK4 jet) [GeV];;", 50, 0., 1000.) ;
    h1_["etaak4jet1_ex_pre"] = pre.make<TH1D>("etaak4jet1_ex_pre", ";#eta(1 st AK4 jet) ;;", 80 ,-4. ,4.) ;
    h1_["phiak4jet1_ex_pre"] = pre.make<TH1D>("phiak4jet1_ex_pre",";#phi(1 st AK4 jet) ;;", 20, -5., 5.) ;
    h1_["massak4jet1_ex_pre"] = pre.make<TH1D>("massak4jet1_ex_pre", ";Mass(1 st AK4 jet) ;;", 50 ,0. ,1000.) ;
    h1_["energyak4jet1_ex_pre"] = pre.make<TH1D>("energyak4jet1_ex_pre", ";Energy(1 st AK4 jet) ;;", 50 ,0. ,1000.) ;

    h1_["ptak4jet2_ex_pre"] = pre.make<TH1D>("ptak4jet2_ex_pre", ";p_{T}(2 nd AK4 jet) [GeV];;", 50, 0., 1000.) ;
    h1_["etaak4jet2_ex_pre"] = pre.make<TH1D>("etaak4jet2_ex_pre", ";#eta(2 nd AK4 jet) ;;", 80 ,-4. ,4.) ;
    h1_["phiak4jet2_ex_pre"] = pre.make<TH1D>("phiak4jet2_ex_pre",";#phi(2 nd AK4 jet) ;;", 20, -5., 5.) ;
    h1_["massak4jet2_ex_pre"] = pre.make<TH1D>("massak4jet2_ex_pre", ";Mass(2 nd AK4 jet) ;;", 50 ,0. ,1000.) ;
    h1_["energyak4jet2_ex_pre"] = pre.make<TH1D>("energyak4jet2_ex_pre", ";Energy(2 nd AK4 jet) ;;", 50 ,0. ,1000.) ;

    h1_["ptak4jet3_ex_pre"] = pre.make<TH1D>("ptak4jet3_ex_pre", ";p_{T}(3 rd AK4 jet) [GeV];;", 50, 0., 1000.) ;
    h1_["etaak4jet3_ex_pre"] = pre.make<TH1D>("etaak4jet3_ex_pre", ";#eta(3 rd AK4 jet) ;;", 80 ,-4. ,4.) ;
    h1_["phiak4jet3_ex_pre"] = pre.make<TH1D>("phiak4jet3_ex_pre",";#phi(3 rd AK4 jet) ;;", 20, -5., 5.) ;
    h1_["massak4jet3_ex_pre"] = pre.make<TH1D>("massak4jet3_ex_pre", ";Mass(3 rd AK4 jet) ;;", 50 ,0. ,1000.) ;
    h1_["energyak4jet3_ex_pre"] = pre.make<TH1D>("energyak4jet3_ex_pre", ";Energy(3 rd AK4 jet) ;;", 50 ,0. ,1000.) ;


    h1_["ptmu1_ex_pre"] = pre.make<TH1D>("ptmu1_ex_pre", ";p_{T}(lep 1) [GeV];;", 50, 0., 500.) ;
    h1_["etamu1_ex_pre"] = pre.make<TH1D>("etamu1_ex_pre", ";#eta(lep 1) ;;", 80 ,-4. ,4.) ;
    h1_["phimu1_ex_pre"] = pre.make<TH1D>("phimu1_ex_pre",";#phi(lep 1) ;;", 20, -5., 5.) ;
    h1_["energymu1_ex_pre"] = pre.make<TH1D>("energymu1_ex_pre", ";Energy(lep 1) ;;", 50 ,0. ,1000.) ;


    h1_["ptmu2_ex_pre"] = pre.make<TH1D>("ptmu2_ex_pre", ";p_{T}(lep 2) [GeV];;", 50, 0., 500.) ;
    h1_["etamu2_ex_pre"] = pre.make<TH1D>("etamu2_ex_pre", ";#eta(lep 2) ;;", 80 ,-4. ,4.) ;
    h1_["phimu2_ex_pre"] = pre.make<TH1D>("phimu2_ex_pre",";#phi(lep 2) ;;", 20, -5., 5.) ;
    h1_["energymu2_ex_pre"] = pre.make<TH1D>("energymu2_ex_pre", ";Energy(lep 2) ;;", 50 ,0. ,1000.) ;


    h1_["dphi_mu1_jet1_pre"] = pre.make<TH1D>("dphi_mu1_jet1_pre", ";#Delta #phi(mu1,jet1);;", 20, -5., 5.) ;
    h1_["dphi_mu1_jet2_pre"] = pre.make<TH1D>("dphi_mu1_jet2_pre", ";#Delta #phi(mu1,jet2);;", 20, -5., 5.) ;
    h1_["dphi_mu1_jet3_pre"] = pre.make<TH1D>("dphi_mu1_jet3_pre", ";#Delta #phi(mu1,jet3);;", 20, -5., 5.) ;

    h1_["dphi_mu2_jet1_pre"] = pre.make<TH1D>("dphi_mu2_jet1_pre", ";#Delta #phi(mu2,jet1);;", 20, -5., 5.) ;
    h1_["dphi_mu2_jet2_pre"] = pre.make<TH1D>("dphi_mu2_jet2_pre", ";#Delta #phi(mu2,jet2);;", 20, -5., 5.) ;
    h1_["dphi_mu2_jet3_pre"] = pre.make<TH1D>("dphi_mu2_jet3_pre", ";#Delta #phi(mu2,jet3);;", 20, -5., 5.) ;


    h1_["ptel1_ex_pre"] = pre.make<TH1D>("ptel1_ex_pre", ";p_{T}(lep 1) [GeV];;", 50, 0., 500.) ;
    h1_["etael1_ex_pre"] = pre.make<TH1D>("etael1_ex_pre", ";#eta(lep 1) ;;", 80 ,-4. ,4.) ;
    h1_["phiel1_ex_pre"] = pre.make<TH1D>("phiel1_ex_pre",";#phi(lep 1) ;;", 20, -5., 5.) ;
    h1_["energyel1_ex_pre"] = pre.make<TH1D>("energyel1_ex_pre", ";Energy(lep 1) ;;", 50 ,0. ,1000.) ;


    h1_["ptel2_ex_pre"] = pre.make<TH1D>("ptel2_ex_pre", ";p_{T}(lep 2) [GeV];;", 50, 0., 500.) ;
    h1_["etael2_ex_pre"] = pre.make<TH1D>("etael2_ex_pre", ";#eta(lep 2) ;;", 80 ,-4. ,4.) ;
    h1_["phiel2_ex_pre"] = pre.make<TH1D>("phiel2_ex_pre",";#phi(lep 2) ;;", 20, -5., 5.) ;
    h1_["energyel2_ex_pre"] = pre.make<TH1D>("energyel2_ex_pre", ";Energy(lep 2) ;;", 50 ,0. ,1000.) ;

    h1_["dphi_el1_jet1_pre"] = pre.make<TH1D>("dphi_el1_jet1_pre", ";#Delta #phi(el1,jet1);;", 20, -5., 5.) ;
    h1_["dphi_el1_jet2_pre"] = pre.make<TH1D>("dphi_el1_jet2_pre", ";#Delta #phi(el1,jet2);;", 20, -5., 5.) ;
    h1_["dphi_el1_jet3_pre"] = pre.make<TH1D>("dphi_el1_jet3_pre", ";#Delta #phi(el1,jet3);;", 20, -5., 5.) ;

    h1_["dphi_el2_jet1_pre"] = pre.make<TH1D>("dphi_el2_jet1_pre", ";#Delta #phi(el2,jet1);;", 20, -5., 5.) ;
    h1_["dphi_el2_jet2_pre"] = pre.make<TH1D>("dphi_el2_jet2_pre", ";#Delta #phi(el2,jet2);;", 20, -5., 5.) ;
    h1_["dphi_el2_jet3_pre"] = pre.make<TH1D>("dphi_el2_jet3_pre", ";#Delta #phi(el2,jet3);;", 20, -5., 5.) ;



    h1_["ak4pt1_cattest1"] = cat.make<TH1D>("ak4pt1_cattest1",";p_{T}(1 leading AK4 jet);;", 50, 0., 1000.) ;
    h1_["ak4eta1_cattest1"] = cat.make<TH1D>("ak4eta1_cattest1",";#eta(1 leading AK4 jet);;", 80, -4., 4.) ;
    h1_["ak4mass1_cattest1"] = cat.make<TH1D>("ak4mass1_cattest1",";M(1 leading AK4 jet);;", 100, 0., 500.) ;

    h1_["ak4pt2_cattest1"] = cat.make<TH1D>("ak4pt2_cattest1",";p_{T}(1 leading AK4 jet);;", 50, 0., 1000.) ;
    h1_["ak4eta2_cattest1"] = cat.make<TH1D>("ak4eta2_cattest1",";#eta(1 leading AK4 jet);;", 80, -4., 4.) ;
    h1_["ak4mass2_cattest1"] = cat.make<TH1D>("ak4mass2_cattest1",";M(1 leading AK4 jet);;", 100, 0., 500.) ;

    h1_["ak4pt1_cattest2"] = cat.make<TH1D>("ak4pt1_cattest2",";p_{T}(1 leading AK4 jet);;", 50, 0., 1000.) ;
    h1_["ak4eta1_cattest2"] = cat.make<TH1D>("ak4eta1_cattest2",";#eta(1 leading AK4 jet);;", 80, -4., 4.) ;
    h1_["ak4mass1_cattest2"] = cat.make<TH1D>("ak4mass1_cattest2",";M(1 leading AK4 jet);;", 100, 0., 500.) ;

    h1_["ak4pt2_cattest2"] = cat.make<TH1D>("ak4pt2_cattest2",";p_{T}(1 leading AK4 jet);;", 50, 0., 1000.) ;
    h1_["ak4eta2_cattest2"] = cat.make<TH1D>("ak4eta2_cattest2",";#eta(1 leading AK4 jet);;", 80, -4., 4.) ;
    h1_["ak4mass2_cattest2"] = cat.make<TH1D>("ak4mass2_cattest2",";M(1 leading AK4 jet);;", 100, 0., 500.) ;


  if(categorize_){
    std::vector<string> suffix1 = { "_cntT1Z1H1b1", "_cntT1Z1H1b2","_cntT1Z1H0b1","_cntT1Z1H0b2","_cntT0Z1H1b2", "_cntT1Z1H1", "_cntT1Z1H0","_cntT0Z1H1","_cntT0Z1H0","_cntT1Z0H1","_cntT1Z0H0","_cntT0Z0H1","_cntT0Z0H0"};
      
      for (int i=0; i<13; i++){
	h1_[("met"+suffix1[i]).c_str()]  =  cat.make<TH1D>( ("met"+suffix1[i]).c_str(), ";MET [GeV]", 100, 0., 1000.);
	h1_[("metPhi"+suffix1[i]).c_str()]  = cat.make<TH1D>( ("metPhi"+suffix1[i]).c_str(), ";#Phi(MET)", 20, -5., 5.);

	//jets                                                                                                                                                 
	for(int j=1; j<4; ++j){
	  string jetPtName1 = Form("ptak4jet%d", j)+suffix1[i]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
	  h1_[jetPtName1.c_str()] = cat.make<TH1D>(jetPtName1.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
	
	  string jetEtaName1 = Form("etaak4jet%d", j)+suffix1[i]; string jetEtaTitle  = Form(";#eta(%d leading AK4 jet) [GeV];;",j);
          h1_[jetEtaName1.c_str()] = cat.make<TH1D>(jetEtaName1.c_str(), jetEtaTitle.c_str(), 80, -4., 4.) ;
	  string jetPhiName1 = Form("phiak4jet%d", j)+suffix1[i]; string jetPhiTitle  = Form(";#phi(%d leading AK4 jet) [GeV];;",j);
          h1_[jetPhiName1.c_str()] = cat.make<TH1D>(jetPhiName1.c_str(), jetPhiTitle.c_str(), 20, -5., 5.) ;


	}
	string pt_Z1 = "pt_z"+lep+lep+suffix1[i];
	h1_[pt_Z1.c_str()] = cat.make<TH1D>(pt_Z1.c_str(), ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
	string eta_Z1 = "eta_z"+lep+lep+suffix1[i];
        h1_[eta_Z1.c_str()] = cat.make<TH1D>(eta_Z1.c_str(), ";#eta (Z#rightarrow l^{+}l^{-}) [GeV]", 80, -4., 4.) ;
	string phi_Z1 = "phi_z"+lep+lep+suffix1[i];
	h1_[phi_Z1.c_str()] = cat.make<TH1D>(phi_Z1.c_str(), ";#phi (Z#rightarrow l^{+}l^{-}) [GeV]", 20, -5., 5.) ;
	string dr_ll = "dr_"+lep+lep+suffix1[i];
	h1_[dr_ll.c_str()] = cat.make<TH1D>(dr_ll.c_str(),";#DeltaR(l^{+}l^{-});;", 40, 0., 4. ) ;


	for(int l=1; l<3; ++l){
	  string lepPtName1 = "pt_"+lep+Form("%d",l)+suffix1[i]; string lepPtTitle = Form(";p_{T}(%d leading lepton) [GeV];;",l);
	  h1_[lepPtName1.c_str()] = cat.make<TH1D>(lepPtName1.c_str(), lepPtTitle.c_str(), 50, 0., 500.) ;
	  string lepEtaName1 = "eta_"+lep+Form("%d",l)+suffix1[i]; string lepEtaTitle  = Form(";#eta(%d leading lepton) ;;",l);
	  h1_[lepEtaName1.c_str()] = cat.make<TH1D>(lepEtaName1.c_str(), lepEtaTitle.c_str(), 80, -4., 4.) ;
	  string lepPhiName1 = "phi_"+lep+Form("%d",l)+suffix1[i]; string lepPhiTitle  = Form(";#phi(%d leading lepton) ;;",l);
          h1_[lepPhiName1.c_str()] = cat.make<TH1D>(lepPhiName1.c_str(), lepPhiTitle.c_str(), 20, -5., 5.) ;



	}
      }
      h1_["ptbjetleading_cntT1Z1H0b1"]  = cat.make<TH1D>("ptbjetleading_cntT1Z1H0b1", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
      h1_["etabjetleading_cntT1Z1H0b1"]  = cat.make<TH1D>("etabjetleading_cntT1Z1H0b1", ";#eta (leading b jet);;" , 80, -4., 4.) ;
      h1_["phibjetleading_cntT1Z1H0b1"]  = cat.make<TH1D>("phibjetleading_cntT1Z1H0b1", ";#phi (leading b jet);;" , 20, -5., 5.) ;
      
      h1_["ptbjetleading_cntT1Z1H1b2"]  = cat.make<TH1D>("ptbjetleading_cntT1Z1H1b2", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
      h1_["etabjetleading_cntT1Z1H1b2"]  = cat.make<TH1D>("etabjetleading_cntT1Z1H1b2", ";#eta (leading b jet);;" , 80, -4., 4.) ;
      h1_["phibjetleading_cntT1Z1H1b2"]  = cat.make<TH1D>("phibjetleading_cntT1Z1H1b2", ";#phi (leading b jet);;" , 20, -5., 5.) ;





      h1_["ptbjetleading_cntT1Z1H1b1"]  = cat.make<TH1D>("ptbjetleading_cntT1Z1H1b1", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
      h1_["etabjetleading_cntT1Z1H1b1"]  = cat.make<TH1D>("etabjetleading_cntT1Z1H1b1", ";#eta (leading b jet);;" , 80, -4., 4.) ;
      h1_["phibjetleading_cntT1Z1H1b1"]  = cat.make<TH1D>("phibjetleading_cntT1Z1H1b1", ";#phi (leading b jet);;" , 20, -5., 5.) ;

      //dr plots 
      // h1_["dr_elel_cntT1Z1H1b1"] = cat.make<TH1D>("dr_elel_cntT1Z1H1b1", ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
      h1_["dr_el1jet1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el1jet1_cntT1Z1H1b1", ";#DeltaR(lep1Jet1);;", 40, 0., 4.) ;  
      h1_["dr_el1jet2_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el1jet2_cntT1Z1H1b1", ";#DeltaR(lep1Jet2);;", 40, 0., 4.) ;
      h1_["dr_el1jet3_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el1jet3_cntT1Z1H1b1", ";#DeltaR(lep1Jet3);;", 40, 0., 4.) ;
      h1_["dr_el2jet1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el2jet1_cntT1Z1H1b1", ";#DeltaR(lep2Jet1);;", 40, 0., 4.) ;
      h1_["dr_el2jet2_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el2jet2_cntT1Z1H1b1", ";#DeltaR(lep2Jet2);;", 40, 0., 4.) ;
      h1_["dr_el2jet3_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el2jet3_cntT1Z1H1b1", ";#DeltaR(lep2Jet3);;", 40, 0., 4.) ;
      h1_["dr_el1bjet1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el1bjet1_cntT1Z1H1b1", ";#DeltaR(lep1bJet1);;", 40, 0., 4.) ;
      h1_["dr_el2bjet1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el2bjet1_cntT1Z1H1b1", ";#DeltaR(lep2bJet1);;", 40, 0., 4.) ;
      h1_["dr_el1Hb1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el1Hb1_cntT1Z1H1b1", ";#DeltaR(lep1HBoosted1);;", 40, 0., 4.) ;
      h1_["dr_el2Hb1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el2Hb1_cntT1Z1H1b1", ";#DeltaR(lep2HBoosted1);;", 40, 0., 4.) ;
      h1_["dr_el1Zb1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el1Zb1_cntT1Z1H1b1", ";#DeltaR(lep1ZBossted1);;", 40, 0., 4.) ;
      h1_["dr_el2Zb1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el2Zb1_cntT1Z1H1b1", ";#DeltaR(lep2ZBossted1);;", 40, 0., 4.) ;

      h1_["dr_el1tb1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el1tb1_cntT1Z1H1b1", ";#DeltaR(lep1topBoosted1);;", 40, 0., 4.) ;
      h1_["dr_el2tb1_cntT1Z1H1b1"] = cat.make<TH1D>("dr_el2tb1_cntT1Z1H1b1", ";#DeltaR(lep2topBoosted1);;", 40, 0., 4.) ;

      

      h1_["dr_el1jet1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el1jet1_cntT1Z1H0b1", ";#DeltaR(lep1Jet1);;", 40, 0., 4.) ;
      h1_["dr_el1jet2_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el1jet2_cntT1Z1H0b1", ";#DeltaR(lep1Jet2);;", 40, 0., 4.) ;
      h1_["dr_el1jet3_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el1jet3_cntT1Z1H0b1", ";#DeltaR(lep1Jet3);;", 40, 0., 4.) ;
      h1_["dr_el2jet1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el2jet1_cntT1Z1H0b1", ";#DeltaR(lep2Jet1);;", 40, 0., 4.) ;
      h1_["dr_el2jet2_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el2jet2_cntT1Z1H0b1", ";#DeltaR(lep2Jet2);;", 40, 0., 4.) ;
      h1_["dr_el2jet3_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el2jet3_cntT1Z1H0b1", ";#DeltaR(lep2Jet3);;", 40, 0., 4.) ;
      h1_["dr_el1bjet1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el1bjet1_cntT1Z1H0b1", ";#DeltaR(lep1bJet1);;", 40, 0., 4.) ;
      h1_["dr_el2bjet1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el2bjet1_cntT1Z1H0b1", ";#DeltaR(lep2bJet1);;", 40, 0., 4.) ;
      h1_["dr_el1Hb1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el1Hb1_cntT1Z1H0b1", ";#DeltaR(lep1HBoosted1);;", 40, 0., 4.) ;
      h1_["dr_el2Hb1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el2Hb1_cntT1Z1H0b1", ";#DeltaR(lep2HBoosted1);;", 40, 0., 4.) ;
      h1_["dr_el1Zb1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el1Zb1_cntT1Z1H0b1", ";#DeltaR(lep1ZBossted1);;", 40, 0., 4.) ;
      h1_["dr_el2Zb1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el2Zb1_cntT1Z1H0b1", ";#DeltaR(lep2ZBossted1);;", 40, 0., 4.) ;

      h1_["dr_el1tb1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el1tb1_cntT1Z1H0b1", ";#DeltaR(lep1topBoosted1);;", 40, 0., 4.) ;
      h1_["dr_el2tb1_cntT1Z1H0b1"] = cat.make<TH1D>("dr_el2tb1_cntT1Z1H0b1", ";#DeltaR(lep2topBoosted1);;", 40, 0., 4.) ;



      h1_["dr_el1jet1_pre"] = pre.make<TH1D>("dr_el1jet1_pre", ";#DeltaR(lep1Jet1);;", 80, 0., 8.) ;
      h1_["dr_el1jet2_pre"] = pre.make<TH1D>("dr_el1jet2_pre", ";#DeltaR(lep1Jet2);;", 80, 0., 8.) ;
      h1_["dr_el1jet3_pre"] = pre.make<TH1D>("dr_el1jet3_pre", ";#DeltaR(lep1Jet3);;", 80, 0., 8.) ;
      h1_["dr_el2jet1_pre"] = pre.make<TH1D>("dr_el2jet1_pre", ";#DeltaR(lep2Jet1);;", 80, 0., 8.) ;
      h1_["dr_el2jet2_pre"] = pre.make<TH1D>("dr_el2jet2_pre", ";#DeltaR(lep2Jet2);;", 80, 0., 8.) ;
      h1_["dr_el2jet3_pre"] = pre.make<TH1D>("dr_el2jet3_pre", ";#DeltaR(lep2Jet3);;", 80, 0., 8.) ;

      h1_["dr_mu1jet1_pre"] = pre.make<TH1D>("dr_mu1jet1_pre", ";#DeltaR(lep1Jet1);;", 80, 0., 8.) ;
      h1_["dr_mu1jet2_pre"] = pre.make<TH1D>("dr_mu1jet2_pre", ";#DeltaR(lep1Jet2);;", 80, 0., 8.) ;
      h1_["dr_mu1jet3_pre"] = pre.make<TH1D>("dr_mu1jet3_pre", ";#DeltaR(lep1Jet3);;", 80, 0., 8.) ;
      h1_["dr_mu2jet1_pre"] = pre.make<TH1D>("dr_mu2jet1_pre", ";#DeltaR(lep2Jet1);;", 80, 0., 8.) ;
      h1_["dr_mu2jet2_pre"] = pre.make<TH1D>("dr_mu2jet2_pre", ";#DeltaR(lep2Jet2);;", 80, 0., 8.) ;
      h1_["dr_mu2jet3_pre"] = pre.make<TH1D>("dr_mu2jet3_pre", ";#DeltaR(lep2Jet3);;", 80, 0., 8.) ;

      h1_["dr_el1minjet_pre"] = pre.make<TH1D>("dr_el1minjet_pre", ";#DeltaR(lep1,closest jet);;", 80, 0., 8.) ;
      h1_["dr_el2minjet_pre"] = pre.make<TH1D>("dr_el2minjet_pre", ";#DeltaR(lep2,closest jet);;", 80, 0., 8.) ;
      h1_["dr_mu1minjet_pre"] = pre.make<TH1D>("dr_mu1minjet_pre", ";#DeltaR(lep1,closest jet);;", 80, 0., 8.) ;
      h1_["dr_mu2minjet_pre"] = pre.make<TH1D>("dr_mu2minjet_pre", ";#DeltaR(lep2,closest jet);;", 80, 0., 8.) ;


      h1_["ptbjetleading_cntT1Z1H0b2"]  = cat.make<TH1D>("ptbjetleading_cntT1Z1H0b2", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
      h1_["etabjetleading_cntT1Z1H0b2"]  = cat.make<TH1D>("etabjetleading_cntT1Z1H0b2", ";#eta (leading b jet);;" , 80, -4., 4.) ;
      h1_["phibjetleading_cntT1Z1H0b2"]  = cat.make<TH1D>("phibjetleading_cntT1Z1H0b2", ";#phi (leading b jet);;" , 20, -5., 5.) ;



      h1_["ptbjetleading_cntT0Z1H1b2"]  = cat.make<TH1D>("ptbjetleading_cntT0Z1H1b2", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.);

      h1_["ptbjetsubleading_cntT1Z1H1b2"]  = cat.make<TH1D>("ptbjetsubleading_cntT1Z1H1b2", ";p_{T}(subleadinleading b jet) [GeV];;" , 50, 0., 1000.) ;
      h1_["etabjetsubleading_cntT1Z1H1b2"]  = cat.make<TH1D>("etabjetsubleading_cntT1Z1H1b2", ";#eta (leading b jet);;" , 80, -4., 4.) ;
      h1_["phibjetsubleading_cntT1Z1H1b2"]  = cat.make<TH1D>("phibjetsubleading_cntT1Z1H1b2", ";#phi (leading b jet);;" , 20, -5., 5.) ;


      h1_["ptbjetsubleading_cntT1Z1H0b2"]  = cat.make<TH1D>("ptbjetsubleading_cntT1Z1H0b2", ";p_{T}(subleadinleading b jet) [GeV];;" , 50, 0., 1000.) ;
      h1_["etabjetsubleading_cntT1Z1H0b2"]  = cat.make<TH1D>("etabjetsubleading_cntT1Z1H0b2", ";#eta (leading b jet);;" , 80, -4., 4.) ;
      h1_["phibjetsubleading_cntT1Z1H0b2"]  = cat.make<TH1D>("phibjetsubleading_cntT1Z1H0b2", ";#phi (leading b jet);;" , 20, -5., 5.) ;


      h1_["ptbjetsubleading_cntT0Z1H1b2"]  = cat.make<TH1D>("ptbjetsubleading_cntT0Z1H1b2", ";p_{T}(subleadinleading b jet) [GeV];;" , 50, 0., 1000.) ;
    }
  
    std::vector<string> suffix2 = { "_st1000_e0b","_st1000_e1b","_st1000_1b","_st1000_2b","_cntT1Z1Hprime1b0","_0b1","_0b2","_0b3"};

    for (int i=0; i<8; i++){
      h1_[("met"+suffix2[i]).c_str()]  =  cat1.make<TH1D>( ("met"+suffix2[i]).c_str(), ";MET [GeV]", 100, 0., 1000.);
      h1_[("st"+suffix2[i]).c_str()]  =  cat1.make<TH1D>( ("st"+suffix2[i]).c_str(), ";ST [GeV]", 100, 0., 1000.);
      h1_[("ht"+suffix2[i]).c_str()]  =  cat1.make<TH1D>( ("ht"+suffix2[i]).c_str(), ";HT [GeV]", 100, 0., 1000.);
      //jets                                                                                                                                                         
      for(int j=1; j<4; ++j){
	string jetPtName2 = Form("ptak4jet%d", j)+suffix2[i]; string jetPtTitle  = Form(";p_{T}(%d leading AK4 jet) [GeV];;",j);
	h1_[jetPtName2.c_str()] = cat1.make<TH1D>(jetPtName2.c_str(), jetPtTitle.c_str(), 50, 0., 1000.) ;
      }
      string pt_Z2 = "pt_z"+lep+lep+suffix2[i];
      h1_[pt_Z2.c_str()] = cat1.make<TH1D>(pt_Z2.c_str(), ";p_{T} (Z#rightarrow l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;

      for(int l=1; l<3; ++l){
	string lepPtName2 = "pt_"+lep+Form("%d",l)+suffix2[i]; string lepPtTitle = Form(";p_{T}(%d leading lepton) [GeV];;",l);
	h1_[lepPtName2.c_str()] = cat1.make<TH1D>(lepPtName2.c_str(), lepPtTitle.c_str(), 50, 0., 500.) ;
	string lepEtaName2 = "eta_"+lep+Form("%d",l)+suffix2[i]; string lepEtaTitle  = Form(";#eta(%d leading lepton) ;;",l);
	h1_[lepEtaName2.c_str()] = cat1.make<TH1D>(lepEtaName2.c_str(), lepEtaTitle.c_str(), 80, -4., 4.) ;
      }
    }
    h1_["ptbjetleading_st1000_e1b"]  = cat1.make<TH1D>("ptbjetleading_st1000_e1b", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["ptbjetsubleading_st1000_e1b"]  = cat1.make<TH1D>("ptbjetsubleading_st1000_e1b", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["ptbjetleading_st1000_1b"]  = cat1.make<TH1D>("ptbjetleading_st1000_1b", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["ptbjetleading_st1000_2b"]  = cat1.make<TH1D>("ptbjetleading_st1000_2b", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["ptbjetsubleading_st1000_2b"]  = cat1.make<TH1D>("ptbjetsubleading_st1000_2b", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
  
  


    h1_["ht1_cnt"]   =  cnt.make<TH1D>("ht1_cnt", ";H_{T} (AK4 jets) [GeV]", 200, 0., 1000.) ;
    h1_["st1_cnt"]   =  cnt.make<TH1D>("st1_cnt",";S_{T} [GeV]", 200, 0., 1000.) ;
    h1_["ht1_cat"]   =  cat.make<TH1D>("ht1_cat", ";H_{T} (AK4 jets) [GeV]", 200, 0., 1000.) ;
    h1_["st1_cat"]   =  cat.make<TH1D>("st1_cat",";S_{T} [GeV]", 200, 0., 1000.) ;  
    h1_["ht1_0cat"]   =  cat.make<TH1D>("ht1_0cat", ";H_{T} (AK4 jets) [GeV]", 200, 0., 1000.) ;
    h1_["st1_0cat"]   =  cat.make<TH1D>("st1_0cat",";S_{T} [GeV]", 200, 0., 1000.) ;
    h1_["nak4_0cat"] = cat.make<TH1D>("nak4_0cat", ";N(AK4 jets);;" , 11, -0.5,10.5) ;
    h1_["dr_mumu_0cat"] = cat.make<TH1D>("dr_mumu_0cat", ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;
    h1_["dr_elel_0cat"] = cat.make<TH1D>("dr_elel_0cat", ";#DeltaR(l^{+}l^{-});;", 40, 0., 4.) ;

    //additional plots
    h1_["nbjets"] = sig.make<TH1D>("nbjets", ";N(b jets);;" , 11, -0.5,10.5) ; 
    h1_["ptbjetleading"]  = sig.make<TH1D>("ptbjetleading", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etabjetleading"] = sig.make<TH1D>("etabjetleading", ";#eta(leading b jet);;" , 80 ,-4. ,4.) ;

    h1_["nak8"] = sig.make<TH1D>("nak8", ";N(AK8 jets);;" , 11, -0.5,10.5) ; 
    h1_["nwjet"] = sig.make<TH1D>("nwjet", ";N(W jets );;" , 6, -0.5,5.5) ; 
    h1_["nhjet"] = sig.make<TH1D>("nhjet", ";N(H jets );;" , 6, -0.5,5.5) ; 
    h1_["ntjet"] = sig.make<TH1D>("ntjet", ";N(top jets);;" , 6, -0.5,5.5) ; 

    h1_["ptak8leading"]  = sig.make<TH1D>("ptak8leading", ";p_{T}(leading AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak8leading"] = sig.make<TH1D>("etaak8leading", ";#eta(leading AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak8leading"] = sig.make<TH1D>("softdropmak8leading", ";M(leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["ptak82nd"]  = sig.make<TH1D>("ptak82nd", ";p_{T}(2nd AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak82nd"] = sig.make<TH1D>("etaak82nd", ";#eta(2nd AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak82nd"] = sig.make<TH1D>("softdropmak82nd", ";M(2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ;


    //AK8 jets pre
    h1_["Wptleading_pre"]  = pre.make<TH1D>("Wptleading_pre", ";p_{T}(leading W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Wetaleading_pre"] = pre.make<TH1D>("Wetaleading_pre", ";#eta(leading W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wprunedleading_pre"] = pre.make<TH1D>("Wprunedleading_pre", ";M(leading W jet) [GeV];;" ,200 ,50., 150.) ;
     
    h1_["Wpt2nd_pre"]  = pre.make<TH1D>("Wpt2nd_pre", ";p_{T}(2nd W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta2nd_pre"] = pre.make<TH1D>("Weta2nd_pre", ";#eta(2nd W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned2nd_pre"] = pre.make<TH1D>("Wpruned2nd_pre", ";M(2nd W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wpt_pre"]  = pre.make<TH1D>("Wpt_pre", ";p_{T}(W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta_pre"] = pre.make<TH1D>("Weta_pre", ";#eta(W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned_pre"] = pre.make<TH1D>("Wpruned_pre", ";Pruned Mass(W jet) [GeV];;" ,200 ,50., 150.) ;


    h1_["Wptleading1_pre"]  = pre.make<TH1D>("Wptleading1_pre", ";p_{T}(leading W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Wetaleading1_pre"] = pre.make<TH1D>("Wetaleading1_pre", ";#eta(leading W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wprunedleading1_pre"] = pre.make<TH1D>("Wprunedleading1_pre", ";M(leading W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wpt2nd1_pre"]  = pre.make<TH1D>("Wpt2nd1_pre", ";p_{T}(2nd W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta2nd1_pre"] = pre.make<TH1D>("Weta2nd1_pre", ";#eta(2nd W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned2nd1_pre"] = pre.make<TH1D>("Wpruned2nd1_pre", ";M(2nd W jet) [GeV];;" ,200 ,50., 150.) ;


    h1_["ptsmear1_pre"]  = pre.make<TH1D>("ptsmear1_pre", ";ptsmear;;" , 200, -0.5, 1.5) ;
    h1_["massCorr1_pre"]  = pre.make<TH1D>("massCorr1_pre", ";massCorr;;" , 200, -0.5, 1.5) ;
    h1_["masssmear1_pre"]  = pre.make<TH1D>("masssmear1_pre", ";masssmear;;" , 200, -0.5, 1.5) ;


    h1_["ptsmear2_pre"]  = pre.make<TH1D>("ptsmear2_pre", ";ptsmear;;" , 200, -0.5, 1.5) ;
    h1_["massCorr2_pre"]  = pre.make<TH1D>("massCorr2_pre", ";massCorr;;" , 200, -0.5, 1.5) ;
    h1_["masssmear2_pre"]  = pre.make<TH1D>("masssmear2_pre", ";masssmear;;" , 200, -0.5, 1.5) ;
    h1_["masssmear21_pre"]  = pre.make<TH1D>("masssmear21_pre", ";masssmear;;" , 400, -0.5, 3.5) ;

    h1_["Hptleading_pre"]  = pre.make<TH1D>("Hptleading_pre", ";p_{T}(leading H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Hetaleading_pre"] = pre.make<TH1D>("Hetaleading_pre", ";#eta(leading H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hprunedleading_pre"] = pre.make<TH1D>("Hprunedleading_pre", ";M(leading H jet) [GeV];;" ,100 ,100., 150.) ;

    h1_["Hpt2nd_pre"]  = pre.make<TH1D>("Hpt2nd_pre", ";p_{T}(2nd H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Heta2nd_pre"] = pre.make<TH1D>("Heta2nd_pre", ";#eta(2nd H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hpruned2nd_pre"] = pre.make<TH1D>("Hpruned2nd_pre", ";M(2nd H jet) [GeV];;" ,100 ,100., 150.) ;

    h1_["Hpt_pre"]  = pre.make<TH1D>("Hpt_pre", ";p_{T}(H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Heta_pre"] = pre.make<TH1D>("Heta_pre", ";#eta(H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hpruned_pre"] = pre.make<TH1D>("Hpruned_pre", ";Pruned Mass(H jet) [GeV];;" ,100 ,100., 150.) ;

    h1_["Topptleading_pre"]  = pre.make<TH1D>("Topptleading_pre", ";p_{T}(leading Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topetaleading_pre"] = pre.make<TH1D>("Topetaleading_pre", ";#eta(leading Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdropleading_pre"] = pre.make<TH1D>("Topsoftdropleading_pre", ";M(leading Top jet) [GeV];;" ,200 ,80., 280.) ;

    h1_["Toppt2nd_pre"]  = pre.make<TH1D>("Toppt2nd_pre", ";p_{T}(2nd Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topeta2nd_pre"] = pre.make<TH1D>("Topeta2nd_pre", ";#eta(2nd Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdrop2nd_pre"] = pre.make<TH1D>("Topsoftdrop2nd_pre", ";M(2nd Top jet) [GeV];;" ,200 ,80., 280.) ;

    h1_["Toppt_pre"]  = pre.make<TH1D>("Toppt_pre", ";p_{T}(Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topeta_pre"] = pre.make<TH1D>("Topeta_pre", ";#eta(Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdrop_pre"] = pre.make<TH1D>("Topsoftdrop_pre", ";SoftDrop Mass(Top jet) [GeV];;" ,200 ,80., 280.) ;

    //AK8 jets cnt                                                                                                                                                  
    h1_["Wptleading_cnt"]  = cnt.make<TH1D>("Wptleading_cnt", ";p_{T}(leading W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Wetaleading_cnt"] = cnt.make<TH1D>("Wetaleading_cnt", ";#eta(leading W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wprunedleading_cnt"] = cnt.make<TH1D>("Wprunedleading_cnt", ";M(leading W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wpt2nd_cnt"]  = cnt.make<TH1D>("Wpt2nd_cnt", ";p_{T}(2nd W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta2nd_cnt"] = cnt.make<TH1D>("Weta2nd_cnt", ";#eta(2nd W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned2nd_cnt"] = cnt.make<TH1D>("Wpruned2nd_cnt", ";M(2nd W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wpt_cnt"]  = cnt.make<TH1D>("Wpt_cnt", ";p_{T}(W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta_cnt"] = cnt.make<TH1D>("Weta_cnt", ";#eta(W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned_cnt"] = cnt.make<TH1D>("Wpruned_cnt", ";Pruned Mass(W jet) [GeV];;" ,200 ,50., 150.) ;


    h1_["Wptleading1_cnt"]  = cnt.make<TH1D>("Wptleading1_cnt", ";p_{T}(leading W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Wetaleading1_cnt"] = cnt.make<TH1D>("Wetaleading1_cnt", ";#eta(leading W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wprunedleading1_cnt"] = cnt.make<TH1D>("Wprunedleading1_cnt", ";M(leading W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wpt2nd1_cnt"]  = cnt.make<TH1D>("Wpt2nd1_cnt", ";p_{T}(2nd W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta2nd1_cnt"] = cnt.make<TH1D>("Weta2nd1_cnt", ";#eta(2nd W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned2nd1_cnt"] = cnt.make<TH1D>("Wpruned2nd1_cnt", ";M(2nd W jet) [GeV];;" ,200 ,50., 150.) ;



    h1_["Hptleading_cnt"]  = cnt.make<TH1D>("Hptleading_cnt", ";p_{T}(leading H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Hetaleading_cnt"] = cnt.make<TH1D>("Hetaleading_cnt", ";#eta(leading H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hprunedleading_cnt"] = cnt.make<TH1D>("Hprunedleading_cnt", ";M(leading H jet) [GeV];;" ,100 ,100., 150.) ;

    h1_["Hpt2nd_cnt"]  = cnt.make<TH1D>("Hpt2nd_cnt", ";p_{T}(2nd H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Heta2nd_cnt"] = cnt.make<TH1D>("Heta2nd_cnt", ";#eta(2nd H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hpruned2nd_cnt"] = cnt.make<TH1D>("Hpruned2nd_cnt", ";M(2nd H jet) [GeV];;" ,100 ,100., 150.) ;

    h1_["Hpt_cnt"]  = cnt.make<TH1D>("Hpt_cnt", ";p_{T}(H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Heta_cnt"] = cnt.make<TH1D>("Heta_cnt", ";#eta(H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hpruned_cnt"] = cnt.make<TH1D>("Hpruned_cnt", ";Pruned Mass(H jet) [GeV];;" ,100 ,100., 150.) ;

    h1_["Topptleading_cnt"]  = cnt.make<TH1D>("Topptleading_cnt", ";p_{T}(leading Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topetaleading_cnt"] = cnt.make<TH1D>("Topetaleading_cnt", ";#eta(leading Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdropleading_cnt"] = cnt.make<TH1D>("Topsoftdropleading_cnt", ";M(leading Top jet) [GeV];;" ,200 ,80., 280.) ;

    h1_["Toppt2nd_cnt"]  = cnt.make<TH1D>("Toppt2nd_cnt", ";p_{T}(2nd Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topeta2nd_cnt"] = cnt.make<TH1D>("Topeta2nd_cnt", ";#eta(2nd Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdrop2nd_cnt"] = cnt.make<TH1D>("Topsoftdrop2nd_cnt", ";M(2nd Top jet) [GeV];;" ,200 ,80., 280.) ;

    h1_["Toppt_cnt"]  = cnt.make<TH1D>("Toppt_cnt", ";p_{T}(Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topeta_cnt"] = cnt.make<TH1D>("Topeta_cnt", ";#eta(Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdrop_cnt"] = cnt.make<TH1D>("Topsoftdrop_cnt", ";SoftDrop Mass(Top jet) [GeV];;" ,200 ,80., 280.) ;


    //AK8 jets cat                                                                                                                                                     
    h1_["Wptleading_cat"]  = cnt.make<TH1D>("Wptleading_cat", ";p_{T}(leading W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Wetaleading_cat"] = cnt.make<TH1D>("Wetaleading_cat", ";#eta(leading W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wprunedleading_cat"] = cnt.make<TH1D>("Wprunedleading_cat", ";M(leading W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wpt2nd_cat"]  = cnt.make<TH1D>("Wpt2nd_cat", ";p_{T}(2nd W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta2nd_cat"] = cnt.make<TH1D>("Weta2nd_cat", ";#eta(2nd W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned2nd_cat"] = cnt.make<TH1D>("Wpruned2nd_cat", ";M(2nd W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wpt_cat"]  = cnt.make<TH1D>("Wpt_cat", ";p_{T}(W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta_cat"] = cnt.make<TH1D>("Weta_cat", ";#eta(W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned_cat"] = cnt.make<TH1D>("Wpruned_cat", ";Pruned Mass(W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wptleading1_cat"]  = cnt.make<TH1D>("Wptleading1_cat", ";p_{T}(leading W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Wetaleading1_cat"] = cnt.make<TH1D>("Wetaleading1_cat", ";#eta(leading W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wprunedleading1_cat"] = cnt.make<TH1D>("Wprunedleading1_cat", ";M(leading W jet) [GeV];;" ,200 ,50., 150.) ;

    h1_["Wpt2nd1_cat"]  = cnt.make<TH1D>("Wpt2nd1_cat", ";p_{T}(2nd W jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Weta2nd1_cat"] = cnt.make<TH1D>("Weta2nd1_cat", ";#eta(2nd W jet);;" , 80 ,-4. ,4.) ;
    h1_["Wpruned2nd1_cat"] = cnt.make<TH1D>("Wpruned2nd1_cat", ";M(2nd W jet) [GeV];;" ,200 ,50., 150.) ;


    h1_["Hptleading_cat"]  = cnt.make<TH1D>("Hptleading_cat", ";p_{T}(leading H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Hetaleading_cat"] = cnt.make<TH1D>("Hetaleading_cat", ";#eta(leading H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hprunedleading_cat"] = cnt.make<TH1D>("Hprunedleading_cat", ";M(leading H jet) [GeV];;" ,100 ,100., 150.) ;

    h1_["Hpt2nd_cat"]  = cnt.make<TH1D>("Hpt2nd_cat", ";p_{T}(2nd H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Heta2nd_cat"] = cnt.make<TH1D>("Heta2nd_cat", ";#eta(2nd H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hpruned2nd_cat"] = cnt.make<TH1D>("Hpruned2nd_cat", ";M(2nd H jet) [GeV];;" ,100 ,100., 150.) ;

    h1_["Hpt_cat"]  = cnt.make<TH1D>("Hpt_cat", ";p_{T}(H jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Heta_cat"] = cnt.make<TH1D>("Heta_cat", ";#eta(H jet);;" , 80 ,-4. ,4.) ;
    h1_["Hpruned_cat"] = cnt.make<TH1D>("Hpruned_cat", ";Pruned Mass(H jet) [GeV];;" ,100 ,100., 150.) ;


    h1_["Topptleading_cat"]  = cnt.make<TH1D>("Topptleading_cat", ";p_{T}(leading Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topetaleading_cat"] = cnt.make<TH1D>("Topetaleading_cat", ";#eta(leading Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdropleading_cat"] = cnt.make<TH1D>("Topsoftdropleading_cat", ";M(leading Top jet) [GeV];;" ,200 ,80., 280.) ;

    h1_["Toppt2nd_cat"]  = cnt.make<TH1D>("Toppt2nd_cat", ";p_{T}(2nd Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topeta2nd_cat"] = cnt.make<TH1D>("Topeta2nd_cat", ";#eta(2nd Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdrop2nd_cat"] = cnt.make<TH1D>("Topsoftdrop2nd_cat", ";M(2nd Top jet) [GeV];;" ,200 ,80., 280.) ;

    h1_["Toppt_cat"]  = cnt.make<TH1D>("Toppt_cat", ";p_{T}(Top jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["Topeta_cat"] = cnt.make<TH1D>("Topeta_cat", ";#eta(Top jet);;" , 80 ,-4. ,4.) ;
    h1_["Topsoftdrop_cat"] = cnt.make<TH1D>("Topsoftdrop_cat", ";SoftDrop Mass(Top jet) [GeV];;" ,200 ,80., 280.) ;


    // h1_["ht1_cat"]   =  cnt.make<TH1D>("ht1_cat", ";H_{T} (AK4 jets) [GeV]", 200, 0., 1000.) ;
    // h1_["st1_cat"]   =  cnt.make<TH1D>("st1_cat",";S_{T} [GeV]", 200, 0., 1000.) ;
    h1_["nak41_cat"] = cnt.make<TH1D>("nak41_cat", ";N(AK4 jets);;" , 11, -0.5,10.5) ;

    h1_["ht2_cat"]   =  cnt.make<TH1D>("ht2_cnt", ";H_{T} (AK4 jets) [GeV]", 200, 0., 1000.) ;
    h1_["st2_cat"]   =  cnt.make<TH1D>("st2_cnt",";S_{T} [GeV]", 200, 0., 1000.) ;
    h1_["nak42_cat"] = cnt.make<TH1D>("nak42_cat", ";N(AK4 jets);;" , 11, -0.5,10.5) ;

    h1_["nak8_cat"] = cnt.make<TH1D>("nak8_cat", ";N(AK8 jets);;" , 11, -0.5,10.5) ;
    h1_["nwjet_cat"] = cnt.make<TH1D>("nwjet_cat", ";N(W jets );;" , 6, -0.5,5.5) ;
    h1_["nhjet_cat"] = cnt.make<TH1D>("nhjet_cat", ";N(H jets );;" , 6, -0.5,5.5) ;
    h1_["ntjet_cat"] = cnt.make<TH1D>("ntjet_cat", ";N(top jets);;" , 6, -0.5,5.5) ;

    h1_["nak8_pre"] = cnt.make<TH1D>("nak8_pre", ";N(AK8 jets);;" , 11, -0.5,10.5) ;
    h1_["nwjet_pre"] = cnt.make<TH1D>("nwjet_pre", ";N(W jets );;" , 6, -0.5,5.5) ;
    h1_["nhjet_pre"] = cnt.make<TH1D>("nhjet_pre", ";N(H jets );;" , 6, -0.5,5.5) ;
    h1_["ntjet_pre"] = cnt.make<TH1D>("ntjet_pre", ";N(top jets);;" , 6, -0.5,5.5) ;

    h1_["nak8_cnt"] = cnt.make<TH1D>("nak8_cnt", ";N(AK8 jets);;" , 11, -0.5,10.5) ;
    h1_["nwjet_cnt"] = cnt.make<TH1D>("nwjet_cnt", ";N(W jets );;" , 6, -0.5,5.5) ;
    h1_["nhjet_cnt"] = cnt.make<TH1D>("nhjet_cnt", ";N(H jets );;" , 6, -0.5,5.5) ;
    h1_["ntjet_cnt"] = cnt.make<TH1D>("ntjet_cnt", ";N(top jets);;" , 6, -0.5,5.5) ;





    //0btag nak4
    h1_["nak4_0b1"] = cat1.make<TH1D>("nak4_0b1", ";N(AK4 jets);;" , 11, -0.5,10.5) ;
    h1_["nak4_0b2"] = cat1.make<TH1D>("nak4_0b2", ";N(AK4 jets);;" , 11, -0.5,10.5) ;
    h1_["nak4_0b3"] = cat1.make<TH1D>("nak4_0b3", ";N(AK4 jets);;" , 11, -0.5,10.5) ;
    h1_["nak4HT_0b3"] = pre.make<TH1D>("nak4HT_0b3", ";N(AK4 jets);;" , 11, -0.5,10.5) ;
    h1_["nak4HT_pre"] = pre.make<TH1D>("nak4HT_pre", ";N(AK4 jets);;" , 11, -0.5,10.5) ;
    h1_["nak4_0b3_l1000"] = cat1.make<TH1D>("nak4_0b3_l1000", ";N(AK4 jets);;" , 11, -0.5,10.5) ;
    h1_["nak4_0b3_g1000"] = cat1.make<TH1D>("nak4_0b3_g1000", ";N(AK4 jets);;" , 11, -0.5,10.5) ;
    h1_["ptak8leading"]  = sig.make<TH1D>("ptak8leading", ";p_{T}(leading AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak8leading"] = sig.make<TH1D>("etaak8leading", ";#eta(leading AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak8leading"] = sig.make<TH1D>("softdropmak8leading", ";M(leading AK8 jet) [GeV];;" ,100 ,0., 200.) ; 
    h1_["ptak82nd"]  = sig.make<TH1D>("ptak82nd", ";p_{T}(2nd AK8 jet) [GeV];;" , 50, 0., 1000.) ; 
    h1_["etaak82nd"] = sig.make<TH1D>("etaak82nd", ";#eta(2nd AK8 jet);;" , 80 ,-4. ,4.) ; 
    h1_["softdropmak82nd"] = sig.make<TH1D>("softdropmak82nd", ";M(2nd AK8 jet) [GeV];;" ,100 ,0., 200.) ;


  if(categorize_){

    //for Z and H categories seperetely along with b
    h1_["cutflow1"] = cat.make<TH1D>("cutflow1", "cut flow", 15, 0.5, 15.5) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(1, "no t,Z,b ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(4, "t1Z1 ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(5, "t0Z1") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(6, "t1Z0 ") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(7, "t0Z0") ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(8, "t1Z1H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(9, "t1Z1H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(10, "t0Z1H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(11, "t0Z1H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(13, "t1Z0H0 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(14, "t0Z0H1 " ) ;
    h1_["cutflow1"] -> GetXaxis() -> SetBinLabel(15, "t0Z0H0 " ) ;
      
    h1_["cutflow2"] = cat.make<TH1D>("cutflow2", "cut flow", 16, 0.5, 16.5) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(1, "t1Z1H1b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(2, "t1Z1H1b2") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(3, "t1Z1H0b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(4, "t1Z1H0b2  ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(5, "t0Z1H1b1 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(6, "t0Z1H1b2 ") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(7, "t0Z1H0b1") ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(8, "t0Z1H0b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(9, "t1Z0H1b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(10, "t1Z0H1b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(11, "t1Z0H0b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H0b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(13, "t0Z0H1b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(14, "t0Z0H1b2 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(15, "t0Z0H0b1 " ) ;
    h1_["cutflow2"] -> GetXaxis() -> SetBinLabel(16, "t0Z0H0b2 " ) ;
      
    //for ZH(combined category)
    h1_["cutflow3"] = cat.make<TH1D>("cutflow3", "cut flow", 15, 0.5, 15.5) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(1, "no t,Z,b ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(4, "t1ZH1 ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(5, "t0ZH1") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(6, "t1ZH0 ") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(7, "t0ZH0") ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(8, "t1ZH1b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(9, "t1ZH1b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(10, "t0ZH1b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(11, "t0ZH0b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(12, "t1ZH0b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(13, "t1ZH0b2 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(14, "t0ZH0b1 " ) ;
    h1_["cutflow3"] -> GetXaxis() -> SetBinLabel(15, "t0ZH0b2 " ) ;
      
      
    h1_["cutflow4"] = cat.make<TH1D>("cutflow4", "cut flow", 16, 0.5, 16.5) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(1, "t1Z1H1b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(2, "t1Z1H1b2") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(3, "t1Z1H0b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(4, "t1Z1H0b2  ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(5, "t0Z1H1b1 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(6, "t0Z1H1b2 ") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(7, "t0Z1H0b1") ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(8, "t0Z1H0b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(9, "t1Z0H1b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(10, "t1Z0H1b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(11, "t1Z0H0b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H0b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(13, "t0Z0H1b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(14, "t0Z0H1b2 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(15, "t0Z0H0b1 " ) ;
    h1_["cutflow4"] -> GetXaxis() -> SetBinLabel(16, "t0Z0H0b2 " ) ;
      
    
    //for ZH(combined category)
    h1_["cutflow5"] = cat.make<TH1D>("cutflow5", "cut flow", 15, 0.5, 15.5) ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(1, "no t,Z,b ") ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(2, "b1  ") ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(3, "b2") ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(4, "t1ZH1 ") ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(5, "t0ZH1") ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(6, "t1ZH0 ") ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(7, "t0ZH0") ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(8, "t1ZH1b1 " ) ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(9, "t1ZH1b2 " ) ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(10, "t0ZH1b1 " ) ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(11, "t0ZH0b2 " ) ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(12, "t1ZH0b1 " ) ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(13, "t1ZH0b2 " ) ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(14, "t0ZH0b1 " ) ;
    h1_["cutflow5"] -> GetXaxis() -> SetBinLabel(15, "t0ZH0b2 " ) ;
      
      
    h1_["cutflow6"] = cat.make<TH1D>("cutflow6", "cut flow", 16, 0.5, 16.5) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(1, "t1Z1H1b1 ") ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(2, "t1Z1H1b2") ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(3, "t1Z1H0b1 ") ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(4, "t1Z1H0b2  ") ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(5, "t0Z1H1b1 ") ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(6, "t0Z1H1b2 ") ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(7, "t0Z1H0b1") ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(8, "t0Z1H0b2 " ) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(9, "t1Z0H1b1 " ) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(10, "t1Z0H1b2 " ) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(11, "t1Z0H0b1 " ) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H0b2 " ) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(13, "t0Z0H1b1 " ) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(14, "t0Z0H1b2 " ) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(15, "t0Z0H0b1 " ) ;
    h1_["cutflow6"] -> GetXaxis() -> SetBinLabel(16, "t0Z0H0b2 " ) ;


    h1_["cutflow10"] = cat.make<TH1D>("cutflow10", "cut flow", 12, 0.5, 12.5) ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(1, "D1ZB1Hb1b1") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(2, "D1ZB1H1b1") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(3, "D1Z1Hb1b1 ") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(4, "D1Z1H1b1") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(5, "BC1ZB1Hb1b1") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(6, "BC1ZB1H1b1") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(7, "BC1Z1Hb1b1 ") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(8, "BC1Z1H1b1") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(9, "t1ZB1Hb1b1") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(10, "t1ZB1H1b1") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(11, "t1Z1Hb1b1 ") ;
    h1_["cutflow10"] -> GetXaxis() -> SetBinLabel(12, "t1Z1H1b1") ;
    

    h1_["cutflow11"] = cat.make<TH1D>("cutflow11", "cut flow", 12, 0.5, 12.5) ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(1, "D1ZB1Hb1b2") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(2, "D1ZB1H1b2") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(3, "D1Z1Hb1b2 ") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(4, "D1Z1H1b2") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(5, "BC1ZB1Hb1b2") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(6, "BC1ZB1H1b2") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(7, "BC1Z1Hb1b2 ") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(8, "BC1Z1H1b2") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(9, "t1ZB1Hb1b2") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(10, "t1ZB1H1b2") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(11, "t1Z1Hb1b2 ") ;
    h1_["cutflow11"] -> GetXaxis() -> SetBinLabel(12, "t1Z1H1b2") ;


    h1_["cutflow12"] = cat.make<TH1D>("cutflow12", "cut flow", 12, 0.5, 12.5) ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(1, "D1ZB1Hb0b1") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(2, "D1ZB1H0b1") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(3, "D1Z1Hb0b1 ") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(4, "D1Z1H0b1") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(5, "BC1ZB1Hb0b1") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(6, "BC1ZB1H0b1") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(7, "BC1Z1Hb0b1 ") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(8, "BC1Z1H0b1") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(9, "t1ZB1Hb0b1") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(10, "t1ZB1H0b1") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(11, "t1Z1Hb0b1 ") ;
    h1_["cutflow12"] -> GetXaxis() -> SetBinLabel(12, "t1Z1H0b1") ;


    h1_["cutflow13"] = cat.make<TH1D>("cutflow13", "cut flow", 12, 0.5, 12.5) ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(1, "D1ZB1Hb0b2") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(2, "D1ZB1H0b2") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(3, "D1Z1Hb0b2 ") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(4, "D1Z1H0b2") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(5, "BC1ZB1Hb0b2") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(6, "BC1ZB1H0b2") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(7, "BC1Z1Hb0b2 ") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(8, "BC1Z1H0b2") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(9, "t1ZB1Hb0b2") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(10, "t1ZB1H0b2") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(11, "t1Z1Hb0b2 ") ;
    h1_["cutflow13"] -> GetXaxis() -> SetBinLabel(12, "t1Z1H0b2") ;

    h1_["cutflow14"] = cat.make<TH1D>("cutflow14", "cut flow", 12, 0.5, 12.5) ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(1, "D1ZB0Hb0b2") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(2, "D1ZB0H0b2") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(3, "D1Z0Hb0b2 ") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(4, "D1Z0H0b2") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(5, "BC1ZB0Hb0b2") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(6, "BC1ZB0H0b2") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(7, "BC1Z0Hb0b2 ") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(8, "BC1Z0H0b2") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(9, "t1ZB0Hb0b2") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(10, "t1ZB0H0b2") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(11, "t1Z0Hb0b2 ") ;
    h1_["cutflow14"] -> GetXaxis() -> SetBinLabel(12, "t1Z0H0b2") ;




    //HPrime Candidates
    h1_["HPrime_mass_b_cnt"]  = cat1.make<TH1D>("HPrimemass-boosted-cnt", ";M(HPrime-boosted) [GeV];;", 100, 0., 400);
    h1_["HPrime_Pt_b_cnt"]  = cat1.make<TH1D>("HPrimePt-boosted-cnt", ";Pt(HPrime-boosted) [GeV];;", 100, 0., 1200);
    h1_["nHPrimecandidatejets_b_cnt"] = cat1.make<TH1D>("nHPrimecandidate-boosted-cnt", ";N(HPrime jets-boosted);;" , 21, -0.5, 20.5);


    h1_["HPrime_mass_nb_cnt"]  = cat1.make<TH1D>("HPrimemassnb-cnt", ";M(HPrime) [GeV];;", 100, 0., 400);
    h1_["HPrime_Pt_nb_cnt"]  = cat1.make<TH1D>("HPrimePtnb-cnt", ";Pt(HPrime) [GeV];;", 100, 0., 1200);
    h1_["nHPrimecandidatejets_nb_cnt"] = cat1.make<TH1D>("nHPrimecandidatesnb-cnt", ";N(HPrime jets);;" , 21, -0.5, 20.5);

    h1_["nHPrimecandidatejets_cnt"] = cat1.make<TH1D>("nHPrimecandidates-tot-cnt", ";N(HPrime jets);;" , 21, -0.5, 20.5);
    h1_["nHPrimecandidatejets1_cnt"] = cat1.make<TH1D>("nHPrimecandidates1-tot-cnt", ";N(HPrime jets);;" , 21, -0.5, 20.5);




    //H candidates                                                                                                           
    h1_["H_mass_b_cnt"]  = cat.make<TH1D>("Hmass-boosted-cnt", ";M(H-boosted) [GeV];;", 100, 0., 400);
    h1_["H_Pt_b_cnt"]  = cat.make<TH1D>("HPt-boosted-cnt", ";Pt(H-boosted) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_b_cnt"] = cat.make<TH1D>("nHcandidate-boosted-cnt", ";N(H jets-boosted);;" , 21, -0.5, 20.5);
      
    h1_["H_mass_nb_cnt"]  = cat.make<TH1D>("Hmassnb-cnt", ";M(H) [GeV];;", 100, 0., 400);
    h1_["H_Pt_nb_cnt"]  = cat.make<TH1D>("HPtnb-cnt", ";Pt(H) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_nb_cnt"] = cat.make<TH1D>("nHcandidatesnb-cnt", ";N(H jets);;" , 21, -0.5, 20.5);
      
    h1_["nHcandidatejets_cnt"] = cat.make<TH1D>("nHcandidates-tot-cnt", ";N(H jets);;" , 21, -0.5, 20.5);
    h1_["nHcandidatejets1_cnt"] = cat.make<TH1D>("nHcandidates1-tot-cnt", ";N(H jets);;" , 21, -0.5, 20.5);
      
    // Z candidates
    h1_["Z_mass_a_cnt"]  = cat.make<TH1D>("Zmass-boosted-cnt", ";M(Z-boosted) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_a_cnt"]  = cat.make<TH1D>("ZPt-boosted-cnt", ";Pt(Z-boosted) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_a_cnt"] = cat.make<TH1D>("nzcandidate-boosted-cnt", ";N(Z jets-boosted);;" , 21, -0.5, 20.5);
      
    h1_["Z_mass_b_cnt"]  = cat.make<TH1D>("Zmass-cnt", ";M(Z) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_b_cnt"]  = cat.make<TH1D>("ZPt-cnt", ";Pt(Z) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_b_cnt"] = cat.make<TH1D>("nzcandidates-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);
      
    h1_["nzcandidatejets_tot_cnt"] = cat.make<TH1D>("nzcandidates-tot-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);
    h1_["nzcandidatejets1_tot_cnt"] = cat.make<TH1D>("nzcandidates1-tot-cnt", ";N(Z jets);;" , 21, -0.5, 20.5);  
    // cat A
    h1_["top_mass_a_cnt"]  = cat.make<TH1D>("topmas-A-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_a_cnt"]  = cat.make<TH1D>("topPt-A-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_a_cnt"] = cat.make<TH1D>("ntopcandidate-A-cnt", ";N(top jets);;" , 21, -0.5, 20.5);

    h1_["top_mass_bc_cnt"]  = cat.make<TH1D>("topmass-Bc-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_bc_cnt"]  = cat.make<TH1D>("topPt-BC-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_bc_cnt"] = cat.make<TH1D>("ntopcandidate-BC-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
      
    // cat D
    h1_["top_mass_d_cnt"]  = cat.make<TH1D>("topmass-D-cnt", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_d_cnt"]  = cat.make<TH1D>("topPt-D-cnt", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_d_cnt"] = cat.make<TH1D>("ntopcandidate-D-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
      
    //W and light jet(BC)
    h1_["W_mass_bc_cnt"]  = cat.make<TH1D>("Wmass-BC-cnt", ";M( W boson) [GeV];;", 100, 0., 400);
    h1_["nWcandidatejets_bc_cnt"] = cat.make<TH1D>("nWcandidate-BC-cnt", ";N(W candidate jets);;" , 21, -0.5, 20.5);
      
    h1_["lightjet_mass_bc_cnt"]  = cat.make<TH1D>("lightjetmass-BC-cnt", ";M( light jet) [GeV];;", 100, 0., 400);
    h1_["nlightjetcandidatejets_bc_cnt"] = cat.make<TH1D>("nlightjetcandidate-cnt", ";N(lightjet candidate jets);;" , 21, -0.5, 20.5);
    //total top ( A+ BC+D)
    h1_["ntopcandidatejets_cnt"] = cat.make<TH1D>("ntopcandidate-tot-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
    h1_["ntopcandidatejets1_cnt"] = cat.make<TH1D>("ntopcandidate1-tot-cnt", ";N(top jets);;" , 21, -0.5, 20.5);
      
    //signal region
      
    //H candidates                                                                                                                                                              
    h1_["H_mass_b_sig"]  = cat.make<TH1D>("Hmass-boosted-sig", ";M(H-boosted) [GeV];;", 100, 0., 400);
    h1_["H_Pt_b_sig"]  = cat.make<TH1D>("HPt-boosted-sig", ";Pt(H-boosted) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_b_sig"] = cat.make<TH1D>("nHcandidate-boosted-sig", ";N(H jets-boosted);;" , 21, -0.5, 20.5);
      
    h1_["H_mass_nb_sig"]  = cat.make<TH1D>("Hmassnb-sig", ";M(H) [GeV];;", 100, 0., 400);
    h1_["H_Pt_nb_sig"]  = cat.make<TH1D>("HPtnb-sig", ";Pt(H) [GeV];;", 100, 0., 1200);
    h1_["nHcandidatejets_nb_sig"] = cat.make<TH1D>("nHcandidatesnb-sig", ";N(H jets);;" , 21, -0.5, 20.5);
      
    h1_["nHcandidatejets_sig"] = cat.make<TH1D>("nHcandidates-tot-sig", ";N(H jets);;" , 21, -0.5, 20.5);
    h1_["nHcandidatejets1_sig"] = cat.make<TH1D>("nHcandidates1-tot-sig", ";N(H jets);;" , 21, -0.5, 20.5);
      
   
    // Z candidates                                                                                                                                                             
    h1_["Z_mass_a_sig"]  = cat.make<TH1D>("Zmass-boosted-sig", ";M(Z-boosted) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_a_sig"]  = cat.make<TH1D>("ZPt-boosted-sig", ";Pt(Z-boosted) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_a_sig"] = cat.make<TH1D>("nzcandidate-boosted-sig", ";N(Z jets-boosted);;" , 21, -0.5, 20.5);
      
    h1_["Z_mass_b_sig"]  = cat.make<TH1D>("Zmass-sig", ";M(Z) [GeV];;", 100, 0., 400);
    h1_["Z_Pt_b_sig"]  = cat.make<TH1D>("ZPt-sig", ";Pt(Z) [GeV];;", 100, 0., 1200);
    h1_["nzcandidatejets_b_sig"] = cat.make<TH1D>("nzcandidates-sig", ";N(Z jets);;" , 21, -0.5, 20.5);
      
    h1_["nzcandidatejets_tot_sig"] = cat.make<TH1D>("nzcandidates-tot-sig", ";N(Z jets);;" , 21, -0.5, 20.5);
    h1_["nzcandidatejets1_tot_sig"] = cat.make<TH1D>("nzcandidates1-tot-sig", ";N(Z jets);;" , 21, -0.5, 20.5);
      
    // cat A                                                                                                                                                                    
    h1_["top_mass_a_sig"]  = cat.make<TH1D>("topmas-A-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_a_sig"]  = cat.make<TH1D>("topPt-A-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_a_sig"] = cat.make<TH1D>("ntopcandidate-A-sig", ";N(top jets);;" , 21, -0.5, 20.5);
      
    h1_["top_mass_bc_sig"]  = cat.make<TH1D>("topmass-Bc-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_bc_sig"]  = cat.make<TH1D>("topPt-BC-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_bc_sig"] = cat.make<TH1D>("ntopcandidate-BC-sig", ";N(top jets);;" , 21, -0.5, 20.5);
      
    // cat D                                                                                                                                                                    
    h1_["top_mass_d_sig"]  = cat.make<TH1D>("topmass-D-sig", ";M( t quark) [GeV];;", 100, 0., 400);
    h1_["top_Pt_d_sig"]  = cat.make<TH1D>("topPt-D-sig", ";Pt( t quark) [GeV];;", 100, 0., 1200);
    h1_["ntopcandidatejets_d_sig"] = cat.make<TH1D>("ntopcandidate-D-sig", ";N(top jets);;" , 21, -0.5, 20.5);
      
    //W and light jet(BC)                                                                                                                                                       
    h1_["W_mass_bc_sig"]  = cat.make<TH1D>("Wmass-BC-sig", ";M( W boson) [GeV];;", 100, 0., 400);
    h1_["nWcandidatejets_bc_sig"] = cat.make<TH1D>("nWcandidate-BC-sig", ";N(W candidate jets);;" , 21, -0.5, 20.5);
      
    h1_["lightjet_mass_bc_sig"]  = cat.make<TH1D>("lightjetmass-BC-sig", ";M( light jet) [GeV];;", 100, 0., 400);
    h1_["nlightjetcandidatejets_bc_sig"] = cat.make<TH1D>("nlightjetcandidate-sig", ";N(lightjet candidate jets);;" , 21, -0.5, 20.5);
    //total top ( A+ BC+D)                                                                                                                                                       
    h1_["ntopcandidatejets_sig"] = cat.make<TH1D>("ntopcandidate-tot-sig", ";N(top jets);;" , 21, -0.5, 20.5);
    h1_["ntopcandidatejets1_sig"] = cat.make<TH1D>("ntopcandidate1-tot-sig", ";N(top jets);;" , 21, -0.5, 20.5);
   
    //ST plots
      
    //top,Z, Higgs
      
    h1_["st_sig"] =cat.make<TH1D>("ST_sig", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_prime"] =cat.make<TH1D>("st_prime", ";S_{T} [Gev];;" , 50, 0.,2500.);
    h1_["st_noprime"] =cat.make<TH1D>("st_noprime", ";S_{T} [Gev];;" , 50, 0.,2500.);  

    h1_["st_prime1"] =cat.make<TH1D>("st_prime1", ";S_{T} [Gev];;" , 50, 0.,5000.);
    h1_["st_noprime1"] =cat.make<TH1D>("st_noprime1", ";S_{T} [Gev];;" , 50, 0.,5000.);




    h1_["st_sig1b"] =cat.make<TH1D>("ST_sig1b", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sig2b"] =cat.make<TH1D>("ST_sig2b", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sigT1Z1"] =cat.make<TH1D>("ST_sigT1Z1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1"] =cat.make<TH1D>("ST_sigT0Z1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0"] =cat.make<TH1D>("ST_sigT1Z0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0"] =cat.make<TH1D>("ST_sigT0Z0", ";S_{T} [Gev];;" , 50, 1000.,2500.);


    h1_["st_sigT1Z1H1"] =cat.make<TH1D>("ST_sigT1Z1H1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0"] =cat.make<TH1D>("ST_sigT1Z1H0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H1"] =cat.make<TH1D>("ST_sigT0Z1H1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0"] =cat.make<TH1D>("ST_sigT0Z1H0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H1"] =cat.make<TH1D>("ST_sigT1Z0H1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0"] =cat.make<TH1D>("ST_sigT1Z0H0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H1"] =cat.make<TH1D>("ST_sigT0Z0H1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0"] =cat.make<TH1D>("ST_sigT0Z0H0", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    
    
    h1_["st_sigT1Z1H1b1"] =cat.make<TH1D>("ST_sigT1Z1H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H1b2"] =cat.make<TH1D>("ST_sigT1Z1H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0b1"] =cat.make<TH1D>("ST_sigT1Z1H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0b2"] =cat.make<TH1D>("ST_sigT1Z1H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    
    h1_["st_sigT0Z1H1b1"] =cat.make<TH1D>("ST_sigT0Z1H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H1b2"] =cat.make<TH1D>("ST_sigT0Z1H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0b1"] =cat.make<TH1D>("ST_sigT0Z1H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0b2"] =cat.make<TH1D>("ST_sigT0Z1H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    
    h1_["st_sigT1Z0H1b1"] =cat.make<TH1D>("ST_sigT1Z0H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H1b2"] =cat.make<TH1D>("ST_sigT1Z0H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0b1"] =cat.make<TH1D>("ST_sigT1Z0H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0b2"] =cat.make<TH1D>("ST_sigT1Z0H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
      
    h1_["st_sigT0Z0H1b1"] =cat.make<TH1D>("ST_sigT0Z0H1b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H1b2"] =cat.make<TH1D>("ST_sigT0Z0H1b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0b1"] =cat.make<TH1D>("ST_sigT0Z0H0b1", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0b2"] =cat.make<TH1D>("ST_sigT0Z0H0b2", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    
   
    /*

    h1_["st_sigT1Z1H1b1_ak8"] =cat.make<TH1D>("ST_sigT1Z1H1b1_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H1b2_ak8"] =cat.make<TH1D>("ST_sigT1Z1H1b2_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0b1_ak8"] =cat.make<TH1D>("ST_sigT1Z1H0b1_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z1H0b2_ak8"] =cat.make<TH1D>("ST_sigT1Z1H0b2_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);

    h1_["st_sigT0Z1H1b1_ak8"] =cat.make<TH1D>("ST_sigT0Z1H1b1_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H1b2_ak8"] =cat.make<TH1D>("ST_sigT0Z1H1b2_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0b1_ak8"] =cat.make<TH1D>("ST_sigT0Z1H0b1_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z1H0b2_ak8"] =cat.make<TH1D>("ST_sigT0Z1H0b2_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);

    h1_["st_sigT1Z0H1b1_ak8"] =cat.make<TH1D>("ST_sigT1Z0H1b1_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H1b2_ak8"] =cat.make<TH1D>("ST_sigT1Z0H1b2_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0b1_ak8"] =cat.make<TH1D>("ST_sigT1Z0H0b1_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT1Z0H0b2_ak8"] =cat.make<TH1D>("ST_sigT1Z0H0b2_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);

    h1_["st_sigT0Z0H1b1_ak8"] =cat.make<TH1D>("ST_sigT0Z0H1b1_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H1b2_ak8"] =cat.make<TH1D>("ST_sigT0Z0H1b2_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0b1_ak8"] =cat.make<TH1D>("ST_sigT0Z0H0b1_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);
    h1_["st_sigT0Z0H0b2_ak8"] =cat.make<TH1D>("ST_sigT0Z0H0b2_ak8", ";S_{T} [Gev];;" , 50, 1000.,2500.);


    */







    h1_["st_cntT1Z1H1b1"] =cat.make<TH1D>("ST_cntT1Z1H1b1", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z1H1b2"] =cat.make<TH1D>("ST_cntT1Z1H1b2", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z1H0b1"] =cat.make<TH1D>("ST_cntT1Z1H0b1", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z1H0b2"] =cat.make<TH1D>("ST_cntT1Z1H0b2", ";S_{T} [Gev];;" ,100,0.,1000.);

  
    h1_["nbjets_cntT1Z1H1b1"] = cat.make<TH1D>("nbjets_cntT1Z1H1b1", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_cntT1Z1H1b2"] = cat.make<TH1D>("nbjets_cntT1Z1H1b2", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_cntT1Z1H0b1"] = cat.make<TH1D>("nbjets_cntT1Z1H0b1", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_cntT1Z1H0b2"] = cat.make<TH1D>("nbjets_cntT1Z1H0b2", ";N(b jets);;" , 11, -0.5,10.5) ;

    h1_["st_cntT1Z1H1"] =cat.make<TH1D>("ST_cntT1Z1H1", ";S_{T} [Gev];;" ,100,0.,4000.);
    h1_["ht_cntT1Z1H1"] =cat.make<TH1D>("HT_cntT1Z1H1", ";H_{T} [Gev];;" ,100,0.,4000.);
    h1_["st_cntT1Z1H0"] =cat.make<TH1D>("ST_cntT1Z1H0", ";S_{T} [Gev];;" ,100,0.,4000.);
    h1_["ht_cntT1Z1H0"] =cat.make<TH1D>("HT_cntT1Z1H0", ";H_{T} [Gev];;" ,100,0.,4000.);


    h1_["st_cntT0Z1H1"] =cat.make<TH1D>("ST_cntT0Z1H1", ";S_{T} [Gev];;" ,100,0.,4000.);
    h1_["ht_cntT0Z1H1"] =cat.make<TH1D>("HT_cntT0Z1H1", ";H_{T} [Gev];;" ,100,0.,4000.);
    h1_["st_cntT0Z1H0"] =cat.make<TH1D>("ST_cntT0Z1H0", ";S_{T} [Gev];;" ,100,0.,4000.);
    h1_["ht_cntT0Z1H0"] =cat.make<TH1D>("HT_cntT0Z1H0", ";H_{T} [Gev];;" ,100,0.,4000.);

    h1_["st_cntT1Z0H1"] =cat.make<TH1D>("ST_cntT1Z0H1", ";S_{T} [Gev];;" ,100,0.,4000.);
    h1_["ht_cntT1Z0H1"] =cat.make<TH1D>("HT_cntT1Z0H1", ";H_{T} [Gev];;" ,100,0.,4000.);
    h1_["st_cntT1Z0H0"] =cat.make<TH1D>("ST_cntT1Z0H0", ";S_{T} [Gev];;" ,100,0.,4000.);
    h1_["ht_cntT1Z0H0"] =cat.make<TH1D>("HT_cntT1Z0H0", ";H_{T} [Gev];;" ,100,0.,4000.);



    h1_["st_cntT0Z0H1"] =cat.make<TH1D>("ST_cntT0Z0H1", ";S_{T} [Gev];;" ,100,0.,4000.);
    h1_["ht_cntT0Z0H1"] =cat.make<TH1D>("HT_cntT0Z0H1", ";H_{T} [Gev];;" ,100,0.,4000.);
    h1_["st_cntT0Z0H0"] =cat.make<TH1D>("ST_cntT0Z0H0", ";S_{T} [Gev];;" ,100,0.,4000.);
    h1_["ht_cntT0Z0H0"] =cat.make<TH1D>("HT_cntT0Z0H0", ";H_{T} [Gev];;" ,100,0.,4000.);




    h1_["st_cntT1Z1H1b1_A"] =cat.make<TH1D>("ST_cntT1Z1H1b1_A", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z1H1b2_A"] =cat.make<TH1D>("ST_cntT1Z1H1b2_A", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z1H0b1_A"] =cat.make<TH1D>("ST_cntT1Z1H0b1_A", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z1H0b2_A"] =cat.make<TH1D>("ST_cntT1Z1H0b2_A", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT0Z1H1b2_A"] =cat.make<TH1D>("ST_cntT0Z1H1b2_A", ";H_{T} [Gev];;" ,100,0.,1000.);

    h1_["ht_cntT1Z1H1b1"] =cat.make<TH1D>("HT_cntT1Z1H1b1", ";H_{T} [Gev];;" ,100,0.,1000.);
    h1_["ht_cntT1Z1H1b2"] =cat.make<TH1D>("HT_cntT1Z1H1b2", ";H_{T} [Gev];;" ,100,0.,1000.);
    h1_["ht_cntT1Z1H0b1"] =cat.make<TH1D>("HT_cntT1Z1H0b1", ";H_{T} [Gev];;" ,100,0.,1000.);
    h1_["ht_cntT1Z1H0b2"] =cat.make<TH1D>("HT_cntT1Z1H0b2", ";H_{T} [Gev];;" ,100,0.,1000.);
    h1_["ht_cntT0Z1H1b2"] =cat.make<TH1D>("HT_cntT0Z1H1b2", ";H_{T} [Gev];;" ,100,0.,1000.);

    h1_["ht_cntT1Z1H1b1_A"] =cat.make<TH1D>("HT_cntT1Z1H1b1_A", ";H_{T} [Gev];;" ,100,0.,1000.);
    h1_["ht_cntT1Z1H1b2_A"] =cat.make<TH1D>("HT_cntT1Z1H1b2_A", ";H_{T} [Gev];;" ,100,0.,1000.);
    h1_["ht_cntT1Z1H0b1_A"] =cat.make<TH1D>("HT_cntT1Z1H0b1_A", ";H_{T} [Gev];;" ,100,0.,1000.);
    h1_["ht_cntT1Z1H0b2_A"] =cat.make<TH1D>("HT_cntT1Z1H0b2_A", ";H_{T} [Gev];;" ,100,0.,1000.);
    h1_["ht_cntT0Z1H1b2_A"] =cat.make<TH1D>("HT_cntT0Z1H1b2_A", ";H_{T} [Gev];;" ,100,0.,1000.);




    h1_["st_cntT1Z1Hprime1b0"] =cat1.make<TH1D>("ST_cntT1Z1Hprime1b0", ";S_{T} [Gev];;" ,100,0.,1000.);

    h1_["st_cntT0Z1H1b1"] =cat.make<TH1D>("ST_cntT0Z1H1b1", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT0Z1H1b2"] =cat.make<TH1D>("ST_cntT0Z1H1b2", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT0Z1H0b1"] =cat.make<TH1D>("ST_cntT0Z1H0b1", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT0Z1H0b2"] =cat.make<TH1D>("ST_cntT0Z1H0b2", ";S_{T} [Gev];;" ,100,0.,1000.);

    h1_["st_cntT1Z0H1b1"] =cat.make<TH1D>("ST_cntT1Z0H1b1", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z0H1b2"] =cat.make<TH1D>("ST_cntT1Z0H1b2", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z0H0b1"] =cat.make<TH1D>("ST_cntT1Z0H0b1", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT1Z0H0b2"] =cat.make<TH1D>("ST_cntT1Z0H0b2", ";S_{T} [Gev];;" ,100,0.,1000.);

    h1_["st_cntT0Z0H1b1"] =cat.make<TH1D>("ST_cntT0Z0H1b1", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT0Z0H1b2"] =cat.make<TH1D>("ST_cntT0Z0H1b2", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT0Z0H0b1"] =cat.make<TH1D>("ST_cntT0Z0H0b1", ";S_{T} [Gev];;" ,100,0.,1000.);
    h1_["st_cntT0Z0H0b2"] =cat.make<TH1D>("ST_cntT0Z0H0b2", ";S_{T} [Gev];;" ,100,0.,1000.); 


    //
    h1_["st_cntD1ZB1Hb1b1"] =cat.make<TH1D>("ST_cntD1ZB1Hb1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1ZB1H1b1"] =cat.make<TH1D>("ST_cntD1ZB1H1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z1Hb1b1"] =cat.make<TH1D>("ST_cntD1Z1Hb1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z1H1b1"] =cat.make<TH1D>("ST_cntD1Z1H1b1", ";S_{T} [Gev];;" ,50,0.,1000.);


    h1_["st_cntBC1ZB1Hb1b1"] =cat.make<TH1D>("ST_cntBC1ZB1Hb1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1ZB1H1b1"] =cat.make<TH1D>("ST_cntBC1ZB1H1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z1Hb1b1"] =cat.make<TH1D>("ST_cntBC1Z1Hb1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z1H1b1"] =cat.make<TH1D>("ST_cntBC1Z1H1b1", ";S_{T} [Gev];;" ,50,0.,1000.);

    h1_["st_cntt1ZB1Hb1b1"] =cat.make<TH1D>("ST_cntt1ZB1Hb1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1ZB1H1b1"] =cat.make<TH1D>("ST_cntt1ZB1H1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z1Hb1b1"] =cat.make<TH1D>("ST_cntt1Z1Hb1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z1H1b1"] =cat.make<TH1D>("ST_cntt1Z1H1b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    //

    h1_["st_cntD1ZB1Hb1b2"] =cat.make<TH1D>("ST_cntD1ZB1Hb1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1ZB1H1b2"] =cat.make<TH1D>("ST_cntD1ZB1H1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z1Hb1b2"] =cat.make<TH1D>("ST_cntD1Z1Hb1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z1H1b2"] =cat.make<TH1D>("ST_cntD1Z1H1b2", ";S_{T} [Gev];;" ,50,0.,1000.);


    h1_["st_cntBC1ZB1Hb1b2"] =cat.make<TH1D>("ST_cntBC1ZB1Hb1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1ZB1H1b2"] =cat.make<TH1D>("ST_cntBC1ZB1H1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z1Hb1b2"] =cat.make<TH1D>("ST_cntBC1Z1Hb1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z1H1b2"] =cat.make<TH1D>("ST_cntBC1Z1H1b2", ";S_{T} [Gev];;" ,50,0.,1000.);

    h1_["st_cntt1ZB1Hb1b2"] =cat.make<TH1D>("ST_cntt1ZB1Hb1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1ZB1H1b2"] =cat.make<TH1D>("ST_cntt1ZB1H1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z1Hb1b2"] =cat.make<TH1D>("ST_cntt1Z1Hb1b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z1H1b2"] =cat.make<TH1D>("ST_cntt1Z1H1b2", ";S_{T} [Gev];;" ,50,0.,1000.);

    //
    h1_["st_cntD1ZB1Hb0b1"] =cat.make<TH1D>("ST_cntD1ZB1Hb0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1ZB1H0b1"] =cat.make<TH1D>("ST_cntD1ZB1H0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z1Hb0b1"] =cat.make<TH1D>("ST_cntD1Z1Hb0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z1H0b1"] =cat.make<TH1D>("ST_cntD1Z1H0b1", ";S_{T} [Gev];;" ,50,0.,1000.);


    h1_["st_cntBC1ZB1Hb0b1"] =cat.make<TH1D>("ST_cntBC1ZB1Hb0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1ZB1H0b1"] =cat.make<TH1D>("ST_cntBC1ZB1H0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z1Hb0b1"] =cat.make<TH1D>("ST_cntBC1Z1Hb0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z1H0b1"] =cat.make<TH1D>("ST_cntBC1Z1H0b1", ";S_{T} [Gev];;" ,50,0.,1000.);

    h1_["st_cntt1ZB1Hb0b1"] =cat.make<TH1D>("ST_cntt1ZB1Hb0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1ZB1H0b1"] =cat.make<TH1D>("ST_cntt1ZB1H0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z1Hb0b1"] =cat.make<TH1D>("ST_cntt1Z1Hb0b1", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z1H0b1"] =cat.make<TH1D>("ST_cntt1Z1H0b1", ";S_{T} [Gev];;" ,50,0.,1000.);


    //

    h1_["st_cntD1ZB1Hb0b2"] =cat.make<TH1D>("ST_cntD1ZB1Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1ZB1H0b2"] =cat.make<TH1D>("ST_cntD1ZB1H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z1Hb0b2"] =cat.make<TH1D>("ST_cntD1Z1Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z1H0b2"] =cat.make<TH1D>("ST_cntD1Z1H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);


    h1_["st_cntBC1ZB1Hb0b2"] =cat.make<TH1D>("ST_cntBC1ZB1Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1ZB1H0b2"] =cat.make<TH1D>("ST_cntBC1ZB1H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z1Hb0b2"] =cat.make<TH1D>("ST_cntBC1Z1Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z1H0b2"] =cat.make<TH1D>("ST_cntBC1Z1H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);

    h1_["st_cntt1ZB1Hb0b2"] =cat.make<TH1D>("ST_cntt1ZB1Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1ZB1H0b2"] =cat.make<TH1D>("ST_cntt1ZB1H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z1Hb0b2"] =cat.make<TH1D>("ST_cntt1Z1Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z1H0b2"] =cat.make<TH1D>("ST_cntt1Z1H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);

    //

    h1_["st_cntD1ZB0Hb0b2"] =cat.make<TH1D>("ST_cntD1ZB0Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1ZB0H0b2"] =cat.make<TH1D>("ST_cntD1ZB0H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z0Hb0b2"] =cat.make<TH1D>("ST_cntD1Z0Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntD1Z0H0b2"] =cat.make<TH1D>("ST_cntD1Z0H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);


    h1_["st_cntBC1ZB0Hb0b2"] =cat.make<TH1D>("ST_cntBC1ZB0Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1ZB0H0b2"] =cat.make<TH1D>("ST_cntBC1ZB0H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z0Hb0b2"] =cat.make<TH1D>("ST_cntBC1Z0Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntBC1Z0H0b2"] =cat.make<TH1D>("ST_cntBC1Z0H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);

    h1_["st_cntt1ZB0Hb0b2"] =cat.make<TH1D>("ST_cntt1ZB0Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1ZB0H0b2"] =cat.make<TH1D>("ST_cntt1ZB0H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z0Hb0b2"] =cat.make<TH1D>("ST_cntt1Z0Hb0b2", ";S_{T} [Gev];;" ,50,0.,1000.);
    h1_["st_cntt1Z0H0b2"] =cat.make<TH1D>("ST_cntt1Z0H0b2", ";S_{T} [Gev];;" ,50,0.,1000.);


  }

  //additional plots
  //o batg region plots
  h1_["nob_ht"]= cnt.make<TH1D>("nob_ht", ";H_{T} [Gev];;", 100, 0., 3000.);
  h1_["nob_st"] = cnt.make<TH1D>("nob_st", ";S_{T} [Gev];;", 50, 0., 4000.) ;

  h1_["nob_1000_ht"]= cnt.make<TH1D>("nob_1000_ht", ";H_{T} [Gev];;", 100, 0., 3000.);
  h1_["nob_1000_st"] = cnt.make<TH1D>("nob_1000_st", ";S_{T} [Gev];;", 50, 0., 4000.) ;


  h1_["ptak8_nob1000"]  = cnt.make<TH1D>("ptak8_nob1000", ";p_{T}(AK8 jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptW_nob1000"]  = cnt.make<TH1D>("ptW_nob1000", ";p_{T}(W tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptH_nob1000"]  = cnt.make<TH1D>("ptH_nob1000", ";p_{T}(H tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptT_nob1000"]  = cnt.make<TH1D>("ptT_nob1000", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);


  h1_["ptak8_st1000"]  = cnt.make<TH1D>("ptak8_st1000", ";p_{T}(AK8 jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptW_st1000"]  = cnt.make<TH1D>("ptW_st1000", ";p_{T}(W tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptH_st1000"]  = cnt.make<TH1D>("ptH_st1000", ";p_{T}(H tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptT_st1000"]  = cnt.make<TH1D>("ptT_st1000", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);

  h1_["ptTHnowwt_st1000"]  = cnt.make<TH1D>("ptTHnowwt_st1000", ";p_{T}(T/H (wt=1.0) tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptTHwt_st1000"]  = cnt.make<TH1D>("ptTHwt_st1000", ";p_{T}(Top/H (wt=1.06) matched jet) [GeV];;" , 100, 0., 2000.);

  h1_["ptTmatched_st1000"]  = cnt.make<TH1D>("ptTmatched_st1000", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptTnonmatched_st1000"]  = cnt.make<TH1D>("ptTnonmatched_st1000", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptHmatched_st1000"]  = cnt.make<TH1D>("ptHmatched_st1000", ";p_{T}(H tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptHnonmatched_st1000"]  = cnt.make<TH1D>("ptHnonmatched_st1000", ";p_{T}(H tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptWmatched_st1000"]  = cnt.make<TH1D>("ptWmatched_st1000", ";p_{T}(W tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptWnonmatched_st1000"]  = cnt.make<TH1D>("ptWnonmatched_st1000", ";p_{T}(W tagged jet) [GeV];;" , 100, 0., 2000.);


  h1_["ptT_st1000_t1Z1H1b1"]  = cnt.make<TH1D>("ptT_st1000_t1Z1H1b1", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptT_st1000_t1Z1H1b2"]  = cnt.make<TH1D>("ptT_st1000_t1Z1H1b2", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptT_st1000_t1Z0H1b2"]  = cnt.make<TH1D>("ptT_st1000_t1Z0H1b2", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);

  h1_["ptTHnowwt_st1000_t1Z1H1b1"]  = cnt.make<TH1D>("ptTHnowwt_st1000_t1Z1H1b1", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptTHwt_st1000_t1Z1H1b1"]  = cnt.make<TH1D>("ptTHwt_st1000_t1Z1H1b1", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptTHnowwt_st1000_t1Z1H1b2"]  = cnt.make<TH1D>("ptTHnowwt_st1000_t1Z1H1b2", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptTHwt_st1000_t1Z1H1b2"]  = cnt.make<TH1D>("ptTHwt_st1000_t1Z1H1b2", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptTHnowwt_st1000_t1Z0H1b2"]  = cnt.make<TH1D>("ptTHnowwt_st1000_t1Z0H1b2", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);
  h1_["ptTHwt_st1000_t1Z0H1b2"]  = cnt.make<TH1D>("ptTHwt_st1000_t1Z0H1b2", ";p_{T}(Top tagged jet) [GeV];;" , 100, 0., 2000.);

  h1_["1b_1000_ht"]= cnt.make<TH1D>("1b_1000_ht", ";H_{T} [Gev];;", 100, 0., 3000.);
  h1_["1b_1000_st"] = cnt.make<TH1D>("1b_1000_st", ";S_{T} [Gev];;", 50, 0., 4000.) ;


  h1_["nob_1000_pt_zelel"] = cnt.make<TH1D>("nob_1000_pt_zelel",";p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV];;", 50, 0., 1000.) ;
  h1_["nob_1000_pt_zmumu"] = cnt.make<TH1D>("nob_1000_pt_zmumu",";p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV];;", 50, 0., 1000.) ;
  h1_["b_pt_z"+lep+lep] = cnt.make<TH1D>(("b_pt_z"+lep+lep).c_str(), "p_{T} (Z#rightarrow  l^{+}l^{-}) [GeV]", 50, 0., 1000.) ;
  h1_["b_st"] = cnt.make<TH1D>("b_st", "ST [GeV]", 50, 0., 4000.) ;
  h1_["pt_zlight_pre"] = fs->make<TH1D>("pt_zlight_pre", "p_{T} (Z + q_{light}) [GeV]", 100, 0., 2000.) ;
  h1_["pt_zb_pre"] = fs->make<TH1D>("pt_zb_pre", "p_{T} (Z + b) [GeV]", 100, 0., 2000.) ;
  h1_["pt_zc_pre"] = fs->make<TH1D>("pt_zc_pre", "p_{T} (Z + c) [GeV]", 100, 0., 2000.) ;
    
  h1_["nmergedtop_bf"] = sig.make<TH1D>("nmergedtop_bf", ";N(top boosted jets);;" , 9, -0.5, 8.5);
  h1_["nmergedZ_bf"] = sig.make<TH1D>("nmergedZ_bf", ";N(Z boosted jets);;" , 9, -0.5, 8.5);
  h1_["nmergedH_bf"] = sig.make<TH1D>("nmergedH_bf", ";N(H boosted jets);;" , 9, -0.5, 8.5);
  h1_["nmergedtotal_bf"] = sig.make<TH1D>("nmergedtotal_bf", ";N(total boosted jets);;" , 9, -0.5, 8.5);
  h1_["nmergedfrac_bf"] = sig.make<TH1D>("nmergedfrac_bf", ";N(fraction boosted jets);;" , 20, 0, 1.0);

  h1_["nmergedtop_af"] = sig.make<TH1D>("nmergedtop_af", ";N(top boosted jets);;" , 9, -0.5, 8.5);
  h1_["nmergedZ_af"] = sig.make<TH1D>("nmergedZ_af", ";N(Z boosted jets);;" , 9, -0.5, 8.5);
  h1_["nmergedH_af"] = sig.make<TH1D>("nmergedH_af", ";N(H boosted jets);;" , 9, -0.5, 8.5);
  h1_["nmergedtotal_af"] = sig.make<TH1D>("nmergedtotal_af", ";N(total boosted jets);;" , 9, -0.5, 8.5);
  h1_["nmergedfrac_af"] = sig.make<TH1D>("nmergedfrac_af", ";N(fraction boosted jets);;" , 20, 0, 1.0);


  //ak4matchedtoak8                                                                                                                                                                      
  h1_["mass_ak4matchedak8"]  = cat.make<TH1D>("mass_ak4matchedak8", ";M(AK4matched to AK8) [GeV];;", 100, 0., 400);
  h1_["pt_ak4matchedak8"]  = cat.make<TH1D>("pt_ak4matchedak8", ";Pt(AK4matched to AK8) [GeV];;", 100, 0., 1200);
  h1_["nak4matchedak8"] = cat.make<TH1D>("nak4matchedak8", ";N(AK4matched to AK8);;" , 21, -0.5, 20.5);


  h1_["dr_Wb_cnt"] = cat.make<TH1D>("dr_Wb_cnt", ";#DeltaR(Wj);;", 40, 0., 4.) ;
  h1_["dr_Wb_sig"] = cat.make<TH1D>("dr_Wb_sig", ";#DeltaR(Wj);;", 40, 0., 4.) ;
  
  h1_["dphi_Wb_cnt"] = cat.make<TH1D>("dphi_Wb_cnt", ";#Delta #phi(Wj);;", 20, -5., 5.) ;
  h1_["dphi_Wb_sig"] = cat.make<TH1D>("dphi_Wb_sig", ";#Delta #phi(Wj);;", 20, -5., 5.) ;
  ////electrons specific varaibles in EE and EB at preselection level
  if (zdecayMode_ == "zelel" && additionalPlots_){
    h1_["Eta_EB_el_pre"] = pre.make<TH1D>("Eta_EB_el_pre", ";Eta (EB);;", 100,-4,4) ;
    h1_["Eta_EE_el_pre"] = pre.make<TH1D>("Eta_EE_el_pre", ";Eta (EE);;", 100,-4,4) ;
    h1_["Iso03_EB_el_pre"] = pre.make<TH1D>("Iso03_EB_el_pre", ";Iso03 (EB);;", 100,0,0.3) ;
    h1_["Iso03_EE_el_pre"] = pre.make<TH1D>("Iso03_EE_el_pre", ";Iso03 (EE);;", 100,0,0.3) ;
    h1_["dEtaInSeed_EB_el_pre"] = pre.make<TH1D>("dEtaInSeed_EB_el_pre", ";dEtaInSeed (EB);;", 200,-0.05,0.05) ;
    h1_["dEtaInSeed_EE_el_pre"] = pre.make<TH1D>("dEtaInSeed_EE_el_pre", ";dEtaInSeed (EE);;", 200,-0.05,0.05) ;
    h1_["dPhiIn_EB_el_pre"] = pre.make<TH1D>("dPhiIn_EB_el_pre", ";dPhiIn (EB);;", 100,-0.2,0.2) ;
    h1_["dPhiIn_EE_el_pre"] = pre.make<TH1D>("dPhiIn_EE_el_pre", ";dPhiIn (EE);;", 100,-0.2,0.2);
    h1_["Dz_EB_el_pre"] = pre.make<TH1D>("Dz_EB_el_pre",";dZ (EB);;", 200,-0.1,0.1) ;
    h1_["Dz_EE_el_pre"] = pre.make<TH1D>("Dz_EE_el_pre", ";dZ (EE);;", 200,-0.4,0.4) ;
    h1_["Dxy_EB_el_pre"] = pre.make<TH1D>("Dxy_EB_el_pre", ";d0 (EB);;", 100,-0.1,0.1) ;
    h1_["Dxy_EE_el_pre"] = pre.make<TH1D>("Dxy_EE_el_pre", ";d0 (EE);;", 100,-0.1,0.1) ;
    h1_["SCETA_EB_el_pre"] = pre.make<TH1D>("scEta_EB_el_pre", ";Eta (EB);;", 100,-4,4) ;
    h1_["SCETA_EE_el_pre"] = pre.make<TH1D>("scEta_EE_el_pre", ";Eta (EE);;", 100,-4,4) ; 
    h1_["Full5x5siee_EB_el_pre"] = pre.make<TH1D>("Full5x5siee_EB_el_pre", ";Full5X5SigmaIEtaIEta (EB);;", 200,0,0.01) ;
    h1_["Full5x5siee_EE_el_pre"] = pre.make<TH1D>("Full5x5siee_EE_el_pre", ";Full5X5SigmaIEtaIEta (EE);;", 100,0,0.03) ;
    h1_["HoE_EB_el_pre"] = pre.make<TH1D>("HoE_EB_el_pre", ";H/E (EB);;", 200,0,0.05) ;
    h1_["HoE_EE_el_pre"] = pre.make<TH1D>("HoE_EE_el_pre", ";H/E (EE);;", 200,0,0.1) ;
    h1_["ooEmooP_EB_el_pre"] = pre.make<TH1D>("ooEmooP_EB_el_pre", ";(1/E - 1/P) (EB);;", 200,0,0.02) ;
    h1_["ooEmooP_EE_el_pre"] = pre.make<TH1D>("ooEmooP_EE_el_pre", ";(1/E - 1/P) (EE);;", 200,0,0.02) ;
    h1_["missHits_EB_el_pre"] = pre.make<TH1D>("missHits_EB_el_pre", ";Expected missing Hits (EB);;", 4,-0.5,3.5) ;
    h1_["missHits_EE_el_pre"] = pre.make<TH1D>("missHits_EE_el_pre", ";Expected missing Hits (EE);;", 4,-0.5,3.5) ;
    h1_["conveto_EB_el_pre"] = pre.make<TH1D>("conveto_EB_El_pre", ";has matched Conveto (EB);;", 4,-0.5,3.5) ;
    h1_["conveto_EE_el_pre"] = pre.make<TH1D>("conveto_EE_el_pre", ";has matched Conveto (EE);;", 4,-0.5,3.5) ;
 
    h1_["ElEB1"] = pre.make<TH1D>("ElEB1", "either gen Electron matched with 1st reco lepton (EE)", 11, 0.5, 11.5) ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(1, "dR matched") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(2, "+ dEtaInseed") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(3, "+ dPhiIn") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(4, "+ full5x5siee") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(5, "+ HoE ") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(6, "+ ooEmooP") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(7, "+ RelIsoEA") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(8, "+ !conv") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(9, "+ missHits") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(10, "+ Dxy") ;
    h1_["ElEB1"] -> GetXaxis() -> SetBinLabel(11, "+ Dz") ;

    h1_["ElEE1"] = pre.make<TH1D>("ElEE1", "either gen Electron matched with 1st reco Electron (EE)", 11, 0.5, 11.5) ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(1, "dR matched") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(2, "+ dEtaInseed") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(3, "+ dPhiIn") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(4, "+ full5x5siee") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(5, "+ HoE ") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(6, "+ ooEmooP") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(7, "+ RelIsoEA") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(8, "+ !conv") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(9, "+ missHits") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(10, "+ Dxy") ;
    h1_["ElEE1"] -> GetXaxis() -> SetBinLabel(11, "+ Dz") ;


    h1_["pt_leadinlepton_drmatched"]   =  pre.make<TH1D>("pt_leadinlepton_drmatched", "Pt [GeV]", 100, 0., 1000.) ;
    h1_["eta_leadinglepton_drmatched"]   =  pre.make<TH1D>("eta_leadinglepton_drmatched", "Eta", 50, -2.5, 2.5) ;

    h1_["pt_leadinlepton_drmatchedEB"]   =  pre.make<TH1D>("pt_leadinlepton_drmatchedEB", "Pt [GeV]", 100, 0., 1000.) ;
    h1_["eta_leadinglepton_drmatchedEB"]   =  pre.make<TH1D>("eta_leadinglepton_drmatchedEB", "Eta", 50, -2.5, 2.5) ;
    h1_["dr_elel_drmatchedEB"]   =  pre.make<TH1D>("dr_elel_drmatchedEB", "Dr", 40, 0., 4.) ;
    h1_["pt_leadinlepton_drmatchedEE"]   =  pre.make<TH1D>("pt_leadinlepton_drmatchedEE", "Pt [GeV]", 100, 0., 1000.) ;
    h1_["eta_leadinglepton_drmatchedEE"]   =  pre.make<TH1D>("eta_leadinglepton_drmatchedEE", "Eta", 50, -2.5, 2.5) ;
    h1_["dr_elel_drmatchedEE"]   =  pre.make<TH1D>("dr_elel_drmatchedEE", "Dr", 40, 0., 4.) ;


    h1_["pt_leadinlepton_idmatched"]   =  pre.make<TH1D>("pt_leadinlepton_idmatched", "Pt [GeV]", 100, 0., 1000.) ;
    h1_["eta_leadinglepton_idmatched"]   =  pre.make<TH1D>("eta_leadinglepton_idmatched", "Eta", 50, -2.5, 2.5) ;

    h1_["pt_leadinlepton_idmatchedEB"]   =  pre.make<TH1D>("pt_leadinlepton_idmatchedEB", "Pt [GeV]", 100, 0., 1000.) ;
    h1_["eta_leadinglepton_idmatchedEB"]   =  pre.make<TH1D>("eta_leadinglepton_idmatchedEB", "Eta", 50, -2.5, 2.5) ;
    h1_["dr_elel_idmatchedEB"]   =  pre.make<TH1D>("dr_elel_idmatchedEB", "Dr", 40, 0., 4.) ;
    h1_["pt_leadinlepton_idmatchedEE"]   =  pre.make<TH1D>("pt_leadinlepton_idmatchedEE", "Pt [GeV]", 100, 0., 1000.) ;
    h1_["eta_leadinglepton_idmatchedEE"]   =  pre.make<TH1D>("eta_leadinglepton_idmatchedEE", "Eta", 50, -2.5, 2.5) ;
    h1_["dr_elel_idmatchedEE"]   =  pre.make<TH1D>("dr_elel_idmatchedEE", "Dr", 40, 0., 4.) ;




 }
 
  if (zdecayMode_ == "zelel" ){
    h1_["Eta_EB_el_ex"] = pre.make<TH1D>("Eta_EB_el_ex", ";Eta (EB);;", 100,-4,4) ;
    h1_["Eta_EE_el_ex"] = pre.make<TH1D>("Eta_EE_el_ex", ";Eta (EE);;", 100,-4,4) ;
    h1_["Iso03_EB_el_ex"] = pre.make<TH1D>("Iso03_EB_el_ex", ";Iso03 (EB);;", 100,0,0.3) ;
    h1_["Iso03_EE_el_ex"] = pre.make<TH1D>("Iso03_EE_el_ex", ";Iso03 (EE);;", 100,0,0.3) ;
    h1_["dEtaInSeed_EB_el_ex"] = pre.make<TH1D>("dEtaInSeed_EB_el_ex", ";dEtaInSeed (EB);;", 200,-0.05,0.05) ;
    h1_["dEtaInSeed_EE_el_ex"] = pre.make<TH1D>("dEtaInSeed_EE_el_ex", ";dEtaInSeed (EE);;", 200,-0.05,0.05) ;
    h1_["dPhiIn_EB_el_ex"] = pre.make<TH1D>("dPhiIn_EB_el_ex", ";dPhiIn (EB);;", 100,-0.2,0.2) ;
    h1_["dPhiIn_EE_el_ex"] = pre.make<TH1D>("dPhiIn_EE_el_ex", ";dPhiIn (EE);;", 100,-0.2,0.2);
    h1_["Dz_EB_el_ex"] = pre.make<TH1D>("Dz_EB_el_ex",";dZ (EB);;", 200,-0.1,0.1) ;
    h1_["Dz_EE_el_ex"] = pre.make<TH1D>("Dz_EE_el_ex", ";dZ (EE);;", 200,-0.4,0.4) ;
    h1_["Dxy_EB_el_ex"] = pre.make<TH1D>("Dxy_EB_el_ex", ";d0 (EB);;", 100,-0.1,0.1) ;
    h1_["Dxy_EE_el_ex"] = pre.make<TH1D>("Dxy_EE_el_ex", ";d0 (EE);;", 100,-0.1,0.1) ;
    h1_["SCETA_EB_el_ex"] = pre.make<TH1D>("scEta_EB_el_ex", ";Eta (EB);;", 100,-4,4) ;
    h1_["SCETA_EE_el_ex"] = pre.make<TH1D>("scEta_EE_el_ex", ";Eta (EE);;", 100,-4,4) ;
    h1_["Full5x5siee_EB_el_ex"] = pre.make<TH1D>("Full5x5siee_EB_el_ex", ";Full5X5SigmaIEtaIEta (EB);;", 200,0,0.01) ;
    h1_["Full5x5siee_EE_el_ex"] = pre.make<TH1D>("Full5x5siee_EE_el_ex", ";Full5X5SigmaIEtaIEta (EE);;", 100,0,0.03) ;
    h1_["HoE_EB_el_ex"] = pre.make<TH1D>("HoE_EB_el_ex", ";H/E (EB);;", 200,0,0.05) ;
    h1_["HoE_EE_el_ex"] = pre.make<TH1D>("HoE_EE_el_ex", ";H/E (EE);;", 200,0,0.1) ;
    h1_["ooEmooP_EB_el_ex"] = pre.make<TH1D>("ooEmooP_EB_el_ex", ";(1/E - 1/P) (EB);;", 200,0,0.02) ;
    h1_["ooEmooP_EE_el_ex"] = pre.make<TH1D>("ooEmooP_EE_el_ex", ";(1/E - 1/P) (EE);;", 200,0,0.02) ;
    h1_["missHits_EB_el_ex"] = pre.make<TH1D>("missHits_EB_el_ex", ";Expected missing Hits (EB);;", 4,-0.5,3.5) ;
    h1_["missHits_EE_el_ex"] = pre.make<TH1D>("missHits_EE_el_ex", ";Expected missing Hits (EE);;", 4,-0.5,3.5) ;
    h1_["conveto_EB_el_ex"] = pre.make<TH1D>("conveto_EB_El_ex", ";has matched Conveto (EB);;", 4,-0.5,3.5) ;
    h1_["conveto_EE_el_ex"] = pre.make<TH1D>("conveto_EE_el_ex", ";has matched Conveto (EE);;", 4,-0.5,3.5) ;
  }





  if (zdecayMode_ == "zmumu"){

    h1_["Iso04_mu_pre"] = pre.make<TH1D>("Iso04_mu_pre", ";Iso04 ;;", 100,0,0.3) ;
    h1_["Iso04_mu1_pre"] = pre.make<TH1D>("Iso04_mu1_pre", ";Iso04 ;;", 100,0,0.3) ;
    h1_["Iso04_mu2_pre"] = pre.make<TH1D>("Iso04_mu2_pre", ";Iso04 ;;", 100,0,0.3) ;

    h1_["Dz_mu_pre"] = pre.make<TH1D>("Dz_mu_pre",";dZ ;;", 200,-0.1,0.1) ;
    h1_["Dxy_mu_pre"] = pre.make<TH1D>("Dxy_mu_pre", ";dxy;;", 100,-0.1,0.1) ;
    h1_["IsGlobalMuon_mu_pre"] = pre.make<TH1D>("IsGlobalMuon_mu_pre", ";IsGlobalMuon;;", 4,-0.5,3.5) ;
    h1_["IsPFMuon_mu_pre"] = pre.make<TH1D>("IsPFMuon_mu_pre", ";IsPFMuon;;", 4,-0.5,3.5) ;
    h1_["GlbTrkNormChi2_mu_pre"] = pre.make<TH1D>("GlbTrkNormChi2_mu_pre", ";GlbTrkNormChi2;;", 200,0.0,10.0) ;
    h1_["NumberValidMuonHits_mu_pre"] = pre.make<TH1D>("NumberValidMuonHits_mu_pre", ";NumberValidMuonHits;;", 20,-0.5,19.5) ;
    h1_["NumberMatchedStations_mu_pre"] = pre.make<TH1D>("NumberMatchedStations_mu_pre", ";NumberMatchedStations;;", 20,-0.5,19.5) ;
    h1_["NumberValidPixelHits_mu_pre"] = pre.make<TH1D>("NumberValidPixelHits_mu_pre", ";NumberValidPixelHits;;", 20,-0.5,19.5) ;
    h1_["NumberTrackerLayers_mu_pre"] = pre.make<TH1D>("NumberTrackerLayers_mu_pre", ";NumberTrackerLayers;;", 20,-0.5,19.5) ;

    
  }
    //additional plots
    h1_["nbjets_cnt"] = cnt.make<TH1D>("nbjets_cnt", ";N(b jets);;" , 11, -0.5,10.5) ; 
    h1_["nbjets_cat"] = cat.make<TH1D>("nbjets_cat", ";N(b jets);;" , 11, -0.5,10.5) ;

    //addition nb plots
    h1_["nbjets_met_sig"] = cnt.make<TH1D>("nbjets_met_sig", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_cnt"] = cnt.make<TH1D>("nbjets_met_cnt", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_0btagcnt"] = cnt.make<TH1D>("nbjets_met_0btagcnt", ";N(b jets);;" , 11, -0.5,10.5) ;
    h1_["nbjets_met_1btagcnt"] = cnt.make<TH1D>("nbjets_met_1btagcnt", ";N(b jets);;" , 11, -0.5,10.5) ;
  
    h1_["ht_met_0btagcnt"]   =  cnt.make<TH1D>( "ht_met_0btagcnt", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["1b_ht"]   =  cnt.make<TH1D>( "1b_ht", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["ht_met_1btagcnt"]   =  cnt.make<TH1D>( "ht_met_1btagcnt", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    h1_["lowmet_ht"]   =  cnt.make<TH1D>("lowmet_ht", ";H_{T} (AK4 jets) [GeV]", 100, 0., 4000.) ;
    


    h1_["st_met_0btagcnt"]   =  cnt.make<TH1D>( "st_met_0btagcnt", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["1b_st"]   =  cnt.make<TH1D>( "1b_st", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["st_met_1btagcnt"]   =  cnt.make<TH1D>( "st_met_1btagcnt", ";S_{T} [GeV]", 100, 0., 4000.) ;
    h1_["lowmet_st"]   =  cnt.make<TH1D>("lowmet_st", ";S_{T} [GeV]", 100, 0., 4000.) ;



    h1_["ptbjetleading_pre"]  = pre.make<TH1D>("ptbjetleading_pre", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["etabjetleading_pre"] = pre.make<TH1D>("etabjetleading_pre", ";#eta(leading b jet);;" , 80 ,-4. ,4.) ;

    h1_["ptbjetsubleading_pre"]  = pre.make<TH1D>("ptbjetsubleading_pre", ";p_{T}(subleading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["etabjetsubleading_pre"] = pre.make<TH1D>("etabjetsubleading_pre", ";#eta(subleading b jet);;" , 80 ,-4. ,4.) ;


    h1_["ptbjetleading_cat"]  = cat.make<TH1D>("ptbjetleading_cat", ";p_{T}(leading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["etabjetleading_cat"] = cat.make<TH1D>("etabjetleading_cat", ";#eta(leading b jet);;" , 80 ,-4. ,4.) ;

    h1_["ptbjetsubleading_cat"]  = cat.make<TH1D>("ptbjetsubleading_cat", ";p_{T}(subleading b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["etabjetsubleading_cat"] = cat.make<TH1D>("etabjetsubleading_cat", ";#eta(subleading b jet);;" , 80 ,-4. ,4.) ;


    h1_["ptbjet_cat"]  = cat.make<TH1D>("ptbjet_cat", ";p_{T}(b jet) [GeV];;" , 50, 0., 1000.) ;
    h1_["etabjet_cat"] = cat.make<TH1D>("etabjet_cat", ";#eta(b jet);;" , 80 ,-4. ,4.) ;




    /*
    if (optimizeReco_){
      //  h1_["genZ"] = sig.make<TH1D>("genZ", ";M (Gen Z Boson) [GeV];;", 20, 0., 200.);
      h1_["hadgentMass"] = sig.make<TH1D>("gentMass", ";M(Gen T quark) [GeV];;", 100, 0., 1200);
      h1_["lepgentMass"] = sig.make<TH1D>("lepgentMass", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
      //  h1_["ZJetMass"] = sig.make<TH1D>("ZJetMass", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
      // h1_["gentMass"] = sig.make<TH1D>("gentMass", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
      h1_["genH"] = sig.make<TH1D>("genH", ";M (Gen H Boson) [GeV];;", 20, 0., 200.); 
      h1_["HJetMass"] = sig.make<TH1D>("HJetMass", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
 
      //W tags (category B) 
      h1_["genH_B"] = sig.make<TH1D>("genH_B", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
      h1_["HJetMass_B"] = sig.make<TH1D>("HJetMass_B", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);     
      h1_["hadgentMass_B"] = sig.make<TH1D>("hadgentMass_B", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
      h1_["lepgentMass_B"] = sig.make<TH1D>("lepgentMass_B", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
     
      //category C
      h1_["genH_C"] = sig.make<TH1D>("genH_C", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
      h1_["HJetMass_C"] = sig.make<TH1D>("HJetMass_C", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
      h1_["hadgentMass_C"] = sig.make<TH1D>("hadgentMass_C", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
      h1_["lepgentMass_C"] = sig.make<TH1D>("lepgentMass_C", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
      
     //category D
      h1_["genH_D"] = sig.make<TH1D>("genH_D", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
      h1_["HJetMass_D"] = sig.make<TH1D>("HJetMass_D", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);      
      h1_["hadgentMass_D"] = sig.make<TH1D>("hadgentMass_D", ";M(Gen t quark) [GeV];;", 100, 0., 1200);
      h1_["lepgentMass_D"] = sig.make<TH1D>("lepgentMass_D", ";M(Gen t quark) [GeV];;", 100, 0., 1200);

      if (bosonMass_ == 91.2){
        //higgs
	//nonboosted
	h1_["hadgenZ_nonb"] = sig.make<TH1D>("hadgenZ_nonb", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadZJetMass_nonb"] = sig.make<TH1D>("hadZJetMass_nonb", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
	h1_["hadZJetPt_nonb"] = sig.make<TH1D>("hadZJetPt_nonb", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
	
	//boosted
	h1_["hadgenZ_b"] = sig.make<TH1D>("hadgenZ_b", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadZJetMass_b"] = sig.make<TH1D>("hadZJetMass_b", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
	h1_["hadZJetPt_b"] = sig.make<TH1D>("hadZJetPt_b", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       
	
	//   h1_["genZ"] = sig.make<TH1D>("genZ", ";M (Gen Z Boson) [GeV];;", 20, 0., 200.);
	//  h1_["ZJetMass"] = sig.make<TH1D>("ZJetMass", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
	//  h1_["ZJetPt"] = sig.make<TH1D>("ZJetPt", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
       
	h1_["hadgenZ"] = sig.make<TH1D>("hadgenZ", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadZJetMass"] = sig.make<TH1D>("hadZJetMass", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
	h1_["hadZJetPt"] = sig.make<TH1D>("hadZJetPt", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
	h1_["lepZ"] = sig.make<TH1D>("lepZ", ";M (lep Z Boson) [GeV];;", 20, 0., 200.);
	h1_["lepZPt"] = sig.make<TH1D>("lepZPt", ";Pt (Hleptonic Z Boson) [GeV];;", 50, 0., 1000.);
	
	//W tags (Cateogory B)
	h1_["hadgenZ_B"] = sig.make<TH1D>("hadgenZ_B", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadZJetMass_B"] = sig.make<TH1D>("hadZJetMass_B", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
	h1_["hadZJetPt_B"] = sig.make<TH1D>("hadZJetPt_B", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
	h1_["lepZ_B"] = sig.make<TH1D>("lepZ_B", ";M (lep Z Boson) [GeV];;", 20, 0., 200.);
	h1_["lepZPt_B"] = sig.make<TH1D>("lepZPt_B", ";Pt (Hleptonic Z Boson) [GeV];;", 50, 0., 1000.);
	
	//category C
	h1_["hadgenZ_C"] = sig.make<TH1D>("hadgenZ_C", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadZJetMass_C"] = sig.make<TH1D>("hadZJetMass_C", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
	h1_["hadZJetPt_C"] = sig.make<TH1D>("hadZJetPt_C", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
	h1_["lepZ_C"] = sig.make<TH1D>("lepZ_C", ";M (lep Z Boson) [GeV];;", 20, 0., 200.);
	h1_["lepZPt_C"] = sig.make<TH1D>("lepZPt_C", ";Pt (Hleptonic Z Boson) [GeV];;", 50, 0., 1000.);
	
	//category D
	h1_["hadgenZ_D"] = sig.make<TH1D>("hadgenZ_D", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadZJetMass_D"] = sig.make<TH1D>("hadZJetMass_D", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
	h1_["hadZJetPt_D"] = sig.make<TH1D>("hadZJetPt_D", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
	h1_["lepZ_D"] = sig.make<TH1D>("lepZ_D", ";M (lep Z Boson) [GeV];;", 20, 0., 200.);
	h1_["lepZPt_D"] = sig.make<TH1D>("lepZPt_D", ";Pt (Hleptonic Z Boson) [GeV];;", 50, 0., 1000.);
	
	//category E(z boson with jet masses equal)

	h1_["hadgenZ_E"] = sig.make<TH1D>("hadgenZ_E", ";M (Gen Z Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadZJetMass_E"] = sig.make<TH1D>("hadZJetMass_E", ";JetM (Hadronic Z Boson) [GeV];;", 20, 0., 200.);
	h1_["hadZJetPt_E"] = sig.make<TH1D>("hadZJetPt_E", ";Pt (Hadronic Z Boson) [GeV];;", 50, 0., 1000.);
	h1_["hadTJetPt_E"] = sig.make<TH1D>("hadTJetPt_E", ";Pt (Hadronic T Prime) [GeV];;", 50, 0., 1000.);
	
      }
      else{
	//higgs
	//nonboosted
	h1_["hadgenH_nonb"] = sig.make<TH1D>("hadgenH_nonb", ";M (Gen H Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadHJetMass_nonb"] = sig.make<TH1D>("hadHJetMass_nonb", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
	h1_["hadHJetPt_nonb"] = sig.make<TH1D>("hadHJetPt_nonb", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
	
	//boosted
	h1_["hadgenH_b"] = sig.make<TH1D>("hadgenH_b", ";M (Gen H Boson -hadronic) [GeV];;", 20, 0., 200.);
	h1_["hadHJetMass_b"] = sig.make<TH1D>("hadHJetMass_b", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
	h1_["hadHJetPt_b"] = sig.make<TH1D>("hadHJetPt_b", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);

	//h1_["genH"] = sig.make<TH1D>("genH", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
	// h1_["HJetMass"] = sig.make<TH1D>("HJetMass", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
	// h1_["HJetPt"] = sig.make<TH1D>("HJetPt", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);

	 h1_["hadgenH"] = sig.make<TH1D>("hadgenH", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetMass"] = sig.make<TH1D>("hadHJetMass", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetPt"] = sig.make<TH1D>("hadHJetPt", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
	 h1_["lepH"] = sig.make<TH1D>("lepH", ";M (leptonic H Boson) [GeV];;", 20, 0., 200.);
         h1_["lepHPt"] = sig.make<TH1D>("lepHPt", ";Pt (leptonic H Boson) [GeV];;", 50, 0., 1000.);

	 //W tags (category B)	 
         h1_["hadgenH_B"] = sig.make<TH1D>("hadgenH_B", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetMass_B"] = sig.make<TH1D>("hadHJetMass_B", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetPt_B"] = sig.make<TH1D>("hadHJetPt_B", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
         h1_["lepH_B"] = sig.make<TH1D>("lepH_B", ";M (leptonic H Boson) [GeV];;", 20, 0., 200.);
         h1_["lepHPt_B"] = sig.make<TH1D>("lepHPt_B", ";Pt (leptonic H Boson) [GeV];;", 50, 0., 1000.);
	 
	 //Category C
         h1_["hadgenH_C"] = sig.make<TH1D>("hadgenH_C", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetMass_C"] = sig.make<TH1D>("hadHJetMass_C", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetPt_C"] = sig.make<TH1D>("hadHJetPt_C", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
         h1_["lepH_C"] = sig.make<TH1D>("lepH_C", ";M (leptonic H Boson) [GeV];;", 20, 0., 200.);
         h1_["lepHPt_C"] = sig.make<TH1D>("lepHPt_C", ";Pt (leptonic H Boson) [GeV];;", 50, 0., 1000.);

	 //category D
         h1_["hadgenH_D"] = sig.make<TH1D>("hadgenH_D", ";M (Gen H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetMass_D"] = sig.make<TH1D>("hadHJetMass_D", ";JetM (Hadronic H Boson) [GeV];;", 20, 0., 200.);
         h1_["hadHJetPt_D"] = sig.make<TH1D>("hadHJetPt_D", ";Pt (Hadronic H Boson) [GeV];;", 50, 0., 1000.);
         h1_["lepH_D"] = sig.make<TH1D>("lepH_D", ";M (leptonic H Boson) [GeV];;", 20, 0., 200.);
         h1_["lepHPt_D"] = sig.make<TH1D>("lepHPt_D", ";Pt (leptonic H Boson) [GeV];;", 50, 0., 1000.);
      }
     
     h1_["hadtJetMass"] = sig.make<TH1D>("tJetMass-hadronic", ";JetM (Hadronic t quark) [GeV];;", 100, 0., 1000.);
     h1_["leptJetMass"] = sig.make<TH1D>("tJetMasslep", ";M (t Jet Leptonic);;", 100, 0., 1000.);
     h1_["hadtJetPt"] = sig.make<TH1D>("tJetPt-hadronic", ";Pt (Hadronic t quark) [GeV];;", 50, 0., 1000.);
     h1_["leptJetPt"] = sig.make<TH1D>("tJetPtlep", ";Pt (t Jet Leptonic);;", 50, 0., 1000.);
     h2_["genhadtJetMasshad"] = sig.make<TH2D>("genthadthad", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
     h2_["genhadtJetMasslep"] = sig.make<TH2D>("genthadtlep", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);

     //W tags (category B)

     h1_["hadtJetMass_B"] = sig.make<TH1D>("tJetMass-hadronic_B", ";JetM (Hadronic t quark) [GeV];;", 100, 0., 1000.);
     h1_["leptJetMass_B"] = sig.make<TH1D>("tJetMasslep_B", ";M (t Jet Leptonic);;", 100, 0., 1000.);
     h1_["hadtJetPt_B"] = sig.make<TH1D>("tJetPt-hadronic_B", ";Pt (Hadronic t quark) [GeV];;", 50, 0., 1000.);
     h1_["leptJetPt_B"] = sig.make<TH1D>("tJetPtlep_B", ";Pt (t Jet Leptonic);;", 50, 0., 1000.);
     h2_["genhadtJetMasshad_B"] = sig.make<TH2D>("genthadthad_B", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
     h2_["genhadtJetMasslep_B"] = sig.make<TH2D>("genthadtlep_B", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);   

     h1_["lepWJetMass_B"]= sig.make<TH1D>("lepWJetMass_B", ";M (leptonic W jet) [GeV];;", 20, 0., 200.);
     h1_["lepWJetPt_B"]= sig.make<TH1D>("lepWJetPt_B", ";Pt (leptonic W jet) [GeV];;", 100, 0., 1200.);
     h1_["hadWJetMass_B"]= sig.make<TH1D>("hadWJetMass_B", ";M (Hadronic W jet) [GeV];;", 20, 0., 200.);
     h1_["hadWJetPt_B"]= sig.make<TH1D>("hadWJetPt_B", ";Pt (hadronic W jet) [GeV];;", 100, 0., 1200.);

     //category C

     h1_["hadtJetMass_C"] = sig.make<TH1D>("tJetMass-hadronic_C", ";JetM (Hadronic t quark) [GeV];;", 100, 0., 1000.);
     h1_["leptJetMass_C"] = sig.make<TH1D>("tJetMasslep_C", ";M (t Jet Leptonic);;", 100, 0., 1000.);
     h1_["hadtJetPt_C"] = sig.make<TH1D>("tJetPt-hadronic_C", ";Pt (Hadronic t quark) [GeV];;", 50, 0., 1000.);
     h1_["leptJetPt_C"] = sig.make<TH1D>("tJetPtlep_C", ";Pt (t Jet Leptonic);;", 50, 0., 1000.);

     h1_["hadmergedJetMass_C"] = sig.make<TH1D>("hadmergedJetMass_C", ";merged jet mass (Hadronic) [GeV];;", 100, 0., 1000.);
     h1_["hadmergedJetPt_C"]   = sig.make<TH1D>("hadmergedJetPt_C", ";merged jet Pt (Hadronic) [GeV];;", 100, 0., 1000.);
     h1_["lepmergedJetMass_C"] = sig.make<TH1D>("lepmergedJetMass_C", ";merged jet mass (leptonic) [GeV];;", 100, 0., 1000.);
     h1_["lepmergedJetPt_C"]   = sig.make<TH1D>("lepmergedJetPt_C", ";merged jet Pt (Leptonic) [GeV];;", 100, 0., 1000.);

     h2_["genhadtJetMasshad_C"] = sig.make<TH2D>("genthadthad_C", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
     h2_["genhadtJetMasslep_C"] = sig.make<TH2D>("genthadtlep_C", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);

     //category D

     h1_["hadtJetMass_D"] = sig.make<TH1D>("tJetMass-hadronic_D", ";JetM (Hadronic t quark) [GeV];;", 100, 0., 1000.);
     h1_["leptJetMass_D"] = sig.make<TH1D>("tJetMasslep_D", ";M (t Jet Leptonic);;", 100, 0., 1000.);
     h1_["hadtJetPt_D"] = sig.make<TH1D>("tJetPt-hadronic_D", ";Pt (Hadronic t quark) [GeV];;", 50, 0., 1000.);
     h1_["leptJetPt_D"] = sig.make<TH1D>("tJetPtlep_D", ";Pt (t Jet Leptonic);;", 50, 0., 1000.);
     h2_["genhadtJetMasshad_D"] = sig.make<TH2D>("genthadthad_D", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
     h2_["genhadtJetMasslep_D"] = sig.make<TH2D>("genthadtlep_D", "; Mass(gen t vs. had t);;",50,0.,500, 50, 0., 1200.);
 
     h2_["hadWjetmassWjetpt_B"] = sig.make<TH2D>("hadWjetmassWjetpt_B", "; Mass Vs.Pt(merged W);;",50,0.,1200., 50, 0., 1200.);
     h2_["lepWjetmassWjetpt_B"] = sig.make<TH2D>("lepWjetmassWjetpt_B", "; Mass Vs.Pt(merged W);;",50,0.,1200., 50, 0., 1200.);
     h2_["hadWjetmassWjetpt_C"] = sig.make<TH2D>("hadWjetmassWjetpt_C", "; Mass Vs.Pt(merged W);;",50,0.,1200., 50, 0., 1200.);
     h2_["lepWjetmassWjetpt_C"] = sig.make<TH2D>("lepWjetmassWjetpt_C", "; Mass Vs.Pt(merged W);;",50,0.,1200., 50, 0., 1200.);
    }
    */
  }


    ////electrons specific varaibles in EE and EB at preselection level
  //    if (zdecayMode_ == "zelel" && additionalPlots_){
  //  h1_["Eta_EB_el_pre"] = pre.make<TH1D>("Eta_EB_el_pre", ";Eta (EB);;", 100,-4,4) ;
  //  h1_["Eta_EE_el_pre"] = pre.make<TH1D>("Eta_EE_el_pre", ";Eta (EE);;", 100,-4,4) ;
  //  h1_["Iso03_EB_el_pre"] = pre.make<TH1D>("Iso03_EB_el_pre", ";Iso03 (EB);;", 100,0,0.3) ;
  //  h1_["Iso03_EE_el_pre"] = pre.make<TH1D>("Iso03_EE_el_pre", ";Iso03 (EE);;", 100,0,0.3) ;
  //  h1_["dEtaInSeed_EB_el_pre"] = pre.make<TH1D>("dEtaInSeed_EB_el_pre", ";dEtaInSeed (EB);;", 200,-0.05,0.05) ;
  //  h1_["dEtaInSeed_EE_el_pre"] = pre.make<TH1D>("dEtaInSeed_EE_el_pre", ";dEtaInSeed (EE);;", 200,-0.05,0.05) ;
  //  h1_["dPhiIn_EB_el_pre"] = pre.make<TH1D>("dPhiIn_EB_el_pre", ";dPhiIn (EB);;", 100,-0.2,0.2) ;
  //  h1_["dPhiIn_EE_el_pre"] = pre.make<TH1D>("dPhiIn_EE_el_pre", ";dPhiIn (EE);;", 100,-0.2,0.2);
  //  h1_["Dz_EB_el_pre"] = pre.make<TH1D>("Dz_EB_el_pre",";dZ (EB);;", 200,-0.1,0.1) ;
  //  h1_["Dz_EE_el_pre"] = pre.make<TH1D>("Dz_EE_el_pre", ";dZ (EE);;", 200,-0.4,0.4) ;
  //  h1_["Dxy_EB_el_pre"] = pre.make<TH1D>("Dxy_EB_el_pre", ";d0 (EB);;", 100,-0.1,0.1) ;
  //  h1_["Dxy_EE_el_pre"] = pre.make<TH1D>("Dxy_EE_el_pre", ";d0 (EE);;", 100,-0.1,0.1) ;
  // }
  // }

  if (maketree_) {
    tree_ = fs->make<TTree>("tree", "HH4b") ; 
    os2ltree_.RegisterTree(tree_) ; 
  } 

}


void OS2LAna::endJob() {

  return ; 
}

bool my_compare(const std::pair<vlq::Jet, double >& firstItem, const std::pair<vlq::Jet, double >& secondItem)
{
  return firstItem.second < secondItem.second;
}


//bool OS2LAna::pairCompare(const std::pair<vlq::JetCollection, double >& firstItem, const std::pair<vlq::JetCollection, double >& secondItem) {
// return firstItem.first > secondItem.first;
//}



DEFINE_FWK_MODULE(OS2LAna);
