// Standard C++ includes
#include <memory>
#include <vector>
#include <iostream>
#include <boost/phoenix.hpp>

// ROOT includes
#include "TLorentzVector.h"
#include "TPRegexp.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"

// CMSSW framework includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

// CMSSW data formats
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/libminifloat.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/EgammaObject.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/fillCovariance.h"
#include "DataFormats/Scouting/interface/Run3ScoutingElectron.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPhoton.h"
#include "DataFormats/Scouting/interface/Run3ScoutingPFJet.h"
#include "DataFormats/Scouting/interface/Run3ScoutingVertex.h"
#include "DataFormats/Scouting/interface/Run3ScoutingTrack.h"
#include "DataFormats/Scouting/interface/Run3ScoutingMuon.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "DataFormats/Common/interface/AssociationMap.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenLumiInfoHeader.h"

// Other relevant CMSSW includes
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

using namespace std;

class EGammaOnly_ScoutingNanoAOD : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns, edm::one::WatchLuminosityBlocks> {
public:
  explicit EGammaOnly_ScoutingNanoAOD(const edm::ParameterSet&);
  ~EGammaOnly_ScoutingNanoAOD();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
    
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  virtual void clearVars();

  const edm::EDGetTokenT<edm::TriggerResults> trgResultsToken;

  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> > pVtxToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> > electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> > photonsToken;
  const edm::EDGetTokenT<double> rhoToken;

  // For L1 trig bits access
  bool doL1;
  edm::InputTag                algInputTag_;       
  edm::InputTag                extInputTag_;       
  edm::EDGetToken              algToken_;
  std::unique_ptr<l1t::L1TGlobalUtil> l1GtUtils_;
  std::vector<std::string>     l1Seeds_;

  // Interesting trigger results
  std::vector<bool> l1Result_;
  Int_t L1_TwoMu, L1_TwoMu_SQOS, L1_HT_ET, L1_OneJet, L1_TwoJet, L1_OneEG, L1_TwoEG;
  Int_t HLT_IsoMu27, HLT_Mu50, DST_Scouting_DoubleMu3, DST_Scouting_EG16EG12, DST_Scouting_EG30, DST_Scouting_JetHT, DST_Run3PFScouting, DST_HLTMuon_Run3PFScouting, HLT_OtherScoutingPFMonitor;

  // Primary Vertex
  UInt_t n_pVtx;
  vector<Float16_t> pVtx_x;
  vector<Float16_t> pVtx_y;
  vector<Float16_t> pVtx_z;
  vector<Float16_t> pVtx_xError;
  vector<Float16_t> pVtx_yError;
  vector<Float16_t> pVtx_zError;
  vector<Int_t> pVtx_trksize;
  vector<Float16_t> pVtx_chi2;
  vector<Int_t> pVtx_ndof;
  vector<Bool_t> pVtx_isvalidvtx;
        
  // Electron
  UInt_t n_ele;
  vector<Float16_t> Electron_pt;
  vector<Float16_t> Electron_eta;
  vector<Float16_t> Electron_phi;
  vector<Float16_t> Electron_m;
  vector< vector<Float16_t> > Electron_trkpt;
  vector< vector<Float16_t> > Electron_trketa;
  vector< vector<Float16_t> > Electron_trkphi;
  vector< vector<Float16_t> > Electron_trkd0;
  vector< vector<Float16_t> > Electron_trkdz;
  vector< vector<Int_t> > Electron_trkcharge;
  vector< vector<Float16_t> > Electron_trkrchi2;
  vector<Float16_t> Electron_detain;
  vector<Float16_t> Electron_dphiin;
  vector<Float16_t> Electron_sigmaietaieta;
  vector<Float16_t> Electron_hoe;
  vector<Float16_t> Electron_ooemoop;
  vector<Int_t>	Electron_missinghits;
  vector<Float16_t> Electron_ecaliso;
  vector<Float16_t> Electron_hcaliso;
  vector<Float16_t> Electron_tkiso;
  vector<Float16_t> Electron_r9;
  vector<Float16_t> Electron_smin;
  vector<Float16_t> Electron_smaj;
  vector<UInt_t> Electron_seedid;
  vector<Bool_t> Electron_rechitzerosuppression;

  //Photon
  UInt_t n_pho;
  vector<Float16_t> Photon_pt;
  vector<Float16_t> Photon_eta;
  vector<Float16_t> Photon_phi;
  vector<Float16_t> Photon_m;
  vector<Float16_t> Photon_sigmaietaieta;
  vector<Float16_t> Photon_hoe;
  vector<Float16_t> Photon_ecaliso;
  vector<Float16_t> Photon_hcaliso;
  vector<Float16_t> Photon_trkiso;
  vector<Float16_t> Photon_r9;
  vector<Float16_t> Photon_smin;
  vector<Float16_t> Photon_smaj;
  vector<Float16_t> Photon_seedid;
  vector<Bool_t> Photon_rechitzerosuppression;

  // Rho
  UInt_t n_rhoval;
  vector<Float16_t> rho;

  // TTree carrying the event weight information
  TTree* tree;

  //Run and lumisection
  Int_t run;
  Int_t event;
  Int_t lumSec;

};

EGammaOnly_ScoutingNanoAOD::EGammaOnly_ScoutingNanoAOD(const edm::ParameterSet& iConfig): 
  trgResultsToken(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerhltres"))),
  pVtxToken(consumes<std::vector<Run3ScoutingVertex> >(iConfig.getParameter<edm::InputTag>("primaryVtx"))), 
  electronsToken(consumes<std::vector<Run3ScoutingElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken(consumes<std::vector<Run3ScoutingPhoton> >(iConfig.getParameter<edm::InputTag>("photons"))),
  rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  doL1(iConfig.existsAs<bool>("doL1")?iConfig.getParameter<bool>("doL1"):false) {

  usesResource("TFileService");
  if (doL1) {
    algInputTag_ = iConfig.getParameter<edm::InputTag>("AlgInputTag");
    extInputTag_ = iConfig.getParameter<edm::InputTag>("l1tExtBlkInputTag");
    algToken_ = consumes<BXVector<GlobalAlgBlk>>(algInputTag_);
    l1Seeds_ = iConfig.getParameter<std::vector<std::string> >("l1Seeds");
    l1GtUtils_ = std::make_unique<l1t::L1TGlobalUtil>
      (iConfig, consumesCollector(), *this, algInputTag_, extInputTag_, l1t::UseEventSetupIn::Event);
  }
  else {
    l1Seeds_ = std::vector<std::string>();
    l1GtUtils_ = 0;
  }

  // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree", "tree");

  // Event details  
  tree->Branch("lumSec", &lumSec, "lumSec/i" );
  tree->Branch("run", &run, "run/i" );
  tree->Branch("event", &event, "event/i" );

  // Trigger Results
  tree->Branch("l1Result", "std::vector<bool>", &l1Result_, 32000, 0);
  tree->Branch("L1_TwoMu", &L1_TwoMu, "L1_TwoMu/i" );
  tree->Branch("L1_TwoMu_SQOS", &L1_TwoMu_SQOS, "L1_TwoMu_SQOS/i" );
  tree->Branch("L1_HT_ET", &L1_HT_ET, "L1_HT_ET/i" );
  tree->Branch("L1_OneJet", &L1_OneJet, "L1_OneJet/i" );
  tree->Branch("L1_TwoJet", &L1_TwoJet, "L1_TwoJet/i" );
  tree->Branch("L1_OneEG", &L1_OneEG, "L1_OneEG/i" );
  tree->Branch("L1_TwoEG", &L1_TwoEG, "L1_TwoEG/i" );
  tree->Branch("HLT_IsoMu27", &HLT_IsoMu27, "HLT_IsoMu27/i" );
  tree->Branch("HLT_Mu50", &HLT_Mu50, "HLT_Mu50/i" );
  tree->Branch("DST_Scouting_DoubleMu3", &DST_Scouting_DoubleMu3, "DST_Scouting_DoubleMu3/i" );
  tree->Branch("DST_Scouting_EG16EG12", &DST_Scouting_EG16EG12, "DST_Scouting_EG16EG12/i" );
  tree->Branch("DST_Scouting_EG30", &DST_Scouting_EG30, "DST_Scouting_EG30/i" );
  tree->Branch("DST_Scouting_JetHT", &DST_Scouting_JetHT, "DST_Scouting_JetHT/i" );
  tree->Branch("DST_Run3PFScouting", &DST_Run3PFScouting, "DST_Run3PFScouting/i" );
  tree->Branch("DST_HLTMuon_Run3PFScouting", &DST_HLTMuon_Run3PFScouting, "DST_HLTMuon_Run3PFScouting/i" );
  tree->Branch("HLT_OtherScoutingPFMonitor", &HLT_OtherScoutingPFMonitor, "HLT_OtherScoutingPFMonitor/i" );

  // Primary vertex info
  tree->Branch("n_pVtx", &n_pVtx, "n_pVtx/i");
  tree->Branch("pVtx_x", &pVtx_x);
  tree->Branch("pVtx_y", &pVtx_y);
  tree->Branch("pVtx_z", &pVtx_z);
  tree->Branch("pVtx_xError", &pVtx_xError);
  tree->Branch("pVtx_yError", &pVtx_yError);
  tree->Branch("pVtx_zError", &pVtx_zError);
  tree->Branch("pVtx_trksize", &pVtx_trksize);
  tree->Branch("pVtx_chi2", &pVtx_chi2);
  tree->Branch("pVtx_ndof", &pVtx_ndof);
  tree->Branch("pVtx_isvalidvtx", &pVtx_isvalidvtx);
  
  // Electrons
  tree->Branch("n_ele", &n_ele, "n_ele/i");
  tree->Branch("Electron_pt", &Electron_pt);
  tree->Branch("Electron_eta", &Electron_eta);
  tree->Branch("Electron_phi", &Electron_phi);
  tree->Branch("Electron_m", &Electron_m);
  tree->Branch("Electron_trkpt", &Electron_trkpt);
  tree->Branch("Electron_trketa", &Electron_trketa);
  tree->Branch("Electron_trkphi", &Electron_trkphi);
  tree->Branch("Electron_trkd0", &Electron_trkd0);
  tree->Branch("Electron_trkdz", &Electron_trkdz);
  tree->Branch("Electron_trkcharge", &Electron_trkcharge);
  tree->Branch("Electron_trkrchi2", &Electron_trkrchi2);
  tree->Branch("Electron_detain", &Electron_detain);
  tree->Branch("Electron_dphiin", &Electron_dphiin);
  tree->Branch("Electron_sigmaietaieta", &Electron_sigmaietaieta);
  tree->Branch("Electron_hoe", &Electron_hoe);
  tree->Branch("Electron_ooemoop", &Electron_ooemoop);
  tree->Branch("Electron_missinghits", &Electron_missinghits);
  tree->Branch("Electron_ecaliso", &Electron_ecaliso);
  tree->Branch("Electron_hcaliso", &Electron_hcaliso);
  tree->Branch("Electron_tkiso", &Electron_tkiso);
  tree->Branch("Electron_r9", &Electron_r9);
  tree->Branch("Electron_smin", &Electron_smaj);
  tree->Branch("Electron_smaj", &Electron_smin);
  tree->Branch("Electron_seedid", &Electron_seedid);
  tree->Branch("Electron_rechitzerosuppression", &Electron_rechitzerosuppression);
  
  // Photons
  tree->Branch("n_pho", &n_pho, "n_pho/i");
  tree->Branch("Photon_pt", &Photon_pt);
  tree->Branch("Photon_eta", &Photon_eta);
  tree->Branch("Photon_phi", &Photon_phi);	
  tree->Branch("Photon_m", &Photon_m);
  tree->Branch("Photon_sigmaietaieta", &Photon_sigmaietaieta);
  tree->Branch("Photon_hoe", &Photon_hoe);
  tree->Branch("Photon_ecaliso", &Photon_ecaliso);
  tree->Branch("Photon_hcaliso", &Photon_hcaliso);
  tree->Branch("Photon_trkiso", &Photon_trkiso);
  tree->Branch("Photon_r9", &Photon_r9);
  tree->Branch("Photon_smin", &Photon_smin);
  tree->Branch("Photon_smaj", &Photon_smaj);
  tree->Branch("Photon_seedid", &Photon_seedid);
  tree->Branch("Photon_rechitzerosuppression", &Photon_rechitzerosuppression);

  // Rho
  tree->Branch("n_rhoval", &n_rhoval, "n_rhoval/i");
  tree->Branch("rho", &rho);
}

EGammaOnly_ScoutingNanoAOD::~EGammaOnly_ScoutingNanoAOD() {
}

void EGammaOnly_ScoutingNanoAOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  using namespace edm;
  using namespace std;
  using namespace reco;
      
  edm::Handle<edm::TriggerResults> trigResults;
  iEvent.getByToken(trgResultsToken, trigResults);

  Handle<vector<Run3ScoutingVertex> > pVtxH;
  iEvent.getByToken(pVtxToken, pVtxH);
  bool pVtxValid = pVtxH.isValid();

  Handle<vector<Run3ScoutingElectron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  bool eleValid = electronsH.isValid();

  Handle<vector<Run3ScoutingPhoton> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  bool phoValid = photonsH.isValid();

  Handle<double > rhoH;
  iEvent.getByToken(rhoToken, rhoH);
  bool rhoValid = rhoH.isValid();

  run = iEvent.eventAuxiliary().run();
  event = iEvent.eventAuxiliary().event();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();

  L1_TwoMu = 0;
  L1_TwoMu_SQOS = 0;
  L1_HT_ET = 0;
  L1_OneJet = 0;
  L1_TwoJet = 0;
  L1_OneEG = 0;
  L1_TwoEG = 0;
  if(doL1) {
    l1GtUtils_->retrieveL1(iEvent,iSetup,algToken_);
    /*
    for( int r = 0; r<512; r++){
      string name ("empty");
      bool algoName_ = false;
      algoName_ = l1GtUtils_->getAlgNameFromBit(r,name);
      cout << "getAlgNameFromBit = " << algoName_  << endl;
      cout << "L1 bit number = " << r << " ; L1 bit name = " << name << endl;
    }
    */
    for( unsigned int iseed = 0; iseed < l1Seeds_.size(); iseed++ ) {
      bool l1htbit = 0;
      
      l1GtUtils_->getFinalDecisionByName(string(l1Seeds_[iseed]), l1htbit);
      //cout<<string(l1Seeds_[iseed])<<"  "<<l1htbit<<endl;
      l1Result_.push_back( l1htbit );
    }
    if(l1Result_[0] || l1Result_[1]) {
      L1_TwoMu = 1;
    }
    if(l1Result_[2] || l1Result_[3] || l1Result_[4] || l1Result_[5]) {
      L1_TwoMu_SQOS = 1;
    }
    if(l1Result_[6] || l1Result_[7] || l1Result_[8] || l1Result_[9] || l1Result_[10] || l1Result_[11] || l1Result_[12] || l1Result_[13]) {
      L1_HT_ET = 1;
    }
    if(l1Result_[14] || l1Result_[15]) {
      L1_OneJet = 1;
    }
    if(l1Result_[16] || l1Result_[17] || l1Result_[18]) {
      L1_TwoJet = 1;
    }
    if(l1Result_[19] || l1Result_[20] || l1Result_[21] || l1Result_[22] || l1Result_[23] || l1Result_[24]) {
      L1_OneEG = 1;
    }
    if(l1Result_[25] || l1Result_[26] || l1Result_[27] || l1Result_[28]) {
      L1_TwoEG = 1;
    }
  }

  // Access the trigger bits
  HLT_IsoMu27 = 0;
  HLT_Mu50 = 0;
  DST_Scouting_DoubleMu3 = 0;
  DST_Scouting_EG16EG12 = 0;
  DST_Scouting_EG30 = 0;
  DST_Scouting_JetHT = 0;
  DST_Run3PFScouting = 0;
  DST_HLTMuon_Run3PFScouting = 0;
  Int_t HLT_Ele115_CaloIdVT_GsfTrkIdT = 0;
  Int_t HLT_Ele35_WPTight_Gsf = 0;
  Int_t HLT_PFHT1050 = 0;
  Int_t HLT_Photon200 = 0;
  if( !trigResults.failedToGet() ) {
    int N_Triggers = trigResults->size();
    //cout<<"Number of triggers: "<<N_Triggers<<endl;
    const edm::TriggerNames & trigName = iEvent.triggerNames(*trigResults);
    for( int i_Trig = 0; i_Trig < N_Triggers; ++i_Trig ) {
      if (trigResults.product()->accept(i_Trig)) {
	//cout << "Path: " <<trigName.triggerName(i_Trig)<<"Results: "<<trigResults.product()->accept(i_Trig)<<endl;
	TString TrigPath =trigName.triggerName(i_Trig);
	if(TrigPath.Index("HLT_IsoMu27_v") >=0) HLT_IsoMu27 = 1; 
	if(TrigPath.Index("HLT_Mu50_v") >=0) HLT_Mu50 = 1; 
	if(TrigPath.Index("DST_Run3_DoubleMu3_PFScoutingPixelTracking_v") >=0) DST_Scouting_DoubleMu3 = 1;
	if(TrigPath.Index("DST_Run3_EG16_EG12_PFScoutingPixelTracking_v") >=0) DST_Scouting_EG16EG12 = 1;
	if(TrigPath.Index("DST_Run3_EG30_PFScoutingPixelTracking_v") >=0) DST_Scouting_EG30 = 1;
	if(TrigPath.Index("DST_Run3_JetHT_PFScoutingPixelTracking_v") >=0) DST_Scouting_JetHT = 1;
	if(TrigPath.Index("DST_Run3_PFScoutingPixelTracking_v") >=0) DST_Run3PFScouting = 1; 
	if(TrigPath.Index("DST_HLTMuon_Run3_PFScoutingPixelTracking_v") >=0) DST_HLTMuon_Run3PFScouting = 1; 
	if(TrigPath.Index("HLT_Ele115_CaloIdVT_GsfTrkIdT_v") >=0) HLT_Ele115_CaloIdVT_GsfTrkIdT = 1; 
	if(TrigPath.Index("HLT_Ele35_WPTight_Gsf_v") >=0) HLT_Ele35_WPTight_Gsf = 1; 
	if(TrigPath.Index("HLT_PFHT1050_v") >=0) HLT_PFHT1050 = 1; 
	if(TrigPath.Index("HLT_Photon200_v") >=0) HLT_Photon200 = 1; 
      }
    }
  } // End of loop for accessing the trigger bits

  if( HLT_Ele115_CaloIdVT_GsfTrkIdT==1 || HLT_Ele35_WPTight_Gsf==1 || HLT_PFHT1050==1 || HLT_Photon200==1 ) {
    HLT_OtherScoutingPFMonitor = 1;
  }
  else {
    HLT_OtherScoutingPFMonitor= 0;
  }  

  // Primary Vertex
  n_pVtx = 0;
  if(pVtxValid) {
    for (auto &pVtx : *pVtxH) {
      auto *pVtx_iter = &pVtx;
      pVtx_x.push_back(pVtx_iter->x());
      pVtx_y.push_back(pVtx_iter->y());
      pVtx_z.push_back(pVtx_iter->z());
      pVtx_xError.push_back(pVtx_iter->xError());
      pVtx_yError.push_back(pVtx_iter->yError());
      pVtx_zError.push_back(pVtx_iter->zError());
      pVtx_trksize.push_back(pVtx_iter->tracksSize());
      pVtx_chi2.push_back(pVtx_iter->chi2());
      pVtx_ndof.push_back(pVtx_iter->ndof());
      pVtx_isvalidvtx.push_back(pVtx_iter->isValidVtx());
      n_pVtx++;
    } // end primary vertex loop
  } // end pVtxValid condition

  // Electrons
  n_ele = 0;
  if(eleValid) {
    for (auto &ele : *electronsH) {
      auto *electrons_iter = &ele;
      Electron_pt.push_back(electrons_iter->pt());
      Electron_eta.push_back(electrons_iter->eta());
      Electron_phi.push_back(electrons_iter->phi());	
      Electron_m.push_back(electrons_iter->m());
      Electron_trkpt.push_back(electrons_iter->trkpt());
      Electron_trketa.push_back(electrons_iter->trketa());
      Electron_trkphi.push_back(electrons_iter->trkphi());	
      Electron_trkd0.push_back(electrons_iter->trkd0());
      Electron_trkdz.push_back(electrons_iter->trkdz());
      Electron_trkcharge.push_back(electrons_iter->trkcharge());
      Electron_trkrchi2.push_back(electrons_iter->trkchi2overndf());
      Electron_detain.push_back(electrons_iter->dEtaIn());
      Electron_dphiin.push_back(electrons_iter->dPhiIn());
      Electron_sigmaietaieta.push_back(electrons_iter->sigmaIetaIeta());
      Electron_hoe.push_back(electrons_iter->hOverE());	
      Electron_ooemoop.push_back(electrons_iter->ooEMOop());
      Electron_missinghits.push_back(electrons_iter->missingHits());
      Electron_ecaliso.push_back(electrons_iter->ecalIso());
      Electron_hcaliso.push_back(electrons_iter->hcalIso());
      Electron_tkiso.push_back(electrons_iter->trackIso());
      Electron_r9.push_back(electrons_iter->r9());
      Electron_smin.push_back(electrons_iter->sMin());
      Electron_smaj.push_back(electrons_iter->sMaj());
      Electron_seedid.push_back(electrons_iter->seedId());
      Electron_rechitzerosuppression.push_back(electrons_iter->rechitZeroSuppression());
      n_ele++;
    } // end electron loop
  } // end eleValid condition

  // Photons
  n_pho = 0;  
  if(phoValid) {

    for (auto &pho : *photonsH) {
      auto *photons_iter = &pho;
      Photon_pt.push_back(photons_iter->pt());
      Photon_eta.push_back(photons_iter->eta());
      Photon_phi.push_back(photons_iter->phi());
      Photon_m.push_back(photons_iter->m());
      Photon_sigmaietaieta.push_back(photons_iter->sigmaIetaIeta());
      Photon_hoe.push_back(photons_iter->hOverE());
      Photon_ecaliso.push_back(photons_iter->ecalIso());
      Photon_hcaliso.push_back(photons_iter->hcalIso());
      Photon_trkiso.push_back(photons_iter->trkIso());
      Photon_r9.push_back(photons_iter->r9());
      Photon_smin.push_back(photons_iter->sMin());
      Photon_smaj.push_back(photons_iter->sMaj());
      Photon_seedid.push_back(photons_iter->seedId());
      n_pho++;

    } // end photon loop
  } // end phoValid condition

  // Rho
  n_rhoval = 0;
  if(rhoValid) {
    rho.push_back(*rhoH);
    n_rhoval++;
  }
    
  tree->Fill();	
  clearVars(); 
}

void EGammaOnly_ScoutingNanoAOD::clearVars(){
  l1Result_.clear();
  pVtx_x.clear();
  pVtx_y.clear();
  pVtx_z.clear();
  pVtx_xError.clear();
  pVtx_yError.clear();
  pVtx_zError.clear();
  pVtx_trksize.clear();
  pVtx_chi2.clear();
  pVtx_ndof.clear();
  pVtx_isvalidvtx.clear();
  Electron_pt.clear();
  Electron_eta.clear();
  Electron_phi.clear();
  Electron_m.clear();
  Electron_trkpt.clear();
  Electron_trketa.clear();
  Electron_trkphi.clear();
  Electron_trkd0.clear();
  Electron_trkdz.clear();
  Electron_trkcharge.clear();
  Electron_trkrchi2.clear();
  Electron_detain.clear();
  Electron_dphiin.clear();
  Electron_sigmaietaieta.clear();
  Electron_hoe.clear();
  Electron_ooemoop.clear();
  Electron_missinghits.clear();
  Electron_ecaliso.clear();
  Electron_hcaliso.clear();
  Electron_tkiso.clear();
  Electron_r9.clear();
  Electron_smin.clear();
  Electron_smaj.clear();
  Electron_seedid.clear();
  Electron_rechitzerosuppression.clear();
  Photon_pt.clear();
  Photon_eta.clear();
  Photon_phi.clear();
  Photon_m.clear();
  Photon_sigmaietaieta.clear();
  Photon_hoe.clear();
  Photon_ecaliso.clear();
  Photon_hcaliso.clear();
  Photon_trkiso.clear();
  Photon_r9.clear();
  Photon_smin.clear();
  Photon_smaj.clear();
  Photon_seedid.clear();
  rho.clear();
}

void EGammaOnly_ScoutingNanoAOD::beginJob() {
  
}

void EGammaOnly_ScoutingNanoAOD::endJob() {
}

void EGammaOnly_ScoutingNanoAOD::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
}

void EGammaOnly_ScoutingNanoAOD::endRun(edm::Run const&, edm::EventSetup const&) {
}

void EGammaOnly_ScoutingNanoAOD::beginLuminosityBlock(edm::LuminosityBlock const& iLumi, edm::EventSetup const&) {
}

void EGammaOnly_ScoutingNanoAOD::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) {
}

void EGammaOnly_ScoutingNanoAOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(EGammaOnly_ScoutingNanoAOD);
