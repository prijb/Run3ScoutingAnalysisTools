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

  const edm::EDGetTokenT<std::vector<reco::Vertex> > oflpVtxToken;
  const edm::EDGetTokenT<std::vector<pat::Electron> > oflelectronsToken;
  const edm::EDGetTokenT<double> oflrhoToken;

  const edm::EDGetTokenT<std::vector<Run3ScoutingVertex> > pVtxToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingElectron> > electronsToken;
  const edm::EDGetTokenT<std::vector<Run3ScoutingPhoton> > photonsToken;
  const edm::EDGetTokenT<double> rhoToken;

  // Interesting trigger results
  Int_t HLT_IsoMu27, HLT_Mu50, DST_Run3PFScouting, DST_HLTMuon_Run3PFScouting, HLT_OtherScoutingPFMonitor;

  // Offline Primary Vertex
  const static int max_oflpv = 1000;
  UInt_t n_oflpVtx;
  vector<Float16_t> oflpVtx_x;
  vector<Float16_t> oflpVtx_y;
  vector<Float16_t> oflpVtx_z;
  vector<Float16_t> oflpVtx_xError;
  vector<Float16_t> oflpVtx_yError;
  vector<Float16_t> oflpVtx_zError;
  vector<Int_t> oflpVtx_trksize;
  vector<Float16_t> oflpVtx_chi2;
  vector<Int_t> oflpVtx_ndof;
  vector<Bool_t> oflpVtx_isvalidvtx;
        
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
        
  // Offline Electron
  UInt_t n_oflele;
  vector<Float16_t> OflElectron_energy;
  vector<Float16_t> OflElectron_pt;
  vector<Float16_t> OflElectron_eta;
  vector<Float16_t> OflElectron_phi;
  vector<Float16_t> OflElectron_d0;
  vector<Float16_t> OflElectron_dz;
  vector<Float16_t> OflElectron_detain;
  vector<Float16_t> OflElectron_dphiin;
  vector<Float16_t> OflElectron_sigmaietaieta;
  vector<Float16_t> OflElectron_hoe;
  vector<Float16_t> OflElectron_ooemoop;
  vector<Int_t> OflElectron_missinghits;
  vector<Int_t> OflElectron_charge;
  vector<Float16_t> OflElectron_photoniso;
  vector<Float16_t> OflElectron_neuthadroniso;
  vector<Float16_t> OflElectron_chrghadroniso;
  vector<Bool_t> OflElectron_conversionveto;

  // Electron
  UInt_t n_ele;
  vector<Float16_t> Electron_pt;
  vector<Float16_t> Electron_eta;
  vector<Float16_t> Electron_phi;
  vector<Float16_t> Electron_m;
  vector<Float16_t> Electron_d0;
  vector<Float16_t> Electron_dz;
  vector<Float16_t> Electron_detain;
  vector<Float16_t> Electron_dphiin;
  vector<Float16_t> Electron_sigmaietaieta;
  vector<Float16_t> Electron_hoe;
  vector<Float16_t> Electron_ooemoop;
  vector<Int_t>	Electron_missinghits;
  vector<Int_t> Electron_charge;
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

  // Offline Rho
  UInt_t n_oflrhoval;
  vector<Float16_t> oflrho;

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
  oflpVtxToken(consumes<std::vector<reco::Vertex> >(iConfig.getParameter<edm::InputTag>("oflprimaryVtx"))), 
  oflelectronsToken(consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("oflelectrons"))), 
  oflrhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("oflrho"))), 
  pVtxToken(consumes<std::vector<Run3ScoutingVertex> >(iConfig.getParameter<edm::InputTag>("primaryVtx"))), 
  electronsToken(consumes<std::vector<Run3ScoutingElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))), 
  photonsToken(consumes<std::vector<Run3ScoutingPhoton> >(iConfig.getParameter<edm::InputTag>("photons"))),
  rhoToken(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))) {

  usesResource("TFileService");

  // Access the TFileService
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("tree", "tree");

  // Event details  
  tree->Branch("lumSec", &lumSec, "lumSec/i" );
  tree->Branch("run", &run, "run/i" );
  tree->Branch("event", &event, "event/i" );

  // Trigger Results
  tree->Branch("HLT_IsoMu27", &HLT_IsoMu27, "HLT_IsoMu27/i" );
  tree->Branch("HLT_Mu50", &HLT_Mu50, "HLT_Mu50/i" );
  tree->Branch("DST_Run3PFScouting", &DST_Run3PFScouting, "DST_Run3PFScouting/i" );
  tree->Branch("DST_HLTMuon_Run3PFScouting", &DST_HLTMuon_Run3PFScouting, "DST_HLTMuon_Run3PFScouting/i" );
  tree->Branch("HLT_OtherScoutingPFMonitor", &HLT_OtherScoutingPFMonitor, "HLT_OtherScoutingPFMonitor/i" );

  // Offline Primary vertex info
  tree->Branch("n_oflpVtx", &n_oflpVtx, "n_oflpVtx/i");
  tree->Branch("oflpVtx_x", &oflpVtx_x);
  tree->Branch("oflpVtx_y", &oflpVtx_y);
  tree->Branch("oflpVtx_z", &oflpVtx_z);
  tree->Branch("oflpVtx_xError", &oflpVtx_xError);
  tree->Branch("oflpVtx_yError", &oflpVtx_yError);
  tree->Branch("oflpVtx_zError", &oflpVtx_zError);
  tree->Branch("oflpVtx_trksize", &oflpVtx_trksize);
  tree->Branch("oflpVtx_chi2", &oflpVtx_chi2);
  tree->Branch("oflpVtx_ndof", &oflpVtx_ndof);
  tree->Branch("oflpVtx_isvalidvtx", &oflpVtx_isvalidvtx);
  
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
  
  // Offline Electrons
  tree->Branch("n_oflele", &n_oflele, "n_oflele/i");
  tree->Branch("OflElectron_energy", &OflElectron_energy);
  tree->Branch("OflElectron_pt", &OflElectron_pt);
  tree->Branch("OflElectron_eta", &OflElectron_eta);
  tree->Branch("OflElectron_phi", &OflElectron_phi);
  tree->Branch("OflElectron_d0", &OflElectron_d0);
  tree->Branch("OflElectron_dz", &OflElectron_dz);
  tree->Branch("OflElectron_detain", &OflElectron_detain);
  tree->Branch("OflElectron_dphiin", &OflElectron_dphiin);
  tree->Branch("OflElectron_sigmaietaieta", &OflElectron_sigmaietaieta);
  tree->Branch("OflElectron_hoe", &OflElectron_hoe);
  tree->Branch("OflElectron_ooemoop", &OflElectron_ooemoop);
  tree->Branch("OflElectron_missinghits", &OflElectron_missinghits);
  tree->Branch("OflElectron_charge", &OflElectron_charge);
  tree->Branch("OflElectron_photoniso", &OflElectron_photoniso);
  tree->Branch("OflElectron_neuthadroniso", &OflElectron_neuthadroniso);
  tree->Branch("OflElectron_chrghadroniso", &OflElectron_chrghadroniso);
  tree->Branch("OflElectron_conversionveto", &OflElectron_conversionveto);

  // Electrons
  tree->Branch("n_ele", &n_ele, "n_ele/i");
  tree->Branch("Electron_pt", &Electron_pt);
  tree->Branch("Electron_eta", &Electron_eta);
  tree->Branch("Electron_phi", &Electron_phi);
  tree->Branch("Electron_m", &Electron_m);
  tree->Branch("Electron_d0", &Electron_d0);
  tree->Branch("Electron_dz", &Electron_dz);
  tree->Branch("Electron_detain", &Electron_detain);
  tree->Branch("Electron_dphiin", &Electron_dphiin);
  tree->Branch("Electron_sigmaietaieta", &Electron_sigmaietaieta);
  tree->Branch("Electron_hoe", &Electron_hoe);
  tree->Branch("Electron_ooemoop", &Electron_ooemoop);
  tree->Branch("Electron_missinghits", &Electron_missinghits);
  tree->Branch("Electron_charge", &Electron_charge);
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

  // Offline Rho
  tree->Branch("n_oflrhoval", &n_oflrhoval, "n_oflrhoval/i");
  tree->Branch("oflrho", &oflrho);

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

  Handle<vector<reco::Vertex> > oflpVtxH;
  iEvent.getByToken(oflpVtxToken, oflpVtxH);
  bool oflpVtxValid = oflpVtxH.isValid();

  Handle<vector<Run3ScoutingVertex> > pVtxH;
  iEvent.getByToken(pVtxToken, pVtxH);
  bool pVtxValid = pVtxH.isValid();

  Handle<vector<pat::Electron> > oflelectronsH;
  iEvent.getByToken(oflelectronsToken, oflelectronsH);
  bool ofleleValid = oflelectronsH.isValid();

  Handle<vector<Run3ScoutingElectron> > electronsH;
  iEvent.getByToken(electronsToken, electronsH);
  bool eleValid = electronsH.isValid();

  Handle<vector<Run3ScoutingPhoton> > photonsH;
  iEvent.getByToken(photonsToken, photonsH);
  bool phoValid = photonsH.isValid();

  Handle<double > oflrhoH;
  iEvent.getByToken(oflrhoToken, oflrhoH);
  bool oflrhoValid = oflrhoH.isValid();

  Handle<double > rhoH;
  iEvent.getByToken(rhoToken, rhoH);
  bool rhoValid = rhoH.isValid();

  run = iEvent.eventAuxiliary().run();
  event = iEvent.eventAuxiliary().event();
  lumSec = iEvent.eventAuxiliary().luminosityBlock();

  // Access the trigger bits
  HLT_IsoMu27 = 0;
  HLT_Mu50 = 0;
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
  // Offline Primary Vertex
  n_oflpVtx = 0;
  if(oflpVtxValid) {
    for (auto &oflpVtx : *oflpVtxH) {
      auto *oflpVtx_iter = &oflpVtx;
      oflpVtx_x.push_back(oflpVtx_iter->x());
      oflpVtx_y.push_back(oflpVtx_iter->y());
      oflpVtx_z.push_back(oflpVtx_iter->z());
      oflpVtx_xError.push_back(oflpVtx_iter->xError());
      oflpVtx_yError.push_back(oflpVtx_iter->yError());
      oflpVtx_zError.push_back(oflpVtx_iter->zError());
      oflpVtx_trksize.push_back(oflpVtx_iter->tracksSize());
      oflpVtx_chi2.push_back(oflpVtx_iter->chi2());
      oflpVtx_ndof.push_back(oflpVtx_iter->ndof());
      oflpVtx_isvalidvtx.push_back(oflpVtx_iter->isValid());
      n_oflpVtx++;
    } // end offline primary vertex loop
  } // end oflpVtxValid condition
  
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

  // Offline Electrons
  n_oflele = 0;
  if(ofleleValid) {
    for (auto &oflele : *oflelectronsH) {
      auto *oflelectrons_iter = &oflele;
      OflElectron_energy.push_back(oflelectrons_iter->energy());
      OflElectron_pt.push_back(oflelectrons_iter->pt());
      OflElectron_eta.push_back(oflelectrons_iter->eta());
      OflElectron_phi.push_back(oflelectrons_iter->phi());	
      OflElectron_d0.push_back(oflelectrons_iter->dB(pat::Electron::PV3D));
      OflElectron_dz.push_back(oflelectrons_iter->dB(pat::Electron::PVDZ));
      OflElectron_detain.push_back(oflelectrons_iter->deltaEtaSuperClusterTrackAtVtx());
      OflElectron_dphiin.push_back(oflelectrons_iter->deltaPhiSuperClusterTrackAtVtx());
      OflElectron_sigmaietaieta.push_back(oflelectrons_iter->full5x5_sigmaIetaIeta());
      OflElectron_hoe.push_back(oflelectrons_iter->hadronicOverEm());
      OflElectron_ooemoop.push_back( (1.0/oflelectrons_iter->correctedEcalEnergy()) - (1.0/oflelectrons_iter->p()) );
      OflElectron_missinghits.push_back(oflelectrons_iter->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      OflElectron_charge.push_back(oflelectrons_iter->charge());
      OflElectron_photoniso.push_back(oflelectrons_iter->photonIso());
      OflElectron_neuthadroniso.push_back(oflelectrons_iter->neutralHadronIso());
      OflElectron_chrghadroniso.push_back(oflelectrons_iter->chargedHadronIso());
      OflElectron_conversionveto.push_back(oflelectrons_iter->passConversionVeto());
      n_oflele++;

    } // end offline electron loop
  } // end ofleleValid condition

  // Electrons
  n_ele = 0;
  if(eleValid) {
    for (auto &ele : *electronsH) {
      auto *electrons_iter = &ele;
      Electron_pt.push_back(electrons_iter->pt());
      Electron_eta.push_back(electrons_iter->eta());
      Electron_phi.push_back(electrons_iter->phi());	
      Electron_m.push_back(electrons_iter->m());
      Electron_d0.push_back(electrons_iter->d0());
      Electron_dz.push_back(electrons_iter->dz());
      Electron_detain.push_back(electrons_iter->dEtaIn());
      Electron_dphiin.push_back(electrons_iter->dPhiIn());
      Electron_sigmaietaieta.push_back(electrons_iter->sigmaIetaIeta());
      Electron_hoe.push_back(electrons_iter->hOverE());	
      Electron_ooemoop.push_back(electrons_iter->ooEMOop());
      Electron_missinghits.push_back(electrons_iter->missingHits());
      Electron_charge.push_back(electrons_iter->charge());
      Electron_ecaliso.push_back(electrons_iter->ecalIso());
      Electron_hcaliso.push_back(electrons_iter->hcalIso());
      Electron_tkiso.push_back(electrons_iter->trackIso());
      Electron_r9.push_back(electrons_iter->r9());
      Electron_smin.push_back(electrons_iter->sMin());
      Electron_smaj.push_back(electrons_iter->sMaj());
      Electron_seedid.push_back(electrons_iter->seedId());
      n_ele++;

    } // end electron loop
  } // end eleValid condition

  // Muons
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

  // Offline Rho
  n_oflrhoval = 0;
  if(oflrhoValid) {
    oflrho.push_back(*oflrhoH);
    n_oflrhoval++;
  }

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
  oflpVtx_x.clear();
  oflpVtx_y.clear();
  oflpVtx_z.clear();
  oflpVtx_xError.clear();
  oflpVtx_yError.clear();
  oflpVtx_zError.clear();
  oflpVtx_trksize.clear();
  oflpVtx_chi2.clear();
  oflpVtx_ndof.clear();
  oflpVtx_isvalidvtx.clear();
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
  OflElectron_energy.clear();
  OflElectron_pt.clear();
  OflElectron_eta.clear();
  OflElectron_phi.clear();
  OflElectron_d0.clear();
  OflElectron_dz.clear();
  OflElectron_detain.clear();
  OflElectron_dphiin.clear();
  OflElectron_sigmaietaieta.clear();
  OflElectron_hoe.clear();
  OflElectron_ooemoop.clear();
  OflElectron_missinghits.clear();
  OflElectron_charge.clear();
  OflElectron_photoniso.clear();
  OflElectron_neuthadroniso.clear();
  OflElectron_chrghadroniso.clear();
  OflElectron_conversionveto.clear();
  Electron_pt.clear();
  Electron_eta.clear();
  Electron_phi.clear();
  Electron_m.clear();
  Electron_d0.clear();
  Electron_dz.clear();
  Electron_detain.clear();
  Electron_dphiin.clear();
  Electron_sigmaietaieta.clear();
  Electron_hoe.clear();
  Electron_ooemoop.clear();
  Electron_missinghits.clear();
  Electron_charge.clear();
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
  oflrho.clear();
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
