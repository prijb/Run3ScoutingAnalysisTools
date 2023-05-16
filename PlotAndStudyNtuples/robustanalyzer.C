#include "robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
robustanalyzer::robustanalyzer(TString filename, TString outfilename, int numCores){

  nC = numCores;

  cout<<"Initializing for file: "<<filename<<endl;
  TChain* chain = new TChain("mmtree/tree");
  chain->Add(filename);

  tree = new TTreeReader(chain);
  l1_2mu = new TTreeReaderValue<unsigned int>((*tree), "L1_TwoMu");
  l1_2musqos = new TTreeReaderValue<unsigned int>((*tree), "L1_TwoMu_SQOS");
  l1_met = new TTreeReaderValue<unsigned int>((*tree), "L1_HT_ET");
  l1_1jet = new TTreeReaderValue<unsigned int>((*tree), "L1_OneJet");
  l1_2jet = new TTreeReaderValue<unsigned int>((*tree), "L1_TwoJet");
  l1_1eg = new TTreeReaderValue<unsigned int>((*tree), "L1_OneEG");
  l1_2eg = new TTreeReaderValue<unsigned int>((*tree), "L1_TwoEG");
  hlt_isomu27 = new TTreeReaderValue<unsigned int>((*tree), "HLT_IsoMu27");
  hlt_mu50 = new TTreeReaderValue<unsigned int>((*tree), "HLT_Mu50");
  hlt_scouting= new TTreeReaderValue<unsigned int>((*tree), "DST_Run3PFScouting");

  n_ofle = new TTreeReaderValue<unsigned int>((*tree), "n_oflele");
  ofle_pt = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_pt");
  ofle_eta = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_eta");
  ofle_phi = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_phi");
  ofle_d0 = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_d0");
  ofle_dz = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_dz");
  ofle_detain = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_detain");
  ofle_dphiin = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_dphiin");
  ofle_sigmaietaieta = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_sigmaietaieta");
  ofle_hoe = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_hoe");
  ofle_ooemoop = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_ooemoop");
  ofle_mhits = new TTreeReaderValue<vector<int>>((*tree), "OflElectron_missinghits");
  ofle_charge = new TTreeReaderValue<vector<int>>((*tree), "OflElectron_charge");
  ofle_phoiso = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_photoniso");
  ofle_nthadiso = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_neuthadroniso");
  ofle_cghadiso = new TTreeReaderValue<vector<float>>((*tree), "OflElectron_chrghadroniso");
  ofle_conveto = new TTreeReaderValue<vector<bool>>((*tree), "OflElectron_conversionveto");

  n_ele = new TTreeReaderValue<unsigned int>((*tree), "n_ele");
  ele_pt = new TTreeReaderValue<vector<float>>((*tree), "Electron_pt");
  ele_eta = new TTreeReaderValue<vector<float>>((*tree), "Electron_eta");
  ele_phi = new TTreeReaderValue<vector<float>>((*tree), "Electron_phi");
  ele_m = new TTreeReaderValue<vector<float>>((*tree), "Electron_m");
  ele_d0 = new TTreeReaderValue<vector<float>>((*tree), "Electron_d0");
  ele_dz = new TTreeReaderValue<vector<float>>((*tree), "Electron_dz");
  ele_detain = new TTreeReaderValue<vector<float>>((*tree), "Electron_detain");
  ele_dphiin = new TTreeReaderValue<vector<float>>((*tree), "Electron_dphiin");
  ele_sigmaietaieta = new TTreeReaderValue<vector<float>>((*tree), "Electron_sigmaietaieta");
  ele_hoe = new TTreeReaderValue<vector<float>>((*tree), "Electron_hoe");
  ele_ooemoop = new TTreeReaderValue<vector<float>>((*tree), "Electron_ooemoop");
  ele_mhits = new TTreeReaderValue<vector<int>>((*tree), "Electron_missinghits");
  ele_charge = new TTreeReaderValue<vector<int>>((*tree), "Electron_charge");
  ele_ecaliso = new TTreeReaderValue<vector<float>>((*tree), "Electron_ecaliso");
  ele_hcaliso = new TTreeReaderValue<vector<float>>((*tree), "Electron_hcaliso");
  ele_tkiso = new TTreeReaderValue<vector<float>>((*tree), "Electron_tkiso");
  ele_r9 = new TTreeReaderValue<vector<float>>((*tree), "Electron_r9");
  ele_smin = new TTreeReaderValue<vector<float>>((*tree), "Electron_smin");
  ele_smaj = new TTreeReaderValue<vector<float>>((*tree), "Electron_smaj");

  n_rho = new TTreeReaderValue<unsigned int>((*tree), "n_rhoval");
  rho = new TTreeReaderValue<vector<float>>((*tree), "rho");
  n_oflrho = new TTreeReaderValue<unsigned int>((*tree), "n_oflrhoval");
  oflrho = new TTreeReaderValue<vector<float>>((*tree), "oflrho");

  n_pho = new TTreeReaderValue<unsigned int>((*tree), "n_pho");
  pho_pt = new TTreeReaderValue<vector<float>>((*tree), "Photon_pt");
  pho_eta = new TTreeReaderValue<vector<float>>((*tree), "Photon_eta");
  pho_phi = new TTreeReaderValue<vector<float>>((*tree), "Photon_phi");
  pho_m = new TTreeReaderValue<vector<float>>((*tree), "Photon_m");
  pho_sigmaietaieta = new TTreeReaderValue<vector<float>>((*tree), "Photon_sigmaietaieta");
  pho_hoe = new TTreeReaderValue<vector<float>>((*tree), "Photon_hoe");
  pho_ecaliso = new TTreeReaderValue<vector<float>>((*tree), "Photon_ecaliso");
  pho_hcaliso = new TTreeReaderValue<vector<float>>((*tree), "Photon_hcaliso");
  pho_r9 = new TTreeReaderValue<vector<float>>((*tree), "Photon_r9");
  pho_smin = new TTreeReaderValue<vector<float>>((*tree), "Photon_smin");
  pho_smaj = new TTreeReaderValue<vector<float>>((*tree), "Photon_smaj");

  outfile = new TFile(outfilename,"RECREATE");
}

// Fill the root file, close the root file, and handle deletions
robustanalyzer::~robustanalyzer() {

  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void robustanalyzer::analyzersinglefile(int splitCnt) { // Assume splitCnt to range from 0 to nCores

  int totEntries = tree->GetEntries(true);
  cout<<"Total number of entries: "<<totEntries<<endl;

  // Verfied that this logic to parallelize works
  int nCores = nC;
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  if(beginevent>=totEntries) return;
  endevent = endevent<totEntries?endevent:totEntries;
  tree->SetEntriesRange(beginevent, endevent);
  cout<<"Processing events in range: [ "<<beginevent<<" , "<<endevent<<" )"<<endl;
  int event = beginevent-1;

  // Count events passing certain selections
  int nosel=0;

  addhist_oflsct_electron("nosel_oflsct_ele");

  addhist_oflsct_electron("ptgt2etalt2p5_oflsct_ele");
  addhist_oflsct_electron("noeg_ptgt2etalt2p5_oflsct_ele");
  addhist_oflsct_electron("eg1_ptgt2etalt2p5_oflsct_ele");
  addhist_oflsct_electron("eg2_ptgt2etalt2p5_oflsct_ele");

  addhist_oflsct_electron("ptgt2etlt1p44_oflsct_ele");
  addhist_oflsct_electron("noeg_ptgt2etlt1p44_oflsct_ele");
  addhist_oflsct_electron("eg1_ptgt2etlt1p44_oflsct_ele");
  addhist_oflsct_electron("eg2_ptgt2etlt1p44_oflsct_ele");

  addhist_oflsct_electron("ptgt2etgt1p57lt2_oflsct_ele");
  addhist_oflsct_electron("noeg_ptgt2etgt1p57lt2_oflsct_ele");
  addhist_oflsct_electron("eg1_ptgt2etgt1p57lt2_oflsct_ele");
  addhist_oflsct_electron("eg2_ptgt2etgt1p57lt2_oflsct_ele");

  addhist_oflsct_electron("ptgt2etgt2lt2p5_oflsct_ele");
  addhist_oflsct_electron("noeg_ptgt2etgt2lt2p5_oflsct_ele");
  addhist_oflsct_electron("eg1_ptgt2etgt2lt2p5_oflsct_ele");
  addhist_oflsct_electron("eg2_ptgt2etgt2lt2p5_oflsct_ele");
  /* For gen matching
  addhist_oflsct_electron("ptgt2lt5etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt5lt8etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt8lt11etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt11lt14etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt14lt18etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt18lt22etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt22lt26etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt26lt30etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt30lt50etalt2p5_oflsct_ele");
  addhist_oflsct_electron("ptgt50lt100etalt2p5_oflsct_ele");
  */
  
  addhist_electron("nosel_ele");
  addhist_electron("mutrigsel_ele");
  addhist_electron("muAscouttrigsel_ele");
  addhist_electron("mutrigsel_noneg_ele");
  addhist_electron("muAscouttrigsel_noneg_ele");
  addhist_electron("mutrigsel_1egonly_ele");
  addhist_electron("muAscouttrigsel_1egonly_ele");
  addhist_electron("mutrigsel_2egmin_ele");
  addhist_electron("muAscouttrigsel_2egmin_ele");
  addhist_electron("mutrig_elptgt50_ele");
  addhist_electron("muAscouttrig_elptgt50_ele");
  
  addhist_photon("nosel_pho");
  addhist_photon("mutrigsel_pho");
  addhist_photon("muAscouttrigsel_pho");
  addhist_photon("mutrigsel_noneg_pho");
  addhist_photon("muAscouttrigsel_noneg_pho");
  addhist_photon("mutrigsel_1egonly_pho");
  addhist_photon("muAscouttrigsel_1egonly_pho");
  addhist_photon("mutrigsel_2egmin_pho");
  addhist_photon("muAscouttrigsel_2egmin_pho");
  addhist_photon("mutrig_elptgt50_pho");
  addhist_photon("muAscouttrig_elptgt50_pho");
  
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%100000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    analyzerelectrons();
    analyzerphotons();

  } // End of loop on events

  cout<<totEntries<<"\t"<<endevent<<"\t"<<beginevent<<"\t"<<nosel<<endl;
}

// Analyzer for electrons
void robustanalyzer::analyzerelectrons() {

  // vector of offline electron indices
  vector<int> noseloflelidx;
  vector<int> ptgt2etalt2p5oflelidx;
  vector<int> ptgt2etlt1p44oflelidx;
  vector<int> ptgt2etgt1p57lt2oflelidx;
  vector<int> ptgt2etgt2lt2p5oflelidx;
  /*
  vector<int> ptgt2lt5etalt2p5oflelidx;
  vector<int> ptgt5lt8etalt2p5oflelidx;
  vector<int> ptgt8lt11etalt2p5oflelidx;
  vector<int> ptgt11lt14etalt2p5oflelidx;
  vector<int> ptgt14lt18etalt2p5oflelidx;
  vector<int> ptgt18lt22etalt2p5oflelidx;
  vector<int> ptgt22lt26etalt2p5oflelidx;
  vector<int> ptgt26lt30etalt2p5oflelidx;
  vector<int> ptgt30lt50etalt2p5oflelidx;
  vector<int> ptgt50lt100etalt2p5oflelidx;
  */
  
  // vector of electron indices
  vector<int> noselelidx;
  vector<int> ptgt2etalt2p5elidx;
  vector<int> ptgt50etalt2p5elidx;
    
  // Sort the electrons based on their pT
  vector<int> sortedoflelidx((*(*n_ofle)));
  iota(begin(sortedoflelidx), end(sortedoflelidx), 0);
  sort(&sortedoflelidx[0], ofle_pt, (*(*n_ofle))); // Verified that the algorithm works fine

  vector<int> sortedelidx((*(*n_ele)));
  iota(begin(sortedelidx), end(sortedelidx), 0);
  sort(&sortedelidx[0], ele_pt, (*(*n_ele))); // Verified that the algorithm works fine
  
  if((*(*n_ofle))<0) throw "Error!!! Wrong technical event processing. Negative number of electrons in event.";
  if((*(*n_ele))<0) throw "Error!!! Wrong technical event processing. Negative number of electrons in event.";
  
  // Loop on offline electrons in the event loop
  for(unsigned int ofle_ctr=0; ofle_ctr<(*(*n_ofle)); ofle_ctr++) {
    
    // Take the sorted index only
    unsigned int ofleidx = sortedoflelidx[ofle_ctr];
    
    noseloflelidx.push_back(ofleidx);
    
    // Get the energy of this electron
    TLorentzVector ele;
    ele.SetPtEtaPhiM((*ofle_pt)->at(ofleidx), (*ofle_eta)->at(ofleidx), (*ofle_phi)->at(ofleidx), 0.0005);
    double ele_energy = ele.Energy();
    double ea = effectivearea((*ofle_eta)->at(ofleidx));
    double neutiso = ((*ofle_nthadiso)->at(ofleidx))+((*ofle_phoiso)->at(ofleidx))-(((*oflrho)->at(0))*ea);
    double reliso = ((*ofle_cghadiso)->at(ofleidx))+(neutiso>0?neutiso:0);
    
    bool ptgt2etalt2p5eldec = true;
    ptgt2etalt2p5eldec *= (*ofle_pt)->at(ofleidx)>2.0;
    ptgt2etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt2etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0117):((*ofle_sigmaietaieta)->at(ofleidx)<0.0298);
    ptgt2etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_detain)->at(ofleidx)<0.0071):((*ofle_detain)->at(ofleidx)<0.0173);
    ptgt2etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_dphiin)->at(ofleidx)<0.208):((*ofle_sigmaietaieta)->at(ofleidx)<0.234);
    ptgt2etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.28/ele_energy)+(0.0422*((*oflrho)->at(0))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.3/ele_energy)+(0.262*((*oflrho)->at(0))/ele_energy)));
    ptgt2etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?reliso<0.406+0.535/((*ofle_pt)->at(ofleidx)):reliso<0.342+0.519/((*ofle_pt)->at(ofleidx));    
    if(ptgt2etalt2p5eldec) ptgt2etalt2p5oflelidx.push_back(ofleidx);

    bool ptgt2etlt1p44eldec = true;
    ptgt2etlt1p44eldec *= (*ofle_pt)->at(ofleidx)>2.0;
    ptgt2etlt1p44eldec *= abs((*ofle_eta)->at(ofleidx))<1.44;
    ptgt2etlt1p44eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0117):((*ofle_sigmaietaieta)->at(ofleidx)<0.0298);
    ptgt2etlt1p44eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_detain)->at(ofleidx)<0.0071):((*ofle_detain)->at(ofleidx)<0.0173);
    ptgt2etlt1p44eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_dphiin)->at(ofleidx)<0.208):((*ofle_sigmaietaieta)->at(ofleidx)<0.234);
    ptgt2etlt1p44eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.28/ele_energy)+(0.0422*((*oflrho)->at(0))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.3/ele_energy)+(0.262*((*oflrho)->at(0))/ele_energy)));
    ptgt2etlt1p44eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?reliso<0.406+0.535/((*ofle_pt)->at(ofleidx)):reliso<0.342+0.519/((*ofle_pt)->at(ofleidx));    
    if(ptgt2etlt1p44eldec) ptgt2etlt1p44oflelidx.push_back(ofleidx);

    bool ptgt2etgt1p57lt2eldec = true;
    ptgt2etgt1p57lt2eldec *= (*ofle_pt)->at(ofleidx)>2.0;
    ptgt2etgt1p57lt2eldec *= ( (abs((*ofle_eta)->at(ofleidx))>1.57) && (abs((*ofle_eta)->at(ofleidx))<2.0) );
    ptgt2etgt1p57lt2eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0117):((*ofle_sigmaietaieta)->at(ofleidx)<0.0298);
    ptgt2etgt1p57lt2eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_detain)->at(ofleidx)<0.0071):((*ofle_detain)->at(ofleidx)<0.0173);
    ptgt2etgt1p57lt2eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_dphiin)->at(ofleidx)<0.208):((*ofle_sigmaietaieta)->at(ofleidx)<0.234);
    ptgt2etgt1p57lt2eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.28/ele_energy)+(0.0422*((*oflrho)->at(0))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.3/ele_energy)+(0.262*((*oflrho)->at(0))/ele_energy)));
    ptgt2etgt1p57lt2eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?reliso<0.406+0.535/((*ofle_pt)->at(ofleidx)):reliso<0.342+0.519/((*ofle_pt)->at(ofleidx));    
    if(ptgt2etgt1p57lt2eldec) ptgt2etgt1p57lt2oflelidx.push_back(ofleidx);

    bool ptgt2etgt2lt2p5eldec = true;
    ptgt2etgt2lt2p5eldec *= (*ofle_pt)->at(ofleidx)>2.0;
    ptgt2etgt2lt2p5eldec *= ( (abs((*ofle_eta)->at(ofleidx))>2.0) && (abs((*ofle_eta)->at(ofleidx))<2.5) );
    ptgt2etgt2lt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0117):((*ofle_sigmaietaieta)->at(ofleidx)<0.0298);
    ptgt2etgt2lt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_detain)->at(ofleidx)<0.0071):((*ofle_detain)->at(ofleidx)<0.0173);
    ptgt2etgt2lt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_dphiin)->at(ofleidx)<0.208):((*ofle_sigmaietaieta)->at(ofleidx)<0.234);
    ptgt2etgt2lt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.28/ele_energy)+(0.0422*((*oflrho)->at(0))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.3/ele_energy)+(0.262*((*oflrho)->at(0))/ele_energy)));
    ptgt2etgt2lt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?reliso<0.406+0.535/((*ofle_pt)->at(ofleidx)):reliso<0.342+0.519/((*ofle_pt)->at(ofleidx));    
    if(ptgt2etgt2lt2p5eldec) ptgt2etgt2lt2p5oflelidx.push_back(ofleidx);

    /* For gen matching
    bool ptgt2lt5etalt2p5eldec = true;
    ptgt2lt5etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>2.0) && ((*ofle_pt)->at(ofleidx)<5.0) );
    ptgt2lt5etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt2lt5etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt2lt5etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt2lt5etalt2p5eldec) ptgt2lt5etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt5lt8etalt2p5eldec = true;
    ptgt5lt8etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>5.0) && ((*ofle_pt)->at(ofleidx)<8.0) );
    ptgt5lt8etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt5lt8etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt5lt8etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt5lt8etalt2p5eldec) ptgt5lt8etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt8lt11etalt2p5eldec = true;
    ptgt8lt11etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>8.0) && ((*ofle_pt)->at(ofleidx)<11.0) );
    ptgt8lt11etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt8lt11etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt8lt11etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt8lt11etalt2p5eldec) ptgt8lt11etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt11lt14etalt2p5eldec = true;
    ptgt11lt14etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>11.0) && ((*ofle_pt)->at(ofleidx)<14.0) );
    ptgt11lt14etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt11lt14etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt11lt14etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt11lt14etalt2p5eldec) ptgt11lt14etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt14lt18etalt2p5eldec = true;
    ptgt14lt18etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>14.0) && ((*ofle_pt)->at(ofleidx)<18.0) );
    ptgt14lt18etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt14lt18etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt14lt18etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt14lt18etalt2p5eldec) ptgt14lt18etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt18lt22etalt2p5eldec = true;
    ptgt18lt22etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>18.0) && ((*ofle_pt)->at(ofleidx)<22.0) );
    ptgt18lt22etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt18lt22etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt18lt22etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt18lt22etalt2p5eldec) ptgt18lt22etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt22lt26etalt2p5eldec = true;
    ptgt22lt26etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>22.0) && ((*ofle_pt)->at(ofleidx)<26.0) );
    ptgt22lt26etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt22lt26etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt22lt26etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt22lt26etalt2p5eldec) ptgt22lt26etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt26lt30etalt2p5eldec = true;
    ptgt26lt30etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>26.0) && ((*ofle_pt)->at(ofleidx)<30.0) );
    ptgt26lt30etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt26lt30etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt26lt30etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt26lt30etalt2p5eldec) ptgt26lt30etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt30lt50etalt2p5eldec = true;
    ptgt30lt50etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>30.0) && ((*ofle_pt)->at(ofleidx)<50.0) );
    ptgt30lt50etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt30lt50etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt30lt50etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt30lt50etalt2p5eldec) ptgt30lt50etalt2p5oflelidx.push_back(ofleidx);
    
    bool ptgt50lt100etalt2p5eldec = true;
    ptgt50lt100etalt2p5eldec *= ( ((*ofle_pt)->at(ofleidx)>50.0) && ((*ofle_pt)->at(ofleidx)<100.0) );
    ptgt50lt100etalt2p5eldec *= abs((*ofle_eta)->at(ofleidx))<2.5;
    ptgt50lt100etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_sigmaietaieta)->at(ofleidx)<0.0126):((*ofle_sigmaietaieta)->at(ofleidx)<0.0457);
    ptgt50lt100etalt2p5eldec *= (abs((*ofle_eta)->at(ofleidx))<1.479)?((*ofle_hoe)->at(ofleidx)<(0.05+(1.16/ele_energy)+(0.0324*(*(*rho))/ele_energy)))
							:((*ofle_hoe)->at(ofleidx)<(0.05+(2.54/ele_energy)+(0.183*(*(*rho))/ele_energy)));
    if(ptgt50lt100etalt2p5eldec) ptgt50lt100etalt2p5oflelidx.push_back(ofleidx);
    */    

  }// End of loop on offline electrons in the event loop

  // Loop on electrons in the event loop
  for(unsigned int ele_ctr=0; ele_ctr<(*(*n_ele)); ele_ctr++) {
    
    // Take the sorted index only
    unsigned int elidx = sortedelidx[ele_ctr];
    
    noselelidx.push_back(elidx);
    
    // Get the energy of this electron
    TLorentzVector ele;
    ele.SetPtEtaPhiM((*ele_pt)->at(elidx), (*ele_eta)->at(elidx), (*ele_phi)->at(elidx), 0.0005);
    double ele_energy = ele.Energy();
    
    bool ptgt2etalt2p5eldec = true;
    ptgt2etalt2p5eldec *= (*ele_pt)->at(elidx)>2.0;
    ptgt2etalt2p5eldec *= abs((*ele_eta)->at(elidx))<2.5;
    if(ptgt2etalt2p5eldec) ptgt2etalt2p5elidx.push_back(elidx);
    
    bool ptgt50etalt2p5eldec = true;
    ptgt50etalt2p5eldec *= (*ele_pt)->at(elidx)>50.0;
    ptgt50etalt2p5eldec *= abs((*ele_eta)->at(elidx))<2.5;
    if(ptgt50etalt2p5eldec) ptgt50etalt2p5elidx.push_back(elidx);
    
  }// End of loop on electrons in the event loop
  
  fillhistinevent_electron("nosel_ele", noselelidx);
  if( ((*(*hlt_isomu27))==1) || ((*(*hlt_mu50))==1) ) {
    fillhistinevent_electron("mutrigsel_ele", ptgt2etalt2p5elidx);
    if( ((*(*l1_2mu))==1) || ((*(*l1_2musqos))==1) || ((*(*l1_met))==1) || ((*(*l1_1jet))==1) || ((*(*l1_2jet))==1) ) fillhistinevent_electron("mutrigsel_noneg_ele", ptgt2etalt2p5elidx);
    if((*(*l1_1eg))==1) fillhistinevent_electron("mutrigsel_1egonly_ele", ptgt2etalt2p5elidx);
    if((*(*l1_2eg))==1) fillhistinevent_electron("mutrigsel_2egmin_ele", ptgt2etalt2p5elidx);
    fillhistinevent_electron("mutrig_elptgt50_ele", ptgt50etalt2p5elidx);
  }
  if( ( ((*(*hlt_isomu27))==1) || ((*(*hlt_mu50))==1) ) && ((*(*hlt_scouting))==1) ) {
    fillhistinevent_electron("muAscouttrigsel_ele", ptgt2etalt2p5elidx);
    if( ((*(*l1_2mu))==1) || ((*(*l1_2musqos))==1) || ((*(*l1_met))==1) || ((*(*l1_1jet))==1) || ((*(*l1_2jet))==1) ) fillhistinevent_electron("muAscouttrigsel_noneg_ele", ptgt2etalt2p5elidx);
    if((*(*l1_1eg))==1) fillhistinevent_electron("muAscouttrigsel_1egonly_ele", ptgt2etalt2p5elidx);
    if((*(*l1_2eg))==1) fillhistinevent_electron("muAscouttrigsel_2egmin_ele", ptgt2etalt2p5elidx);
    fillhistinevent_electron("muAscouttrig_elptgt50_ele", ptgt50etalt2p5elidx);
  }

  bool l1_noeg = ( ((*(*l1_2mu))==1) || ((*(*l1_2musqos))==1) || ((*(*l1_met))==1) || ((*(*l1_1jet))==1) || ((*(*l1_2jet))==1) );
  fillhistinevent_oflsct_electron("nosel_oflsct_ele", noseloflelidx, noselelidx);
  fillhistinevent_oflsct_electron("ptgt2etalt2p5_oflsct_ele", ptgt2etalt2p5oflelidx, ptgt2etalt2p5elidx);
  if(l1_noeg) fillhistinevent_oflsct_electron("noeg_ptgt2etalt2p5_oflsct_ele", ptgt2etalt2p5oflelidx, ptgt2etalt2p5elidx);
  if((*(*l1_1eg))==1) fillhistinevent_oflsct_electron("eg1_ptgt2etalt2p5_oflsct_ele", ptgt2etalt2p5oflelidx, ptgt2etalt2p5elidx);
  if((*(*l1_2eg))==1) fillhistinevent_oflsct_electron("eg2_ptgt2etalt2p5_oflsct_ele", ptgt2etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt2etlt1p44_oflsct_ele", ptgt2etlt1p44oflelidx, ptgt2etalt2p5elidx);
  if(l1_noeg) fillhistinevent_oflsct_electron("noeg_ptgt2etlt1p44_oflsct_ele", ptgt2etlt1p44oflelidx, ptgt2etalt2p5elidx);
  if((*(*l1_1eg))==1) fillhistinevent_oflsct_electron("eg1_ptgt2etlt1p44_oflsct_ele", ptgt2etlt1p44oflelidx, ptgt2etalt2p5elidx);
  if((*(*l1_2eg))==1) fillhistinevent_oflsct_electron("eg2_ptgt2etlt1p44_oflsct_ele", ptgt2etlt1p44oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt2etgt1p57lt2_oflsct_ele", ptgt2etgt1p57lt2oflelidx, ptgt2etalt2p5elidx);
  if(l1_noeg) fillhistinevent_oflsct_electron("noeg_ptgt2etgt1p57lt2_oflsct_ele", ptgt2etgt1p57lt2oflelidx, ptgt2etalt2p5elidx);
  if((*(*l1_1eg))==1) fillhistinevent_oflsct_electron("eg1_ptgt2etgt1p57lt2_oflsct_ele", ptgt2etgt1p57lt2oflelidx, ptgt2etalt2p5elidx);
  if((*(*l1_2eg))==1) fillhistinevent_oflsct_electron("eg2_ptgt2etgt1p57lt2_oflsct_ele", ptgt2etgt1p57lt2oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt2etgt2lt2p5_oflsct_ele", ptgt2etgt2lt2p5oflelidx, ptgt2etalt2p5elidx);
  if(l1_noeg) fillhistinevent_oflsct_electron("noeg_ptgt2etgt2lt2p5_oflsct_ele", ptgt2etgt2lt2p5oflelidx, ptgt2etalt2p5elidx);
  if((*(*l1_1eg))==1) fillhistinevent_oflsct_electron("eg1_ptgt2etgt2lt2p5_oflsct_ele", ptgt2etgt2lt2p5oflelidx, ptgt2etalt2p5elidx);
  if((*(*l1_2eg))==1) fillhistinevent_oflsct_electron("eg2_ptgt2etgt2lt2p5_oflsct_ele", ptgt2etgt2lt2p5oflelidx, ptgt2etalt2p5elidx);
  /* For gen matching
  fillhistinevent_oflsct_electron("ptgt2lt5etalt2p5_oflsct_ele", ptgt2lt5etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt5lt8etalt2p5_oflsct_ele", ptgt5lt8etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt8lt11etalt2p5_oflsct_ele", ptgt8lt11etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt11lt14etalt2p5_oflsct_ele", ptgt11lt14etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt14lt18etalt2p5_oflsct_ele", ptgt14lt18etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt18lt22etalt2p5_oflsct_ele", ptgt18lt22etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt22lt26etalt2p5_oflsct_ele", ptgt22lt26etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt26lt30etalt2p5_oflsct_ele", ptgt26lt30etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt30lt50etalt2p5_oflsct_ele", ptgt30lt50etalt2p5oflelidx, ptgt2etalt2p5elidx);
  fillhistinevent_oflsct_electron("ptgt50lt100etalt2p5_oflsct_ele", ptgt50lt100etalt2p5oflelidx, ptgt2etalt2p5elidx);
  */
  
  // Clear all vector
  noseloflelidx.clear();
  ptgt2etalt2p5oflelidx.clear();
  ptgt2etlt1p44oflelidx.clear();
  ptgt2etgt1p57lt2oflelidx.clear();
  ptgt2etgt2lt2p5oflelidx.clear();
  /* For gen matching
  ptgt2lt5etalt2p5oflelidx.clear();
  ptgt5lt8etalt2p5oflelidx.clear();
  ptgt8lt11etalt2p5oflelidx.clear();
  ptgt11lt14etalt2p5oflelidx.clear();
  ptgt14lt18etalt2p5oflelidx.clear();
  ptgt18lt22etalt2p5oflelidx.clear();
  ptgt22lt26etalt2p5oflelidx.clear();
  ptgt26lt30etalt2p5oflelidx.clear();
  ptgt30lt50etalt2p5oflelidx.clear();
  ptgt50lt100etalt2p5oflelidx.clear();
  */
  
  noselelidx.clear();
  ptgt2etalt2p5elidx.clear();
  ptgt50etalt2p5elidx.clear();
    
}

// Analyzer for photons
void robustanalyzer::analyzerphotons() {

  // vector of photon indices
  vector<int> noselphidx;
  vector<int> ptgt2etalt2p5phidx;
  vector<int> ptgt50etalt2p5phidx;
  
  // Sort the photons based on their pT
  vector<int> sortedphidx((*(*n_pho)));
  iota(begin(sortedphidx), end(sortedphidx), 0);
  sort(&sortedphidx[0], pho_pt, (*(*n_pho))); // Verified that the algorithm works fine
  
  if((*(*n_pho))<0) throw "Error!!! Wrong technical event processing. Negative number of photons in event.";
  
  // Loop on photons in the event loop
  for(unsigned int pho_ctr=0; pho_ctr<(*(*n_pho)); pho_ctr++) {
    
    // Take the sorted index only
    unsigned int phidx = sortedphidx[pho_ctr];
    
    noselphidx.push_back(phidx);
    
    // Get the energy of this electron
    TLorentzVector pho;
    pho.SetPtEtaPhiM((*pho_pt)->at(phidx), (*pho_eta)->at(phidx), (*pho_phi)->at(phidx), 0.0005);
    double pho_energy = pho.Energy();
    
    bool ptgt2etalt2p5phdec = true;
    ptgt2etalt2p5phdec *= (*pho_pt)->at(phidx)>2.0;
    ptgt2etalt2p5phdec *= abs((*pho_eta)->at(phidx))<2.5;
    if(ptgt2etalt2p5phdec) ptgt2etalt2p5phidx.push_back(phidx);
    
    bool ptgt50etalt2p5phdec = true;
    ptgt50etalt2p5phdec *= (*pho_pt)->at(phidx)>50.0;
    ptgt50etalt2p5phdec *= abs((*pho_eta)->at(phidx))<2.5;
    if(ptgt50etalt2p5phdec) ptgt50etalt2p5phidx.push_back(phidx);

  }// End of loop on photons in the event loop
  
  fillhistinevent_photon("nosel_pho", noselphidx);
  if( ((*(*hlt_isomu27))==1) || ((*(*hlt_mu50))==1) ) {
    fillhistinevent_photon("mutrigsel_pho", ptgt2etalt2p5phidx);
    if( ((*(*l1_2mu))==1) || ((*(*l1_2musqos))==1) || ((*(*l1_met))==1) || ((*(*l1_1jet))==1) || ((*(*l1_2jet))==1) ) fillhistinevent_photon("mutrigsel_noneg_pho", ptgt2etalt2p5phidx);
    if((*(*l1_1eg))==1) fillhistinevent_photon("mutrigsel_1egonly_pho", ptgt2etalt2p5phidx);
    if((*(*l1_2eg))==1) fillhistinevent_photon("mutrigsel_2egmin_pho", ptgt2etalt2p5phidx);
    fillhistinevent_photon("mutrig_elptgt50_pho", ptgt50etalt2p5phidx);
  }
  if( ( ((*(*hlt_isomu27))==1) || ((*(*hlt_mu50))==1) ) && ((*(*hlt_scouting))==1) ) {
    fillhistinevent_photon("muAscouttrigsel_pho", ptgt2etalt2p5phidx);
    if( ((*(*l1_2mu))==1) || ((*(*l1_2musqos))==1) || ((*(*l1_met))==1) || ((*(*l1_1jet))==1) || ((*(*l1_2jet))==1) ) fillhistinevent_photon("muAscouttrigsel_noneg_pho", ptgt2etalt2p5phidx);
    if((*(*l1_1eg))==1) fillhistinevent_photon("muAscouttrigsel_1egonly_pho", ptgt2etalt2p5phidx);
    if((*(*l1_2eg))==1) fillhistinevent_photon("muAscouttrigsel_2egmin_pho", ptgt2etalt2p5phidx);
    fillhistinevent_photon("muAscouttrig_elptgt50_pho", ptgt50etalt2p5phidx);
  }
  
  // Clear all vector
  noselphidx.clear();
  ptgt2etalt2p5phidx.clear();
  ptgt50etalt2p5phidx.clear();
    
}

// Function to add a set of histograms for scouting electrons
void robustanalyzer::addhist_electron(TString selection) {

  // Event Wide
  all1dhists.push_back(new TH1F(selection+"sct_rho","rho",10000,-10,90));
  all1dhists.push_back(new TH1F(selection+"sct_isomu27","HLT_IsoMu27",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_mu50","HLT_Mu50",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_scouting","HLT_Run3PixelOnlyPFScouting",10,-5,5));

  // Scouting Electron
  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_dielM","all M(e,e)",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadsublead_dielM","M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadsublead_chargeprod","q_{e1}*q_{e2}",3,-1,1));
  all1dhists.push_back(new TH1F(selection+"sct_leadbarsubleadbar_dielM","EB M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadecsubleadec_dielM","EE M(e_{1},e_{2})",100000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sctbar_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctbar_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctbar_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldetain","#Delta#eta_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldphiin","#Delta#phi_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sctbar_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_elcharge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sctbar_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sctec_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctec_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctec_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctec_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctec_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctec_eldetain","#Delta#eta_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_eldphiin","#Delta#phi_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_elhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctec_elooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sctec_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctec_elcharge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sctec_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctec_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctec_elsmaj","smaj",1500,-0.1,1.4));
}

// Function to fill a set of histograms for scouting electrons
void robustanalyzer::fillhistinevent_electron(TString selection, vector<int> elidx) {

  if(elidx.size()==0) return;

  TH1F* rhohist = (TH1F*) outfile->Get(selection+"sct_rho");
  TH1F* isomu27 = (TH1F*) outfile->Get(selection+"sct_isomu27");
  TH1F* mu50 = (TH1F*) outfile->Get(selection+"sct_mu50");
  TH1F* scouting = (TH1F*) outfile->Get(selection+"sct_scouting");

  TH1F* elmult = (TH1F*) outfile->Get(selection+"sct_elmult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"sct_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"sct_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"sct_elphi");
  TH1F* dielM = (TH1F*) outfile->Get(selection+"sct_dielM");
  TH1F* leadsubleaddielM = (TH1F*) outfile->Get(selection+"sct_leadsublead_dielM");
  TH1F* leadsubleadqprod = (TH1F*) outfile->Get(selection+"sct_leadsublead_chargeprod");
  TH1F* leadbarsubleadbardielM = (TH1F*) outfile->Get(selection+"sct_leadbarsubleadbar_dielM");
  TH1F* leadecsubleadecdielM = (TH1F*) outfile->Get(selection+"sct_leadecsubleadec_dielM");
  
  TH1F* barelpt = (TH1F*) outfile->Get(selection+"sctbar_elpt");
  TH1F* bareleta = (TH1F*) outfile->Get(selection+"sctbar_eleta");
  TH1F* barelphi = (TH1F*) outfile->Get(selection+"sctbar_elphi");
  TH1F* barelm = (TH1F*) outfile->Get(selection+"sctbar_elm");
  TH1F* bareld0 = (TH1F*) outfile->Get(selection+"sctbar_eld0");
  TH1F* barellog10d0 = (TH1F*) outfile->Get(selection+"sctbar_ellog10d0");
  TH1F* bareldz = (TH1F*) outfile->Get(selection+"sctbar_eldz");
  TH1F* bareldetain = (TH1F*) outfile->Get(selection+"sctbar_eldetain");
  TH1F* bareldphiin = (TH1F*) outfile->Get(selection+"sctbar_eldphiin");
  TH1F* barelsigmaietaieta = (TH1F*) outfile->Get(selection+"sctbar_elsigmaietaieta");
  TH1F* barelhoe = (TH1F*) outfile->Get(selection+"sctbar_elhoe");
  TH1F* barelooemoop = (TH1F*) outfile->Get(selection+"sctbar_elooemoop");
  TH1F* barelmhits = (TH1F*) outfile->Get(selection+"sctbar_elmhits");
  TH1F* barelcharge = (TH1F*) outfile->Get(selection+"sctbar_elcharge");
  TH1F* barelecaliso = (TH1F*) outfile->Get(selection+"sctbar_elecaliso");
  TH1F* barelhcaliso = (TH1F*) outfile->Get(selection+"sctbar_elhcaliso");
  TH1F* bareltkiso = (TH1F*) outfile->Get(selection+"sctbar_eltkiso");
  TH1F* barelr9 = (TH1F*) outfile->Get(selection+"sctbar_elr9");
  TH1F* barelsmin = (TH1F*) outfile->Get(selection+"sctbar_elsmin");
  TH1F* barelsmaj = (TH1F*) outfile->Get(selection+"sctbar_elsmaj");
  
  TH1F* ecelpt = (TH1F*) outfile->Get(selection+"sctec_elpt");
  TH1F* eceleta = (TH1F*) outfile->Get(selection+"sctec_eleta");
  TH1F* ecelphi = (TH1F*) outfile->Get(selection+"sctec_elphi");
  TH1F* ecelm = (TH1F*) outfile->Get(selection+"sctec_elm");
  TH1F* eceld0 = (TH1F*) outfile->Get(selection+"sctec_eld0");
  TH1F* ecellog10d0 = (TH1F*) outfile->Get(selection+"sctec_ellog10d0");
  TH1F* eceldz = (TH1F*) outfile->Get(selection+"sctec_eldz");
  TH1F* eceldetain = (TH1F*) outfile->Get(selection+"sctec_eldetain");
  TH1F* eceldphiin = (TH1F*) outfile->Get(selection+"sctec_eldphiin");
  TH1F* ecelsigmaietaieta = (TH1F*) outfile->Get(selection+"sctec_elsigmaietaieta");
  TH1F* ecelhoe = (TH1F*) outfile->Get(selection+"sctec_elhoe");
  TH1F* ecelooemoop = (TH1F*) outfile->Get(selection+"sctec_elooemoop");
  TH1F* ecelmhits = (TH1F*) outfile->Get(selection+"sctec_elmhits");
  TH1F* ecelcharge = (TH1F*) outfile->Get(selection+"sctec_elcharge");
  TH1F* ecelecaliso = (TH1F*) outfile->Get(selection+"sctec_elecaliso");
  TH1F* ecelhcaliso = (TH1F*) outfile->Get(selection+"sctec_elhcaliso");
  TH1F* eceltkiso = (TH1F*) outfile->Get(selection+"sctec_eltkiso");
  TH1F* ecelr9 = (TH1F*) outfile->Get(selection+"sctec_elr9");
  TH1F* ecelsmin = (TH1F*) outfile->Get(selection+"sctec_elsmin");
  TH1F* ecelsmaj = (TH1F*) outfile->Get(selection+"sctec_elsmaj");

  rhohist->Fill((*rho)->at(0));
  isomu27->Fill((*(*hlt_isomu27)));
  mu50->Fill((*(*hlt_mu50)));
  scouting->Fill((*(*hlt_scouting)));
  elmult->Fill(elidx.size());

  if(elidx.size()>=2) {
    TLorentzVector leadel, subleadel;
    leadel.SetPtEtaPhiM((*ele_pt)->at(elidx[0]),(*ele_eta)->at(elidx[0]),(*ele_phi)->at(elidx[0]),0.0005);
    subleadel.SetPtEtaPhiM((*ele_pt)->at(elidx[1]),(*ele_eta)->at(elidx[1]),(*ele_phi)->at(elidx[1]),0.0005);
    leadsubleaddielM->Fill((leadel+subleadel).M());
    leadsubleadqprod->Fill((*ele_charge)->at(elidx[0])*(*ele_charge)->at(elidx[1]));
    if(abs((*ele_eta)->at(elidx[0]))<1.479 && abs((*ele_eta)->at(elidx[1]))<1.479)  {
      leadbarsubleadbardielM->Fill((leadel+subleadel).M());
    }
    if(abs((*ele_eta)->at(elidx[0]))>1.479 && abs((*ele_eta)->at(elidx[1]))>1.479){
      leadecsubleadecdielM->Fill((leadel+subleadel).M());
    }
  }

  for(unsigned int ctr=0; ctr<elidx.size(); ctr++) {
    elpt->Fill((*ele_pt)->at(elidx[ctr]));
    eleta->Fill((*ele_eta)->at(elidx[ctr]));
    elphi->Fill((*ele_phi)->at(elidx[ctr]));
    for(unsigned int ctr2=ctr+1; ctr2<elidx.size(); ctr2++) {
      TLorentzVector el, el2;
      el.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr]),(*ele_eta)->at(elidx[ctr]),(*ele_phi)->at(elidx[ctr]),0.0005);
      el2.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr2]),(*ele_eta)->at(elidx[ctr2]),(*ele_phi)->at(elidx[ctr2]),0.0005);
      dielM->Fill((el+el2).M());
    }

    if(TMath::Abs((*ele_eta)->at(elidx[ctr]))<1.479) {
      barelpt->Fill((*ele_pt)->at(elidx[ctr]));
      bareleta->Fill((*ele_eta)->at(elidx[ctr]));
      barelphi->Fill((*ele_phi)->at(elidx[ctr]));
      barelm->Fill((*ele_m)->at(elidx[ctr]));
      bareld0->Fill((*ele_d0)->at(elidx[ctr]));
      barellog10d0->Fill(TMath::Log10(TMath::Abs((*ele_d0)->at(elidx[ctr]))));
      bareldz->Fill((*ele_dz)->at(elidx[ctr]));
      bareldetain->Fill((*ele_detain)->at(elidx[ctr]));
      bareldphiin->Fill((*ele_dphiin)->at(elidx[ctr]));
      barelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]));
      barelhoe->Fill((*ele_hoe)->at(elidx[ctr]));
      barelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]));
      barelmhits->Fill((*ele_mhits)->at(elidx[ctr]));
      barelcharge->Fill((*ele_charge)->at(elidx[ctr]));
      barelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]));
      barelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]));
      bareltkiso->Fill((*ele_tkiso)->at(elidx[ctr]));
      barelr9->Fill((*ele_r9)->at(elidx[ctr]));
      barelsmin->Fill((*ele_smin)->at(elidx[ctr]));
      barelsmaj->Fill((*ele_smaj)->at(elidx[ctr]));
    }

    else {
      ecelpt->Fill((*ele_pt)->at(elidx[ctr]));
      eceleta->Fill((*ele_eta)->at(elidx[ctr]));
      ecelphi->Fill((*ele_phi)->at(elidx[ctr]));
      ecelm->Fill((*ele_m)->at(elidx[ctr]));
      eceld0->Fill((*ele_d0)->at(elidx[ctr]));
      ecellog10d0->Fill(TMath::Log10(TMath::Abs((*ele_d0)->at(elidx[ctr]))));
      eceldz->Fill((*ele_dz)->at(elidx[ctr]));
      eceldetain->Fill((*ele_detain)->at(elidx[ctr]));
      eceldphiin->Fill((*ele_dphiin)->at(elidx[ctr]));
      ecelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]));
      ecelhoe->Fill((*ele_hoe)->at(elidx[ctr]));
      ecelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]));
      ecelmhits->Fill((*ele_mhits)->at(elidx[ctr]));
      ecelcharge->Fill((*ele_charge)->at(elidx[ctr]));
      ecelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]));
      ecelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]));
      eceltkiso->Fill((*ele_tkiso)->at(elidx[ctr]));
      ecelr9->Fill((*ele_r9)->at(elidx[ctr]));
      ecelsmin->Fill((*ele_smin)->at(elidx[ctr]));
      ecelsmaj->Fill((*ele_smaj)->at(elidx[ctr]));
    }
  } // End of main electron for loop

}  

// Function to add a set of histograms for comapring offline and scouting electrons
void robustanalyzer::addhist_oflsct_electron(TString selection) {

  // Event Wide
  all1dhists.push_back(new TH1F(selection+"evt_rho","rho",10000,-10,90));
  all1dhists.push_back(new TH1F(selection+"evt_isomu27","HLT_IsoMu27",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_mu50","HLT_Mu50",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_scouting","HLT_Run3PixelOnlyPFScouting",10,-5,5));

  // Offline Electron  
  all1dhists.push_back(new TH1F(selection+"ofl_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"ofl_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"ofl_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"ofl_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"ofl_dielM","all M(e,e)",100000,-10,990));
  
  all1dhists.push_back(new TH1F(selection+"ofleb_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"ofleb_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"ofleb_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"ofleb_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"ofleb_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"ofleb_lead_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"ofleb_lead_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"ofleb_lead_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"ofleb_lead_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"ofleb_lead_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"ofleb_sublead_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"ofleb_sublead_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"ofleb_sublead_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"ofleb_sublead_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"ofleb_sublead_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"oflee_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"oflee_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"oflee_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"oflee_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"oflee_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"oflee_lead_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"oflee_lead_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"oflee_lead_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"oflee_lead_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"oflee_lead_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"oflee_sublead_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"oflee_sublead_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"oflee_sublead_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"oflee_sublead_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"oflee_sublead_elhoe","H/E",100000,0,1));

  // Scouting Electron  
  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_dielM","all M(e,e)",100000,-10,990));
  
  all1dhists.push_back(new TH1F(selection+"scteb_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"scteb_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"scteb_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"scteb_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"scteb_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"sctee_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctee_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctee_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctee_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctee_elhoe","H/E",100000,0,1));

  // Angular matching condition
  all1dhists.push_back(new TH1F(selection+"sct_ofl_eb_deta","#Delta#eta(ofl, sct)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_eb_dphi","#Delta#phi(ofl, sct)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_eb_dr","#DeltaR(ofl, sct)",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_ee_deta","#Delta#eta(ofl, sct)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_ee_dphi","#Delta#phi(ofl, sct)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_ee_dr","#DeltaR(ofl, sct)",10000,0,10));

  // After angular matching
  all1dhists.push_back(new TH1F(selection+"sct_ofl_aftmch_eb_deta","#Delta#eta(ofl, sct)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_aftmch_eb_dphi","#Delta#phi(ofl, sct)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_aftmch_eb_dr","#DeltaR(ofl, sct)",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_aftmch_ee_deta","#Delta#eta(ofl, sct)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_aftmch_ee_dphi","#Delta#phi(ofl, sct)",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ofl_aftmch_ee_dr","#DeltaR(ofl, sct)",10000,0,10));

  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_lead_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_lead_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_lead_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_lead_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_lead_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_sublead_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_sublead_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_sublead_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_sublead_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"ofleb_aftmch_sublead_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_lead_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_lead_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_lead_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_lead_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_lead_elhoe","H/E",100000,0,1));

  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_sublead_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_sublead_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_sublead_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_sublead_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"oflee_aftmch_sublead_elhoe","H/E",100000,0,1));
}

// Function to fill a set of histograms for scouting electrons
void robustanalyzer::fillhistinevent_oflsct_electron(TString selection, vector<int> oflelidx, vector<int> sctelidx) {

  TH1F* rhohist = (TH1F*) outfile->Get(selection+"evt_rho");
  TH1F* isomu27 = (TH1F*) outfile->Get(selection+"evt_isomu27");
  TH1F* mu50 = (TH1F*) outfile->Get(selection+"evt_mu50");
  TH1F* scouting = (TH1F*) outfile->Get(selection+"evt_scouting");

  TH1F* oflelmult = (TH1F*) outfile->Get(selection+"ofl_elmult");
  TH1F* oflelpt = (TH1F*) outfile->Get(selection+"ofl_elpt");
  TH1F* ofleleta = (TH1F*) outfile->Get(selection+"ofl_eleta");
  TH1F* oflelphi = (TH1F*) outfile->Get(selection+"ofl_elphi");
  TH1F* dioflelM = (TH1F*) outfile->Get(selection+"ofl_dielM");
  
  TH1F* eboflelpt = (TH1F*) outfile->Get(selection+"ofleb_elpt");
  TH1F* ebofleleta = (TH1F*) outfile->Get(selection+"ofleb_eleta");
  TH1F* eboflelphi = (TH1F*) outfile->Get(selection+"ofleb_elphi");
  TH1F* eboflelsigmaietaieta = (TH1F*) outfile->Get(selection+"ofleb_elsigmaietaieta");
  TH1F* eboflelhoe = (TH1F*) outfile->Get(selection+"ofleb_elhoe");
  
  TH1F* eboflleadelpt = (TH1F*) outfile->Get(selection+"ofleb_lead_elpt");
  TH1F* eboflleadeleta = (TH1F*) outfile->Get(selection+"ofleb_lead_eleta");
  TH1F* eboflleadelphi = (TH1F*) outfile->Get(selection+"ofleb_lead_elphi");
  TH1F* eboflleadelsigmaietaieta = (TH1F*) outfile->Get(selection+"ofleb_lead_elsigmaietaieta");
  TH1F* eboflleadelhoe = (TH1F*) outfile->Get(selection+"ofleb_lead_elhoe");
  
  TH1F* eboflsubleadelpt = (TH1F*) outfile->Get(selection+"ofleb_sublead_elpt");
  TH1F* eboflsubleadeleta = (TH1F*) outfile->Get(selection+"ofleb_sublead_eleta");
  TH1F* eboflsubleadelphi = (TH1F*) outfile->Get(selection+"ofleb_sublead_elphi");
  TH1F* eboflsubleadelsigmaietaieta = (TH1F*) outfile->Get(selection+"ofleb_sublead_elsigmaietaieta");
  TH1F* eboflsubleadelhoe = (TH1F*) outfile->Get(selection+"ofleb_sublead_elhoe");
  
  TH1F* eeoflelpt = (TH1F*) outfile->Get(selection+"oflee_elpt");
  TH1F* eeofleleta = (TH1F*) outfile->Get(selection+"oflee_eleta");
  TH1F* eeoflelphi = (TH1F*) outfile->Get(selection+"oflee_elphi");
  TH1F* eeoflelsigmaietaieta = (TH1F*) outfile->Get(selection+"oflee_elsigmaietaieta");
  TH1F* eeoflelhoe = (TH1F*) outfile->Get(selection+"oflee_elhoe");

  TH1F* eeoflleadelpt = (TH1F*) outfile->Get(selection+"oflee_lead_elpt");
  TH1F* eeoflleadeleta = (TH1F*) outfile->Get(selection+"oflee_lead_eleta");
  TH1F* eeoflleadelphi = (TH1F*) outfile->Get(selection+"oflee_lead_elphi");
  TH1F* eeoflleadelsigmaietaieta = (TH1F*) outfile->Get(selection+"oflee_lead_elsigmaietaieta");
  TH1F* eeoflleadelhoe = (TH1F*) outfile->Get(selection+"oflee_lead_elhoe");

  TH1F* eeoflsubleadelpt = (TH1F*) outfile->Get(selection+"oflee_sublead_elpt");
  TH1F* eeoflsubleadeleta = (TH1F*) outfile->Get(selection+"oflee_sublead_eleta");
  TH1F* eeoflsubleadelphi = (TH1F*) outfile->Get(selection+"oflee_sublead_elphi");
  TH1F* eeoflsubleadelsigmaietaieta = (TH1F*) outfile->Get(selection+"oflee_sublead_elsigmaietaieta");
  TH1F* eeoflsubleadelhoe = (TH1F*) outfile->Get(selection+"oflee_sublead_elhoe");

  TH1F* sctelmult = (TH1F*) outfile->Get(selection+"sct_elmult");
  TH1F* sctelpt = (TH1F*) outfile->Get(selection+"sct_elpt");
  TH1F* scteleta = (TH1F*) outfile->Get(selection+"sct_eleta");
  TH1F* sctelphi = (TH1F*) outfile->Get(selection+"sct_elphi");
  TH1F* disctelM = (TH1F*) outfile->Get(selection+"sct_dielM");
  
  TH1F* ebsctelpt = (TH1F*) outfile->Get(selection+"scteb_elpt");
  TH1F* ebscteleta = (TH1F*) outfile->Get(selection+"scteb_eleta");
  TH1F* ebsctelphi = (TH1F*) outfile->Get(selection+"scteb_elphi");
  TH1F* ebsctelsigmaietaieta = (TH1F*) outfile->Get(selection+"scteb_elsigmaietaieta");
  TH1F* ebsctelhoe = (TH1F*) outfile->Get(selection+"scteb_elhoe");
  
  TH1F* eesctelpt = (TH1F*) outfile->Get(selection+"sctee_elpt");
  TH1F* eescteleta = (TH1F*) outfile->Get(selection+"sctee_eleta");
  TH1F* eesctelphi = (TH1F*) outfile->Get(selection+"sctee_elphi");
  TH1F* eesctelsigmaietaieta = (TH1F*) outfile->Get(selection+"sctee_elsigmaietaieta");
  TH1F* eesctelhoe = (TH1F*) outfile->Get(selection+"sctee_elhoe");

  TH1F* ebdeta = (TH1F*) outfile->Get(selection+"sct_ofl_eb_deta");
  TH1F* ebdphi = (TH1F*) outfile->Get(selection+"sct_ofl_eb_dphi");
  TH1F* ebdr = (TH1F*) outfile->Get(selection+"sct_ofl_eb_dr");
  TH1F* eedeta = (TH1F*) outfile->Get(selection+"sct_ofl_ee_deta");
  TH1F* eedphi = (TH1F*) outfile->Get(selection+"sct_ofl_ee_dphi");
  TH1F* eedr = (TH1F*) outfile->Get(selection+"sct_ofl_ee_dr");

  TH1F* amebdeta = (TH1F*) outfile->Get(selection+"sct_ofl_aftmch_eb_deta");
  TH1F* amebdphi = (TH1F*) outfile->Get(selection+"sct_ofl_aftmch_eb_dphi");
  TH1F* amebdr = (TH1F*) outfile->Get(selection+"sct_ofl_aftmch_eb_dr");
  TH1F* ameedeta = (TH1F*) outfile->Get(selection+"sct_ofl_aftmch_ee_deta");
  TH1F* ameedphi = (TH1F*) outfile->Get(selection+"sct_ofl_aftmch_ee_dphi");
  TH1F* ameedr = (TH1F*) outfile->Get(selection+"sct_ofl_aftmch_ee_dr");

  TH1F* ameboflelpt = (TH1F*) outfile->Get(selection+"ofleb_aftmch_elpt");
  TH1F* amebofleleta = (TH1F*) outfile->Get(selection+"ofleb_aftmch_eleta");
  TH1F* ameboflelphi = (TH1F*) outfile->Get(selection+"ofleb_aftmch_elphi");
  TH1F* ameboflelsigmaietaieta = (TH1F*) outfile->Get(selection+"ofleb_aftmch_elsigmaietaieta");
  TH1F* ameboflelhoe = (TH1F*) outfile->Get(selection+"ofleb_aftmch_elhoe");

  TH1F* amebleadoflelpt = (TH1F*) outfile->Get(selection+"ofleb_aftmch_lead_elpt");
  TH1F* amebleadofleleta = (TH1F*) outfile->Get(selection+"ofleb_aftmch_lead_eleta");
  TH1F* amebleadoflelphi = (TH1F*) outfile->Get(selection+"ofleb_aftmch_lead_elphi");
  TH1F* amebleadoflelsigmaietaieta = (TH1F*) outfile->Get(selection+"ofleb_aftmch_lead_elsigmaietaieta");
  TH1F* amebleadoflelhoe = (TH1F*) outfile->Get(selection+"ofleb_aftmch_lead_elhoe");
  
  TH1F* amebsubleadoflelpt = (TH1F*) outfile->Get(selection+"ofleb_aftmch_sublead_elpt");
  TH1F* amebsubleadofleleta = (TH1F*) outfile->Get(selection+"ofleb_aftmch_sublead_eleta");
  TH1F* amebsubleadoflelphi = (TH1F*) outfile->Get(selection+"ofleb_aftmch_sublead_elphi");
  TH1F* amebsubleadoflelsigmaietaieta = (TH1F*) outfile->Get(selection+"ofleb_aftmch_sublead_elsigmaietaieta");
  TH1F* amebsubleadoflelhoe = (TH1F*) outfile->Get(selection+"ofleb_aftmch_sublead_elhoe");
  
  TH1F* ameeoflelpt = (TH1F*) outfile->Get(selection+"oflee_aftmch_elpt");
  TH1F* ameeofleleta = (TH1F*) outfile->Get(selection+"oflee_aftmch_eleta");
  TH1F* ameeoflelphi = (TH1F*) outfile->Get(selection+"oflee_aftmch_elphi");
  TH1F* ameeoflelsigmaietaieta = (TH1F*) outfile->Get(selection+"oflee_aftmch_elsigmaietaieta");
  TH1F* ameeoflelhoe = (TH1F*) outfile->Get(selection+"oflee_aftmch_elhoe");
  
  TH1F* ameeleadoflelpt = (TH1F*) outfile->Get(selection+"oflee_aftmch_lead_elpt");
  TH1F* ameeleadofleleta = (TH1F*) outfile->Get(selection+"oflee_aftmch_lead_eleta");
  TH1F* ameeleadoflelphi = (TH1F*) outfile->Get(selection+"oflee_aftmch_lead_elphi");
  TH1F* ameeleadoflelsigmaietaieta = (TH1F*) outfile->Get(selection+"oflee_aftmch_lead_elsigmaietaieta");
  TH1F* ameeleadoflelhoe = (TH1F*) outfile->Get(selection+"oflee_aftmch_lead_elhoe");
  
  TH1F* ameesubleadoflelpt = (TH1F*) outfile->Get(selection+"oflee_aftmch_sublead_elpt");
  TH1F* ameesubleadofleleta = (TH1F*) outfile->Get(selection+"oflee_aftmch_sublead_eleta");
  TH1F* ameesubleadoflelphi = (TH1F*) outfile->Get(selection+"oflee_aftmch_sublead_elphi");
  TH1F* ameesubleadoflelsigmaietaieta = (TH1F*) outfile->Get(selection+"oflee_aftmch_sublead_elsigmaietaieta");
  TH1F* ameesubleadoflelhoe = (TH1F*) outfile->Get(selection+"oflee_aftmch_sublead_elhoe");
 
  if(oflelidx.size()>=1) {
    for(unsigned int ctr=0; ctr<oflelidx.size(); ctr++) {
      int oflidx = oflelidx[ctr];
      oflelpt->Fill((*ofle_pt)->at(oflidx));
      ofleleta->Fill((*ofle_eta)->at(oflidx));
      oflelphi->Fill((*ofle_phi)->at(oflidx));
      for(unsigned int ctr2=ctr+1; ctr2<oflelidx.size(); ctr2++) {
	int oflidx2 = oflelidx[ctr2];
	TLorentzVector el, el2;
	el.SetPtEtaPhiM((*ofle_pt)->at(oflidx),(*ofle_eta)->at(oflidx),(*ofle_phi)->at(oflidx),0.0005);
	el2.SetPtEtaPhiM((*ofle_pt)->at(oflidx2),(*ofle_eta)->at(oflidx2),(*ofle_phi)->at(oflidx2),0.0005);
	dioflelM->Fill((el+el2).M());
      }
      
      if(TMath::Abs((*ofle_eta)->at(oflidx))<1.479) {
	eboflelpt->Fill((*ofle_pt)->at(oflidx));
	ebofleleta->Fill((*ofle_eta)->at(oflidx));
	eboflelphi->Fill((*ofle_phi)->at(oflidx));
	eboflelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(oflidx));
	eboflelhoe->Fill((*ofle_hoe)->at(oflidx));
      }
      
      else {
	eeoflelpt->Fill((*ofle_pt)->at(oflidx));
	eeofleleta->Fill((*ofle_eta)->at(oflidx));
	eeoflelphi->Fill((*ofle_phi)->at(oflidx));
	eeoflelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(oflidx));
	eeoflelhoe->Fill((*ofle_hoe)->at(oflidx));
      }
    } // End of offline electron for loop

    if(TMath::Abs((*ofle_eta)->at(oflelidx[0]))<1.479) {
      eboflleadelpt->Fill((*ofle_pt)->at(oflelidx[0]));
      eboflleadeleta->Fill((*ofle_eta)->at(oflelidx[0]));
      eboflleadelphi->Fill((*ofle_phi)->at(oflelidx[0]));
      eboflleadelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(oflelidx[0]));
      eboflleadelhoe->Fill((*ofle_hoe)->at(oflelidx[0]));
    }
    else {
      eeoflleadelpt->Fill((*ofle_pt)->at(oflelidx[0]));
      eeoflleadeleta->Fill((*ofle_eta)->at(oflelidx[0]));
      eeoflleadelphi->Fill((*ofle_phi)->at(oflelidx[0]));
      eeoflleadelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(oflelidx[0]));
      eeoflleadelhoe->Fill((*ofle_hoe)->at(oflelidx[0]));
    }
  }

  if(oflelidx.size()>=2) {
    if(TMath::Abs((*ofle_eta)->at(oflelidx[1]))<1.479) {
      eboflsubleadelpt->Fill((*ofle_pt)->at(oflelidx[1]));
      eboflsubleadeleta->Fill((*ofle_eta)->at(oflelidx[1]));
      eboflsubleadelphi->Fill((*ofle_phi)->at(oflelidx[1]));
      eboflsubleadelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(oflelidx[1]));
      eboflsubleadelhoe->Fill((*ofle_hoe)->at(oflelidx[1]));
    }
    else {
      eeoflsubleadelpt->Fill((*ofle_pt)->at(oflelidx[1]));
      eeoflsubleadeleta->Fill((*ofle_eta)->at(oflelidx[1]));
      eeoflsubleadelphi->Fill((*ofle_phi)->at(oflelidx[1]));
      eeoflsubleadelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(oflelidx[1]));
      eeoflsubleadelhoe->Fill((*ofle_hoe)->at(oflelidx[1]));
    }
  }

  if(sctelidx.size()>=1) {
    for(unsigned int ctr=0; ctr<sctelidx.size(); ctr++) {
      int sctidx = sctelidx[ctr];
      sctelpt->Fill((*ele_pt)->at(sctidx));
      scteleta->Fill((*ele_eta)->at(sctidx));
      sctelphi->Fill((*ele_phi)->at(sctidx));
      for(unsigned int ctr2=ctr+1; ctr2<sctelidx.size(); ctr2++) {
	int sctidx2 = sctelidx[ctr2];
	TLorentzVector el, el2;
	el.SetPtEtaPhiM((*ele_pt)->at(sctidx),(*ele_eta)->at(sctidx),(*ele_phi)->at(sctidx),0.0005);
	el2.SetPtEtaPhiM((*ele_pt)->at(sctidx2),(*ele_eta)->at(sctidx2),(*ele_phi)->at(sctidx2),0.0005);
	disctelM->Fill((el+el2).M());
      }
      
      if(TMath::Abs((*ele_eta)->at(sctidx))<1.479) {
	ebsctelpt->Fill((*ele_pt)->at(sctidx));
	ebscteleta->Fill((*ele_eta)->at(sctidx));
	ebsctelphi->Fill((*ele_phi)->at(sctidx));
	ebsctelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(sctidx));
	ebsctelhoe->Fill((*ele_hoe)->at(sctidx));
      }
      
      else {
	eesctelpt->Fill((*ele_pt)->at(sctidx));
	eescteleta->Fill((*ele_eta)->at(sctidx));
	eesctelphi->Fill((*ele_phi)->at(sctidx));
	eesctelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(sctidx));
	eesctelhoe->Fill((*ele_hoe)->at(sctidx));
      }
    } // End of scouting electron for loop
  }

  if(oflelidx.size()>=1 && sctelidx.size()>=1) {
    vector<int> amoflidx;
    for(unsigned int oflctr=0; oflctr<oflelidx.size(); oflctr++) {
      bool pass_matching = false;
      
      int oflidx = oflelidx[oflctr];
      TLorentzVector ofl_el;
      ofl_el.SetPtEtaPhiM((*ofle_pt)->at(oflidx),(*ofle_eta)->at(oflidx),(*ofle_phi)->at(oflidx),0.0005);
      
      for(unsigned int sctctr=0; sctctr<sctelidx.size(); sctctr++) {
	int sctidx = sctelidx[sctctr];
	TLorentzVector sct_el;
	sct_el.SetPtEtaPhiM((*ele_pt)->at(sctidx),(*ele_eta)->at(sctidx),(*ele_phi)->at(sctidx),0.0005);
	
	pass_matching = offline_scouting_angmch((*ofle_eta)->at(oflidx), ofl_el.DeltaR(sct_el));
	
	if(TMath::Abs((*ofle_eta)->at(oflidx))<1.479) {
	  ebdeta->Fill((*ofle_eta)->at(oflidx)-(*ele_eta)->at(sctidx));
	  ebdphi->Fill(ofl_el.DeltaPhi(sct_el));
	  ebdr->Fill(ofl_el.DeltaR(sct_el));
	  if(! pass_matching) {
	    amebdeta->Fill((*ofle_eta)->at(oflidx)-(*ele_eta)->at(sctidx));
	    amebdphi->Fill(ofl_el.DeltaPhi(sct_el));
	    amebdr->Fill(ofl_el.DeltaR(sct_el));
	  }
	}
	else {
	  eedeta->Fill((*ofle_eta)->at(oflidx)-(*ele_eta)->at(sctidx));
	  eedphi->Fill(ofl_el.DeltaPhi(sct_el));
	  eedr->Fill(ofl_el.DeltaR(sct_el));
	  if(! pass_matching) {
	    ameedeta->Fill((*ofle_eta)->at(oflidx)-(*ele_eta)->at(sctidx));
	    ameedphi->Fill(ofl_el.DeltaPhi(sct_el));
	    ameedr->Fill(ofl_el.DeltaR(sct_el));
	  }
	}
	if(pass_matching) break;
      } // End of socuting electron loop
      if(pass_matching) amoflidx.push_back(oflidx);
    } // End of offline electron loop

    // Fill the offline electrons with matching angular scouting electrons
    for(unsigned int amctr=0; amctr<amoflidx.size(); amctr++) {
      int amidx = amoflidx[amctr];
      if(amctr>0 && ((*ofle_pt)->at(amidx)>(*ofle_pt)->at(amoflidx[amctr-1]))) throw "Error!!! Pt ordering violated by offline electrons after scout matching.";
      if(TMath::Abs((*ofle_eta)->at(amidx))<1.479) {
	ameboflelpt->Fill((*ofle_pt)->at(amidx));
	amebofleleta->Fill((*ofle_eta)->at(amidx));
	ameboflelphi->Fill((*ofle_phi)->at(amidx));
	ameboflelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(amidx));
	ameboflelhoe->Fill((*ofle_hoe)->at(amidx));
      }
      else {
	ameeoflelpt->Fill((*ofle_pt)->at(amidx));
	ameeofleleta->Fill((*ofle_eta)->at(amidx));
	ameeoflelphi->Fill((*ofle_phi)->at(amidx));
	ameeoflelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(amidx));
	ameeoflelhoe->Fill((*ofle_hoe)->at(amidx));
      }
    } // End of loop on atfer matching offline electrons
    if(amoflidx.size()>=1) {
      if(TMath::Abs((*ofle_eta)->at(amoflidx[0]))<1.479) {
	amebleadoflelpt->Fill((*ofle_pt)->at(amoflidx[0]));
	amebleadofleleta->Fill((*ofle_eta)->at(amoflidx[0]));
	amebleadoflelphi->Fill((*ofle_phi)->at(amoflidx[0]));
	amebleadoflelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(amoflidx[0]));
	amebleadoflelhoe->Fill((*ofle_hoe)->at(amoflidx[0]));
      }
      else {
	ameeleadoflelpt->Fill((*ofle_pt)->at(amoflidx[0]));
	ameeleadofleleta->Fill((*ofle_eta)->at(amoflidx[0]));
	ameeleadoflelphi->Fill((*ofle_phi)->at(amoflidx[0]));
	ameeleadoflelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(amoflidx[0]));
	ameeleadoflelhoe->Fill((*ofle_hoe)->at(amoflidx[0]));
      }
    }
    if(amoflidx.size()>=2) {
      if(TMath::Abs((*ofle_eta)->at(amoflidx[1]))<1.479) {
	amebsubleadoflelpt->Fill((*ofle_pt)->at(amoflidx[1]));
	amebsubleadofleleta->Fill((*ofle_eta)->at(amoflidx[1]));
	amebsubleadoflelphi->Fill((*ofle_phi)->at(amoflidx[1]));
	amebsubleadoflelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(amoflidx[1]));
	amebsubleadoflelhoe->Fill((*ofle_hoe)->at(amoflidx[1]));
      }
      else {
	ameesubleadoflelpt->Fill((*ofle_pt)->at(amoflidx[1]));
	ameesubleadofleleta->Fill((*ofle_eta)->at(amoflidx[1]));
	ameesubleadoflelphi->Fill((*ofle_phi)->at(amoflidx[1]));
	ameesubleadoflelsigmaietaieta->Fill((*ofle_sigmaietaieta)->at(amoflidx[1]));
	ameesubleadoflelhoe->Fill((*ofle_hoe)->at(amoflidx[1]));
      }
    }
    amoflidx.clear();
  } // End of if for one offline and one socuting electron
}

// Function to add a set of histograms for scouting photons
void robustanalyzer::addhist_photon(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_rho","rho",10000,-10,90));
  all1dhists.push_back(new TH1F(selection+"sct_isomu27","HLT_IsoMu27",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_mu50","HLT_Mu50",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_scouting","HLT_Run3PixelOnlyPFScouting",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_phmult","N #gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_phpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_pheta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_phphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_diphM","all M(#gamma,#gamma)",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadsublead_diphM","M(#gamma_{1},#gamma_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadbarsubleadbar_diphM","EB M(#gamma_{1},#gamma_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadecsubleadec_diphM","EE M(#gamma_{1},#gamma_{2})",100000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sctbar_phpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctbar_pheta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_phphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctbar_phm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctbar_phsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_phhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_phooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sctbar_phecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_phhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_phr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctbar_phsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctbar_phsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sctec_phpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctec_pheta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_phphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctec_phm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctec_phsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_phhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctec_phecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_phhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_phr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctec_phsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctec_phsmaj","smaj",1500,-0.1,1.4));
}

// Function to fill a set of histograms for scouting photons
void robustanalyzer::fillhistinevent_photon(TString selection, vector<int> phidx) {

  if(phidx.size()==0) return;

  TH1F* rhohist = (TH1F*) outfile->Get(selection+"sct_rho");
  TH1F* isomu27 = (TH1F*) outfile->Get(selection+"sct_isomu27");
  TH1F* mu50 = (TH1F*) outfile->Get(selection+"sct_mu50");
  TH1F* scouting = (TH1F*) outfile->Get(selection+"sct_scouting");
  TH1F* phmult = (TH1F*) outfile->Get(selection+"sct_phmult");
  TH1F* phpt = (TH1F*) outfile->Get(selection+"sct_phpt");
  TH1F* pheta = (TH1F*) outfile->Get(selection+"sct_pheta");
  TH1F* phphi = (TH1F*) outfile->Get(selection+"sct_phphi");
  TH1F* diphM = (TH1F*) outfile->Get(selection+"sct_diphM");
  TH1F* leadsubleaddiphM = (TH1F*) outfile->Get(selection+"sct_leadsublead_diphM");
  TH1F* leadbarsubleadbardiphM = (TH1F*) outfile->Get(selection+"sct_leadbarsubleadbar_diphM");
  TH1F* leadecsubleadecdiphM = (TH1F*) outfile->Get(selection+"sct_leadecsubleadec_diphM");
  
  TH1F* barphpt = (TH1F*) outfile->Get(selection+"sctbar_phpt");
  TH1F* barpheta = (TH1F*) outfile->Get(selection+"sctbar_pheta");
  TH1F* barphphi = (TH1F*) outfile->Get(selection+"sctbar_phphi");
  TH1F* barphm = (TH1F*) outfile->Get(selection+"sctbar_phm");
  TH1F* barphsigmaietaieta = (TH1F*) outfile->Get(selection+"sctbar_phsigmaietaieta");
  TH1F* barphhoe = (TH1F*) outfile->Get(selection+"sctbar_phhoe");
  TH1F* barphecaliso = (TH1F*) outfile->Get(selection+"sctbar_phecaliso");
  TH1F* barphhcaliso = (TH1F*) outfile->Get(selection+"sctbar_phhcaliso");
  TH1F* barphr9 = (TH1F*) outfile->Get(selection+"sctbar_phr9");
  TH1F* barphsmin = (TH1F*) outfile->Get(selection+"sctbar_phsmin");
  TH1F* barphsmaj = (TH1F*) outfile->Get(selection+"sctbar_phsmaj");
  
  TH1F* ecphpt = (TH1F*) outfile->Get(selection+"sctec_phpt");
  TH1F* ecpheta = (TH1F*) outfile->Get(selection+"sctec_pheta");
  TH1F* ecphphi = (TH1F*) outfile->Get(selection+"sctec_phphi");
  TH1F* ecphm = (TH1F*) outfile->Get(selection+"sctec_phm");
  TH1F* ecphsigmaietaieta = (TH1F*) outfile->Get(selection+"sctec_phsigmaietaieta");
  TH1F* ecphhoe = (TH1F*) outfile->Get(selection+"sctec_phhoe");
  TH1F* ecphecaliso = (TH1F*) outfile->Get(selection+"sctec_phecaliso");
  TH1F* ecphhcaliso = (TH1F*) outfile->Get(selection+"sctec_phhcaliso");
  TH1F* ecphr9 = (TH1F*) outfile->Get(selection+"sctec_phr9");
  TH1F* ecphsmin = (TH1F*) outfile->Get(selection+"sctec_phsmin");
  TH1F* ecphsmaj = (TH1F*) outfile->Get(selection+"sctec_phsmaj");

  rhohist->Fill((*rho)->at(0));
  isomu27->Fill((*(*hlt_isomu27)));
  mu50->Fill((*(*hlt_mu50)));
  scouting->Fill((*(*hlt_scouting)));
  phmult->Fill(phidx.size());

  if(phidx.size()>=2) {
    TLorentzVector leadph, subleadph;
    leadph.SetPtEtaPhiM((*pho_pt)->at(phidx[0]),(*pho_eta)->at(phidx[0]),(*pho_phi)->at(phidx[0]),0.0005);
    subleadph.SetPtEtaPhiM((*pho_pt)->at(phidx[1]),(*pho_eta)->at(phidx[1]),(*pho_phi)->at(phidx[1]),0.0005);
    leadsubleaddiphM->Fill((leadph+subleadph).M());
    if(abs((*pho_eta)->at(phidx[0]))<1.479 && abs((*pho_eta)->at(phidx[1]))<1.479)  {
      leadbarsubleadbardiphM->Fill((leadph+subleadph).M());
    }
    if(abs((*pho_eta)->at(phidx[0]))>1.479 && abs((*pho_eta)->at(phidx[1]))>1.479){
      leadecsubleadecdiphM->Fill((leadph+subleadph).M());
    }
  }

  for(unsigned int ctr=0; ctr<phidx.size(); ctr++) {
    phpt->Fill((*pho_pt)->at(phidx[ctr]));
    pheta->Fill((*pho_eta)->at(phidx[ctr]));
    phphi->Fill((*pho_phi)->at(phidx[ctr]));
    for(unsigned int ctr2=ctr+1; ctr2<phidx.size(); ctr2++) {
      TLorentzVector ph, ph2;
      ph.SetPtEtaPhiM((*pho_pt)->at(phidx[ctr]),(*pho_eta)->at(phidx[ctr]),(*pho_phi)->at(phidx[ctr]),0.0005);
      ph2.SetPtEtaPhiM((*pho_pt)->at(phidx[ctr2]),(*pho_eta)->at(phidx[ctr2]),(*pho_phi)->at(phidx[ctr2]),0.0005);
      diphM->Fill((ph+ph2).M());
    }

    if(TMath::Abs((*pho_eta)->at(phidx[ctr]))<1.479) {
      barphpt->Fill((*pho_pt)->at(phidx[ctr]));
      barpheta->Fill((*pho_eta)->at(phidx[ctr]));
      barphphi->Fill((*pho_phi)->at(phidx[ctr]));
      barphm->Fill((*pho_m)->at(phidx[ctr]));
      barphsigmaietaieta->Fill((*pho_sigmaietaieta)->at(phidx[ctr]));
      barphhoe->Fill((*pho_hoe)->at(phidx[ctr]));
      barphecaliso->Fill((*pho_ecaliso)->at(phidx[ctr]));
      barphhcaliso->Fill((*pho_hcaliso)->at(phidx[ctr]));
      barphr9->Fill((*pho_r9)->at(phidx[ctr]));
      barphsmin->Fill((*pho_smin)->at(phidx[ctr]));
      barphsmaj->Fill((*pho_smaj)->at(phidx[ctr]));
    }

    else {
      ecphpt->Fill((*pho_pt)->at(phidx[ctr]));
      ecpheta->Fill((*pho_eta)->at(phidx[ctr]));
      ecphphi->Fill((*pho_phi)->at(phidx[ctr]));
      ecphm->Fill((*pho_m)->at(phidx[ctr]));
      ecphsigmaietaieta->Fill((*pho_sigmaietaieta)->at(phidx[ctr]));
      ecphhoe->Fill((*pho_hoe)->at(phidx[ctr]));
      ecphecaliso->Fill((*pho_ecaliso)->at(phidx[ctr]));
      ecphhcaliso->Fill((*pho_hcaliso)->at(phidx[ctr]));
      ecphr9->Fill((*pho_r9)->at(phidx[ctr]));
      ecphsmin->Fill((*pho_smin)->at(phidx[ctr]));
      ecphsmaj->Fill((*pho_smaj)->at(phidx[ctr]));
    }
  } // End of main photon for loop

}  

// Function to perform angular matching between offline and scouting
bool robustanalyzer::offline_scouting_angmch(double eta, double dr) {

  bool eval = false;

  if(abs(eta)<1.479) {
    if(dr<0.15) eval = true;
  }
  else {
    if(dr<0.075) eval = true;
  }
  
  return eval;
}

// Function to sort the indices based on a factor (Usually pT)
void robustanalyzer::sort(int* idx, TTreeReaderValue<std::vector<float>> *factor, int n) {
  for(unsigned int i=0; i<n; i++) {
    for(unsigned int j=i+1; j<n; j++) {
      if((*factor)->at((*(idx+j)))>(*factor)->at((*(idx+i)))) { // Sort in decreasing value of factor
	double temp = *(idx+i);
	*(idx+i) = *(idx+j);
	*(idx+j) = temp;
      }
    }
  }
}

double robustanalyzer::effectivearea(double eta) {

  double ea = -1.0;
  double abseta = abs(eta);

  if(abseta < 1.0) {
    ea = 0.1243;
  }
  else if(abseta >= 1.0 && abseta < 1.479) {
    ea = 0.1458;
  }
  else if(abseta >= 1.479 && abseta < 2.0) {
    ea = 0.0992;
  }
  else if(abseta >= 2.0 && abseta < 2.2) {
    ea = 0.0794;
  }
  else if(abseta >= 2.2 && abseta < 2.3) {
    ea = 0.0762;
  }
  else if(abseta >= 2.3 && abseta < 2.4) {
    ea = 0.0766;
  }
  else if(abseta >= 2.4 && abseta < 2.5) {
    ea = 0.1003;
  }
  else if(abseta >= 2.5 && abseta < 2.6) {
    ea = 0.2322;
  }
  else if(abseta >= 2.6 && abseta < 2.65) {
    ea = 0.2537;
  }
  else if(abseta >= 2.65 && abseta < 2.7) {
    ea = 0.2529;
  }
  else if(abseta >= 2.7 && abseta < 2.8) {
    ea = 0.2563;
  }
  else if(abseta >= 2.8 && abseta < 3.0) {
    ea = 0.2423;
  }
  else {
    ea = 0.0;
  }
  
  return ea;
}
