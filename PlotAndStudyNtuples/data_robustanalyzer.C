#include "data_robustanalyzer.hh"
#include <iostream>
#include <numeric>
#include <boost/range/combine.hpp>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
data_robustanalyzer::data_robustanalyzer(TString filename, TString outfilename, bool isDoubleElectron, bool isDYToLL, bool isSimulation){

  isDiEl = isDoubleElectron;
  isDY = isDYToLL;
  isMC = isSimulation;
  
  TFile *inpfile = TFile::Open(filename,"READ");
  cout<<"Initializing for file: "<<filename<<endl;

  tree = new TTreeReader("mmtree/tree",inpfile);
  l1Result = new TTreeReaderValue<vector<bool>>((*tree), "l1Result");
  bsx = new TTreeReaderArray<float>((*tree), "beamspot_x");
  bsy = new TTreeReaderArray<float>((*tree), "beamspot_y");
  bsz = new TTreeReaderArray<float>((*tree), "beamspot_z");
  if(isMC) {
    n_gen = new TTreeReaderValue<unsigned int>((*tree), "n_genpart");
    gen_pdg = new TTreeReaderValue<vector<int>>((*tree), "genpart_pdg");
    gen_pt = new TTreeReaderValue<vector<float>>((*tree), "genpart_pt");
    gen_eta = new TTreeReaderValue<vector<float>>((*tree), "genpart_eta");
    gen_phi = new TTreeReaderValue<vector<float>>((*tree), "genpart_phi");
    gen_vx = new TTreeReaderValue<vector<float>>((*tree), "genpart_vx");
    gen_vy = new TTreeReaderValue<vector<float>>((*tree), "genpart_vy");
    gen_vz = new TTreeReaderValue<vector<float>>((*tree), "genpart_vz");
    gen_nmoms = new TTreeReaderValue<vector<int>>((*tree), "genpart_nmoms");
    gen_mompdg = new TTreeReaderValue<vector<int>>((*tree), "genpart_mompdg");
    gen_isPFS = new TTreeReaderValue<vector<bool>>((*tree), "genpart_isPromptFS");
    gen_isLastC = new TTreeReaderValue<vector<bool>>((*tree), "genpart_isLastCopy");
  }
  n_ele = new TTreeReaderValue<UInt_t>((*tree), "n_ele");
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
  ele_seedid = new TTreeReaderValue<vector<unsigned int>>((*tree), "Electron_seedid");
  ele_enemat = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_energymatrix");
  ele_timmat = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_timingmatrix");
  
  outfile = new TFile(outfilename,"RECREATE");
}

// Fill the root file, close the root file, and handle deletions
data_robustanalyzer::~data_robustanalyzer() {

  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void data_robustanalyzer::analyzersinglefile(int splitCnt) { // Assume splitCnt to range from 0 to nCores

  int totEntries = tree->GetEntries();
  cout<<"Total number of entries: "<<totEntries<<endl;

  // Verfied that this logic to parallelize works
  int nCores = 6; // Assume parallel processing over 7 cores where
  // there is a lesser no.of events in the last core
  int beginevent = splitCnt*(totEntries/nCores);
  int endevent = (splitCnt+1)*(totEntries/nCores);
  if(beginevent>=totEntries) return;
  endevent = endevent<totEntries?endevent:totEntries;
  tree->SetEntriesRange(beginevent, endevent);
  int event = beginevent-1;

  // Count events passing certain selections
  double nosel=0, oldL1sel=0, unpreEGL1sel=0, EG_25_12_L1sel=0, EG_25_14_L1sel=0, EG_27_14_L1sel=0, EGLIso_22_12_L1sel=0, EGLIso_25_12_L1sel=0, EG_18_17_8_L1sel=0, EG_18_18_12_L1sel=0, EG_16_16_16_L1sel=0, EG8_HTT300_L1sel=0, preEGL1sel=0, zeroEGL1sel=0, oldscoutsel=0, oldscoutvetoidsel=0;

  // Define the histograms
  addgenhist("noselgen");
  addgenhist("basicselgen");
  addhist("nosel");
  addhist("oldL1sel");
  addhist("unpreEGL1sel");
  addhist("EG_25_12_L1sel");
  addhist("EG_25_14_L1sel");
  addhist("EG_27_14_L1sel");
  addhist("EGLIso_22_12_L1sel");
  addhist("EGLIso_25_12_L1sel");
  addhist("EG_18_17_8_L1sel");
  addhist("EG_18_18_12_L1sel");
  addhist("EG_16_16_16_L1sel");
  addhist("EG8_HTT300_L1sel");
  addhist("preEGL1sel");
  addhist("zeroEGL1sel");
  addhist("oldscoutsel");
  addhist("oldscoutvetosel");
  addhist("randsel");
  addgenmchhist("noselgenAnosel");
  addgenmchhist("basicselgenAnosel");
  addgenmchhist("noselgenAoldL1sel");
  addgenmchhist("basicselgenAoldL1sel");
  addgenmchhist("noselgenAunpreEGL1sel");
  addgenmchhist("noselgenAEG_25_12_L1sel");
  addgenmchhist("noselgenAEG_25_14_L1sel");
  addgenmchhist("noselgenAEG_27_14_L1sel");
  addgenmchhist("noselgenAEGLIso_22_12_L1sel");
  addgenmchhist("noselgenAEGLIso_25_12_L1sel");
  addgenmchhist("noselgenAEG_18_17_8_L1sel");
  addgenmchhist("noselgenAEG_18_18_12_L1sel");
  addgenmchhist("noselgenAEG_16_16_16_L1sel");
  addgenmchhist("noselgenAEG8_HTT300_L1sel");
  addgenmchhist("noselgenApreEGL1sel");
  addgenmchhist("noselgenAzeroEGL1sel");
  addgenmchhist("noselgenAoldscoutsel");
  addgenmchhist("basicselgenAoldscoutsel");
  addgenmchhist("noselgenAoldscoutvetosel");
  addgenmchhist("basicselgenAoldscoutvetosel");
  addgenmchhist("noselgenArandsel");
  addgenmchhist("basicselgenArandsel");
  
  // vector of electron indices
  vector<int> noselgenidx;
  vector<int> basicselgenidx;
  vector<int> noselelidx;
  vector<int> oldL1selelidx;
  vector<int> unpreEGL1selelidx;
  vector<int> EG_25_12_L1selelidx;
  vector<int> EG_25_14_L1selelidx;
  vector<int> EG_27_14_L1selelidx;
  vector<int> EGLIso_22_12_L1selelidx;
  vector<int> EGLIso_25_12_L1selelidx;
  vector<int> EG_18_17_8_L1selelidx;
  vector<int> EG_18_18_12_L1selelidx;
  vector<int> EG_16_16_16_L1selelidx;
  vector<int> EG8_HTT300_L1selelidx;
  vector<int> preEGL1selelidx;
  vector<int> zeroEGL1selelidx;
  vector<int> oldscoutselelidx;
  vector<int> oldscoutvetoselelidx;
  vector<int> randselelidx;
  
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>1000) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // L1 selection decision
    unsigned int l1size = l1name.size();
    double l1_nosel = 0.0; 
    l1_nosel = getL1decision({0},{l1size-1});
    double l1_oldsel = 0.0; 
    l1_oldsel = getL1decision({0},{18});
    double l1_unpreEGsel = 0.0; 
    l1_unpreEGsel = getL1decision({0,33},{18,41});
    double l1_EG_25_12_sel = 0.0; 
    l1_EG_25_12_sel = getL1decision({0,33},{18,35});
    double l1_EG_25_14_sel = 0.0; 
    l1_EG_25_14_sel = getL1decision({0,34},{18,34});
    double l1_EG_27_14_sel = 0.0; 
    l1_EG_27_14_sel = getL1decision({0,35},{18,35});
    double l1_EGLIso_22_12_sel = 0.0; 
    l1_EGLIso_22_12_sel = getL1decision({0,36},{18,37});
    double l1_EGLIso_25_12_sel = 0.0; 
    l1_EGLIso_25_12_sel = getL1decision({0,37},{18,37});
    double l1_EG_18_17_8_sel = 0.0; 
    l1_EG_18_17_8_sel = getL1decision({0,38},{18,40});
    double l1_EG_18_18_12_sel = 0.0; 
    l1_EG_18_18_12_sel = getL1decision({0,39},{18,39});
    double l1_EG_16_16_16_sel = 0.0; 
    l1_EG_16_16_16_sel = getL1decision({0,40},{18,40});
    double l1_EG8_HTT300_sel = 0.0; 
    l1_EG8_HTT300_sel = getL1decision({0,41},{18,41});
    double l1_preEGsel = 0.0; 
    l1_preEGsel = getL1decision({0,36},{18,37});
    double l1_zeroEGsel = 0.0;
    l1_zeroEGsel = getL1decision({0,24},{18,32});

    if(l1name.size()!=42 || l1prescale.size()!=42 || (*l1Result)->size()!=42) throw "Error!! Inconsistent l1 bit size.";
    
    if(isMC) {
      // In DY MC, choose the DiEl topology
      bool isDYDiEl = false;
      if(isDY) {
	int countHPelectron = 0;
	for(unsigned int gen_ctr=0; gen_ctr<(*(*n_gen)); gen_ctr++) {
	  if( TMath::Abs((*gen_pdg)->at(gen_ctr))==11 && (*gen_isPFS)->at(gen_ctr) ) {
	    countHPelectron++;
	  }
	}
	if(countHPelectron==2) isDYDiEl = true;
      }
      
      // Loop on Gen particles
      for(unsigned int gen_ctr=0; gen_ctr<(*(*n_gen)); gen_ctr++) {
	
	bool noselcond = true;
	noselcond *= TMath::Abs((*gen_pdg)->at(gen_ctr))==11;
	noselcond *= isDiEl?true:( isDY?(*gen_isPFS)->at(gen_ctr):(*gen_isLastC)->at(gen_ctr) );
	noselcond *= isDY?isDYDiEl:true;
	if(noselcond) noselgenidx.push_back(gen_ctr);
	
	bool basicselcond = true;
	basicselcond *= TMath::Abs((*gen_pdg)->at(gen_ctr))==11;
	basicselcond *= isDiEl?true:( isDY?(*gen_isPFS)->at(gen_ctr):(*gen_isLastC)->at(gen_ctr) );
	basicselcond *= isDY?isDYDiEl:true;
	basicselcond *= TMath::Abs((*gen_eta)->at(gen_ctr))<=2.7;
	if(basicselcond) basicselgenidx.push_back(gen_ctr);

      } // End of loop on gen electrons
      
      fillgenhistinevent("noselgen", noselgenidx);
      fillgenhistinevent("basicselgen", basicselgenidx);
    }
    
    // Sort the electrons based on their pT
    vector<int> sortedelidx((*(*n_ele)));
    iota(begin(sortedelidx), end(sortedelidx), 0);
    sort(&sortedelidx[0], ele_pt, (*(*n_ele))); // Verified that the algorithm works fine
    
    if((*(*n_ele))<0) cout<<"Error!!! Wrong technical event processing. Negative number of electrons in event."<<endl;;
    
    // Loop on scouting electrons in the event loop
    for(unsigned int ele_ctr=0; ele_ctr<(*(*n_ele)); ele_ctr++) {

      // Take the sorted index only
      unsigned int elidx = sortedelidx[ele_ctr];

      bool noselpass = true;
      noselpass *= l1_nosel>0;
      noselelidx.push_back(elidx);
      /*
      bool cut1elpass = true;
      cut1elpass *= TMath::Abs((*ele_eta)->at(elidx))<2.5;
      cut1elpass *= TMath::Abs((*ele_eta)->at(elidx))<1.479?(*ele_sigmaietaieta)->at(elidx)<0.0126:(*ele_sigmaietaieta)->at(elidx)<0.0457;
      if(cut1elpass) cut1elidx.push_back(elidx);
      */

      bool oldL1elpass = true;
      oldL1elpass *= l1_oldsel>0;
      if(oldL1elpass) oldL1selelidx.push_back(elidx);
	
      bool unpreEGL1elpass = true;
      unpreEGL1elpass *= l1_unpreEGsel>0;
      if(unpreEGL1elpass) unpreEGL1selelidx.push_back(elidx);

      bool EG_25_12_L1elpass = true;
      EG_25_12_L1elpass *= l1_EG_25_12_sel>0;
      if(EG_25_12_L1elpass) EG_25_12_L1selelidx.push_back(elidx);
	
      bool EG_25_14_L1elpass = true;
      EG_25_14_L1elpass *= l1_EG_25_14_sel>0;
      if(EG_25_14_L1elpass) EG_25_14_L1selelidx.push_back(elidx);
	
      bool EG_27_14_L1elpass = true;
      EG_27_14_L1elpass *= l1_EG_27_14_sel>0;
      if(EG_27_14_L1elpass) EG_27_14_L1selelidx.push_back(elidx);
	
      bool EGLIso_22_12_L1elpass = true;
      EGLIso_22_12_L1elpass *= l1_EGLIso_22_12_sel>0;
      if(EGLIso_22_12_L1elpass) EGLIso_22_12_L1selelidx.push_back(elidx);
	
      bool EGLIso_25_12_L1elpass = true;
      EGLIso_25_12_L1elpass *= l1_EGLIso_25_12_sel>0;
      if(EGLIso_25_12_L1elpass) EGLIso_25_12_L1selelidx.push_back(elidx);
	
      bool EG_18_17_8_L1elpass = true;
      EG_18_17_8_L1elpass *= l1_EG_18_17_8_sel>0;
      if(EG_18_17_8_L1elpass) EG_18_17_8_L1selelidx.push_back(elidx);
	
      bool EG_18_18_12_L1elpass = true;
      EG_18_18_12_L1elpass *= l1_EG_18_18_12_sel>0;
      if(EG_18_18_12_L1elpass) EG_18_18_12_L1selelidx.push_back(elidx);
	
      bool EG_16_16_16_L1elpass = true;
      EG_16_16_16_L1elpass *= l1_EG_16_16_16_sel>0;
      if(EG_16_16_16_L1elpass) EG_16_16_16_L1selelidx.push_back(elidx);
	
      bool EG8_HTT300_L1elpass = true;
      EG8_HTT300_L1elpass *= l1_EG8_HTT300_sel>0;
      if(EG8_HTT300_L1elpass) EG8_HTT300_L1selelidx.push_back(elidx);
	
      bool preEGL1elpass = true;
      preEGL1elpass *= l1_preEGsel>0;
      if(preEGL1elpass) preEGL1selelidx.push_back(elidx);
	
      bool zeroEGL1elpass = true;
      zeroEGL1elpass *= l1_zeroEGsel>0;
      if(zeroEGL1elpass) zeroEGL1selelidx.push_back(elidx);
	
      bool oldscoutselpass = true;
      oldscoutselpass *= l1_oldsel>0;
      oldscoutselpass *= TMath::Abs((*ele_eta)->at(elidx))<2.5;
      oldscoutselpass *= (*ele_pt)->at(elidx)>4;
      oldscoutselpass *= (*ele_hoe)->at(elidx)<0.2;
      if(oldscoutselpass) oldscoutselelidx.push_back(elidx);
	
      bool oldscoutvetoselpass = true;
      oldscoutvetoselpass *= l1_oldsel>0;
      oldscoutvetoselpass *= TMath::Abs((*ele_eta)->at(elidx))<2.5;
      oldscoutvetoselpass *= (*ele_pt)->at(elidx)>4;
      oldscoutvetoselpass *= (*ele_hoe)->at(elidx)<0.2;
      oldscoutvetoselpass *= TMath::Abs((*ele_eta)->at(elidx))<1.479?(*ele_sigmaietaieta)->at(elidx)<0.0126:(*ele_sigmaietaieta)->at(elidx)<0.0457;
      if(oldscoutvetoselpass) oldscoutvetoselelidx.push_back(elidx);
	
      bool randselpass = true;
      randselpass *= l1_nosel>0;
      randselpass *= TMath::Abs((*ele_eta)->at(elidx))<2.5;
      randselpass *= (*ele_pt)->at(elidx)>4;
      randselpass *= (*ele_hoe)->at(elidx)<0.2;
      randselpass *= TMath::Abs((*ele_eta)->at(elidx))<1.479?(*ele_sigmaietaieta)->at(elidx)<0.0126:(*ele_sigmaietaieta)->at(elidx)<0.0457;
      if(randselpass) randselelidx.push_back(elidx);
	
    }// End of loop on electrons in the event loop

    /*
    // Event level selection on the electrons
    pair<int,int> noselZwindels = inZwindow(noselelidx);
    if(noselZwindels.first!=-1) {
      noselZwindelidx.push_back(noselZwindels.first);
      noselZwindelidx.push_back(noselZwindels.second);
    }
    pair<int,int> cut1Zwindels = inZwindow(cut1elidx);
    if(cut1Zwindels.first!=-1) {
      cut1Zwindelidx.push_back(cut1Zwindels.first);
      cut1Zwindelidx.push_back(cut1Zwindels.second);
    }
    */
    
    if(noselelidx.size()>0) nosel += l1_nosel;
    if(oldL1selelidx.size()>0) oldL1sel += l1_oldsel;
    if(unpreEGL1selelidx.size()>0) unpreEGL1sel += l1_unpreEGsel;
    if(EG_25_12_L1selelidx.size()>0) EG_25_12_L1sel += l1_EG_25_12_sel;
    if(EG_25_14_L1selelidx.size()>0) EG_25_14_L1sel += l1_EG_25_14_sel;
    if(EG_27_14_L1selelidx.size()>0) EG_27_14_L1sel += l1_EG_27_14_sel;
    if(EGLIso_22_12_L1selelidx.size()>0) EGLIso_22_12_L1sel += l1_EGLIso_22_12_sel;
    if(EGLIso_25_12_L1selelidx.size()>0) EGLIso_25_12_L1sel += l1_EGLIso_25_12_sel;
    if(EG_18_17_8_L1selelidx.size()>0) EG_18_17_8_L1sel += l1_EG_18_17_8_sel;
    if(EG_18_18_12_L1selelidx.size()>0) EG_18_18_12_L1sel += l1_EG_18_18_12_sel;
    if(EG_16_16_16_L1selelidx.size()>0) EG_16_16_16_L1sel += l1_EG_16_16_16_sel;
    if(EG8_HTT300_L1selelidx.size()>0) EG8_HTT300_L1sel += l1_EG8_HTT300_sel;
    if(preEGL1selelidx.size()>0) preEGL1sel += l1_preEGsel;
    if(zeroEGL1selelidx.size()>0) zeroEGL1sel += l1_zeroEGsel;
    if(oldscoutselelidx.size()>0) oldscoutsel += l1_oldsel;
    fillhistinevent("nosel", noselelidx, l1_nosel);
    fillhistinevent("oldL1sel", oldL1selelidx, l1_oldsel);
    fillhistinevent("unpreEGL1sel", unpreEGL1selelidx, l1_unpreEGsel);
    fillhistinevent("EG_25_12_L1sel", EG_25_12_L1selelidx, l1_EG_25_12_sel);
    fillhistinevent("EG_25_14_L1sel", EG_25_14_L1selelidx, l1_EG_25_14_sel);
    fillhistinevent("EG_27_14_L1sel", EG_27_14_L1selelidx, l1_EG_27_14_sel);
    fillhistinevent("EGLIso_22_12_L1sel", EGLIso_22_12_L1selelidx, l1_EGLIso_22_12_sel);
    fillhistinevent("EGLIso_25_12_L1sel", EGLIso_25_12_L1selelidx, l1_EGLIso_25_12_sel);
    fillhistinevent("EG_18_17_8_L1sel", EG_18_17_8_L1selelidx, l1_EG_18_17_8_sel);
    fillhistinevent("EG_18_18_12_L1sel", EG_18_18_12_L1selelidx, l1_EG_18_18_12_sel);
    fillhistinevent("EG_16_16_16_L1sel", EG_16_16_16_L1selelidx, l1_EG_16_16_16_sel);
    fillhistinevent("EG8_HTT300_L1sel", EG8_HTT300_L1selelidx, l1_EG8_HTT300_sel);
    fillhistinevent("preEGL1sel", preEGL1selelidx, l1_preEGsel);
    fillhistinevent("zeroEGL1sel", zeroEGL1selelidx, l1_zeroEGsel);
    fillhistinevent("oldscoutsel", oldscoutselelidx, l1_oldsel);
    fillhistinevent("oldscoutvetosel", oldscoutvetoselelidx, l1_oldsel);
    fillhistinevent("randsel", randselelidx, l1_nosel);
    if(noselelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAnosel", noselgenidx, noselelidx, l1_nosel);
    if(oldL1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAoldL1sel", noselgenidx, oldL1selelidx, l1_oldsel);
    if(noselelidx.size()>0 && basicselgenidx.size()>0) fillgenmchhistinevent("basicselgenAnosel", basicselgenidx, noselelidx, l1_nosel);
    if(oldL1selelidx.size()>0 && basicselgenidx.size()>0) fillgenmchhistinevent("basicselgenAoldL1sel", basicselgenidx, oldL1selelidx, l1_oldsel);
    if(unpreEGL1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAunpreEGL1sel", noselgenidx, unpreEGL1selelidx, l1_unpreEGsel);
    if(EG_25_12_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEG_25_12_L1sel", noselgenidx, EG_25_12_L1selelidx, l1_EG_25_12_sel);
    if(EG_25_14_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEG_25_14_L1sel", noselgenidx, EG_25_14_L1selelidx, l1_EG_25_14_sel);
    if(EG_27_14_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEG_27_14_L1sel", noselgenidx, EG_27_14_L1selelidx, l1_EG_27_14_sel);
    if(EGLIso_22_12_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEGLIso_22_12_L1sel", noselgenidx, EGLIso_22_12_L1selelidx, l1_EGLIso_22_12_sel);
    if(EGLIso_25_12_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEGLIso_25_12_L1sel", noselgenidx, EGLIso_25_12_L1selelidx, l1_EGLIso_25_12_sel);
    if(EG_18_17_8_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEG_18_17_8_L1sel", noselgenidx, EG_18_17_8_L1selelidx, l1_EG_18_17_8_sel);
    if(EG_18_18_12_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEG_18_18_12_L1sel", noselgenidx, EG_18_18_12_L1selelidx, l1_EG_18_18_12_sel);
    if(EG_16_16_16_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEG_16_16_16_L1sel", noselgenidx, EG_16_16_16_L1selelidx, l1_EG_16_16_16_sel);
    if(EG8_HTT300_L1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAEG8_HTT300_L1sel", noselgenidx, EG8_HTT300_L1selelidx, l1_EG8_HTT300_sel);
    if(preEGL1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenApreEGL1sel", noselgenidx, preEGL1selelidx, l1_preEGsel);
    if(zeroEGL1selelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAzeroEGL1sel", noselgenidx, zeroEGL1selelidx, l1_zeroEGsel);
    if(oldscoutselelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAoldscoutsel", noselgenidx, oldscoutselelidx, l1_oldsel);
    if(oldscoutselelidx.size()>0 && basicselgenidx.size()>0) fillgenmchhistinevent("basicselgenAoldscoutsel", basicselgenidx, oldscoutselelidx, l1_oldsel);
    if(oldscoutvetoselelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAoldscoutvetosel", noselgenidx, oldscoutvetoselelidx, l1_oldsel);
    if(oldscoutvetoselelidx.size()>0 && basicselgenidx.size()>0) fillgenmchhistinevent("basicselgenAoldscoutvetosel", basicselgenidx, oldscoutvetoselelidx, l1_oldsel);
    if(randselelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenArandsel", noselgenidx, randselelidx, l1_nosel);
    if(randselelidx.size()>0 && basicselgenidx.size()>0) fillgenmchhistinevent("basicselgenArandsel", basicselgenidx, randselelidx, l1_nosel);

    // Clear all vector
    noselgenidx.clear();
    basicselgenidx.clear();
    noselelidx.clear();
    oldL1selelidx.clear();
    unpreEGL1selelidx.clear();
    EG_25_12_L1selelidx.clear();
    EG_25_14_L1selelidx.clear();
    EG_27_14_L1selelidx.clear();
    EGLIso_22_12_L1selelidx.clear();
    EGLIso_25_12_L1selelidx.clear();
    EG_18_17_8_L1selelidx.clear();
    EG_18_18_12_L1selelidx.clear();
    EG_16_16_16_L1selelidx.clear();
    EG8_HTT300_L1selelidx.clear();
    preEGL1selelidx.clear();
    zeroEGL1selelidx.clear();
    oldscoutselelidx.clear();
    oldscoutvetoselelidx.clear();
    randselelidx.clear();
    
  } // End of loop on events

  cout<<totEntries<<"\t"<<endevent-beginevent<<"\t"<<nosel<<"\t"<<oldL1sel<<"\t"<<unpreEGL1sel<<"\t"<<preEGL1sel<<"\t"<<zeroEGL1sel<<"\t"<<oldscoutsel<<endl;
}

// Fill histograms for gen matching
void data_robustanalyzer::fillgenmchhistinevent(TString selection, vector<int> genidx, vector<int> elidx, double w) {

  // Variables before gen match
  TH1F* bardEta = (TH1F*) outfile->Get(selection+"genelsctbar_dEta");
  TH1F* barqdPhi = (TH1F*) outfile->Get(selection+"genelsctbar_qdPhi");
  TH1F* eedEta = (TH1F*) outfile->Get(selection+"genelsctee_dEta");
  TH1F* eeqdPhi = (TH1F*) outfile->Get(selection+"genelsctee_qdPhi");
  
  // Variables after gen match
  TH1F* mchbardEta = (TH1F*) outfile->Get(selection+"genelsctmchbar_dEta");
  TH1F* mchbarqdPhi = (TH1F*) outfile->Get(selection+"genelsctmchbar_qdPhi");
  TH1F* mcheedEta = (TH1F*) outfile->Get(selection+"genelsctmchee_dEta");
  TH1F* mcheeqdPhi = (TH1F*) outfile->Get(selection+"genelsctmchee_qdPhi");

  TH1F* genelmult = (TH1F*) outfile->Get(selection+"sctmchgenel_elmult");
  TH1F* genelpt = (TH1F*) outfile->Get(selection+"sctmchgenel_elpt");
  TH1F* geneleta = (TH1F*) outfile->Get(selection+"sctmchgenel_eleta");
  TH1F* genelphi = (TH1F*) outfile->Get(selection+"sctmchgenel_elphi");
  TH1F* gendielM = (TH1F*) outfile->Get(selection+"sctmchgenel_dielM");
  TH1F* genleadelpt = (TH1F*) outfile->Get(selection+"sctmchgenel_lead_elpt");
  TH1F* gensubleadelpt = (TH1F*) outfile->Get(selection+"sctmchgenel_sublead_elpt");

  TH1F* genmchsctelmult = (TH1F*) outfile->Get(selection+"genmchsct_elmult");
  TH1F* genmchsctelpt = (TH1F*) outfile->Get(selection+"genmchsct_elpt");
  TH1F* genmchscteleta = (TH1F*) outfile->Get(selection+"genmchsct_eleta");
  TH1F* genmchsctelphi = (TH1F*) outfile->Get(selection+"genmchsct_elphi");
  TH1F* genmchsctdielM = (TH1F*) outfile->Get(selection+"genmchsct_dielM");

  // Fill variables before gen match
  for(unsigned int sct=0; sct<elidx.size(); sct++) {
    for(unsigned int gen=0; gen<genidx.size(); gen++) {
      double charge = ((*gen_pdg)->at(genidx[gen]))/TMath::Abs((*gen_pdg)->at(genidx[gen]));
      if(TMath::Abs((*ele_eta)->at(elidx[sct]))<1.479) {
	bardEta->Fill((*ele_eta)->at(elidx[sct])-(*gen_eta)->at(genidx[gen]));
	barqdPhi->Fill(charge*((*ele_phi)->at(elidx[sct])-(*gen_phi)->at(genidx[gen])));
      }
      else {
	eedEta->Fill((*ele_eta)->at(elidx[sct])-(*gen_eta)->at(genidx[gen]));
	eeqdPhi->Fill(charge*((*ele_phi)->at(elidx[sct])-(*gen_phi)->at(genidx[gen])));
      }
    }
  }

  // Fill variables after gen match
  vector< pair<int,int> > sctelgenmch = diElecGenMatching(genidx, elidx);

  int gen1 = sctelgenmch[0].first;
  int el1 = sctelgenmch[0].second;
  int gen2 = sctelgenmch[1].first;
  int el2 = sctelgenmch[1].second;
  
  unsigned int leadptpos=-1, subleadptpos=-1;
  if(gen1!=-1 && gen2!=-1) {
    if((*gen_pt)->at(gen1) >= (*gen_pt)->at(gen2)) {
      leadptpos = gen1;
      subleadptpos = gen2;
    }
    else {
      leadptpos = gen2;
      subleadptpos = gen1;
    }
  }
  else if(gen1==-1) {
    leadptpos = gen1;
  }
  else {
    if(gen2==-1) leadptpos = gen2;
  }
  
  int genmchfound = 0;
  if(el1!=-1) {
    genmchfound++;
    genelpt->Fill((*gen_pt)->at(gen1), w);
    geneleta->Fill((*gen_eta)->at(gen1), w);
    genelphi->Fill((*gen_phi)->at(gen1), w);
  }
  if(el2!=-1) {
    genmchfound++;
    genelpt->Fill((*gen_pt)->at(gen2), w);
    geneleta->Fill((*gen_eta)->at(gen2), w);
    genelphi->Fill((*gen_phi)->at(gen2), w);
  } 
  genelmult->Fill(genmchfound, w);
  
  if(leadptpos!=-1) {
    genleadelpt->Fill((*gen_pt)->at(leadptpos), w);
  }
  if(subleadptpos!=-1) {
    gensubleadelpt->Fill((*gen_pt)->at(subleadptpos), w);
  }
  
  if(genmchfound==2) {
    TLorentzVector elec1, elec2;
    elec1.SetPtEtaPhiM((*gen_pt)->at(gen1),(*gen_eta)->at(gen1),(*gen_phi)->at(gen1),0.0005);
    elec2.SetPtEtaPhiM((*gen_pt)->at(gen2),(*gen_eta)->at(gen2),(*gen_phi)->at(gen2),0.0005);
    gendielM->Fill((elec1+elec2).M(), w);
  }
  
  if(el1!=-1) {
    double charge = ((*gen_pdg)->at(gen1))/TMath::Abs((*gen_pdg)->at(gen1));
    if(TMath::Abs((*ele_eta)->at(el1))<1.479) {
      mchbardEta->Fill((*ele_eta)->at(el1)-(*gen_eta)->at(gen1));
      mchbarqdPhi->Fill(charge*((*ele_phi)->at(el1)-(*gen_phi)->at(gen1)));
    }
    else {
      mcheedEta->Fill((*ele_eta)->at(el1)-(*gen_eta)->at(gen1));
      mcheeqdPhi->Fill(charge*((*ele_phi)->at(el1)-(*gen_phi)->at(gen1)));
    }
    genmchsctelpt->Fill((*ele_pt)->at(el1), w);
    genmchscteleta->Fill((*ele_eta)->at(el1), w);
    genmchsctelphi->Fill((*ele_phi)->at(el1), w);
  }
  
  if(el2!=-1) {
    double charge = ((*gen_pdg)->at(gen2))/TMath::Abs((*gen_pdg)->at(gen2));
    if(TMath::Abs((*ele_eta)->at(el2))<1.479) {
      mchbardEta->Fill((*ele_eta)->at(el2)-(*gen_eta)->at(gen2));
      mchbarqdPhi->Fill(charge*((*ele_phi)->at(el2)-(*gen_phi)->at(gen2)));
    }
    else {
      mcheedEta->Fill((*ele_eta)->at(el2)-(*gen_eta)->at(gen2));
      mcheeqdPhi->Fill(charge*((*ele_phi)->at(el2)-(*gen_phi)->at(gen2)));
    }
    genmchsctelpt->Fill((*ele_pt)->at(el2), w);
    genmchscteleta->Fill((*ele_eta)->at(el2), w);
    genmchsctelphi->Fill((*ele_phi)->at(el2), w);
  }
  genmchsctelmult->Fill(genmchfound, w);
  if(genmchfound==2){
    TLorentzVector elec1, elec2;
    elec1.SetPtEtaPhiM((*ele_pt)->at(el1),(*ele_eta)->at(el1),(*ele_phi)->at(el1),0.0005);
    elec2.SetPtEtaPhiM((*ele_pt)->at(el2),(*ele_eta)->at(el2),(*ele_phi)->at(el2),0.0005);
    genmchsctdielM->Fill((elec1+elec2).M(), w);
  }

}

// Function to fill a set of histograms for gen particles
void data_robustanalyzer::fillgenhistinevent(TString selection, vector<int> genidx) {

  if(genidx.size()==0) return;

  TH1F* elmult = (TH1F*) outfile->Get(selection+"gen_elmult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"gen_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"gen_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"gen_elphi");

  TH1F* leadelpt = (TH1F*) outfile->Get(selection+"gen_lead_elpt");
  TH1F* leadeleta = (TH1F*) outfile->Get(selection+"gen_lead_eleta");
  TH1F* leadelphi = (TH1F*) outfile->Get(selection+"gen_lead_elphi");

  TH1F* subleadelpt = (TH1F*) outfile->Get(selection+"gen_sublead_elpt");
  TH1F* subleadeleta = (TH1F*) outfile->Get(selection+"gen_sublead_eleta");
  TH1F* subleadelphi = (TH1F*) outfile->Get(selection+"gen_sublead_elphi");

  TH1F* dielM = (TH1F*) outfile->Get(selection+"gen_diel_M");
  TH1F* dieldeta = (TH1F*) outfile->Get(selection+"gen_diel_deta");
  TH1F* dieldphi = (TH1F*) outfile->Get(selection+"gen_diel_dphi");
  TH1F* dieldR = (TH1F*) outfile->Get(selection+"gen_diel_dR");

  unsigned int leadptpos=-1, subleadptpos=-1;
  
  elmult->Fill(genidx.size());
  for(unsigned int ctr=0; ctr<genidx.size(); ctr++) {
    elpt->Fill((*gen_pt)->at(genidx[ctr]));
    eleta->Fill((*gen_eta)->at(genidx[ctr]));
    elphi->Fill((*gen_phi)->at(genidx[ctr]));

    if(leadptpos==-1 || ( (*gen_pt)->at(genidx[ctr])>(*gen_pt)->at(leadptpos) )) {
      subleadptpos = leadptpos;
      leadptpos = genidx[ctr];
    }
    if(leadptpos!=genidx[ctr] &&
       (subleadptpos==-1 || ( (*gen_pt)->at(genidx[ctr])>(*gen_pt)->at(subleadptpos) ) )
       ) {
      subleadptpos = genidx[ctr];
    }
  }
  if(leadptpos!=-1) {
    leadelpt->Fill((*gen_pt)->at(leadptpos));
    leadeleta->Fill((*gen_eta)->at(leadptpos));
    leadelphi->Fill((*gen_phi)->at(leadptpos));
  }
  if(subleadptpos!=-1) {
    subleadelpt->Fill((*gen_pt)->at(subleadptpos));
    subleadeleta->Fill((*gen_eta)->at(subleadptpos));
    subleadelphi->Fill((*gen_phi)->at(subleadptpos));
  }
  if(leadptpos!=-1 && subleadptpos!=-1) {
    TLorentzVector el1, el2;
    el1.SetPtEtaPhiM((*gen_pt)->at(leadptpos),(*gen_eta)->at(leadptpos),(*gen_phi)->at(leadptpos),0.0005);
    el2.SetPtEtaPhiM((*gen_pt)->at(subleadptpos),(*gen_eta)->at(subleadptpos),(*gen_phi)->at(subleadptpos),0.0005);
    dielM->Fill((el1+el2).M());
    dieldeta->Fill((*gen_eta)->at(leadptpos)-(*gen_eta)->at(subleadptpos));
    dieldphi->Fill(el1.DeltaPhi(el2));
    dieldR->Fill(el1.DeltaR(el2));
  }

}

// Function to fill a set of histograms for scouting electrons
void data_robustanalyzer::fillhistinevent(TString selection, vector<int> elidx, double w) {

  if(elidx.size()==0) return;

  TH1F* elmult = (TH1F*) outfile->Get(selection+"sct_elmult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"sct_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"sct_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"sct_elphi");
  TH1F* dielM = (TH1F*) outfile->Get(selection+"sct_dielM");
  TH1F* elminpt = (TH1F*) outfile->Get(selection+"sct_elminpt");
  TH1F* elmaxpt = (TH1F*) outfile->Get(selection+"sct_elmaxpt");
  
  TH1F* barelpt = (TH1F*) outfile->Get(selection+"sctbar_elpt");
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
  
  elmult->Fill(elidx.size(), w);
  for(unsigned int ctr=0; ctr<elidx.size(); ctr++) {
    elpt->Fill((*ele_pt)->at(elidx[ctr]), w);
    eleta->Fill((*ele_eta)->at(elidx[ctr]), w);
    elphi->Fill((*ele_phi)->at(elidx[ctr]), w);
    for(unsigned int ctr2=ctr+1; ctr2<elidx.size(); ctr2++) {
      TLorentzVector el, el2;
      el.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr]),(*ele_eta)->at(elidx[ctr]),(*ele_phi)->at(elidx[ctr]),0.0005);
      el2.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr2]),(*ele_eta)->at(elidx[ctr2]),(*ele_phi)->at(elidx[ctr2]),0.0005);
      dielM->Fill((el+el2).M(), w);
    }

    if(TMath::Abs((*ele_eta)->at(elidx[ctr]))<1.479) {
      barelpt->Fill((*ele_pt)->at(elidx[ctr]), w);
      barelm->Fill((*ele_m)->at(elidx[ctr]), w);
      bareld0->Fill((*ele_d0)->at(elidx[ctr]), w);
      barellog10d0->Fill(TMath::Log10(TMath::Abs((*ele_d0)->at(elidx[ctr]))), w);
      bareldz->Fill((*ele_dz)->at(elidx[ctr]), w);
      bareldetain->Fill((*ele_detain)->at(elidx[ctr]), w);
      bareldphiin->Fill((*ele_dphiin)->at(elidx[ctr]), w);
      barelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]), w);
      barelhoe->Fill((*ele_hoe)->at(elidx[ctr]), w);
      barelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]), w);
      barelmhits->Fill((*ele_mhits)->at(elidx[ctr]), w);
      barelcharge->Fill((*ele_charge)->at(elidx[ctr]), w);
      barelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]), w);
      barelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]), w);
      bareltkiso->Fill((*ele_tkiso)->at(elidx[ctr]), w);
      barelr9->Fill((*ele_r9)->at(elidx[ctr]), w);
      barelsmin->Fill((*ele_smin)->at(elidx[ctr]), w);
      barelsmaj->Fill((*ele_smaj)->at(elidx[ctr]), w);
    }

    else {
      ecelpt->Fill((*ele_pt)->at(elidx[ctr]), w);
      ecelm->Fill((*ele_m)->at(elidx[ctr]), w);
      eceld0->Fill((*ele_d0)->at(elidx[ctr]), w);
      ecellog10d0->Fill(TMath::Log10(TMath::Abs((*ele_d0)->at(elidx[ctr]))), w);
      eceldz->Fill((*ele_dz)->at(elidx[ctr]), w);
      eceldetain->Fill((*ele_detain)->at(elidx[ctr]), w);
      eceldphiin->Fill((*ele_dphiin)->at(elidx[ctr]), w);
      ecelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]), w);
      ecelhoe->Fill((*ele_hoe)->at(elidx[ctr]), w);
      ecelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]), w);
      ecelmhits->Fill((*ele_mhits)->at(elidx[ctr]), w);
      ecelcharge->Fill((*ele_charge)->at(elidx[ctr]), w);
      ecelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]), w);
      ecelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]), w);
      eceltkiso->Fill((*ele_tkiso)->at(elidx[ctr]), w);
      ecelr9->Fill((*ele_r9)->at(elidx[ctr]), w);
      ecelsmin->Fill((*ele_smin)->at(elidx[ctr]), w);
      ecelsmaj->Fill((*ele_smaj)->at(elidx[ctr]), w);
    }
  } // End of main electron for loop

  elminpt->Fill((*ele_pt)->at(elidx[elidx.size()-1]), w);
  elmaxpt->Fill((*ele_pt)->at(elidx[0]), w);

}  

// Add histograms for gen matching
void data_robustanalyzer::addgenmchhist(TString selection) {
  
  // Variables before gen match
  all1dhists.push_back(new TH1F(selection+"genelsctbar_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctbar_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctee_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctee_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));
  
  // Variables after gen match
  all1dhists.push_back(new TH1F(selection+"genelsctmchbar_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmchbar_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmchee_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmchee_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));

  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elmult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elpt","gen e p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_eleta","gen e #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elphi","gen e #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_dielM","M(e,e)",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_lead_elpt","gen e_{1} p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_sublead_elpt","gen e_{2} p_{T} / GeV",1000,-10,990));

  all1dhists.push_back(new TH1F(selection+"genmchsct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"genmchsct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"genmchsct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genmchsct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"genmchsct_dielM","all M(e,e)",1000,-10,990));
}

// Function to add a set of histograms for gen electrons
void data_robustanalyzer::addgenhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"gen_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"gen_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_elphi","#phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_lead_elpt","e_{1} p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_lead_eleta","e_{1} #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_lead_elphi","e_{1} #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_sublead_elpt","e_{1} p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_sublead_eleta","e_{1} #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_sublead_elphi","e_{1} #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_diel_M","all M(e,e)",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_diel_deta","#Delta#eta(e_{1}, e_{2})",1400,-7,7));
  all1dhists.push_back(new TH1F(selection+"gen_diel_dphi","#Delta#phi(e_{1}, e_{2})",1400,-7,7));
  all1dhists.push_back(new TH1F(selection+"gen_diel_dR","#Delta R(e_{1}, e_{2})",1000,0,10));

}

// Function to add a set of histograms for scouting electrons
void data_robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_dielM","all M(e,e)",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_elminpt","min. p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_elmaxpt","max. p_{T} / GeV",1000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sctbar_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctbar_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldetain","#Delta#eta_{in}",10000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldphiin","#Delta#phi_{in}",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhoe","H/E",20000,0,20));
  all1dhists.push_back(new TH1F(selection+"sctbar_elooemoop","1/E-1/p",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_elcharge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sctbar_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sctec_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctec_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctec_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctec_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctec_eldetain","#Delta#eta_{in}",10000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_eldphiin","#Delta#phi_{in}",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctec_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_elhoe","H/E",20000,0,20));
  all1dhists.push_back(new TH1F(selection+"sctec_elooemoop","1/E-1/p",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctec_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctec_elcharge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sctec_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctec_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctec_elsmaj","smaj",1500,-0.1,1.4));
}

// Function to sort the indices based on a factor (Usually pT)
void data_robustanalyzer::sort(int* idx, TTreeReaderValue<std::vector<float>> *factor, int n) {
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

// Function to find an electron pair in the Z mass window 80<M<100
pair<int,int> data_robustanalyzer::inZwindow(vector<int> elidx) {
  
  if(elidx.size()<2) { // Require atleast two electrons
    return make_pair(-1,-1);
  }

  pair<int,int> foundZels = make_pair(-1,-1);
  TLorentzVector el1, el2;
  for(unsigned int ctr1=0; ctr1<elidx.size(); ctr1++) {
    for(unsigned int ctr2=ctr1+1; ctr2<elidx.size(); ctr2++) {
      el1.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr1]),(*ele_eta)->at(elidx[ctr1]),(*ele_phi)->at(elidx[ctr1]),0.0005);
      el2.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr2]),(*ele_eta)->at(elidx[ctr2]),(*ele_phi)->at(elidx[ctr2]),0.0005);
      if((el1+el2).M()>75 && (el1+el2).M()<100) {
	foundZels = make_pair(elidx[ctr1],elidx[ctr2]);
	break;
      }
    }
  }

  return foundZels;
}

vector< pair<int,int> > data_robustanalyzer::diElecGenMatching(vector<int> genidx, vector<int> sctelidx) {

  if(genidx.size()==0) {
    return { make_pair(-1,-1), make_pair(-1,-1) };
  }
  
  if(genidx.size()>2) {
    throw "Error! Analysis code only suitable for <=2 gen electrons";
  }
  
  if(!isMC) {
    throw "Error! Cannot do gen matching. Not MC file.";
  }
  
  vector< pair<int,int> > *gensctelmch = new vector< pair<int,int> >;
  for(unsigned int ctr=0; ctr<2; ctr++) {
    if(ctr<genidx.size()) {
      gensctelmch->push_back(make_pair(genidx[ctr],-1));
    }
    else {
      gensctelmch->push_back(make_pair(-1,-1));
    }
  }
  
  // Find the sctelidx with the best angular match
  for(int elidx : sctelidx) {
    // Loop over the gen particles
    for(auto it=gensctelmch->begin(); it!=gensctelmch->end(); it++) {

      int geneidx = (*it).first;
      int mchsctelidx = (*it).second;
      if(geneidx==-1 || mchsctelidx!=-1) continue;

      double diffeta = abs((*ele_eta)->at(elidx)-(*gen_eta)->at(geneidx)); 
      TLorentzVector vecsctel, vecgen;
      vecgen.SetPtEtaPhiM((*gen_pt)->at(geneidx),(*gen_eta)->at(geneidx),(*gen_phi)->at(geneidx),0.0005);
      vecsctel.SetPtEtaPhiM((*ele_pt)->at(elidx),(*ele_eta)->at(elidx),(*ele_phi)->at(elidx),0.0005);
      double qdiffphi = ( (*gen_pdg)->at(geneidx)/abs((*gen_pdg)->at(geneidx)) )*( vecgen.DeltaPhi(vecsctel) );
      // Condition for gen matching
      if(abs((*ele_eta)->at(elidx))<1.479) {
	if(diffeta<0.2 && qdiffphi<0.05 && qdiffphi>-0.3) {
	  (*it).first = geneidx;
	  (*it).second = elidx;
	}
      }
      else {
	if(diffeta<0.1 && qdiffphi<0.05 && qdiffphi>-0.2) {
	  (*it).first = geneidx;
	  (*it).second = elidx;
	}
      }
    } // End of gen loop
    
  } // End of unseeded egamma object loop
  
  return (*gensctelmch);
}

// Get the L1 decision based on their position in the l1name variable
double data_robustanalyzer::getL1decision(vector<unsigned int> startPos, vector<unsigned int> endPos) {
  
  if(startPos.size() != endPos.size()) throw "Error!!! Not matching position vector size for l1 trigger decision";
  
  double prescale = 0;
  
  for(auto tup: boost::combine(startPos, endPos)) {
    unsigned int start, end;
    boost::tie(start, end) = tup;

    for(; start<=end && start<(*l1Result)->size(); start++) {
      if( (*l1Result)->at(start) == true ) {
	if(prescale==0) prescale = l1prescale[start]*1.0;
	if(l1prescale[start]!=0) prescale = l1prescale[start]<prescale?l1prescale[start]*1.0:prescale*1.0;
      }
    }
  } // End of check on all sections of L1 triggers

  if(prescale==0.0) return 0.0;
  else return 1.0/prescale;
}
