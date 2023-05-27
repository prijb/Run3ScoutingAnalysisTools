#include "robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
robustanalyzer::robustanalyzer(TString filename, TString outfilename, int numCores, bool isMonteCarlo){

  nC = numCores;
  isMC = isMonteCarlo; 

  cout<<"Initializing for file: "<<filename<<endl;
  TChain* chain = new TChain("mmtree/tree");
  chain->Add(filename);

  tree = new TTreeReader(chain);
  hlt_isomu27 = new TTreeReaderValue<unsigned int>((*tree), "HLT_IsoMu27");
  hlt_mu50 = new TTreeReaderValue<unsigned int>((*tree), "HLT_Mu50");
  hlt_scouting = new TTreeReaderValue<unsigned int>((*tree), "DST_Run3PFScouting");
  hlt_sct_dimu3 = new TTreeReaderValue<unsigned int>((*tree), "DST_Scouting_DoubleMu3");
  hlt_sct_dieg = new TTreeReaderValue<unsigned int>((*tree), "DST_Scouting_EG16EG12");
  hlt_sct_eg30 = new TTreeReaderValue<unsigned int>((*tree), "DST_Scouting_EG30");
  hlt_sct_jet = new TTreeReaderValue<unsigned int>((*tree), "DST_Scouting_JetHT");
  hltmu_sct = new TTreeReaderValue<unsigned int>((*tree), "DST_HLTMuon_Run3PFScouting");
  hltoth_sctmnt = new TTreeReaderValue<unsigned int>((*tree), "HLT_OtherScoutingPFMonitor");
  n_ele = new TTreeReaderValue<unsigned int>((*tree), "n_ele");
  ele_pt = new TTreeReaderValue<vector<float>>((*tree), "Electron_pt");
  ele_eta = new TTreeReaderValue<vector<float>>((*tree), "Electron_eta");
  ele_phi = new TTreeReaderValue<vector<float>>((*tree), "Electron_phi");
  ele_trkpt = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_trkpt");
  ele_trketa = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_trketa");
  ele_trkphi = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_trkphi");
  ele_trkd0 = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_trkd0");
  ele_trkdz = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_trkdz");
  ele_trkcharge = new TTreeReaderValue<vector<vector<int>>>((*tree), "Electron_trkcharge");
  ele_trkchi2 = new TTreeReaderValue<vector<vector<float>>>((*tree), "Electron_trkrchi2");
  ele_m = new TTreeReaderValue<vector<float>>((*tree), "Electron_m");
  ele_detain = new TTreeReaderValue<vector<float>>((*tree), "Electron_detain");
  ele_dphiin = new TTreeReaderValue<vector<float>>((*tree), "Electron_dphiin");
  ele_sigmaietaieta = new TTreeReaderValue<vector<float>>((*tree), "Electron_sigmaietaieta");
  ele_hoe = new TTreeReaderValue<vector<float>>((*tree), "Electron_hoe");
  ele_ooemoop = new TTreeReaderValue<vector<float>>((*tree), "Electron_ooemoop");
  ele_mhits = new TTreeReaderValue<vector<int>>((*tree), "Electron_missinghits");
  ele_ecaliso = new TTreeReaderValue<vector<float>>((*tree), "Electron_ecaliso");
  ele_hcaliso = new TTreeReaderValue<vector<float>>((*tree), "Electron_hcaliso");
  ele_tkiso = new TTreeReaderValue<vector<float>>((*tree), "Electron_tkiso");
  ele_r9 = new TTreeReaderValue<vector<float>>((*tree), "Electron_r9");
  ele_smin = new TTreeReaderValue<vector<float>>((*tree), "Electron_smin");
  ele_smaj = new TTreeReaderValue<vector<float>>((*tree), "Electron_smaj");
  n_rho = new TTreeReaderValue<unsigned int>((*tree), "n_rhoval");
  rho = new TTreeReaderValue<vector<float>>((*tree), "rho");

  if(isMC) {
    n_genlep = new TTreeReaderValue<unsigned int>((*tree), "n_genpartlep");
    genlep_pdg = new TTreeReaderValue<vector<int>>((*tree), "genpartlep_pdg");
    genlep_pt = new TTreeReaderValue<vector<float>>((*tree), "genpartlep_pt");
    genlep_eta = new TTreeReaderValue<vector<float>>((*tree), "genpartlep_eta");
    genlep_phi = new TTreeReaderValue<vector<float>>((*tree), "genpartlep_phi");
    genlep_m = new TTreeReaderValue<vector<float>>((*tree), "genpartlep_m");
    genlep_vx = new TTreeReaderValue<vector<float>>((*tree), "genpartlep_vx");
    genlep_vy = new TTreeReaderValue<vector<float>>((*tree), "genpartlep_vy");
    genlep_vz = new TTreeReaderValue<vector<float>>((*tree), "genpartlep_vz");
    n_genjpsi = new TTreeReaderValue<unsigned int>((*tree), "n_genpartjpsi");
    genjpsi_pdg = new TTreeReaderValue<vector<int>>((*tree), "genpartjpsi_pdg");
    genjpsi_pt = new TTreeReaderValue<vector<float>>((*tree), "genpartjpsi_pt");
    genjpsi_eta = new TTreeReaderValue<vector<float>>((*tree), "genpartjpsi_eta");
    genjpsi_phi = new TTreeReaderValue<vector<float>>((*tree), "genpartjpsi_phi");
    genjpsi_m = new TTreeReaderValue<vector<float>>((*tree), "genpartjpsi_m");
    genjpsi_vx = new TTreeReaderValue<vector<float>>((*tree), "genpartjpsi_vx");
    genjpsi_vy = new TTreeReaderValue<vector<float>>((*tree), "genpartjpsi_vy");
    genjpsi_vz = new TTreeReaderValue<vector<float>>((*tree), "genpartjpsi_vz");
  }

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
   
  // Define the gen histograms
  if(isMC) {
    addgenhist("noselgen_");
    addgenhist("noselgen_dieg_");
    addgenhist("noselgen_eg30_");
    addgenhist("noselgen_muht_");
    addgenmchhist("noselAnosel_genmch_");

    addhist("noselAnosel_genmchSC_dieg_");
    addhist("noselAnosel_genmchSC_eg30_");
    addhist("noselAnosel_genmchTRK_dieg_");
    addhist("noselAnosel_genmchTRK_eg30_");
    addhist("noselAnosel_genmchTRK_diegOeg30_");
  }
  
  // Count events passing certain selections
  int nosel=0;
  
  addhist("nosel_");
  addhist("nosel_dieg_");
  addhist("nosel_eg30_");
  addhist("nosel_muht_");
  addhist("v1selSC_diegOeg30_");
  addhist("v1selTK_diegOeg30_");
  addhist("v1sel_diegOeg30_");

  // vector of gen indices
  vector<int> noselgenjpidx;
  vector<int> noselgenelidx;

  // vector of electron indices
  vector<int> noselelidx;
  vector<vector<int>> noseleltrkidx;
  vector<int> v1selscelidx;
  vector<vector<int>> v1selsceltrkidx;
  vector<int> v1seltkelidx;
  vector<vector<int>> v1seltkeltrkidx;
  vector<int> v1selelidx;
  vector<vector<int>> v1seleltrkidx;
    
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%1==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // Loop on Gen particles to select good gen electrons
    if(isMC) {

      // Check for 2 daughters per JPsi
      if((*(*n_genlep))/2!=(*(*n_genjpsi))) throw "Expecting two daughters for every JPsi";

      for(unsigned int gen_ctr=0; gen_ctr<(*(*n_genjpsi)); gen_ctr++) {
	
	noselgenjpidx.push_back(gen_ctr);
	noselgenelidx.push_back(2*gen_ctr);
	noselgenelidx.push_back(2*gen_ctr+1);
	
      } // End of loop on gen electrons
      if(noselgenjpidx.size()==1) fillgenhistinevent("noselgen_", noselgenelidx, noselgenjpidx);
      if( (noselgenjpidx.size()==1) && ((*(*hlt_sct_dieg))==1) ) fillgenhistinevent("noselgen_dieg_", noselgenelidx, noselgenjpidx);
      if( (noselgenjpidx.size()==1) && ((*(*hlt_sct_eg30))==1) ) fillgenhistinevent("noselgen_eg30_", noselgenelidx, noselgenjpidx);
      if( (noselgenjpidx.size()==1) && ((*(*hlt_sct_dimu3))==1) && ((*(*hlt_sct_jet))==1) ) fillgenhistinevent("noselgen_muht_", noselgenelidx, noselgenjpidx);
    }

    // Sort the electrons based on their pT
    vector<int> sortedelidx((*(*n_ele)));
    iota(begin(sortedelidx), end(sortedelidx), 0);
    sort(&sortedelidx[0], ele_pt, (*(*n_ele))); // Verified that the algorithm works fine

    if((*(*n_ele))<0) throw "Error!!! Wrong technical event processing. Negative number of electrons in event.";

    // Loop on electrons in the event loop
    for(unsigned int ele_ctr=0; ele_ctr<(*(*n_ele)); ele_ctr++) {
      
      // Take the sorted index only
      unsigned int elidx = sortedelidx[ele_ctr];
      
      noselelidx.push_back(elidx);
      vector<int> trki;
      for(unsigned int tkc=0; tkc<((*ele_trkpt)->at(elidx)).size(); tkc++) {
	trki.push_back(tkc);
      }
      noseleltrkidx.push_back(trki);
      trki.clear();

      bool v1selscsc = true;
      v1selscsc *= v1selfunc(-1, elidx, true);
      if(v1selscsc) {
	for(unsigned int tkc=0; tkc<((*ele_trkpt)->at(elidx)).size(); tkc++) {
	  trki.push_back(tkc);
	}
	if(trki.size()>0) {
	  v1selscelidx.push_back(elidx);
	  v1selsceltrkidx.push_back(trki);
	}
	trki.clear();
      }
      
      bool v1seltksc = true;
      if(v1seltksc) {
	for(unsigned int tkc=0; tkc<((*ele_trkpt)->at(elidx)).size(); tkc++) {
	  bool v1seltktk = true;
	  v1seltktk *= v1selfunc(tkc, elidx, false);
	  if(v1seltktk) trki.push_back(tkc);
	}
	if(trki.size()>0) {
	  v1seltkelidx.push_back(elidx);
	  v1seltkeltrkidx.push_back(trki);
	}
	trki.clear();
      }
      
      bool v1selsc = true;
      v1selsc *= v1selfunc(-1, elidx, true);
      if(v1selsc) {
	for(unsigned int tkc=0; tkc<((*ele_trkpt)->at(elidx)).size(); tkc++) {
	  bool v1seltk = true;
	  v1seltk *= v1selfunc(tkc, elidx, false);
	  if(v1seltk) trki.push_back(tkc);
	}
	if(trki.size()>0) {
	  v1selelidx.push_back(elidx);
	  v1seleltrkidx.push_back(trki);
	}
	trki.clear();
      }
      
    }// End of loop on electrons in the event loop
    
    if(noselelidx.size()>0) nosel++;
    fillhistinevent("nosel_", noselelidx, noseleltrkidx);
    if( ((*(*hlt_sct_dieg))==1) ) fillhistinevent("nosel_dieg_", noselelidx, noseleltrkidx);
    if( ((*(*hlt_sct_eg30))==1) ) fillhistinevent("nosel_eg30_", noselelidx, noseleltrkidx);
    if( ((*(*hlt_sct_dimu3))==1) && ((*(*hlt_sct_jet))==1) ) fillhistinevent("nosel_muht_", noselelidx, noseleltrkidx);
    if( ((*(*hlt_sct_eg30))==1) || ((*(*hlt_sct_dieg))==1) ) fillhistinevent("v1selSC_diegOeg30_", v1selscelidx, v1selsceltrkidx);
    if( ((*(*hlt_sct_eg30))==1) || ((*(*hlt_sct_dieg))==1) ) fillhistinevent("v1selTK_diegOeg30_", v1seltkelidx, v1seltkeltrkidx);
    if( ((*(*hlt_sct_eg30))==1) || ((*(*hlt_sct_dieg))==1) ) fillhistinevent("v1sel_diegOeg30_", v1selelidx, v1seleltrkidx);
  
    if(isMC) {

      fillgenmchhistinevent("noselAnosel_genmch_", noselgenelidx, noselelidx);
      
      std::pair< vector<int>,vector<int> > noselAnosel_gensctpair = getGenMatched(noselgenelidx, noselelidx);
      vector<int> noselAnosel_matchedgenidx = noselAnosel_gensctpair.first;
      vector<int> noselAnosel_matchedsctidx = noselAnosel_gensctpair.second;
      vector<vector<int>> noselAnosel_matchedscttrkidx;
      for(unsigned int si=0; si<noselAnosel_matchedsctidx.size(); si++) {
	vector<int> trki;
	for(unsigned int tkc=0; tkc<((*ele_trkpt)->at(noselAnosel_matchedsctidx[si])).size(); tkc++) {
	  trki.push_back(tkc);
	}
	noselAnosel_matchedscttrkidx.push_back(trki);
	trki.clear();
      }
      if( ((*(*hlt_sct_dieg))==1) ) fillhistinevent("noselAnosel_genmchSC_dieg_", noselAnosel_matchedsctidx, noselAnosel_matchedscttrkidx);
      if( ((*(*hlt_sct_eg30))==1) ) fillhistinevent("noselAnosel_genmchSC_eg30_", noselAnosel_matchedsctidx, noselAnosel_matchedscttrkidx);

      std::pair< vector<int>,vector< std::pair<int,int> > > noselAnosel_gentrkpair = getGenMatchedTrk(noselgenelidx, noselelidx);
      vector<int> noselAnosel_trkmatchedgenidx = noselAnosel_gentrkpair.first;
      vector< std::pair<int,int> > noselAnosel_trkmatchedidx = noselAnosel_gentrkpair.second;
      vector<int> noselAnosel_trkmatchedsctidx;
      vector<vector<int>> noselAnosel_trkmatchedtrkidx;
      for(unsigned int si=0; si<noselAnosel_trkmatchedidx.size(); si++) {
	noselAnosel_trkmatchedsctidx.push_back(noselAnosel_trkmatchedidx[si].first);
	vector<int> trki;
	trki.push_back(noselAnosel_trkmatchedidx[si].second);
        noselAnosel_trkmatchedtrkidx.push_back(trki);
	trki.clear();
      }
      if( ((*(*hlt_sct_dieg))==1) ) fillhistinevent("noselAnosel_genmchTRK_dieg_", noselAnosel_trkmatchedsctidx, noselAnosel_trkmatchedtrkidx);
      if( ((*(*hlt_sct_eg30))==1) ) fillhistinevent("noselAnosel_genmchTRK_eg30_", noselAnosel_trkmatchedsctidx, noselAnosel_trkmatchedtrkidx);
      if( ((*(*hlt_sct_dieg))==1) || ((*(*hlt_sct_eg30))==1) ) fillhistinevent("noselAnosel_genmchTRK_diegOeg30_", noselAnosel_trkmatchedsctidx, noselAnosel_trkmatchedtrkidx);
    }
    
    // Clear all vector
    noselgenjpidx.clear();
    noselgenelidx.clear();

    noselelidx.clear();
    noseleltrkidx.clear();
    v1selscelidx.clear();
    v1selsceltrkidx.clear();
    v1seltkelidx.clear();
    v1seltkeltrkidx.clear();
    v1selelidx.clear();
    v1seleltrkidx.clear();
    
  } // End of loop on events

  cout<<totEntries<<"\t"<<endevent-beginevent<<"\t"<<nosel<<endl;
}

// Function to perform gen matching for the super-cluster
std::pair< vector<int>,vector<int> > robustanalyzer::getGenMatched(vector<int> geni, vector<int> ei) {

  vector<int> genfoundi, sctfoundi;
  
  for(int gid: geni) {

    /////////////////////////
    for(unsigned int sct=0; sct<ei.size(); sct++) {
      double deta = ((*genlep_eta)->at(gid)) - ((*ele_eta)->at(ei[sct]));
      double dphi = ((*genlep_phi)->at(gid)) - ((*ele_phi)->at(ei[sct]));
      double qdph = (charge((*genlep_pdg)->at(gid))) * dphi;
      double eta = (*ele_eta)->at(ei[sct]);
      if(abs(eta)<1.479 && abs(deta)<0.05 && qdph>-0.06 && qdph<0) {
	genfoundi.push_back(gid);
	sctfoundi.push_back(ei[sct]);
	ei.erase(ei.begin()+sct);
	break;
      }
      else {
	if(abs(deta)<0.03 && qdph>-0.06 && qdph<0.01) {
	  genfoundi.push_back(gid);
	  sctfoundi.push_back(ei[sct]);
	  ei.erase(ei.begin()+sct);
	  break;
	}
      }
    }
    ////////////////////////
    
  }

  return std::make_pair(genfoundi, sctfoundi);
}

// Function to perform gen matching for the track
std::pair< vector<int>, vector< std::pair<int,int> > > robustanalyzer::getGenMatchedTrk(vector<int> geni, vector<int> ei) {

  vector<int> genfoundi;
  vector< pair<int,int> > trkfoundi;

  for(int gid: geni) {

    for(unsigned int sct=0; sct<ei.size(); sct++) {
      int sci = ei[sct];
      unsigned int ntk = ((*ele_trkpt)->at(sci)).size();
      double eta = (*ele_eta)->at(sci);

      /////////////////////////
      bool trkfound = false;
      for(unsigned int tki=0; tki<ntk; tki++) {
	double tdeta = ((*genlep_eta)->at(gid)) - ((*ele_trketa)->at(sci))[tki];
	double tdphi = ((*genlep_phi)->at(gid)) - ((*ele_trkphi)->at(sci))[tki];
	double tqdph = ((*ele_trkcharge)->at(sci))[tki] * tdphi;
	if(abs(eta)<1.479 && abs(tdeta)<0.0015 && tqdph>-0.004 && tqdph<0.002) {
	  genfoundi.push_back(gid);
	  trkfoundi.push_back(std::make_pair(sci, tki));
	  trkfound = true;
	  ei.erase(ei.begin()+sct);
	  break;
	}
	else {
	  if(abs(tdeta)<0.003 && tqdph>-0.006 && tqdph<0.003) {
	    genfoundi.push_back(gid);
	    trkfoundi.push_back(std::make_pair(sci, tki));
	    trkfound = true;
	    ei.erase(ei.begin()+sct);
	    break;
	  }
	}
      }
      if(trkfound) break;
      ////////////////////////
    }
  }

  return std::make_pair(genfoundi, trkfoundi);
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

// Function to obtain charge from pdg (Open to overloading)
int robustanalyzer::charge(int pdg) {
  return pdg/abs(pdg);
}

bool robustanalyzer::v1selfunc(int tkidx, int scidx, bool isSC) {

  bool sel = true;

  // Get the energy of this electron
  TLorentzVector se;
  se.SetPtEtaPhiM((*ele_pt)->at(scidx), (*ele_eta)->at(scidx), (*ele_phi)->at(scidx), 0.0005);
  double sE = se.Energy();
      
  if(isSC) {
    sel *= abs((*ele_eta)->at(scidx))<1.479 ? ((*ele_sigmaietaieta)->at(scidx))<0.0105 : ((*ele_sigmaietaieta)->at(scidx))<0.035 ;
    sel *= ((*ele_hoe)->at(scidx))<0.1;
    sel *= abs((*ele_eta)->at(scidx))<1.479 ? ((*ele_hcaliso)->at(scidx))<0.15*sE : ((*ele_hcaliso)->at(scidx))<0.2*sE;
    sel *= abs((*ele_eta)->at(scidx))<1.479 ? ((*ele_ooemoop)->at(scidx))<0.015 : ((*ele_ooemoop)->at(scidx))<0.01;
    sel *= abs((*ele_eta)->at(scidx))<1.479 ? ((*ele_detain)->at(scidx))<0.005 : ((*ele_detain)->at(scidx))<0.01;
    sel *= abs((*ele_eta)->at(scidx))<1.479 ? ((*ele_dphiin)->at(scidx))<0.03 : ((*ele_dphiin)->at(scidx))<0.05;
  }
  else {
    if(!isSC && tkidx<0) throw "Error!!! Tracker Index <0 in v1 selection.";
    TLorentzVector te;
    te.SetPtEtaPhiM(((*ele_trkpt)->at(scidx))[tkidx], ((*ele_trketa)->at(scidx))[tkidx], ((*ele_trkphi)->at(scidx))[tkidx], 0.0005);
    double tE = te.Energy();
    double dEoE = abs(sE - tE)/sE;
    double deta = abs((*ele_eta)->at(scidx) - ((*ele_trketa)->at(scidx))[tkidx]);
    double dphi = abs((*ele_phi)->at(scidx) - ((*ele_trkphi)->at(scidx))[tkidx]);
    sel *= abs(((*ele_trkd0)->at(scidx))[tkidx])<0.25;
    sel *= abs(((*ele_trkpt)->at(scidx))[tkidx])>12;
    sel *= dEoE<1;
    //sel *= abs((*ele_eta)->at(scidx))<1.479 ? deta<0.04 : deta<0.02;
    sel *= dphi<0.06;
    sel *= abs((*ele_eta)->at(scidx))<1.479 ? (((*ele_trkchi2)->at(scidx))[tkidx])<3 : (((*ele_trkchi2)->at(scidx))[tkidx])<2;
  }
  
  return sel;
}

// Function to add a set of histograms for scouting electrons
void robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_rho","rho",10000,-10,90));
  all1dhists.push_back(new TH1F(selection+"sct_isomu27","HLT_IsoMu27",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_mu50","HLT_Mu50",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_scouting","HLT_Run3PixelOnlyPFScouting",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_dielM","all M(e,e)",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadsublead_dielM","M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadsublead_chargeprod","q_{e1}*q_{e2}",3,-1,1));
  all1dhists.push_back(new TH1F(selection+"sct_leadbarsubleadbar_dielM","bar.e_{1} M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadecsubleadec_dielM","ec. e_{1} M(e_{1},e_{2})",100000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sct_trk_dielM","all trk M(e,e)",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_trk_elmult","N trk e",50,-5,45));

  all1dhists.push_back(new TH1F(selection+"sct_bar_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_bar_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_eldetain","#Delta#eta_{in}",100000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_eldphiin","#Delta#phi_{in}",100000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_bar_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elecalisoovere","ecal. iso. / E",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elhcalisoovere","hcal. iso. / E",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_bar_eltkisoovere","tk. iso. / E",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sct_bar_elsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_elpt","trk p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_eleta","trk #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_elphi","trk #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_eld0","trk d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_ellog10d0","trk log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_eldz","trk d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_elcharge","trk charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_elrchi2","trk reduced chi sq.",10000,0,100));

  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_el_min_dE","min #DeltaE(SC, trk)",1000,-50,50));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_el_min_dEoverE","min #DeltaE(SC, trk) / SC E",1000,-50,50));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_el_min_deta","min #Delta#eta(SC, trk)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_trk_el_min_dphi","min #Delta#phi(SC, trk)",10000,-5,5));
  
  all1dhists.push_back(new TH1F(selection+"sct_ec_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_ec_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_eldetain","#Delta#eta_{in}",100000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_eldphiin","#Delta#phi_{in}",100000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_ec_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elecalisoovere","ecal. iso. / E",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elhcalisoovere","hcal. iso. / E",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_ec_eltkisoovere","tk. iso. / E",10000,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sct_ec_elsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_elpt","trk p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_eleta","trk #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_elphi","trk #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_eld0","trk d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_ellog10d0","trk log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_eldz","trk d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_elcharge","trk charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_elrchi2","trk reduced chi sq.",10000,0,100));

  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_el_min_dE","min #DeltaE(SC, trk)",1000,-50,50));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_el_min_dEoverE","min #DeltaE(SC, trk) / SC E",1000,-50,50));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_el_min_deta","min #Delta#eta(SC, trk)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_trk_el_min_dphi","min #Delta#phi(SC, trk)",10000,-5,5));
}

// Function to fill a set of histograms for scouting electrons
void robustanalyzer::fillhistinevent(TString selection, vector<int> elidx, vector<vector<int>> tkidx) {

  if(tkidx.size()!=elidx.size()) throw "SC index count should match TRK index count";
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

  TH1F* trkdielM = (TH1F*) outfile->Get(selection+"sct_trk_dielM");
  TH1F* trkmult = (TH1F*) outfile->Get(selection+"sct_trk_elmult");
  
  TH1F* barelpt = (TH1F*) outfile->Get(selection+"sct_bar_elpt");
  TH1F* bareleta = (TH1F*) outfile->Get(selection+"sct_bar_eleta");
  TH1F* barelphi = (TH1F*) outfile->Get(selection+"sct_bar_elphi");
  TH1F* barelm = (TH1F*) outfile->Get(selection+"sct_bar_elm");
  TH1F* bareldetain = (TH1F*) outfile->Get(selection+"sct_bar_eldetain");
  TH1F* bareldphiin = (TH1F*) outfile->Get(selection+"sct_bar_eldphiin");
  TH1F* barelsigmaietaieta = (TH1F*) outfile->Get(selection+"sct_bar_elsigmaietaieta");
  TH1F* barelhoe = (TH1F*) outfile->Get(selection+"sct_bar_elhoe");
  TH1F* barelooemoop = (TH1F*) outfile->Get(selection+"sct_bar_elooemoop");
  TH1F* barelmhits = (TH1F*) outfile->Get(selection+"sct_bar_elmhits");
  TH1F* barelecaliso = (TH1F*) outfile->Get(selection+"sct_bar_elecaliso");
  TH1F* barelhcaliso = (TH1F*) outfile->Get(selection+"sct_bar_elhcaliso");
  TH1F* bareltkiso = (TH1F*) outfile->Get(selection+"sct_bar_eltkiso");
  TH1F* barelecalisoovere = (TH1F*) outfile->Get(selection+"sct_bar_elecalisoovere");
  TH1F* barelhcalisoovere = (TH1F*) outfile->Get(selection+"sct_bar_elhcalisoovere");
  TH1F* bareltkisoovere = (TH1F*) outfile->Get(selection+"sct_bar_eltkisoovere");
  TH1F* barelr9 = (TH1F*) outfile->Get(selection+"sct_bar_elr9");
  TH1F* barelsmin = (TH1F*) outfile->Get(selection+"sct_bar_elsmin");
  TH1F* barelsmaj = (TH1F*) outfile->Get(selection+"sct_bar_elsmaj");
  
  TH1F* trkbarelpt = (TH1F*) outfile->Get(selection+"sct_bar_trk_elpt");
  TH1F* trkbareleta = (TH1F*) outfile->Get(selection+"sct_bar_trk_eleta");
  TH1F* trkbarelphi = (TH1F*) outfile->Get(selection+"sct_bar_trk_elphi");
  TH1F* trkbareld0 = (TH1F*) outfile->Get(selection+"sct_bar_trk_eld0");
  TH1F* trkbarellog10d0 = (TH1F*) outfile->Get(selection+"sct_bar_trk_ellog10d0");
  TH1F* trkbareldz = (TH1F*) outfile->Get(selection+"sct_bar_trk_eldz");
  TH1F* trkbarelcharge = (TH1F*) outfile->Get(selection+"sct_bar_trk_elcharge");
  TH1F* trkbarelrchi2 = (TH1F*) outfile->Get(selection+"sct_bar_trk_elrchi2");

  TH1F* trkbarmindE = (TH1F*) outfile->Get(selection+"sct_bar_trk_el_min_dE");
  TH1F* trkbarmindEoverE = (TH1F*) outfile->Get(selection+"sct_bar_trk_el_min_dEoverE");
  TH1F* trkbarmindeta = (TH1F*) outfile->Get(selection+"sct_bar_trk_el_min_deta");
  TH1F* trkbarmindphi = (TH1F*) outfile->Get(selection+"sct_bar_trk_el_min_dphi");

  TH1F* ecelpt = (TH1F*) outfile->Get(selection+"sct_ec_elpt");
  TH1F* eceleta = (TH1F*) outfile->Get(selection+"sct_ec_eleta");
  TH1F* ecelphi = (TH1F*) outfile->Get(selection+"sct_ec_elphi");
  TH1F* ecelm = (TH1F*) outfile->Get(selection+"sct_ec_elm");
  TH1F* eceldetain = (TH1F*) outfile->Get(selection+"sct_ec_eldetain");
  TH1F* eceldphiin = (TH1F*) outfile->Get(selection+"sct_ec_eldphiin");
  TH1F* ecelsigmaietaieta = (TH1F*) outfile->Get(selection+"sct_ec_elsigmaietaieta");
  TH1F* ecelhoe = (TH1F*) outfile->Get(selection+"sct_ec_elhoe");
  TH1F* ecelooemoop = (TH1F*) outfile->Get(selection+"sct_ec_elooemoop");
  TH1F* ecelmhits = (TH1F*) outfile->Get(selection+"sct_ec_elmhits");
  TH1F* ecelecaliso = (TH1F*) outfile->Get(selection+"sct_ec_elecaliso");
  TH1F* ecelhcaliso = (TH1F*) outfile->Get(selection+"sct_ec_elhcaliso");
  TH1F* eceltkiso = (TH1F*) outfile->Get(selection+"sct_ec_eltkiso");
  TH1F* ecelecalisoovere = (TH1F*) outfile->Get(selection+"sct_ec_elecalisoovere");
  TH1F* ecelhcalisoovere = (TH1F*) outfile->Get(selection+"sct_ec_elhcalisoovere");
  TH1F* eceltkisoovere = (TH1F*) outfile->Get(selection+"sct_ec_eltkisoovere");
  TH1F* ecelr9 = (TH1F*) outfile->Get(selection+"sct_ec_elr9");
  TH1F* ecelsmin = (TH1F*) outfile->Get(selection+"sct_ec_elsmin");
  TH1F* ecelsmaj = (TH1F*) outfile->Get(selection+"sct_ec_elsmaj");

  TH1F* trkecelpt = (TH1F*) outfile->Get(selection+"sct_ec_trk_elpt");
  TH1F* trkeceleta = (TH1F*) outfile->Get(selection+"sct_ec_trk_eleta");
  TH1F* trkecelphi = (TH1F*) outfile->Get(selection+"sct_ec_trk_elphi");
  TH1F* trkeceld0 = (TH1F*) outfile->Get(selection+"sct_ec_trk_eld0");
  TH1F* trkecellog10d0 = (TH1F*) outfile->Get(selection+"sct_ec_trk_ellog10d0");
  TH1F* trkeceldz = (TH1F*) outfile->Get(selection+"sct_ec_trk_eldz");
  TH1F* trkecelcharge = (TH1F*) outfile->Get(selection+"sct_ec_trk_elcharge");
  TH1F* trkecelrchi2 = (TH1F*) outfile->Get(selection+"sct_ec_trk_elrchi2");

  TH1F* trkecmindE = (TH1F*) outfile->Get(selection+"sct_ec_trk_el_min_dE");
  TH1F* trkecmindEoverE = (TH1F*) outfile->Get(selection+"sct_ec_trk_el_min_dEoverE");
  TH1F* trkecmindeta = (TH1F*) outfile->Get(selection+"sct_ec_trk_el_min_deta");
  TH1F* trkecmindphi = (TH1F*) outfile->Get(selection+"sct_ec_trk_el_min_dphi");

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
    //leadsubleadqprod->Fill((*ele_charge)->at(elidx[0])*(*ele_charge)->at(elidx[1]));
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

    TLorentzVector el;
    el.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr]),(*ele_eta)->at(elidx[ctr]),(*ele_phi)->at(elidx[ctr]),0.0005);
    
    for(unsigned int ctr2=ctr+1; ctr2<elidx.size(); ctr2++) {
      TLorentzVector el2;
      el2.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr2]),(*ele_eta)->at(elidx[ctr2]),(*ele_phi)->at(elidx[ctr2]),0.0005);
      dielM->Fill((el+el2).M());
    }

    if(TMath::Abs((*ele_eta)->at(elidx[ctr]))<1.479) {
      barelpt->Fill((*ele_pt)->at(elidx[ctr]));
      bareleta->Fill((*ele_eta)->at(elidx[ctr]));
      barelphi->Fill((*ele_phi)->at(elidx[ctr]));
      barelm->Fill((*ele_m)->at(elidx[ctr]));
      bareldetain->Fill((*ele_detain)->at(elidx[ctr]));
      bareldphiin->Fill((*ele_dphiin)->at(elidx[ctr]));
      barelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]));
      barelhoe->Fill((*ele_hoe)->at(elidx[ctr]));
      barelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]));
      barelmhits->Fill((*ele_mhits)->at(elidx[ctr]));
      barelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]));
      barelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]));
      bareltkiso->Fill((*ele_tkiso)->at(elidx[ctr]));
      barelecalisoovere->Fill( (*ele_ecaliso)->at(elidx[ctr])/el.E() );
      barelhcalisoovere->Fill( (*ele_hcaliso)->at(elidx[ctr])/el.E() );
      bareltkisoovere->Fill( (*ele_tkiso)->at(elidx[ctr])/el.E() );
      barelr9->Fill((*ele_r9)->at(elidx[ctr]));
      barelsmin->Fill((*ele_smin)->at(elidx[ctr]));
      barelsmaj->Fill((*ele_smaj)->at(elidx[ctr]));
    }

    else {
      ecelpt->Fill((*ele_pt)->at(elidx[ctr]));
      eceleta->Fill((*ele_eta)->at(elidx[ctr]));
      ecelphi->Fill((*ele_phi)->at(elidx[ctr]));
      ecelm->Fill((*ele_m)->at(elidx[ctr]));
      eceldetain->Fill((*ele_detain)->at(elidx[ctr]));
      eceldphiin->Fill((*ele_dphiin)->at(elidx[ctr]));
      ecelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]));
      ecelhoe->Fill((*ele_hoe)->at(elidx[ctr]));
      ecelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]));
      ecelmhits->Fill((*ele_mhits)->at(elidx[ctr]));
      ecelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]));
      ecelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]));
      eceltkiso->Fill((*ele_tkiso)->at(elidx[ctr]));
      ecelecalisoovere->Fill( (*ele_ecaliso)->at(elidx[ctr])/el.E() );
      ecelhcalisoovere->Fill( (*ele_hcaliso)->at(elidx[ctr])/el.E() );
      eceltkisoovere->Fill( (*ele_tkiso)->at(elidx[ctr])/el.E() );
      ecelr9->Fill((*ele_r9)->at(elidx[ctr]));
      ecelsmin->Fill((*ele_smin)->at(elidx[ctr]));
      ecelsmaj->Fill((*ele_smaj)->at(elidx[ctr]));
    }

    trkmult->Fill( tkidx[ctr].size() );
    double minde, mindeta, mindphi;
    TLorentzVector trk;
    trk.SetPtEtaPhiM(((*ele_trkpt)->at(elidx[ctr]))[tkidx[ctr][0]], ((*ele_trketa)->at(elidx[ctr]))[tkidx[ctr][0]], ((*ele_trkphi)->at(elidx[ctr]))[tkidx[ctr][0]], 0.0005);
    minde = el.E() - trk.E();
    mindeta = (*ele_eta)->at(elidx[ctr]) - ((*ele_trketa)->at(elidx[ctr]))[tkidx[ctr][0]];
    mindphi = (*ele_phi)->at(elidx[ctr]) - ((*ele_trkphi)->at(elidx[ctr]))[tkidx[ctr][0]];
    // Fill the track variables
    for(unsigned int trkc=0; trkc<tkidx[ctr].size(); trkc++) {
      int trki = tkidx[ctr][trkc];
      if(TMath::Abs((*ele_eta)->at(elidx[ctr]))<1.479) {
	trkbarelpt->Fill(((*ele_trkpt)->at(elidx[ctr]))[trki]);
	trkbareleta->Fill(((*ele_trketa)->at(elidx[ctr]))[trki]);
	trkbarelphi->Fill(((*ele_trkphi)->at(elidx[ctr]))[trki]);
	trkbareld0->Fill(((*ele_trkd0)->at(elidx[ctr]))[trki]);
	trkbarellog10d0->Fill(TMath::Log10(TMath::Abs(((*ele_trkd0)->at(elidx[ctr]))[trki])));
	trkbareldz->Fill(((*ele_trkdz)->at(elidx[ctr]))[trki]);
	trkbarelcharge->Fill(((*ele_trkcharge)->at(elidx[ctr]))[trki]);
	trkbarelrchi2->Fill(((*ele_trkchi2)->at(elidx[ctr]))[trki]);
      }
      else {
	trkecelpt->Fill(((*ele_trkpt)->at(elidx[ctr]))[trki]);
	trkeceleta->Fill(((*ele_trketa)->at(elidx[ctr]))[trki]);
	trkecelphi->Fill(((*ele_trkphi)->at(elidx[ctr]))[trki]);
	trkeceld0->Fill(((*ele_trkd0)->at(elidx[ctr]))[trki]);
	trkecellog10d0->Fill(TMath::Log10(TMath::Abs(((*ele_trkd0)->at(elidx[ctr]))[trki])));
	trkeceldz->Fill(((*ele_trkdz)->at(elidx[ctr]))[trki]);
	trkecelcharge->Fill(((*ele_trkcharge)->at(elidx[ctr]))[trki]);
	trkecelrchi2->Fill(((*ele_trkchi2)->at(elidx[ctr]))[trki]);
      }
      trk.SetPtEtaPhiM(((*ele_trkpt)->at(elidx[ctr]))[trki], ((*ele_trketa)->at(elidx[ctr]))[trki], ((*ele_trkphi)->at(elidx[ctr]))[trki], 0.0005);
      double de = el.E() - trk.E();
      double deta = (*ele_eta)->at(elidx[ctr]) - ((*ele_trketa)->at(elidx[ctr]))[trki];
      double dphi = (*ele_phi)->at(elidx[ctr]) - ((*ele_trkphi)->at(elidx[ctr]))[trki];
      if(abs(de)<minde) {
	minde = de;
      }
      if(abs(deta)<mindeta) {
	mindeta = deta;
      }
      if(abs(dphi)<mindphi) {
	mindphi = dphi;
      }
    } // End of track loop
    if(TMath::Abs((*ele_eta)->at(elidx[ctr]))<1.479) {
      trkbarmindE->Fill(minde);
      trkbarmindEoverE->Fill(minde/el.E());
      trkbarmindeta->Fill(mindeta);
      trkbarmindphi->Fill(mindphi);
    }
    else {
      trkecmindE->Fill(minde);
      trkecmindEoverE->Fill(minde/el.E());
      trkecmindeta->Fill(mindeta);
      trkecmindphi->Fill(mindphi);
    }
    
  } // End of main electron for loop

  // For all track invariant mass
  for(unsigned int sc1=0; sc1<elidx.size(); sc1++) {
    for(unsigned int tk1=0; tk1<tkidx[sc1].size(); tk1++) {
      TLorentzVector trk1;
      int tki1 = tkidx[sc1][tk1];
      trk1.SetPtEtaPhiM(((*ele_trkpt)->at(elidx[sc1]))[tki1], ((*ele_trketa)->at(elidx[sc1]))[tki1], ((*ele_trkphi)->at(elidx[sc1]))[tki1], 0.0005);
      for(unsigned int sc2=sc1+1; sc2<elidx.size(); sc2++) {
	for(unsigned int tk2=0; tk2<tkidx[sc2].size(); tk2++) {
	  TLorentzVector trk2;
	  int tki2 = tkidx[sc2][tk2];
	  trk2.SetPtEtaPhiM(((*ele_trkpt)->at(elidx[sc2]))[tki2], ((*ele_trketa)->at(elidx[sc2]))[tki2], ((*ele_trkphi)->at(elidx[sc2]))[tki2], 0.0005);
	  trkdielM->Fill((trk1+trk2).M());
	}
      }
    }
  }
    
}  

// Function to add a set of histograms for gen electrons
void robustanalyzer::addgenhist(TString selection) {

  if(!isMC) throw "Error!!! Trying to define gen hists for a non MC file.";

  all1dhists.push_back(new TH1F(selection+"evt_isomu27","HLT_IsoMu27",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_mu50","HLT_Mu50",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_scouting","DST_Run3PixelOnlyPFScouting",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_sctdimu3","DST_Scouting_DoubleMu3",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_sctdieg","DST_Scouting_EG16EG12",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_scteg30","DST_Scouting_EG30",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_sctjet","DST_Scouting_JetHT",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_musct","DST_HLTMuon_Run3PFScouting",10,-5,5));
  all1dhists.push_back(new TH1F(selection+"evt_othsctmnt","HLT_OtherScoutingPFMonitor",10,-5,5));

  all1dhists.push_back(new TH1F(selection+"gen_jpsi_mult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"gen_jpsi_m","M / GeV",1500,-1,14));
  all1dhists.push_back(new TH1F(selection+"gen_jpsi_pt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_jpsi_eta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_jpsi_phi","#phi",66,-3.3,3.3));
  
  all1dhists.push_back(new TH1F(selection+"gen_el_mult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"gen_el_pt","e p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_el_eta","e #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_el_phi","e #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_leadel_pt","e_{1} p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_leadel_eta","e_{1} #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_leadel_phi","e_{1} #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_subleadel_pt","e_{2} p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_subleadel_eta","e_{2} #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"gen_subleadel_phi","e_{2} #phi",66,-3.3,3.3));

  all1dhists.push_back(new TH1F(selection+"gen_diel_M","all M(e_{1},e_{2})",1000,-10,99));
  all1dhists.push_back(new TH1F(selection+"gen_diel_deta","#Delta#eta(e_{1}, e_{2})",1400,-7,7));
  all1dhists.push_back(new TH1F(selection+"gen_diel_dphi","#Delta#phi(e_{1}, e_{2})",1400,-7,7));
  all1dhists.push_back(new TH1F(selection+"gen_diel_dR","#Delta R(e_{1}, e_{2})",1000,0,10));

}

// Function to fill a set of histograms for gen particles
void robustanalyzer::fillgenhistinevent(TString selection, vector<int> genelidx, vector<int> genjpidx) {

  if(!isMC) throw "Error!!! Trying to fill gen info. for a non MC file.";
    
  if(genjpidx.size()==0) return;
  
  TH1F* isomu27 = (TH1F*) outfile->Get(selection+"evt_isomu27");
  TH1F* mu50 = (TH1F*) outfile->Get(selection+"evt_mu50");
  TH1F* scouting = (TH1F*) outfile->Get(selection+"evt_scouting");
  TH1F* dimu3 = (TH1F*) outfile->Get(selection+"evt_sctdimu3");
  TH1F* dieg = (TH1F*) outfile->Get(selection+"evt_sctdieg");
  TH1F* eg30 = (TH1F*) outfile->Get(selection+"evt_scteg30");
  TH1F* jet = (TH1F*) outfile->Get(selection+"evt_sctjet");
  TH1F* musct = (TH1F*) outfile->Get(selection+"evt_musct");
  TH1F* othsctmnt = (TH1F*) outfile->Get(selection+"evt_othsctmnt");

  TH1F* jpsimult = (TH1F*) outfile->Get(selection+"gen_jpsi_mult");
  TH1F* jpsim = (TH1F*) outfile->Get(selection+"gen_jpsi_m");
  TH1F* jpsipt = (TH1F*) outfile->Get(selection+"gen_jpsi_pt");
  TH1F* jpsieta = (TH1F*) outfile->Get(selection+"gen_jpsi_eta");
  TH1F* jpsiphi = (TH1F*) outfile->Get(selection+"gen_jpsi_phi");
  
  TH1F* elmult = (TH1F*) outfile->Get(selection+"gen_el_mult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"gen_el_pt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"gen_el_eta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"gen_el_phi");

  TH1F* leadelpt = (TH1F*) outfile->Get(selection+"gen_leadel_pt");
  TH1F* leadeleta = (TH1F*) outfile->Get(selection+"gen_leadel_eta");
  TH1F* leadelphi = (TH1F*) outfile->Get(selection+"gen_leadel_phi");

  TH1F* subleadelpt = (TH1F*) outfile->Get(selection+"gen_subleadel_pt");
  TH1F* subleadeleta = (TH1F*) outfile->Get(selection+"gen_subleadel_eta");
  TH1F* subleadelphi = (TH1F*) outfile->Get(selection+"gen_subleadel_phi");

  TH1F* dielM = (TH1F*) outfile->Get(selection+"gen_diel_M");
  TH1F* dieldeta = (TH1F*) outfile->Get(selection+"gen_diel_deta");
  TH1F* dieldphi = (TH1F*) outfile->Get(selection+"gen_diel_dphi");
  TH1F* dieldR = (TH1F*) outfile->Get(selection+"gen_diel_dR");

  isomu27->Fill( (*(*hlt_isomu27)) );
  mu50->Fill( (*(*hlt_mu50)) );
  scouting->Fill( (*(*hlt_scouting)) );
  dimu3->Fill( (*(*hlt_sct_dimu3)) );
  dieg->Fill( (*(*hlt_sct_dieg)) );
  eg30->Fill( (*(*hlt_sct_eg30)) );
  jet->Fill( (*(*hlt_sct_jet)) );
  musct->Fill( (*(*hltmu_sct)) );
  othsctmnt->Fill( (*(*hltoth_sctmnt)) );
  
  jpsimult->Fill(genjpidx.size());
  for(int jp : genjpidx) {
    jpsim->Fill((*genjpsi_m)->at(jp));
    jpsipt->Fill((*genjpsi_pt)->at(jp));
    jpsieta->Fill((*genjpsi_eta)->at(jp));
    jpsiphi->Fill((*genjpsi_phi)->at(jp));
  }
  
  // Assume the daughters of JPsi are stored sequentially.
  
  // Throw error if uneven number of electron indices
  if(genelidx.size()%2 != 0) throw "Error!!! Uneven number of daughter gen electrons. Assuming two gen electrons per JPsi.";
  
  elmult->Fill(genelidx.size());
  for(unsigned int ctr=0; ctr<genelidx.size(); ctr=ctr+2) {
    
    // Two electrons of opposite charge
    if((*genlep_pdg)->at(genelidx[ctr])*(*genlep_pdg)->at(genelidx[ctr+1])!=-121) throw "Error!!! Assumed sequential storage of JPsi daughters. Should be oppositely charged.";
    
    unsigned int leadelptpos=-1, subleadelptpos=-1;

    elpt->Fill((*genlep_pt)->at(genelidx[ctr]));
    eleta->Fill((*genlep_eta)->at(genelidx[ctr]));
    elphi->Fill((*genlep_phi)->at(genelidx[ctr]));
    elpt->Fill((*genlep_pt)->at(genelidx[ctr+1]));
    eleta->Fill((*genlep_eta)->at(genelidx[ctr+1]));
    elphi->Fill((*genlep_phi)->at(genelidx[ctr+1]));

    TLorentzVector el1, el2;
    el1.SetPtEtaPhiM((*genlep_pt)->at(genelidx[ctr]),(*genlep_eta)->at(genelidx[ctr]),(*genlep_phi)->at(genelidx[ctr]),0.0005);
    el2.SetPtEtaPhiM((*genlep_pt)->at(genelidx[ctr+1]),(*genlep_eta)->at(genelidx[ctr+1]),(*genlep_phi)->at(genelidx[ctr+1]),0.0005);
    dielM->Fill((el1+el2).M());
    dieldeta->Fill((*genlep_eta)->at(genelidx[ctr])-(*genlep_eta)->at(genelidx[ctr]));
    dieldphi->Fill(el1.DeltaPhi(el2));
    dieldR->Fill(el1.DeltaR(el2));

    if( ((*genlep_pt)->at(genelidx[ctr])) > ((*genlep_pt)->at(genelidx[ctr+1])) ) {
      leadelptpos = genelidx[ctr];
      subleadelptpos = genelidx[ctr+1];
    }      
    else {
      leadelptpos = genelidx[ctr+1];
      subleadelptpos = genelidx[ctr];
    }
    leadelpt->Fill((*genlep_pt)->at(leadelptpos));
    leadeleta->Fill((*genlep_eta)->at(leadelptpos));
    leadelphi->Fill((*genlep_phi)->at(leadelptpos));
    subleadelpt->Fill((*genlep_pt)->at(subleadelptpos));
    subleadeleta->Fill((*genlep_eta)->at(subleadelptpos));
    subleadelphi->Fill((*genlep_phi)->at(subleadelptpos));
    
  }
}

// Function to add a set of histograms for gen matching scouting electrons
void robustanalyzer::addgenmchhist(TString selection) {

  // Variables before gen match
  all1dhists.push_back(new TH1F(selection+"genelsct_eb_dE","#Delta E(gen e, sct. e)",10000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"genelsct_eb_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsct_eb_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsct_ee_dE","#Delta E(gen e, sct. e)",10000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"genelsct_ee_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsct_ee_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));

  // Variables before gen match with tracks
  all1dhists.push_back(new TH1F(selection+"genelsct_eb_trk_dE","trk. #Delta E(gen e, sct. e)",10000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"genelsct_eb_trk_dEta","trk. #Delta#eta(gen e, sct. e)",10000,-0.5,0.5));
  all1dhists.push_back(new TH1F(selection+"genelsct_eb_trk_qdPhi","trk. q#Delta#phi(gen e, sct. e)",10000,-0.5,0.5));
  all1dhists.push_back(new TH1F(selection+"genelsct_ee_trk_dE","trk. #Delta E(gen e, sct. e)",10000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"genelsct_ee_trk_dEta","trk. #Delta#eta(gen e, sct. e)",10000,-0.5,0.5));
  all1dhists.push_back(new TH1F(selection+"genelsct_ee_trk_qdPhi","trk. q#Delta#phi(gen e, sct. e)",10000,-0.5,0.5));

  // Variables after gen match
  all1dhists.push_back(new TH1F(selection+"genelsctmch_eb_dE","#Delta E(gen e, sct. e)",10000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_eb_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_eb_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_ee_dE","#Delta E(gen e, sct. e)",10000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_ee_dEta","#Delta#eta(gen e, sct. e)",10000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_ee_qdPhi","q#Delta#phi(gen e, sct. e)",10000,-5,5));

  // Variables after gen match with tracks
  all1dhists.push_back(new TH1F(selection+"genelsctmch_eb_trk_dE","trk. #Delta E(gen e, sct. e)",10000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_eb_trk_dEta","trk. #Delta#eta(gen e, sct. e)",10000,-0.5,0.5));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_eb_trk_qdPhi","trk. q#Delta#phi(gen e, sct. e)",10000,-0.5,0.5));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_ee_trk_dE","trk. #Delta E(gen e, sct. e)",10000,-1000,1000));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_ee_trk_dEta","trk. #Delta#eta(gen e, sct. e)",10000,-0.5,0.5));
  all1dhists.push_back(new TH1F(selection+"genelsctmch_ee_trk_qdPhi","trk. q#Delta#phi(gen e, sct. e)",10000,-0.5,0.5));

}

// Function to fill a set of histograms for genmatching particles
void robustanalyzer::fillgenmchhistinevent(TString selection, vector<int> genidx, vector<int> sctidx) {

  if(!isMC) throw "Error!!! Trying to obtain gen info. for a non MC file.";
    
  if(genidx.size()==0 || sctidx.size()==0) return;
  
  // Variables before gen match
  TH1F* eb_dE = (TH1F*) outfile->Get(selection+"genelsct_eb_dE");
  TH1F* eb_dEta = (TH1F*) outfile->Get(selection+"genelsct_eb_dEta");
  TH1F* eb_qdPhi = (TH1F*) outfile->Get(selection+"genelsct_eb_qdPhi");
  TH1F* ee_dE = (TH1F*) outfile->Get(selection+"genelsct_ee_dE");
  TH1F* ee_dEta = (TH1F*) outfile->Get(selection+"genelsct_ee_dEta");
  TH1F* ee_qdPhi = (TH1F*) outfile->Get(selection+"genelsct_ee_qdPhi");
  
  // Variables before gen match with tracks
  TH1F* trkeb_dE = (TH1F*) outfile->Get(selection+"genelsct_eb_trk_dE");
  TH1F* trkeb_dEta = (TH1F*) outfile->Get(selection+"genelsct_eb_trk_dEta");
  TH1F* trkeb_qdPhi = (TH1F*) outfile->Get(selection+"genelsct_eb_trk_qdPhi");
  TH1F* trkee_dE = (TH1F*) outfile->Get(selection+"genelsct_ee_trk_dE");
  TH1F* trkee_dEta = (TH1F*) outfile->Get(selection+"genelsct_ee_trk_dEta");
  TH1F* trkee_qdPhi = (TH1F*) outfile->Get(selection+"genelsct_ee_trk_qdPhi");
  
  // Variables after gen match
  TH1F* mcheb_dE = (TH1F*) outfile->Get(selection+"genelsctmch_eb_dE");
  TH1F* mcheb_dEta = (TH1F*) outfile->Get(selection+"genelsctmch_eb_dEta");
  TH1F* mcheb_qdPhi = (TH1F*) outfile->Get(selection+"genelsctmch_eb_qdPhi");
  TH1F* mchee_dE = (TH1F*) outfile->Get(selection+"genelsctmch_ee_dE");
  TH1F* mchee_dEta = (TH1F*) outfile->Get(selection+"genelsctmch_ee_dEta");
  TH1F* mchee_qdPhi = (TH1F*) outfile->Get(selection+"genelsctmch_ee_qdPhi");

  // Variables after gen match with tracks
  TH1F* trkmcheb_dE = (TH1F*) outfile->Get(selection+"genelsctmch_eb_trk_dE");
  TH1F* trkmcheb_dEta = (TH1F*) outfile->Get(selection+"genelsctmch_eb_trk_dEta");
  TH1F* trkmcheb_qdPhi = (TH1F*) outfile->Get(selection+"genelsctmch_eb_trk_qdPhi");
  TH1F* trkmchee_dE = (TH1F*) outfile->Get(selection+"genelsctmch_ee_trk_dE");
  TH1F* trkmchee_dEta = (TH1F*) outfile->Get(selection+"genelsctmch_ee_trk_dEta");
  TH1F* trkmchee_qdPhi = (TH1F*) outfile->Get(selection+"genelsctmch_ee_trk_qdPhi");

  for(int ge : genidx) {
    TLorentzVector gene, scte, trke;
    gene.SetPtEtaPhiM((*genlep_pt)->at(ge), (*genlep_eta)->at(ge), (*genlep_phi)->at(ge), 0.0005);
    scte.SetPtEtaPhiM((*ele_pt)->at(sctidx[0]), (*ele_eta)->at(sctidx[0]), (*ele_phi)->at(sctidx[0]), 0.0005);
    trke.SetPtEtaPhiM((*ele_trkpt)->at(sctidx[0])[0], (*ele_trketa)->at(sctidx[0])[0], (*ele_trkphi)->at(sctidx[0])[0], 0.0005);

    double mdene = gene.E() - scte.E();
    double mdeta = ((*genlep_eta)->at(ge)) - ((*ele_eta)->at(sctidx[0]));
    double mdphi = ((*genlep_phi)->at(ge)) - ((*ele_phi)->at(sctidx[0]));
    double mqdph = (charge((*genlep_pdg)->at(ge))) * mdphi;
    double eta = (*ele_eta)->at(sctidx[0]);

    double mtrkdene = gene.E() - trke.E();
    double mtrkdeta = ((*genlep_eta)->at(ge)) - ((*ele_trketa)->at(sctidx[0]))[0];
    double mtrkdphi = ((*genlep_phi)->at(ge)) - ((*ele_trkphi)->at(sctidx[0]))[0];
    double mtrkqdph = ((*ele_trkcharge)->at(sctidx[0]))[0] * mtrkdphi;

    for(int se : sctidx) {
      scte.SetPtEtaPhiM((*ele_pt)->at(se), (*ele_eta)->at(se), (*ele_phi)->at(se), 0.0005);
      double dene = gene.E() - scte.E();
      double deta = ((*genlep_eta)->at(ge)) - ((*ele_eta)->at(se));
      double dphi = ((*genlep_phi)->at(ge)) - ((*ele_phi)->at(se));
      double qdph = (charge((*genlep_pdg)->at(ge))) * dphi;      
      if(abs(dene)<mdene) {
	mdene = dene;
      }
      if(abs(deta)<mdeta) {
	mdeta = deta;
	eta = (*ele_eta)->at(se);
      }
      if(abs(qdph)<mqdph) {
	mqdph = qdph;
      }

      int ntrk = ((*ele_trkpt)->at(se)).size();
      for(unsigned int tc=0; tc<ntrk; tc++) {
	trke.SetPtEtaPhiM((*ele_trkpt)->at(se)[tc], (*ele_trketa)->at(se)[tc], (*ele_trkphi)->at(se)[tc], 0.0005);
	double tdene = gene.E() - trke.E();
	double tdeta = ((*genlep_eta)->at(ge)) - ((*ele_trketa)->at(sctidx[se]))[tc];
	double tdphi = ((*genlep_phi)->at(ge)) - ((*ele_trkphi)->at(sctidx[se]))[tc];
	double tqdph = ((*ele_trkcharge)->at(sctidx[se]))[tc] * tdphi;

	if(abs(tdene)<mtrkdene) {
	  mtrkdene = tdene;
	}
	if(abs(tdeta)<mtrkdeta) {
	  mtrkdeta = tdeta;
	}
	if(abs(tqdph)<mtrkqdph) {
	  mtrkqdph = tqdph;
	}
      } // End of loop on track
      
    } // End of loop on SC
    if(abs(eta)<1.479) {
      eb_dE->Fill(mdene);
      eb_dEta->Fill(mdeta);
      eb_qdPhi->Fill(mqdph);
      trkeb_dE->Fill(mtrkdene);
      trkeb_dEta->Fill(mtrkdeta);
      trkeb_qdPhi->Fill(mtrkqdph);
    }
    else {
      ee_dE->Fill(mdene);
      ee_dEta->Fill(mdeta);
      ee_qdPhi->Fill(mqdph);
      trkee_dE->Fill(mtrkdene);
      trkee_dEta->Fill(mtrkdeta);
      trkee_qdPhi->Fill(mtrkqdph);
    }
  }

  // Get the gen-matched indices for SC
  std::pair< vector<int>,vector<int> > gensctpair = getGenMatched(genidx, sctidx);
  vector<int> matchedgenidx = gensctpair.first;
  vector<int> matchedsctidx = gensctpair.second;

  // Check for duplications
  for(unsigned int ctr1=0; ctr1<matchedsctidx.size(); ctr1++) {
    for(unsigned int ctr2=ctr1+1; ctr2<matchedsctidx.size(); ctr2++) {
      if(matchedsctidx[ctr1]==matchedsctidx[ctr2]) throw "Error!!! Faulty gen matching. Duplicated Scouting Index.";
    }
  }  
  
  for(unsigned int ctr=0; ctr<matchedsctidx.size(); ctr++) {
    TLorentzVector gene, scte;
    gene.SetPtEtaPhiM((*genlep_pt)->at(matchedgenidx[ctr]), (*genlep_eta)->at(matchedgenidx[ctr]), (*genlep_phi)->at(matchedgenidx[ctr]), 0.0005);
    scte.SetPtEtaPhiM((*ele_pt)->at(matchedsctidx[ctr]), (*ele_eta)->at(matchedsctidx[ctr]), (*ele_phi)->at(matchedsctidx[ctr]), 0.0005);
    double dene = gene.E() - scte.E();
    double deta = ((*genlep_eta)->at(matchedgenidx[ctr])) - ((*ele_eta)->at(matchedsctidx[ctr]));
    double dphi = ((*genlep_phi)->at(matchedgenidx[ctr])) - ((*ele_phi)->at(matchedsctidx[ctr]));
    double qdph = (charge((*genlep_pdg)->at(matchedgenidx[ctr]))) * dphi;
    double eta = (*ele_eta)->at(matchedsctidx[ctr]);
    if(abs(eta)<1.479) {
      mcheb_dE->Fill(dene);
      mcheb_dEta->Fill(deta);
      mcheb_qdPhi->Fill(qdph);
    }
    else {
      mchee_dE->Fill(dene);
      mchee_dEta->Fill(deta);
      mchee_qdPhi->Fill(qdph);
    }
  }

  // Get the gen-matched indices for trk
  std::pair< vector<int>,vector< std::pair<int,int> > > gentrkpair = getGenMatchedTrk(genidx, sctidx);
  vector<int> matchedgentrkidx = gentrkpair.first;
  vector< pair<int,int> > matchedtrkidx = gentrkpair.second;
  
  // Check for duplications
  for(unsigned int ctr1=0; ctr1<matchedtrkidx.size(); ctr1++) {
    for(unsigned int ctr2=ctr1+1; ctr2<matchedtrkidx.size(); ctr2++) {
      if(matchedtrkidx[ctr1].first==matchedtrkidx[ctr2].first) throw "Error!!! Faulty gen matching. Duplicated Scouting Index.";
    }
  }  
  
  for(unsigned int ctr=0; ctr<matchedtrkidx.size(); ctr++) {
    int ti = matchedtrkidx[ctr].second;
    int si = matchedtrkidx[ctr].first;
    TLorentzVector gene, trke;
    gene.SetPtEtaPhiM((*genlep_pt)->at(matchedgentrkidx[ctr]), (*genlep_eta)->at(matchedgentrkidx[ctr]), (*genlep_phi)->at(matchedgentrkidx[ctr]), 0.0005);
    trke.SetPtEtaPhiM(((*ele_trkpt)->at(si))[ti], ((*ele_trketa)->at(si))[ti], ((*ele_trkphi)->at(si))[ti], 0.0005);
    double tdene = gene.E() - trke.E();
    double tdeta = ((*genlep_eta)->at(matchedgentrkidx[ctr])) - ((*ele_trketa)->at(si))[ti];
    double tdphi = ((*genlep_phi)->at(matchedgentrkidx[ctr])) - ((*ele_trkphi)->at(si))[ti];
    double tqdph = ((*ele_trkcharge)->at(si))[ti] * tdphi;
    double eta = (*ele_eta)->at(si);
    if(abs(eta)<1.479) {
      trkmcheb_dE->Fill(tdene);
      trkmcheb_dEta->Fill(tdeta);
      trkmcheb_qdPhi->Fill(tqdph);
    }
    else {
      trkmchee_dE->Fill(tdene);
      trkmchee_dEta->Fill(tdeta);
      trkmchee_qdPhi->Fill(tqdph);
    }
  }

}
