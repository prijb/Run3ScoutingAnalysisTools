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
    addgenhist("noselgen_othr_");
  }
  
  // Count events passing certain selections
  int nosel=0;
  
  addhist("nosel");
  addhist("mutrigsel");
  addhist("muAscouttrigsel");

  // vector of gen indices
  vector<int> noselgenjpidx;
  vector<int> noselgenelidx;

  // vector of electron indices
  vector<int> noselelidx;
  vector<int> mutrigselelidx;
  vector<int> muAscouttrigselelidx;
    
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%100000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

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
      if( (noselgenjpidx.size()==1) && ((*(*hlt_sct_dimu3))==1) && ((*(*hlt_sct_jet))==1) ) fillgenhistinevent("noselgen_othr_", noselgenelidx, noselgenjpidx);
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
      
      // Get the energy of this electron
      TLorentzVector ele;
      ele.SetPtEtaPhiM((*ele_pt)->at(elidx), (*ele_eta)->at(elidx), (*ele_phi)->at(elidx), 0.0005);
      double ele_energy = ele.Energy();

      bool mutrigdec = true;
      mutrigdec *= ( ((*(*hlt_isomu27))==1) || ((*(*hlt_mu50))==1) );
      mutrigdec *= (*ele_pt)->at(elidx)>2.0;
      mutrigdec *= abs((*ele_eta)->at(elidx))<2.5;
      //mutrigdec *= 
      
      if( ( ((*(*hlt_isomu27))==1) || ((*(*hlt_mu50))==1) ) && ((*(*hlt_scouting))==1) ) muAscouttrigselelidx.push_back(elidx);
      
    }// End of loop on electrons in the event loop
    
    if(noselelidx.size()>0) nosel++;
    fillhistinevent("nosel", noselelidx);
    fillhistinevent("mutrigsel", mutrigselelidx);
    fillhistinevent("muAscouttrigsel", muAscouttrigselelidx);

    // Clear all vector
    noselelidx.clear();
    mutrigselelidx.clear();
    muAscouttrigselelidx.clear();
    noselgenjpidx.clear();
    noselgenelidx.clear();
    
  } // End of loop on events

  cout<<totEntries<<"\t"<<endevent-beginevent<<"\t"<<nosel<<endl;
}

// Function to fill a set of histograms for scouting electrons
void robustanalyzer::fillhistinevent(TString selection, vector<int> elidx) {

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
  TH1F* bareldetain = (TH1F*) outfile->Get(selection+"sctbar_eldetain");
  TH1F* bareldphiin = (TH1F*) outfile->Get(selection+"sctbar_eldphiin");
  TH1F* barelsigmaietaieta = (TH1F*) outfile->Get(selection+"sctbar_elsigmaietaieta");
  TH1F* barelhoe = (TH1F*) outfile->Get(selection+"sctbar_elhoe");
  TH1F* barelooemoop = (TH1F*) outfile->Get(selection+"sctbar_elooemoop");
  TH1F* barelmhits = (TH1F*) outfile->Get(selection+"sctbar_elmhits");
  TH1F* barelecaliso = (TH1F*) outfile->Get(selection+"sctbar_elecaliso");
  TH1F* barelhcaliso = (TH1F*) outfile->Get(selection+"sctbar_elhcaliso");
  TH1F* bareltkiso = (TH1F*) outfile->Get(selection+"sctbar_eltkiso");
  TH1F* barelr9 = (TH1F*) outfile->Get(selection+"sctbar_elr9");
  TH1F* barelsmin = (TH1F*) outfile->Get(selection+"sctbar_elsmin");
  TH1F* barelsmaj = (TH1F*) outfile->Get(selection+"sctbar_elsmaj");
  
  TH1F* bareld0 = (TH1F*) outfile->Get(selection+"sctbar_eld0");
  TH1F* barellog10d0 = (TH1F*) outfile->Get(selection+"sctbar_ellog10d0");
  TH1F* bareldz = (TH1F*) outfile->Get(selection+"sctbar_eldz");
  TH1F* barelcharge = (TH1F*) outfile->Get(selection+"sctbar_elcharge");

  TH1F* ecelpt = (TH1F*) outfile->Get(selection+"sctec_elpt");
  TH1F* eceleta = (TH1F*) outfile->Get(selection+"sctec_eleta");
  TH1F* ecelphi = (TH1F*) outfile->Get(selection+"sctec_elphi");
  TH1F* ecelm = (TH1F*) outfile->Get(selection+"sctec_elm");
  TH1F* eceldetain = (TH1F*) outfile->Get(selection+"sctec_eldetain");
  TH1F* eceldphiin = (TH1F*) outfile->Get(selection+"sctec_eldphiin");
  TH1F* ecelsigmaietaieta = (TH1F*) outfile->Get(selection+"sctec_elsigmaietaieta");
  TH1F* ecelhoe = (TH1F*) outfile->Get(selection+"sctec_elhoe");
  TH1F* ecelooemoop = (TH1F*) outfile->Get(selection+"sctec_elooemoop");
  TH1F* ecelmhits = (TH1F*) outfile->Get(selection+"sctec_elmhits");
  TH1F* ecelecaliso = (TH1F*) outfile->Get(selection+"sctec_elecaliso");
  TH1F* ecelhcaliso = (TH1F*) outfile->Get(selection+"sctec_elhcaliso");
  TH1F* eceltkiso = (TH1F*) outfile->Get(selection+"sctec_eltkiso");
  TH1F* ecelr9 = (TH1F*) outfile->Get(selection+"sctec_elr9");
  TH1F* ecelsmin = (TH1F*) outfile->Get(selection+"sctec_elsmin");
  TH1F* ecelsmaj = (TH1F*) outfile->Get(selection+"sctec_elsmaj");

  TH1F* eceld0 = (TH1F*) outfile->Get(selection+"sctec_eld0");
  TH1F* ecellog10d0 = (TH1F*) outfile->Get(selection+"sctec_ellog10d0");
  TH1F* eceldz = (TH1F*) outfile->Get(selection+"sctec_eldz");
  TH1F* ecelcharge = (TH1F*) outfile->Get(selection+"sctec_elcharge");

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
      //bareld0->Fill((*ele_d0)->at(elidx[ctr]));
      //barellog10d0->Fill(TMath::Log10(TMath::Abs((*ele_d0)->at(elidx[ctr]))));
      //bareldz->Fill((*ele_dz)->at(elidx[ctr]));
      bareldetain->Fill((*ele_detain)->at(elidx[ctr]));
      bareldphiin->Fill((*ele_dphiin)->at(elidx[ctr]));
      barelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]));
      barelhoe->Fill((*ele_hoe)->at(elidx[ctr]));
      barelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]));
      barelmhits->Fill((*ele_mhits)->at(elidx[ctr]));
      //barelcharge->Fill((*ele_charge)->at(elidx[ctr]));
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
      //eceld0->Fill((*ele_d0)->at(elidx[ctr]));
      //ecellog10d0->Fill(TMath::Log10(TMath::Abs((*ele_d0)->at(elidx[ctr]))));
      //eceldz->Fill((*ele_dz)->at(elidx[ctr]));
      eceldetain->Fill((*ele_detain)->at(elidx[ctr]));
      eceldphiin->Fill((*ele_dphiin)->at(elidx[ctr]));
      ecelsigmaietaieta->Fill((*ele_sigmaietaieta)->at(elidx[ctr]));
      ecelhoe->Fill((*ele_hoe)->at(elidx[ctr]));
      ecelooemoop->Fill((*ele_ooemoop)->at(elidx[ctr]));
      ecelmhits->Fill((*ele_mhits)->at(elidx[ctr]));
      //ecelcharge->Fill((*ele_charge)->at(elidx[ctr]));
      ecelecaliso->Fill((*ele_ecaliso)->at(elidx[ctr]));
      ecelhcaliso->Fill((*ele_hcaliso)->at(elidx[ctr]));
      eceltkiso->Fill((*ele_tkiso)->at(elidx[ctr]));
      ecelr9->Fill((*ele_r9)->at(elidx[ctr]));
      ecelsmin->Fill((*ele_smin)->at(elidx[ctr]));
      ecelsmaj->Fill((*ele_smaj)->at(elidx[ctr]));
    }
  } // End of main electron for loop

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

  all1dhists.push_back(new TH1F(selection+"sctbar_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctbar_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctbar_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldetain","#Delta#eta_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldphiin","#Delta#phi_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sctbar_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctbar_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sctbar_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctbar_elcharge","charge",5,-2,3));

  all1dhists.push_back(new TH1F(selection+"sctec_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctec_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctec_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctec_eldetain","#Delta#eta_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_eldphiin","#Delta#phi_{in}",200000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctec_elhoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctec_elooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sctec_elmhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sctec_elecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_elhcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_eltkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sctec_elr9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sctec_elsmin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sctec_elsmaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sctec_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctec_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctec_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctec_elcharge","charge",5,-2,3));
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
