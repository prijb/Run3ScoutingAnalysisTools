#include "robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
robustanalyzer::robustanalyzer(TString filename, TString outfilename, int numCores, bool isDrellYan, bool isMonteCarlo, bool isJPsiMC){

  nC = numCores;
  isDY = isDrellYan;
  isMC = isMonteCarlo;
  isJP = isJPsiMC;

  cout<<"Initializing for file: "<<filename<<endl;
  TChain* chain = new TChain("mmtree/tree");
  chain->Add(filename);

  tree = new TTreeReader(chain);
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
  n_ele = new TTreeReaderValue<UInt_t>((*tree), "n_ele");
  ele_pt = new TTreeReaderValue<vector<float>>((*tree), "Electron_pt");
  ele_eta = new TTreeReaderValue<vector<float>>((*tree), "Electron_eta");
  ele_phi = new TTreeReaderValue<vector<float>>((*tree), "Electron_phi");
  ele_m = new TTreeReaderValue<vector<float>>((*tree), "Electron_m");
  //ele_trkpt = new TTreeReaderValue<vector<float>>((*tree), "Electron_trkpt");
  //ele_trketa = new TTreeReaderValue<vector<float>>((*tree), "Electron_trketa");
  //ele_trkphi = new TTreeReaderValue<vector<float>>((*tree), "Electron_trkphi");
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

  rho = new TTreeReaderValue<vector<float>>((*tree), "rho");

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
  int nosel=0, vetoidsel=0, looseidsel=0, mediumidsel=0, tightidsel=0, tightnoneidsel=0, nonetightidsel=0, tightlooseidsel=0, loosetightidsel=0, tightmediumidesl=0, mediumtightidsel=0;

  // Define the histograms
  if(isMC) {
    addgenhist("noselgen_");
    addgenmchhist("noselgenAnosel_");
  }
  addhist("nosel_");
  addhist("loosesel_");
  addhist("mediumsel_");

  // vector of electron indices
  vector<int> noselgenjpidx;
  vector<int> noselgenelidx;
  vector<int> noselelidx;
  vector<int> looseselelidx;
  vector<int> mediumselelidx;

  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%100000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    double evtrho = (*rho)->at(0);

    bool noselcond = false;
    // Loop on Gen particles to select good gen electrons
    if(isMC) {

      // Check for 2 daughters per JPsi
      if((*(*n_genlep))/2!=(*(*n_genjpsi))) throw "Expecting two daughters for every JPsi";
      
      for(unsigned int gen_ctr=0; gen_ctr<(*(*n_genjpsi)); gen_ctr++) {

	noselcond = true;
	if(noselcond) {
	  noselgenjpidx.push_back(gen_ctr);
	  noselgenelidx.push_back(2*gen_ctr);
	  noselgenelidx.push_back(2*gen_ctr+1);
	}
	
      } // End of loop on gen electrons
      if(noselgenjpidx.size()==1) fillgenhistinevent("noselgen_", noselgenelidx, noselgenjpidx);
    }

    bool looseselcond_ = false;
    bool mediumselcond_ = false;
    
    // Loop on electrons in the event loop
    for(unsigned int ele_ctr=0; ele_ctr<(*(*n_ele)); ele_ctr++) {

      // Get the energy of the electron
      TLorentzVector el;
      el.SetPtEtaPhiM((*ele_pt)->at(ele_ctr),(*ele_eta)->at(ele_ctr),(*ele_phi)->at(ele_ctr),0.0005);
      double ele_energy = el.E();

      // Get the isolation with effective area of the electron
      
      noselelidx.push_back(ele_ctr);

      looseselcond_ = true;
      looseselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_sigmaietaieta)->at(ele_ctr)<0.0107):((*ele_sigmaietaieta)->at(ele_ctr)<0.0275);
      looseselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_detain)->at(ele_ctr)<0.00691):((*ele_detain)->at(ele_ctr)<0.0121);
      looseselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_dphiin)->at(ele_ctr)<0.175):((*ele_dphiin)->at(ele_ctr)<0.228);
      looseselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_hoe)->at(ele_ctr)<(0.05+(1.28/ele_energy)+(0.0422*evtrho/ele_energy)))
	:((*ele_hoe)->at(ele_ctr)<(0.05+(2.3/ele_energy)+(0.262*evtrho/ele_energy)));
      //looseselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_tkiso)->at(ele_ctr)<(0.0478+(0.506/(*ele_pt)->at(ele_ctr))))
      //:((*ele_tkiso)->at(ele_ctr)<(0.0658+(0.963/(*ele_pt)->at(ele_ctr))));
      looseselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_ooemoop)->at(ele_ctr)<0.138):((*ele_ooemoop)->at(ele_ctr)<0.127);
      //looseselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_mhits)->at(ele_ctr)<=1):((*ele_mhits)->at(ele_ctr)<=1);
      if(looseselcond_) looseselelidx.push_back(ele_ctr);
      
      mediumselcond_ = true;
      mediumselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_sigmaietaieta)->at(ele_ctr)<0.0103):((*ele_sigmaietaieta)->at(ele_ctr)<0.0272);
      mediumselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_detain)->at(ele_ctr)<0.00481):((*ele_detain)->at(ele_ctr)<0.00951);
      mediumselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_dphiin)->at(ele_ctr)<0.127):((*ele_dphiin)->at(ele_ctr)<0.221);
      mediumselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_hoe)->at(ele_ctr)<(0.0241+(1.28/ele_energy)+(0.0422*evtrho/ele_energy)))
	:((*ele_hoe)->at(ele_ctr)<(0.05+(2.3/ele_energy)+(0.262*evtrho/ele_energy)));
      //mediumselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_tkiso)->at(ele_ctr)<(0.0478+(0.506/(*ele_pt)->at(ele_ctr))))
      //:((*ele_tkiso)->at(ele_ctr)<(0.0658+(0.963/(*ele_pt)->at(ele_ctr))));
      mediumselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_ooemoop)->at(ele_ctr)<0.0966):((*ele_ooemoop)->at(ele_ctr)<0.111);
      //mediumselcond_ *= (abs((*ele_eta)->at(ele_ctr))<1.479)?((*ele_mhits)->at(ele_ctr)<=1):((*ele_mhits)->at(ele_ctr)<=1);
      if(mediumselcond_) mediumselelidx.push_back(ele_ctr);
      
    } // End of loop on electrons in the event loop
    
    if(noselelidx.size()>=1) fillhistinevent("nosel_", noselelidx);
    if(looseselelidx.size()>=1) fillhistinevent("loosesel_", looseselelidx);
    if(mediumselelidx.size()>=1) fillhistinevent("mediumsel_", mediumselelidx);

    if(isMC) {
      if(noselelidx.size()>=1 && noselgenelidx.size()==2) fillgenmchhistinevent("noselgenAnosel_", noselgenelidx, noselelidx, 1);
    }
    
    // Clear all vector
    noselgenjpidx.clear();
    noselgenelidx.clear();
    noselelidx.clear();
    looseselelidx.clear();
    mediumselelidx.clear();
    
  } // End of a single event loop

  cout<<totEntries<<"\t"<<endevent<<"\t"<<beginevent<<"\t"<<endevent-beginevent<<endl;

} // End of analyzer for a single file

// Function to add a set of histograms for gen electrons
void robustanalyzer::addgenhist(TString selection) {

  if(!isMC) throw "Error!!! Trying to define gen hists for a non MC file.";

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

// Function to add a set of histograms for scouting electrons
void robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_el_mult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_el_pt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_el_eta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_el_phi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_el_trkpt","trk p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_el_trketa","trk #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_el_trkphi","trk #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_diel_M","all M(e,e)",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_diel_trkM","all trk M(e,e)",100000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sct_leadsubleaddiel_chargeprod","q_{e1}*q_{e2}",3,-1,1));
  all1dhists.push_back(new TH1F(selection+"sct_leadsubleaddiel_M","M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadbarsubleadbardiel_M","bar.e_{1} M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadecsubleadecdiel_M","ec. e_{1} M(e_{1},e_{2})",100000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sct_leadsubleaddiel_trkchargeprod","trk q_{e1}*q_{e2}",3,-1,1));
  all1dhists.push_back(new TH1F(selection+"sct_leadsubleaddiel_trkM","trk M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadbarsubleadbardiel_trkM","bar.e_{1} trk M(e_{1},e_{2})",100000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_leadecsubleadecdiel_trkM","ec. e_{1} trk M(e_{1},e_{2})",100000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sct_bar_el_pt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_eta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_phi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_trkpt","trk p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_trketa","trk #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_trkphi","trk #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_trksc_detain","#Delta#eta{trk, SC}",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_trksc_dphiin","#Delta#phi{trk, SC}",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_m","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_d0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_log10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_dz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_detain","#Delta#eta_{in}",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_dphiin","#Delta#phi_{in}",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_sigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_hoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_ooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_mhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_charge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_ecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_hcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_tkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_r9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_smin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sct_bar_el_smaj","smaj",1500,-0.1,1.4));

  all1dhists.push_back(new TH1F(selection+"sct_ec_el_pt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_eta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_phi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_trkpt","trk p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_trketa","trk #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_trkphi","trk #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_trksc_detain","#Delta#eta{trk, SC}",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_trksc_dphiin","#Delta#phi{trk, SC}",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_m","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_d0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_log10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_dz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_detain","#Delta#eta_{in}",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_dphiin","#Delta#phi_{in}",2000,-0.1,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_sigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_hoe","H/E",100000,0,1));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_ooemoop","1/E-1/p",20000,0,0.2));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_mhits","missing hits",10,0,10));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_charge","charge",5,-2,3));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_ecaliso","ecal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_hcaliso","hcal. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_tkiso","tk. iso.",10000,0,50));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_r9","r9",1500,0,15));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_smin","smin",3500,-0.1,3.4));
  all1dhists.push_back(new TH1F(selection+"sct_ec_el_smaj","smaj",1500,-0.1,1.4));
}

// Function to fill a set of histograms for scouting electrons
void robustanalyzer::fillhistinevent(TString selection, vector<int> elidx) {

  if(elidx.size()==0) return;

  TH1F* elmult = (TH1F*) outfile->Get(selection+"sct_el_mult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"sct_el_pt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"sct_el_eta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"sct_el_phi");
  TH1F* eltrkpt = (TH1F*) outfile->Get(selection+"sct_el_trkpt");
  TH1F* eltrketa = (TH1F*) outfile->Get(selection+"sct_el_trketa");
  TH1F* eltrkphi = (TH1F*) outfile->Get(selection+"sct_el_trkphi");
  TH1F* dielM = (TH1F*) outfile->Get(selection+"sct_diel_M");
  TH1F* dieltrkM = (TH1F*) outfile->Get(selection+"sct_diel_trkM");

  TH1F* leadsubleadqprod = (TH1F*) outfile->Get(selection+"sct_leadsubleaddiel_chargeprod");
  TH1F* leadsubleaddielM = (TH1F*) outfile->Get(selection+"sct_leadsubleaddiel_M");
  TH1F* leadbarsubleadbardielM = (TH1F*) outfile->Get(selection+"sct_leadbarsubleadbardiel_M");
  TH1F* leadecsubleadecdielM = (TH1F*) outfile->Get(selection+"sct_leadecsubleadecdiel_M");

  TH1F* leadsubleadtrkqprod = (TH1F*) outfile->Get(selection+"sct_leadsubleaddiel_trkchargeprod");
  TH1F* leadsubleaddieltrkM = (TH1F*) outfile->Get(selection+"sct_leadsubleaddiel_trkM");
  TH1F* leadbarsubleadbardieltrkM = (TH1F*) outfile->Get(selection+"sct_leadbarsubleadbardiel_trkM");
  TH1F* leadecsubleadecdieltrkM = (TH1F*) outfile->Get(selection+"sct_leadecsubleadecdiel_trkM");
  
  TH1F* barelpt = (TH1F*) outfile->Get(selection+"sct_bar_el_pt");
  TH1F* bareleta = (TH1F*) outfile->Get(selection+"sct_bar_el_eta");
  TH1F* barelphi = (TH1F*) outfile->Get(selection+"sct_bar_el_phi");
  TH1F* bareltrkpt = (TH1F*) outfile->Get(selection+"sct_bar_el_trkpt");
  TH1F* bareltrketa = (TH1F*) outfile->Get(selection+"sct_bar_el_trketa");
  TH1F* bareltrkphi = (TH1F*) outfile->Get(selection+"sct_bar_el_trkphi");
  TH1F* bareltrkscdeltaeta = (TH1F*) outfile->Get(selection+"sct_bar_el_trksc_detain");
  TH1F* bareltrkscdeltaphi = (TH1F*) outfile->Get(selection+"sct_bar_el_trksc_dphiin");
  TH1F* barelm = (TH1F*) outfile->Get(selection+"sct_bar_el_m");
  TH1F* bareld0 = (TH1F*) outfile->Get(selection+"sct_bar_el_d0");
  TH1F* barellog10d0 = (TH1F*) outfile->Get(selection+"sct_bar_el_log10d0");
  TH1F* bareldz = (TH1F*) outfile->Get(selection+"sct_bar_el_dz");
  TH1F* bareldetain = (TH1F*) outfile->Get(selection+"sct_bar_el_detain");
  TH1F* bareldphiin = (TH1F*) outfile->Get(selection+"sct_bar_el_dphiin");
  TH1F* barelsigmaietaieta = (TH1F*) outfile->Get(selection+"sct_bar_el_sigmaietaieta");
  TH1F* barelhoe = (TH1F*) outfile->Get(selection+"sct_bar_el_hoe");
  TH1F* barelooemoop = (TH1F*) outfile->Get(selection+"sct_bar_el_ooemoop");
  TH1F* barelmhits = (TH1F*) outfile->Get(selection+"sct_bar_el_mhits");
  TH1F* barelcharge = (TH1F*) outfile->Get(selection+"sct_bar_el_charge");
  TH1F* barelecaliso = (TH1F*) outfile->Get(selection+"sct_bar_el_ecaliso");
  TH1F* barelhcaliso = (TH1F*) outfile->Get(selection+"sct_bar_el_hcaliso");
  TH1F* bareltkiso = (TH1F*) outfile->Get(selection+"sct_bar_el_tkiso");
  TH1F* barelr9 = (TH1F*) outfile->Get(selection+"sct_bar_el_r9");
  TH1F* barelsmin = (TH1F*) outfile->Get(selection+"sct_bar_el_smin");
  TH1F* barelsmaj = (TH1F*) outfile->Get(selection+"sct_bar_el_smaj");
  
  TH1F* ecelpt = (TH1F*) outfile->Get(selection+"sct_ec_el_pt");
  TH1F* eceleta = (TH1F*) outfile->Get(selection+"sct_ec_el_eta");
  TH1F* ecelphi = (TH1F*) outfile->Get(selection+"sct_ec_el_phi");
  TH1F* eceltrkpt = (TH1F*) outfile->Get(selection+"sct_ec_el_trkpt");
  TH1F* eceltrketa = (TH1F*) outfile->Get(selection+"sct_ec_el_trketa");
  TH1F* eceltrkphi = (TH1F*) outfile->Get(selection+"sct_ec_el_trkphi");
  TH1F* eceltrkscdeltaeta = (TH1F*) outfile->Get(selection+"sct_ec_el_trksc_detain");
  TH1F* eceltrkscdeltaphi = (TH1F*) outfile->Get(selection+"sct_ec_el_trksc_dphiin");
  TH1F* ecelm = (TH1F*) outfile->Get(selection+"sct_ec_el_m");
  TH1F* eceld0 = (TH1F*) outfile->Get(selection+"sct_ec_el_d0");
  TH1F* ecellog10d0 = (TH1F*) outfile->Get(selection+"sct_ec_el_log10d0");
  TH1F* eceldz = (TH1F*) outfile->Get(selection+"sct_ec_el_dz");
  TH1F* eceldetain = (TH1F*) outfile->Get(selection+"sct_ec_el_detain");
  TH1F* eceldphiin = (TH1F*) outfile->Get(selection+"sct_ec_el_dphiin");
  TH1F* ecelsigmaietaieta = (TH1F*) outfile->Get(selection+"sct_ec_el_sigmaietaieta");
  TH1F* ecelhoe = (TH1F*) outfile->Get(selection+"sct_ec_el_hoe");
  TH1F* ecelooemoop = (TH1F*) outfile->Get(selection+"sct_ec_el_ooemoop");
  TH1F* ecelmhits = (TH1F*) outfile->Get(selection+"sct_ec_el_mhits");
  TH1F* ecelcharge = (TH1F*) outfile->Get(selection+"sct_ec_el_charge");
  TH1F* ecelecaliso = (TH1F*) outfile->Get(selection+"sct_ec_el_ecaliso");
  TH1F* ecelhcaliso = (TH1F*) outfile->Get(selection+"sct_ec_el_hcaliso");
  TH1F* eceltkiso = (TH1F*) outfile->Get(selection+"sct_ec_el_tkiso");
  TH1F* ecelr9 = (TH1F*) outfile->Get(selection+"sct_ec_el_r9");
  TH1F* ecelsmin = (TH1F*) outfile->Get(selection+"sct_ec_el_smin");
  TH1F* ecelsmaj = (TH1F*) outfile->Get(selection+"sct_ec_el_smaj");

  elmult->Fill(elidx.size());
  
  int leadelpos=-1.0, subleadelpos=-1.0, leadeltrkpos=-1.0,subleadeltrkpos=-1.0;
  for(unsigned int ctr=0; ctr<elidx.size(); ctr++) {
    elpt->Fill((*ele_pt)->at(elidx[ctr]));
    eleta->Fill((*ele_eta)->at(elidx[ctr]));
    elphi->Fill((*ele_phi)->at(elidx[ctr]));
    //eltrkpt->Fill((*ele_trkpt)->at(elidx[ctr]));
    //eltrketa->Fill((*ele_trketa)->at(elidx[ctr]));
    //eltrkphi->Fill((*ele_trkphi)->at(elidx[ctr]));
    for(unsigned int ctr2=ctr+1; ctr2<elidx.size(); ctr2++) {
      TLorentzVector el, el2;
      el.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr]),(*ele_eta)->at(elidx[ctr]),(*ele_phi)->at(elidx[ctr]),0.0005);
      el2.SetPtEtaPhiM((*ele_pt)->at(elidx[ctr2]),(*ele_eta)->at(elidx[ctr2]),(*ele_phi)->at(elidx[ctr2]),0.0005);
      dielM->Fill((el+el2).M());
      //TLorentzVector trkel, trkel2;
      //trkel.SetPtEtaPhiM((*ele_trkpt)->at(elidx[ctr]),(*ele_trketa)->at(elidx[ctr]),(*ele_trkphi)->at(elidx[ctr]),0.0005);
      //trkel2.SetPtEtaPhiM((*ele_trkpt)->at(elidx[ctr2]),(*ele_trketa)->at(elidx[ctr2]),(*ele_trkphi)->at(elidx[ctr2]),0.0005);
      //dieltrkM->Fill((trkel+trkel2).M());
    }

    // Segregation in barrel and endcap is based on the calo angular parameters
    if(TMath::Abs((*ele_eta)->at(elidx[ctr]))<1.479) {
      barelpt->Fill((*ele_pt)->at(elidx[ctr]));
      bareleta->Fill((*ele_eta)->at(elidx[ctr]));
      barelphi->Fill((*ele_phi)->at(elidx[ctr]));
      //bareltrkpt->Fill((*ele_trkpt)->at(elidx[ctr]));
      //bareltrketa->Fill((*ele_trketa)->at(elidx[ctr]));
      //bareltrkphi->Fill((*ele_trkphi)->at(elidx[ctr]));
      //bareltrkscdeltaeta->Fill( (*ele_trketa)->at(elidx[ctr])-(*ele_eta)->at(elidx[ctr]) );
      //bareltrkscdeltaphi->Fill( (*ele_trkphi)->at(elidx[ctr])-(*ele_phi)->at(elidx[ctr]) );
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
      //eceltrkpt->Fill((*ele_trkpt)->at(elidx[ctr]));
      //eceltrketa->Fill((*ele_trketa)->at(elidx[ctr]));
      //eceltrkphi->Fill((*ele_trkphi)->at(elidx[ctr]));
      //eceltrkscdeltaeta->Fill( (*ele_trketa)->at(elidx[ctr])-(*ele_eta)->at(elidx[ctr]) );
      //eceltrkscdeltaphi->Fill( (*ele_trkphi)->at(elidx[ctr])-(*ele_phi)->at(elidx[ctr]) );
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
    
    // Find the lead and the sublead electron in the event
    if(leadelpos == -1) {
      leadelpos = elidx[ctr];
    }
    else {
      int charge1 = (*ele_charge)->at(elidx[ctr]);
      int charge2 = (*ele_charge)->at(leadelpos);
      if( (charge1*charge2)==-1 ) {
	if( ((*ele_pt)->at(elidx[ctr]))>((*ele_pt)->at(leadelpos)) ) {
	  subleadelpos = leadelpos;
	  leadelpos = elidx[ctr];
	}
	else {
	  if(subleadelpos==-1) {
	    subleadelpos = elidx[ctr];	    
	  }
	  else {
	    if( ((*ele_pt)->at(elidx[ctr]))>((*ele_pt)->at(subleadelpos)) ) {
	      subleadelpos = elidx[ctr];	    
	    }
	  }
	}
      }
    }

    if(leadeltrkpos == -1) {
      leadeltrkpos = elidx[ctr];
    }
    else {
      int charge1 = (*ele_charge)->at(elidx[ctr]);
      int charge2 = (*ele_charge)->at(leadeltrkpos);
      if( (charge1*charge2)==-1 ) {
	if( ((*ele_pt)->at(elidx[ctr]))>((*ele_pt)->at(leadeltrkpos)) ) {
	  subleadeltrkpos = leadeltrkpos;
	  leadeltrkpos = elidx[ctr];
	}
	else {
	  if(subleadeltrkpos==-1) {
	    subleadeltrkpos = elidx[ctr];	    
	  }
	  else {
	    if( ((*ele_pt)->at(elidx[ctr]))>((*ele_pt)->at(subleadeltrkpos)) ) {
	      subleadeltrkpos = elidx[ctr];	    
	    }
	  }
	}
      }
    }
    
  } // End of main electron for loop
  
  if(leadelpos!=-1 && subleadelpos!=-1) {
    int charge1 = (*ele_charge)->at(leadelpos);
    int charge2 = (*ele_charge)->at(subleadelpos);
    leadsubleadqprod->Fill(charge1*charge2);
    TLorentzVector leadel, subleadel;
    leadel.SetPtEtaPhiM((*ele_pt)->at(leadelpos),(*ele_eta)->at(leadelpos),(*ele_phi)->at(leadelpos),0.0005);
    subleadel.SetPtEtaPhiM((*ele_pt)->at(subleadelpos),(*ele_eta)->at(subleadelpos),(*ele_phi)->at(subleadelpos),0.0005);
    leadsubleaddielM->Fill((leadel+subleadel).M());
    if(abs((*ele_eta)->at(leadelpos))<1.479 && abs((*ele_eta)->at(subleadelpos))<1.479)  {
    leadbarsubleadbardielM->Fill((leadel+subleadel).M());
    }
    if(abs((*ele_eta)->at(leadelpos))>1.479 && abs((*ele_eta)->at(subleadelpos))>1.479){
    leadecsubleadecdielM->Fill((leadel+subleadel).M());
    }
  }
  if(leadeltrkpos!=-1 && subleadeltrkpos!=-1) {
    leadsubleadtrkqprod->Fill(((*ele_charge)->at(leadeltrkpos))*((*ele_charge)->at(subleadeltrkpos)));
    //TLorentzVector leadel, subleadel;
    //leadel.SetPtEtaPhiM((*ele_trkpt)->at(leadeltrkpos),(*ele_trketa)->at(leadeltrkpos),(*ele_trkphi)->at(leadeltrkpos),0.0005);
    //subleadel.SetPtEtaPhiM((*ele_trkpt)->at(subleadeltrkpos),(*ele_trketa)->at(subleadeltrkpos),(*ele_trkphi)->at(subleadeltrkpos),0.0005);
    //leadsubleaddieltrkM->Fill((leadel+subleadel).M());
    //if(abs((*ele_eta)->at(leadeltrkpos))<1.479 && abs((*ele_eta)->at(subleadeltrkpos))<1.479)  {
    //leadbarsubleadbardieltrkM->Fill((leadel+subleadel).M());
    //}
    //if(abs((*ele_eta)->at(leadeltrkpos))>1.479 && abs((*ele_eta)->at(subleadeltrkpos))>1.479){
    //leadecsubleadecdieltrkM->Fill((leadel+subleadel).M());
    //}
  }
  
}  

// Add histograms for gen matching
void robustanalyzer::addgenmchhist(TString selection) {
  
  if(!isMC) throw "Error!!! Trying to define gen hists for a non MC file.";
    
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

  all1dhists.push_back(new TH1F(selection+"sctmchgenel_el_mult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_el_pt","gen e p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_el_eta","gen e #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_el_phi","gen e #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_diel_M","M(e,e)",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_lead_el_pt","gen e p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_sublead_el_pt","gen e p_{T} / GeV",1000,-10,990));

  all1dhists.push_back(new TH1F(selection+"genmchsct_el_mult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"genmchsct_el_pt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"genmchsct_el_eta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genmchsct_el_phi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"genmchsct_diel_M","all M(e,e)",1000,-10,990));
}

// Fill histograms for gen matching
void robustanalyzer::fillgenmchhistinevent(TString selection, vector<int> genidx, vector<int> elidx, double w) {

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

  TH1F* genelmult = (TH1F*) outfile->Get(selection+"sctmchgenel_el_mult");
  TH1F* genelpt = (TH1F*) outfile->Get(selection+"sctmchgenel_el_pt");
  TH1F* geneleta = (TH1F*) outfile->Get(selection+"sctmchgenel_el_eta");
  TH1F* genelphi = (TH1F*) outfile->Get(selection+"sctmchgenel_el_phi");
  TH1F* gendielM = (TH1F*) outfile->Get(selection+"sctmchgenel_diel_M");
  TH1F* genleadelpt = (TH1F*) outfile->Get(selection+"sctmchgenel_lead_el_pt");
  TH1F* gensubleadelpt = (TH1F*) outfile->Get(selection+"sctmchgenel_sublead_el_pt");

  TH1F* genmchsctelmult = (TH1F*) outfile->Get(selection+"genmchsct_el_mult");
  TH1F* genmchsctelpt = (TH1F*) outfile->Get(selection+"genmchsct_el_pt");
  TH1F* genmchscteleta = (TH1F*) outfile->Get(selection+"genmchsct_el_eta");
  TH1F* genmchsctelphi = (TH1F*) outfile->Get(selection+"genmchsct_el_phi");
  TH1F* genmchsctdielM = (TH1F*) outfile->Get(selection+"genmchsct_diel_M");

  // Fill variables before gen match
  for(unsigned int sct=0; sct<elidx.size(); sct++) {
    for(unsigned int gen=0; gen<genidx.size(); gen++) {
      double charge = ((*genlep_pdg)->at(genidx[gen]))/TMath::Abs((*genlep_pdg)->at(genidx[gen]));
      if(TMath::Abs((*ele_eta)->at(elidx[sct]))<1.479) {
	bardEta->Fill((*ele_eta)->at(elidx[sct])-(*genlep_eta)->at(genidx[gen]));
	barqdPhi->Fill(charge*((*ele_phi)->at(elidx[sct])-(*genlep_phi)->at(genidx[gen])));
      }
      else {
	eedEta->Fill((*ele_eta)->at(elidx[sct])-(*genlep_eta)->at(genidx[gen]));
	eeqdPhi->Fill(charge*((*ele_phi)->at(elidx[sct])-(*genlep_phi)->at(genidx[gen])));
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
    if((*genlep_pt)->at(gen1) >= (*genlep_pt)->at(gen2)) {
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
    genelpt->Fill((*genlep_pt)->at(gen1), w);
    geneleta->Fill((*genlep_eta)->at(gen1), w);
    genelphi->Fill((*genlep_phi)->at(gen1), w);
  }
  if(el2!=-1) {
    genmchfound++;
    genelpt->Fill((*genlep_pt)->at(gen2), w);
    geneleta->Fill((*genlep_eta)->at(gen2), w);
    genelphi->Fill((*genlep_phi)->at(gen2), w);
  } 
  genelmult->Fill(genmchfound, w);
  
  if(leadptpos!=-1) {
    genleadelpt->Fill((*genlep_pt)->at(leadptpos), w);
  }
  if(subleadptpos!=-1) {
    gensubleadelpt->Fill((*genlep_pt)->at(subleadptpos), w);
  }
  
  if(genmchfound==2) {
    TLorentzVector elec1, elec2;
    elec1.SetPtEtaPhiM((*genlep_pt)->at(gen1),(*genlep_eta)->at(gen1),(*genlep_phi)->at(gen1),0.0005);
    elec2.SetPtEtaPhiM((*genlep_pt)->at(gen2),(*genlep_eta)->at(gen2),(*genlep_phi)->at(gen2),0.0005);
    gendielM->Fill((elec1+elec2).M(), w);
  }
  
  if(el1!=-1) {
    double charge = ((*genlep_pdg)->at(gen1))/TMath::Abs((*genlep_pdg)->at(gen1));
    if(TMath::Abs((*ele_eta)->at(el1))<1.479) {
      mchbardEta->Fill((*ele_eta)->at(el1)-(*genlep_eta)->at(gen1));
      mchbarqdPhi->Fill(charge*((*ele_phi)->at(el1)-(*genlep_phi)->at(gen1)));
    }
    else {
      mcheedEta->Fill((*ele_eta)->at(el1)-(*genlep_eta)->at(gen1));
      mcheeqdPhi->Fill(charge*((*ele_phi)->at(el1)-(*genlep_phi)->at(gen1)));
    }
    genmchsctelpt->Fill((*ele_pt)->at(el1), w);
    genmchscteleta->Fill((*ele_eta)->at(el1), w);
    genmchsctelphi->Fill((*ele_phi)->at(el1), w);
  }
  
  if(el2!=-1) {
    double charge = ((*genlep_pdg)->at(gen2))/TMath::Abs((*genlep_pdg)->at(gen2));
    if(TMath::Abs((*ele_eta)->at(el2))<1.479) {
      mchbardEta->Fill((*ele_eta)->at(el2)-(*genlep_eta)->at(gen2));
      mchbarqdPhi->Fill(charge*((*ele_phi)->at(el2)-(*genlep_phi)->at(gen2)));
    }
    else {
      mcheedEta->Fill((*ele_eta)->at(el2)-(*genlep_eta)->at(gen2));
      mcheeqdPhi->Fill(charge*((*ele_phi)->at(el2)-(*genlep_phi)->at(gen2)));
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

vector< pair<int,int> > robustanalyzer::diElecGenMatching(vector<int> genidx, vector<int> sctelidx) {

  if(genidx.size()==0) {
    return { make_pair(-1,-1), make_pair(-1,-1) };
  }
  
  if(genidx.size()>2) {
    cout<<"Found "<<genidx.size()<<" gen electrons."<<endl;
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

      double diffeta = abs((*ele_eta)->at(elidx)-(*genlep_eta)->at(geneidx)); 
      TLorentzVector vecsctel, vecgen;
      vecgen.SetPtEtaPhiM((*genlep_pt)->at(geneidx),(*genlep_eta)->at(geneidx),(*genlep_phi)->at(geneidx),0.0005);
      vecsctel.SetPtEtaPhiM((*ele_pt)->at(elidx),(*ele_eta)->at(elidx),(*ele_phi)->at(elidx),0.0005);
      double qdiffphi = ( (*genlep_pdg)->at(geneidx)/abs((*genlep_pdg)->at(geneidx)) )*( vecgen.DeltaPhi(vecsctel) );
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
/*
  vector<int> noselelidx;
  vector<int> noselZwindelidx;
  vector<int> leadepemselelidx;
  vector<int> leadepemonlyselelidx;
  vector<int> leadepemptgt20selelidx;
  vector<int> leadepemselZwindelidx;
  vector<int> leadepemZwindptgt20selelidx;
  vector<int> vetoselelidx;
  vector<int> looseselelidx;
  vector<int> looseselZwindelidx;
  vector<int> mediumselelidx;
  vector<int> mediumselZwindelidx;
  vector<int> tightselelidx;
  vector<int> tightselZwindelidx;
  vector<int> tightselSideBand1elidx;
  vector<int> tightselSideBand2elidx;
  vector<int> tightselSideBand3elidx;
  vector<int> tightnoneselelidx;
  vector<int> tightnoneselZwindelidx;
  vector<int> tightnoneselSideBand1elidx;
  vector<int> tightnoneselSideBand2elidx;
  vector<int> tightnoneselSideBand3elidx;
  vector<int> nonetightselelidx;
  vector<int> nonetightselZwindelidx;
  vector<int> nonetightselSideBand1elidx;
  vector<int> nonetightselSideBand2elidx;
  vector<int> nonetightselSideBand3elidx;
  vector<int> tightlooseselelidx;
  vector<int> tightlooseselZwindelidx;
  vector<int> tightlooseselSideBand1elidx;
  vector<int> tightlooseselSideBand2elidx;
  vector<int> tightlooseselSideBand3elidx;
  vector<int> loosetightselelidx;
  vector<int> loosetightselZwindelidx;
  vector<int> loosetightselSideBand1elidx;
  vector<int> loosetightselSideBand2elidx;
  vector<int> loosetightselSideBand3elidx;
  vector<int> tightmediumselelidx;
  vector<int> tightmediumselZwindelidx;
  vector<int> tightmediumselSideBand1elidx;
  vector<int> tightmediumselSideBand2elidx;
  vector<int> tightmediumselSideBand3elidx;
  vector<int> mediumtightselelidx;
  vector<int> mediumtightselZwindelidx;
  vector<int> mediumtightselSideBand1elidx;
  vector<int> mediumtightselSideBand2elidx;
  vector<int> mediumtightselSideBand3elidx;
    
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%100000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    
    if(isMC && isDY) {
      if(noselgenidx.size()<2) continue;
    }
    
    // Sort the electrons based on their pT
    vector<int> sortedelidx((*(*n_ele)));
    iota(begin(sortedelidx), end(sortedelidx), 0);
    sort(&sortedelidx[0], ele_pt, (*(*n_ele))); // Verified that the algorithm works fine

    if((*(*n_ele))<0) throw "Error!!! Wrong technical event processing. Negative number of electrons in event.";

    bool vetoselcond = false;
    bool looseselcond = false;
    bool mediumselcond = false;
    bool tightselcond = false;
    
    // Loop on electrons in the event loop
    for(unsigned int ele_ctr=0; ele_ctr<(*(*n_ele)); ele_ctr++) {

      // Take the sorted index only
      unsigned int elidx = sortedelidx[ele_ctr];

      noselelidx.push_back(elidx);

      // Get the energy of this electron
      double ele_energy=0;
      for(unsigned int ctr=0; ctr<((*ele_enemat)->at(elidx)).size(); ctr++) {
	//if(((*ele_enemat)->at(elidx))[ctr]<0) throw "Error!!! Negative energy value in the electron energy matrix." ;
	if(((*ele_enemat)->at(elidx))[ctr]>0) {
	  ele_energy += ((*ele_enemat)->at(elidx))[ctr];
	}
      }
      
   }// End of loop on electrons in the event loop

    // At least 2 opp. charged electrons in event
    if(noselelidx.size()>=2) {
      if((*ele_charge)->at(noselelidx[0])*(*ele_charge)->at(noselelidx[1])<0) {

        leadepemselelidx = noselelidx;
	leadepemonlyselelidx.push_back(leadepemselelidx[0]);
	leadepemonlyselelidx.push_back(leadepemselelidx[1]);
	if((*ele_pt)->at(leadepemselelidx[0])>20 && (*ele_pt)->at(leadepemselelidx[1])>20) {
	  leadepemptgt20selelidx.push_back(leadepemselelidx[0]);
	  leadepemptgt20selelidx.push_back(leadepemselelidx[1]);
	}
	
	bool leadloosecond=false, leadmediumcond=false, leadtightcond=false;
	bool subleadloosecond=false, subleadmediumcond=false, subleadtightcond=false;
	
	double lead_ele_energy=0;
	for(unsigned int ctr=0; ctr<((*ele_enemat)->at(noselelidx[0])).size(); ctr++) {
	  //if(((*ele_enemat)->at(noselelidx[0]))[ctr]<0) throw "Error!!! Negative energy value in the electron energy matrix." ;
	  lead_ele_energy += ((*ele_enemat)->at(noselelidx[0]))[ctr];
	}
	
	double sublead_ele_energy=0;
	for(unsigned int ctr=0; ctr<((*ele_enemat)->at(noselelidx[1])).size(); ctr++) {
	  //if(((*ele_enemat)->at(noselelidx[1]))[ctr]<0) throw "Error!!! Negative energy value in the electron energy matrix." ;
	  sublead_ele_energy += ((*ele_enemat)->at(noselelidx[1]))[ctr];
	}
	

	if(leadtightcond) {
	  tightnoneselelidx.push_back(noselelidx[0]);
	  tightnoneselelidx.push_back(noselelidx[1]);
	  if(subleadloosecond) {
	    tightlooseselelidx.push_back(noselelidx[0]);
	    tightlooseselelidx.push_back(noselelidx[1]);
	  }
	  if(subleadmediumcond) {
	    tightmediumselelidx.push_back(noselelidx[0]);
	    tightmediumselelidx.push_back(noselelidx[1]);
	  }
	}

	if(subleadtightcond) {
	  nonetightselelidx.push_back(noselelidx[0]);
	  nonetightselelidx.push_back(noselelidx[1]);
	  if(leadloosecond) {
	    loosetightselelidx.push_back(noselelidx[0]);
	    loosetightselelidx.push_back(noselelidx[1]);
	  }
	  if(leadmediumcond) {
	    mediumtightselelidx.push_back(noselelidx[0]);
	    mediumtightselelidx.push_back(noselelidx[1]);
	  }
	}
      }	
    } // Atleast 2 opp. charged e in event

    // At least 2 opp. charged electrons with no sel in event
    if(noselelidx.size()>=2) {
      if((*ele_charge)->at(noselelidx[0])*(*ele_charge)->at(noselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(noselelidx[0]),(*ele_eta)->at(noselelidx[0]),(*ele_phi)->at(noselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(noselelidx[1]),(*ele_eta)->at(noselelidx[1]),(*ele_phi)->at(noselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  noselZwindelidx.push_back(noselelidx[0]);
	  noselZwindelidx.push_back(noselelidx[1]);
	}
      }
    } // Atleast 2 opp. charged e with nosel in event

    if(leadepemonlyselelidx.size()>=2) {
      TLorentzVector leadel, subleadel;
      leadel.SetPtEtaPhiM((*ele_pt)->at(leadepemonlyselelidx[0]),(*ele_eta)->at(leadepemonlyselelidx[0]),(*ele_phi)->at(leadepemonlyselelidx[0]),0.0005);
      subleadel.SetPtEtaPhiM((*ele_pt)->at(leadepemonlyselelidx[1]),(*ele_eta)->at(leadepemonlyselelidx[1]),(*ele_phi)->at(leadepemonlyselelidx[1]),0.0005);
      if(inZwind(leadel, subleadel)) {
	leadepemselZwindelidx.push_back(noselelidx[0]);
	leadepemselZwindelidx.push_back(noselelidx[1]);
	if((*ele_pt)->at(leadepemonlyselelidx[0])>20 && (*ele_pt)->at(leadepemonlyselelidx[1])>20) {
	  leadepemZwindptgt20selelidx.push_back(noselelidx[0]);
	  leadepemZwindptgt20selelidx.push_back(noselelidx[1]);
	}
      }
    }
    
    // At least 2 opp. charged electrons with loose ID in event
    if(looseselelidx.size()>=2) {
      if((*ele_charge)->at(looseselelidx[0])*(*ele_charge)->at(looseselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(looseselelidx[0]),(*ele_eta)->at(looseselelidx[0]),(*ele_phi)->at(looseselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(looseselelidx[1]),(*ele_eta)->at(looseselelidx[1]),(*ele_phi)->at(looseselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  looseselZwindelidx.push_back(looseselelidx[0]);
	  looseselZwindelidx.push_back(looseselelidx[1]);
	}
      }
    } // Atleast 2 opp. charged e with loose ID in event
    
    // At least 2 opp. charged electrons with medium ID in event
    if(mediumselelidx.size()>=2) {
      if((*ele_charge)->at(mediumselelidx[0])*(*ele_charge)->at(mediumselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(mediumselelidx[0]),(*ele_eta)->at(mediumselelidx[0]),(*ele_phi)->at(mediumselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(mediumselelidx[1]),(*ele_eta)->at(mediumselelidx[1]),(*ele_phi)->at(mediumselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  mediumselZwindelidx.push_back(mediumselelidx[0]);
	  mediumselZwindelidx.push_back(mediumselelidx[1]);
	}
      }
    } // Atleast 2 opp. charged e with medium ID in event
    
    // At least 2 opp. charged electrons with tight ID in event
    if(tightselelidx.size()>=2) {
      if((*ele_charge)->at(tightselelidx[0])*(*ele_charge)->at(tightselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(tightselelidx[0]),(*ele_eta)->at(tightselelidx[0]),(*ele_phi)->at(tightselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(tightselelidx[1]),(*ele_eta)->at(tightselelidx[1]),(*ele_phi)->at(tightselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  tightselZwindelidx.push_back(tightselelidx[0]);
	  tightselZwindelidx.push_back(tightselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  tightselSideBand1elidx.push_back(tightselelidx[0]);
	  tightselSideBand1elidx.push_back(tightselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  tightselSideBand2elidx.push_back(tightselelidx[0]);
	  tightselSideBand2elidx.push_back(tightselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  tightselSideBand3elidx.push_back(tightselelidx[0]);
	  tightselSideBand3elidx.push_back(tightselelidx[1]);
	}
      }
    } // Atleast 2 opp. charged e with tight ID in event
    
    // Lead tight Sub-lead none in event
    if(tightnoneselelidx.size()>=2) {
      if((*ele_charge)->at(tightnoneselelidx[0])*(*ele_charge)->at(tightnoneselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(tightnoneselelidx[0]),(*ele_eta)->at(tightnoneselelidx[0]),(*ele_phi)->at(tightnoneselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(tightnoneselelidx[1]),(*ele_eta)->at(tightnoneselelidx[1]),(*ele_phi)->at(tightnoneselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  tightnoneselZwindelidx.push_back(tightnoneselelidx[0]);
	  tightnoneselZwindelidx.push_back(tightnoneselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  tightnoneselSideBand1elidx.push_back(tightnoneselelidx[0]);
	  tightnoneselSideBand1elidx.push_back(tightnoneselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  tightnoneselSideBand2elidx.push_back(tightnoneselelidx[0]);
	  tightnoneselSideBand2elidx.push_back(tightnoneselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  tightnoneselSideBand3elidx.push_back(tightnoneselelidx[0]);
	  tightnoneselSideBand3elidx.push_back(tightnoneselelidx[1]);
	}
      }
    } // Lead tight Sub-lead none in event

    // Lead none Sub-lead tight in event
    if(nonetightselelidx.size()>=2) {
      if((*ele_charge)->at(nonetightselelidx[0])*(*ele_charge)->at(nonetightselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(nonetightselelidx[0]),(*ele_eta)->at(nonetightselelidx[0]),(*ele_phi)->at(nonetightselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(nonetightselelidx[1]),(*ele_eta)->at(nonetightselelidx[1]),(*ele_phi)->at(nonetightselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  nonetightselZwindelidx.push_back(nonetightselelidx[0]);
	  nonetightselZwindelidx.push_back(nonetightselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  nonetightselSideBand1elidx.push_back(nonetightselelidx[0]);
	  nonetightselSideBand1elidx.push_back(nonetightselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  nonetightselSideBand2elidx.push_back(nonetightselelidx[0]);
	  nonetightselSideBand2elidx.push_back(nonetightselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  nonetightselSideBand3elidx.push_back(nonetightselelidx[0]);
	  nonetightselSideBand3elidx.push_back(nonetightselelidx[1]);
	}
      }
    } // Lead none Sub-lead tight in event

    // Lead tight Sub-lead loose in event
    if(tightlooseselelidx.size()>=2) {
      if((*ele_charge)->at(tightlooseselelidx[0])*(*ele_charge)->at(tightlooseselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(tightlooseselelidx[0]),(*ele_eta)->at(tightlooseselelidx[0]),(*ele_phi)->at(tightlooseselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(tightlooseselelidx[1]),(*ele_eta)->at(tightlooseselelidx[1]),(*ele_phi)->at(tightlooseselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  tightlooseselZwindelidx.push_back(tightlooseselelidx[0]);
	  tightlooseselZwindelidx.push_back(tightlooseselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  tightlooseselSideBand1elidx.push_back(tightlooseselelidx[0]);
	  tightlooseselSideBand1elidx.push_back(tightlooseselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  tightlooseselSideBand2elidx.push_back(tightlooseselelidx[0]);
	  tightlooseselSideBand2elidx.push_back(tightlooseselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  tightlooseselSideBand3elidx.push_back(tightlooseselelidx[0]);
	  tightlooseselSideBand3elidx.push_back(tightlooseselelidx[1]);
	}
      }
    } // Lead tight Sub-lead loose in event

    // Lead loose Sub-lead tight in event
    if(loosetightselelidx.size()>=2) {
      if((*ele_charge)->at(loosetightselelidx[0])*(*ele_charge)->at(loosetightselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(loosetightselelidx[0]),(*ele_eta)->at(loosetightselelidx[0]),(*ele_phi)->at(loosetightselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(loosetightselelidx[1]),(*ele_eta)->at(loosetightselelidx[1]),(*ele_phi)->at(loosetightselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  loosetightselZwindelidx.push_back(loosetightselelidx[0]);
	  loosetightselZwindelidx.push_back(loosetightselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  loosetightselSideBand1elidx.push_back(loosetightselelidx[0]);
	  loosetightselSideBand1elidx.push_back(loosetightselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  loosetightselSideBand2elidx.push_back(loosetightselelidx[0]);
	  loosetightselSideBand2elidx.push_back(loosetightselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  loosetightselSideBand3elidx.push_back(loosetightselelidx[0]);
	  loosetightselSideBand3elidx.push_back(loosetightselelidx[1]);
	}
      }
    } // Lead loose Sub-lead tight in event

    // Lead tight Sub-lead medium in event
    if(tightmediumselelidx.size()>=2) {
      if((*ele_charge)->at(tightmediumselelidx[0])*(*ele_charge)->at(tightmediumselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(tightmediumselelidx[0]),(*ele_eta)->at(tightmediumselelidx[0]),(*ele_phi)->at(tightmediumselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(tightmediumselelidx[1]),(*ele_eta)->at(tightmediumselelidx[1]),(*ele_phi)->at(tightmediumselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  tightmediumselZwindelidx.push_back(tightmediumselelidx[0]);
	  tightmediumselZwindelidx.push_back(tightmediumselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  tightmediumselSideBand1elidx.push_back(tightmediumselelidx[0]);
	  tightmediumselSideBand1elidx.push_back(tightmediumselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  tightmediumselSideBand2elidx.push_back(tightmediumselelidx[0]);
	  tightmediumselSideBand2elidx.push_back(tightmediumselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  tightmediumselSideBand3elidx.push_back(tightmediumselelidx[0]);
	  tightmediumselSideBand3elidx.push_back(tightmediumselelidx[1]);
	}
      }
    } // Lead tight Sub-lead medium in event

    // Lead medium Sub-lead tight in event
    if(mediumtightselelidx.size()>=2) {
      if((*ele_charge)->at(mediumtightselelidx[0])*(*ele_charge)->at(mediumtightselelidx[1])<0) {
	TLorentzVector leadel, subleadel;
	leadel.SetPtEtaPhiM((*ele_pt)->at(mediumtightselelidx[0]),(*ele_eta)->at(mediumtightselelidx[0]),(*ele_phi)->at(mediumtightselelidx[0]),0.0005);
	subleadel.SetPtEtaPhiM((*ele_pt)->at(mediumtightselelidx[1]),(*ele_eta)->at(mediumtightselelidx[1]),(*ele_phi)->at(mediumtightselelidx[1]),0.0005);
	if(inZwind(leadel, subleadel)) {
	  mediumtightselZwindelidx.push_back(mediumtightselelidx[0]);
	  mediumtightselZwindelidx.push_back(mediumtightselelidx[1]);
	}
	if(inSideBand1(leadel, subleadel)) {
	  mediumtightselSideBand1elidx.push_back(mediumtightselelidx[0]);
	  mediumtightselSideBand1elidx.push_back(mediumtightselelidx[1]);
	}
	if(inSideBand2(leadel, subleadel)) {
	  mediumtightselSideBand2elidx.push_back(mediumtightselelidx[0]);
	  mediumtightselSideBand2elidx.push_back(mediumtightselelidx[1]);
	}
	if(inSideBand3(leadel, subleadel)) {
	  mediumtightselSideBand3elidx.push_back(mediumtightselelidx[0]);
	  mediumtightselSideBand3elidx.push_back(mediumtightselelidx[1]);
	}
      }
    } // Lead medium Sub-lead tight in event

    if(noselelidx.size()>0) nosel++;
    fillhistinevent("nosel", noselelidx);
    fillhistinevent("nosel_Zwind_", noselZwindelidx);
    fillhistinevent("leadepemsel", leadepemselelidx);
    fillhistinevent("leadepem_only_sel", leadepemonlyselelidx);
    fillhistinevent("leadepem_ptgt20_sel", leadepemptgt20selelidx);
    fillhistinevent("leadepem_Zwind_sel", leadepemselZwindelidx);
    fillhistinevent("leadepem_Zwindptgt20_sel", leadepemZwindptgt20selelidx);
    fillhistinevent("vetosel", vetoselelidx);
    fillhistinevent("loosesel", looseselelidx);
    fillhistinevent("loosesel_Zwind_", looseselZwindelidx);
    fillhistinevent("mediumsel", mediumselelidx);
    fillhistinevent("mediumsel_Zwind_", mediumselZwindelidx);
    fillhistinevent("tightsel", tightselelidx);
    fillhistinevent("tightsel_Zwind_", tightselZwindelidx);
    fillhistinevent("tightsel_SideBand1_", tightselSideBand1elidx);
    fillhistinevent("tightsel_SideBand2_", tightselSideBand2elidx);
    fillhistinevent("tightsel_SideBand3_", tightselSideBand3elidx);
    fillhistinevent("tightnonesel", tightnoneselelidx);
    fillhistinevent("tightnonesel_Zwind_", tightnoneselZwindelidx);
    fillhistinevent("tightnonesel_SideBand1_", tightnoneselSideBand1elidx);
    fillhistinevent("tightnonesel_SideBand2_", tightnoneselSideBand2elidx);
    fillhistinevent("tightnonesel_SideBand3_", tightnoneselSideBand3elidx);
    fillhistinevent("nonetightsel", nonetightselelidx);
    fillhistinevent("nonetightsel_Zwind_", nonetightselZwindelidx);
    fillhistinevent("nonetightsel_SideBand1_", nonetightselSideBand1elidx);
    fillhistinevent("nonetightsel_SideBand2_", nonetightselSideBand2elidx);
    fillhistinevent("nonetightsel_SideBand3_", nonetightselSideBand3elidx);
    fillhistinevent("tightloosesel", tightlooseselelidx);
    fillhistinevent("tightloosesel_Zwind_", tightlooseselZwindelidx);
    fillhistinevent("tightloosesel_SideBand1_", tightlooseselSideBand1elidx);
    fillhistinevent("tightloosesel_SideBand2_", tightlooseselSideBand2elidx);
    fillhistinevent("tightloosesel_SideBand3_", tightlooseselSideBand3elidx);
    fillhistinevent("loosetightsel", loosetightselelidx);
    fillhistinevent("loosetightsel_Zwind_", loosetightselZwindelidx);
    fillhistinevent("loosetightsel_SideBand1_", loosetightselSideBand1elidx);
    fillhistinevent("loosetightsel_SideBand2_", loosetightselSideBand2elidx);
    fillhistinevent("loosetightsel_SideBand3_", loosetightselSideBand3elidx);
    fillhistinevent("tightmediumsel", tightmediumselelidx);
    fillhistinevent("tightmediumsel_Zwind_", tightmediumselZwindelidx);
    fillhistinevent("tightmediumsel_SideBand1_", tightmediumselSideBand1elidx);
    fillhistinevent("tightmediumsel_SideBand2_", tightmediumselSideBand2elidx);
    fillhistinevent("tightmediumsel_SideBand3_", tightmediumselSideBand3elidx);
    fillhistinevent("mediumtightsel", mediumtightselelidx);
    fillhistinevent("mediumtightsel_Zwind_", mediumtightselZwindelidx);
    fillhistinevent("mediumtightsel_SideBand1_", mediumtightselSideBand1elidx);
    fillhistinevent("mediumtightsel_SideBand2_", mediumtightselSideBand2elidx);
    fillhistinevent("mediumtightsel_SideBand3_", mediumtightselSideBand3elidx);
    //if(isMC && noselelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAnosel", noselgenidx, noselelidx);

    // Clear all vector
    noselgenidx.clear();
    noselelidx.clear();
    noselZwindelidx.clear();
    leadepemselelidx.clear();
    leadepemonlyselelidx.clear();
    leadepemptgt20selelidx.clear();
    leadepemselZwindelidx.clear();
    leadepemZwindptgt20selelidx.clear();
    vetoselelidx.clear();
    looseselelidx.clear();
    looseselZwindelidx.clear();
    mediumselelidx.clear();
    mediumselZwindelidx.clear();
    tightselelidx.clear();
    tightselZwindelidx.clear();
    tightselSideBand1elidx.clear();
    tightselSideBand2elidx.clear();
    tightselSideBand3elidx.clear();
    tightnoneselelidx.clear();
    tightnoneselZwindelidx.clear();
    tightnoneselSideBand1elidx.clear();
    tightnoneselSideBand2elidx.clear();
    tightnoneselSideBand3elidx.clear();
    nonetightselelidx.clear();
    nonetightselZwindelidx.clear();
    nonetightselSideBand1elidx.clear();
    nonetightselSideBand2elidx.clear();
    nonetightselSideBand3elidx.clear();
    tightlooseselelidx.clear();
    tightlooseselZwindelidx.clear();
    tightlooseselSideBand1elidx.clear();
    tightlooseselSideBand2elidx.clear();
    tightlooseselSideBand3elidx.clear();
    loosetightselelidx.clear();
    loosetightselZwindelidx.clear();
    loosetightselSideBand1elidx.clear();
    loosetightselSideBand2elidx.clear();
    loosetightselSideBand3elidx.clear();
    tightmediumselelidx.clear();
    tightmediumselZwindelidx.clear();
    tightmediumselSideBand1elidx.clear();
    tightmediumselSideBand2elidx.clear();
    tightmediumselSideBand3elidx.clear();
    mediumtightselelidx.clear();
    mediumtightselZwindelidx.clear();
    mediumtightselSideBand1elidx.clear();
    mediumtightselSideBand2elidx.clear();
    mediumtightselSideBand3elidx.clear();
    
  } // End of loop on events

// Fill histograms for gen matching
void robustanalyzer::fillgenmchhistinevent(TString selection, vector<int> genidx, vector<int> elidx) {

  if(!isMC) throw "Error!!! Trying to fill gen info. for a non MC file.";
    
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
  
  int genmchfound = 0;
  if(el1!=-1) {
    genmchfound++;
    genelpt->Fill((*gen_pt)->at(gen1));
    geneleta->Fill((*gen_eta)->at(gen1));
    genelphi->Fill((*gen_phi)->at(gen1));
  }
  if(el2!=-1) {
    genmchfound++;
    genelpt->Fill((*gen_pt)->at(gen2));
    geneleta->Fill((*gen_eta)->at(gen2));
    genelphi->Fill((*gen_phi)->at(gen2));
  }
  genelmult->Fill(genmchfound);
  if(genmchfound==2) {
    TLorentzVector elec1, elec2;
    elec1.SetPtEtaPhiM((*gen_pt)->at(gen1),(*gen_eta)->at(gen1),(*gen_phi)->at(gen1),0.0005);
    elec2.SetPtEtaPhiM((*gen_pt)->at(gen2),(*gen_eta)->at(gen2),(*gen_phi)->at(gen2),0.0005);
    gendielM->Fill((elec1+elec2).M());
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
    genmchsctelpt->Fill((*ele_pt)->at(el1));
    genmchscteleta->Fill((*ele_eta)->at(el1));
    genmchsctelphi->Fill((*ele_phi)->at(el1));
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
    genmchsctelpt->Fill((*ele_pt)->at(el2));
    genmchscteleta->Fill((*ele_eta)->at(el2));
    genmchsctelphi->Fill((*ele_phi)->at(el2));
  }
  genmchsctelmult->Fill(genmchfound);
  if(genmchfound==2){
    TLorentzVector elec1, elec2;
    elec1.SetPtEtaPhiM((*ele_pt)->at(el1),(*ele_eta)->at(el1),(*ele_phi)->at(el1),0.0005);
    elec2.SetPtEtaPhiM((*ele_pt)->at(el2),(*ele_eta)->at(el2),(*ele_phi)->at(el2),0.0005);
    genmchsctdielM->Fill((elec1+elec2).M());
  }
  
}


// Function to do gen matching for to electrons
vector< pair<int,int> > robustanalyzer::diElecGenMatching(vector<int> genidx, vector<int> sctelidx) {

  if(genidx.size()!=2) {
    throw "Error! Analysis code only suitable for 2 gen electrons";
  }
  
  if(!isMC) {
    throw "Error! Cannot do gen matching. Not MC file.";
  }
  
  // Counter for the total no.of gen matches
  // and 1 history variable for the first gen matched position
  // Logic is not made to work with more than 2 gen el.
  vector< pair<int,int> > *gensctelmch = new vector< pair<int,int> >;
  gensctelmch->push_back(make_pair(genidx[0],-1));
  gensctelmch->push_back(make_pair(genidx[1],-1));
  
  // Find the sctelidx with the best angular match
  for(int elidx : sctelidx) {
    // Loop over the gen particles
    for(auto it=gensctelmch->begin(); it!=gensctelmch->end(); it++) {

      int genidx = (*it).first;
      int mchsctelidx = (*it).second;
      if(genidx==-1 || mchsctelidx!=-1) continue;

      double diffeta = abs((*ele_eta)->at(elidx)-(*gen_eta)->at(genidx)); 
      TLorentzVector vecsctel, vecgen;
      vecgen.SetPtEtaPhiM((*gen_pt)->at(genidx),(*gen_eta)->at(genidx),(*gen_phi)->at(genidx),0.0005);
      vecsctel.SetPtEtaPhiM((*ele_pt)->at(elidx),(*ele_eta)->at(elidx),(*ele_phi)->at(elidx),0.0005);
      double qdiffphi = ( (*gen_pdg)->at(genidx)/abs((*gen_pdg)->at(genidx)) )*( vecgen.DeltaPhi(vecsctel) );
      // Condition for gen matching
      if(abs((*ele_eta)->at(elidx))<1.479) {
	if(diffeta<0.2 && qdiffphi<0.05 && qdiffphi>-0.3) {
	  (*it).first = genidx;
	  (*it).second = elidx;
	}
      }
      else {
	if(diffeta<0.1 && qdiffphi<0.05 && qdiffphi>-0.2) {
	  (*it).first = genidx;
	  (*it).second = elidx;
	}
      }
    } // End of gen loop
    
  } // End of unseeded egamma object loop
  
  return (*gensctelmch);
}

// Return if selection in Z window
bool robustanalyzer::inZwind(TLorentzVector lead, TLorentzVector sublead) {

  double invm = (lead+sublead).M();
  if(abs(lead.Eta())<1.479 || abs(sublead.Eta())<1.479) {
    if(invm>86 && invm<97) {
      return true;
    }
    else {
      return false;
    }
  }
  else {
    if(invm>79 && invm<98) {
      return true;
    }
    else {
      return false;
    }
  }
}

// Return if selection in side-band
bool robustanalyzer::inSideBand1(TLorentzVector lead, TLorentzVector sublead) {

  double invm = (lead+sublead).M();
  if((invm>40 && invm<60) || invm>120) {
    return true;
  }
  else {
    return false;
  }
}

// Return if selection in side-band
bool robustanalyzer::inSideBand2(TLorentzVector lead, TLorentzVector sublead) {

  double invm = (lead+sublead).M();
  if((invm>45 && invm<60) || invm>120) {
    return true;
  }
  else {
    return false;
  }
}

// Return if selection in side-band
bool robustanalyzer::inSideBand3(TLorentzVector lead, TLorentzVector sublead) {

  double invm = (lead+sublead).M();
  if((invm>50 && invm<60) || invm>120) {
    return true;
  }
  else {
    return false;
  }
}
*/
