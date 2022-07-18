#include "robustanalyzer.hh"
#include <iostream>
#include <numeric>
#include <boost/range/combine.hpp>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
robustanalyzer::robustanalyzer(TString filename, TString outfilename, bool isDoubleElectron, bool isMC){

  isDiEl = isDoubleElectron;
  isMC = isMC;
  
  TFile *inpfile = TFile::Open(filename,"READ");
  cout<<"Initializing for file: "<<filename<<endl;

  tree = new TTreeReader("mmtree/tree",inpfile);
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
robustanalyzer::~robustanalyzer() {

  outfile->Write();
  outfile->Close();
}

// Analyzer for a single file
void robustanalyzer::analyzersinglefile(int splitCnt) { // Assume splitCnt to range from 0 to nCores

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
  int nosel=0, noselZwind=0, cut1=0, cut1Zwind=0;

  // Define the histograms
  if(isMC) {
    addgenhist("noselgen");
    addgenmchhist("noselgenAnosel");
  }
  addhist("nosel");
  
  // vector of electron indices
  vector<int> noselgenidx;
  vector<int> noselelidx;
  
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>1000) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    // Loop on Gen particles to select good gen electrons
    if(isMC) {
      for(unsigned int gen_ctr=0; gen_ctr<(*(*n_gen)); gen_ctr++) {
	noselgenidx.push_back(gen_ctr);
      } // End of loop on gen electrons
      
      fillgenhistinevent("noselgen", noselgenidx);
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
    }// End of loop on electrons in the event loop

    if(noselelidx.size()>0) nosel++;
    fillhistinevent("nosel", noselelidx);
    if(isMC && noselelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselgenAnosel", noselgenidx, noselelidx);

    // Clear all vector
    noselgenidx.clear();
    noselelidx.clear();
    
  } // End of loop on events

  cout<<totEntries<<"\t"<<endevent-beginevent<<"\t"<<nosel<<"\t"<<noselZwind<<"\t"<<cut1<<"\t"<<cut1Zwind<<endl;
}

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

// Function to fill a set of histograms for gen particles
void robustanalyzer::fillgenhistinevent(TString selection, vector<int> genidx) {

  if(!isMC) throw "Error!!! Trying to fill gen info. for a non MC file.";
    
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
    leadelpt->Fill((*gen_pt)->at(genidx[leadptpos]));
    leadeleta->Fill((*gen_eta)->at(genidx[leadptpos]));
    leadelphi->Fill((*gen_phi)->at(genidx[leadptpos]));
  }
  if(subleadptpos!=-1) {
    subleadelpt->Fill((*gen_pt)->at(genidx[subleadptpos]));
    subleadeleta->Fill((*gen_eta)->at(genidx[subleadptpos]));
    subleadelphi->Fill((*gen_phi)->at(genidx[subleadptpos]));
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
  else {
    cout<<"Warning!!! Unequal to two electron found in event."<<endl;
  }

}

// Function to fill a set of histograms for scouting electrons
void robustanalyzer::fillhistinevent(TString selection, vector<int> elidx) {

  if(elidx.size()==0) return;

  TH1F* elmult = (TH1F*) outfile->Get(selection+"sct_elmult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"sct_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"sct_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"sct_elphi");
  TH1F* dielM = (TH1F*) outfile->Get(selection+"sct_dielM");
  
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

  elmult->Fill(elidx.size());
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

  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elmult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elpt","gen e p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_eleta","gen e #eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elphi","gen e #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_dielM","M(e,e)",1000,-10,990));

  all1dhists.push_back(new TH1F(selection+"genmchsct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"genmchsct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"genmchsct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"genmchsct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"genmchsct_dielM","all M(e,e)",1000,-10,990));
}

// Function to add a set of histograms for gen electrons
void robustanalyzer::addgenhist(TString selection) {

  if(!isMC) throw "Error!!! Trying to define gen hists for a non MC file.";
    
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
void robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_eleta","#eta",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sct_dielM","all M(e,e)",1000,-10,990));

  all1dhists.push_back(new TH1F(selection+"sctbar_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctbar_elm","m / GeV",1000,-1e-5,1e-5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eld0","d_{0} / cm",20000,-10,10));
  all1dhists.push_back(new TH1F(selection+"sctbar_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldz","d_{z} / cm",4000,-20,20));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldetain","#Delta#eta_{in}",10000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_eldphiin","#Delta#phi_{in}",10000,0,1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elsigmaietaieta","#sigmai#etai#eta",1000,0,0.1));
  all1dhists.push_back(new TH1F(selection+"sctbar_elhoe","H/E",20000,0,0.2));
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
  all1dhists.push_back(new TH1F(selection+"sctec_elhoe","H/E",20000,0,0.2));
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

// Function to find an electron pair in the Z mass window 80<M<100
pair<int,int> robustanalyzer::inZwindow(vector<int> elidx) {
  
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
