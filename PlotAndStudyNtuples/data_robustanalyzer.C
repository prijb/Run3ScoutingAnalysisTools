#include "data_robustanalyzer.hh"
#include <iostream>
#include <numeric>

#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

// Initialize and open the root file in the constructor
data_robustanalyzer::data_robustanalyzer(TString filename, TString outfilename, bool isDoubleElectron){

  isDiEl = isDoubleElectron;
  
  TFile *inpfile = TFile::Open(filename,"READ");
  cout<<"Initializing for file: "<<filename<<endl;

  tree = new TTreeReader("mmtree/tree",inpfile);
  bsx = new TTreeReaderArray<float>((*tree), "beamspot_x");
  bsy = new TTreeReaderArray<float>((*tree), "beamspot_y");
  bsz = new TTreeReaderArray<float>((*tree), "beamspot_z");
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
  int nosel=0, noselZwind=0, cut1=0, cut1Zwind=0;

  // Define the histograms
  addgenhist("nosel");
  addgenhist("ptetaminsel");
  addgenhist("ptminetabarsel");
  addgenhist("narZwindsel");
  addhist("nosel");
  addhist("noselZwind");
  addgenmchhist("noselAnosel");
  addgenmchhist("noselAptetaminsel");
  addgenmchhist("noselZwindAptetaminsel");
  addgenmchhist("noselZwindAnarZwindsel");
  addhist("cut1");
  addhist("cut1Zwind");

  // vector of electron indices
  vector<int> noselgenidx;
  vector<int> ptetaminselgenidx;
  vector<int> ptminetabarselgenidx;
  vector<int> narZwindselgenidx;
  vector<int> noselelidx;
  vector<int> noselZwindelidx;
  vector<int> cut1elidx;
  vector<int> cut1Zwindelidx;
  
  // Loop beginning on events
  while(tree->Next()) {

    event++;
    //if(event>100) break;
    //if(event!=283991 && event!=326114) continue;
    if(event%10000==0) std::cout<<"Processed event: "<<event+1<<std::endl;

    if((*(*n_ele))<0) cout<<"Error!!! Wrong technical event processing. Negative number of electrons in event."<<endl;;
    if((*(*n_ele))==0) continue;
    
    // Loop on Gen particles to select good gen electrons from Z
    for(unsigned int gen_ctr=0; gen_ctr<(*(*n_gen)); gen_ctr++) {

      // Check that the gen particle is electron coming from a mother Z
      bool mincond = true;
      mincond *= TMath::Abs((*gen_pdg)->at(gen_ctr))==11;
      mincond *= (*gen_nmoms)->at(gen_ctr)==1;
      mincond *= (*gen_mompdg)->at(gen_ctr)==23;
      if(!mincond) continue;

      noselgenidx.push_back(gen_ctr);

      double ptetaminselcond = true;
      ptetaminselcond *= (*gen_pt)->at(gen_ctr)>5;
      ptetaminselcond *= TMath::Abs((*gen_eta)->at(gen_ctr))<2.3;
      if(ptetaminselcond) ptetaminselgenidx.push_back(gen_ctr);

      double ptminetabarselcond = true;
      ptminetabarselcond *= (*gen_pt)->at(gen_ctr)>5;
      ptminetabarselcond *= TMath::Abs((*gen_eta)->at(gen_ctr))<1.4 || TMath::Abs((*gen_eta)->at(gen_ctr))>1.6;
      ptminetabarselcond *= TMath::Abs((*gen_eta)->at(gen_ctr))<2.3;
      if(ptminetabarselcond) ptminetabarselgenidx.push_back(gen_ctr);
    }
    
    // Sort the electrons based on their pT
    vector<int> sortedelidx((*(*n_ele)));
    iota(begin(sortedelidx), end(sortedelidx), 0);
    sort(&sortedelidx[0], ele_pt, (*(*n_ele))); // Verified that the algorithm works fine

    // Loop on electrons in the event loop
    for(unsigned int ele_ctr=0; ele_ctr<(*(*n_ele)); ele_ctr++) {

      // Take the sorted index only
      unsigned int elidx = sortedelidx[ele_ctr];

      noselelidx.push_back(elidx);

      bool cut1elpass = true;
      cut1elpass *= TMath::Abs((*ele_eta)->at(elidx))<2.5;
      cut1elpass *= TMath::Abs((*ele_eta)->at(elidx))<1.479?(*ele_sigmaietaieta)->at(elidx)<0.0126:(*ele_sigmaietaieta)->at(elidx)<0.0457;
      if(cut1elpass) cut1elidx.push_back(elidx);
      
    }// End of loop on electrons in the event loop

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
        
    fillgenhistinevent("nosel", noselgenidx);
    if(ptetaminselgenidx.size()==2) fillgenhistinevent("ptetaminsel", ptetaminselgenidx);
    if(ptminetabarselgenidx.size()==2) {
      fillgenhistinevent("ptminetabarsel", ptminetabarselgenidx);
      TLorentzVector el1, el2;
      el1.SetPtEtaPhiM((*gen_pt)->at(ptminetabarselgenidx[0]),(*gen_eta)->at(ptminetabarselgenidx[0]),(*gen_phi)->at(ptminetabarselgenidx[0]),0.0005);
      el2.SetPtEtaPhiM((*gen_pt)->at(ptminetabarselgenidx[1]),(*gen_eta)->at(ptminetabarselgenidx[1]),(*gen_phi)->at(ptminetabarselgenidx[1]),0.0005);
      double invM = (el1+el2).M();
      if(invM>85 && invM<100) {
	narZwindselgenidx.push_back(ptminetabarselgenidx[0]);
	narZwindselgenidx.push_back(ptminetabarselgenidx[1]);
      }
    }
    if(narZwindselgenidx.size()==2) fillgenhistinevent("narZwindsel", narZwindselgenidx);
    if(noselelidx.size()>0) nosel++;
    fillhistinevent("nosel", noselelidx);
    if(noselZwindelidx.size()>0) noselZwind++;
    fillhistinevent("noselZwind", noselZwindelidx);
    if(noselelidx.size()>0 && noselgenidx.size()>0) fillgenmchhistinevent("noselAnosel", noselgenidx, noselelidx);
    if(noselelidx.size()>0 && ptetaminselgenidx.size()==2) fillgenmchhistinevent("noselAptetaminsel", ptetaminselgenidx, noselelidx);
    if(noselZwindelidx.size()>=2 && ptetaminselgenidx.size()==2) fillgenmchhistinevent("noselZwindAptetaminsel", ptetaminselgenidx, noselZwindelidx);
    if(noselZwindelidx.size()>=2 && narZwindselgenidx.size()==2) fillgenmchhistinevent("noselZwindAnarZwindsel", narZwindselgenidx, noselZwindelidx);
    if(cut1elidx.size()>0) cut1++;
    fillhistinevent("cut1", cut1elidx);
    if(cut1Zwindelidx.size()>0) cut1Zwind++;
    fillhistinevent("cut1Zwind", cut1Zwindelidx);

    // Clear all vector
    noselgenidx.clear();
    ptetaminselgenidx.clear();
    ptminetabarselgenidx.clear();
    narZwindselgenidx.clear();
    noselelidx.clear();
    noselZwindelidx.clear();
    cut1elidx.clear();
    cut1Zwindelidx.clear();
    
  } // End of loop on events

  cout<<totEntries<<"\t"<<endevent-beginevent<<"\t"<<nosel<<"\t"<<noselZwind<<"\t"<<cut1<<"\t"<<cut1Zwind<<endl;
}

// Fill histograms for gen matching
void data_robustanalyzer::fillgenmchhistinevent(TString selection, vector<int> genidx, vector<int> elidx) {

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
  pair< pair<int,int>, pair<int,int> > sctelgenmch = diElecGenMatching(genidx, elidx);

  int el1 = sctelgenmch.first.first;
  int gen1 = sctelgenmch.first.second;
  int el2 = sctelgenmch.second.first;
  int gen2 = sctelgenmch.second.second;

  int genmchfound = 0;
  if(gen1!=-1) {
    genmchfound++;
    genelpt->Fill((*gen_pt)->at(gen1));
    geneleta->Fill((*gen_eta)->at(gen1));
    genelphi->Fill((*gen_phi)->at(gen1));
  }
  if(gen2!=-1) {
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
void data_robustanalyzer::fillgenhistinevent(TString selection, vector<int> genidx) {

  if(genidx.size()==0) return;

  TH1F* elmult = (TH1F*) outfile->Get(selection+"gen_elmult");
  TH1F* elpt = (TH1F*) outfile->Get(selection+"gen_elpt");
  TH1F* eleta = (TH1F*) outfile->Get(selection+"gen_eleta");
  TH1F* elphi = (TH1F*) outfile->Get(selection+"gen_elphi");
  TH1F* dielM = (TH1F*) outfile->Get(selection+"gen_dielM");


  elmult->Fill(genidx.size());
  for(unsigned int ctr=0; ctr<genidx.size(); ctr++) {
    elpt->Fill((*gen_pt)->at(genidx[ctr]));
    eleta->Fill((*gen_eta)->at(genidx[ctr]));
    elphi->Fill((*gen_phi)->at(genidx[ctr]));
  }
  if(genidx.size()==2) {
    TLorentzVector el1, el2;
    el1.SetPtEtaPhiM((*gen_pt)->at(genidx[0]),(*gen_eta)->at(genidx[0]),(*gen_phi)->at(genidx[0]),0.0005);
    el2.SetPtEtaPhiM((*gen_pt)->at(genidx[1]),(*gen_eta)->at(genidx[1]),(*gen_phi)->at(genidx[1]),0.0005);
    dielM->Fill((el1+el2).M());
  }
  else {
    cout<<"Warning!!! Unequal to two electron found in event."<<endl;
  }
  
}

// Function to fill a set of histograms for scouting electrons
void data_robustanalyzer::fillhistinevent(TString selection, vector<int> elidx) {

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
void data_robustanalyzer::addgenmchhist(TString selection) {
  
  // Variables before gen match
  all1dhists.push_back(new TH1F(selection+"genelsctbar_dEta","#Delta#eta(gen e, sct. e)",4000,-2,2));
  all1dhists.push_back(new TH1F(selection+"genelsctbar_qdPhi","q#Delta#phi(gen e, sct. e)",4000,-2,2));
  all1dhists.push_back(new TH1F(selection+"genelsctee_dEta","#Delta#eta(gen e, sct. e)",4000,-2,2));
  all1dhists.push_back(new TH1F(selection+"genelsctee_qdPhi","q#Delta#phi(gen e, sct. e)",4000,-2,2));
  
  // Variables after gen match
  all1dhists.push_back(new TH1F(selection+"genelsctmchbar_dEta","#Delta#eta(gen e, sct. e)",4000,-2,2));
  all1dhists.push_back(new TH1F(selection+"genelsctmchbar_qdPhi","q#Delta#phi(gen e, sct. e)",4000,-2,2));
  all1dhists.push_back(new TH1F(selection+"genelsctmchee_dEta","#Delta#eta(gen e, sct. e)",4000,-2,2));
  all1dhists.push_back(new TH1F(selection+"genelsctmchee_qdPhi","q#Delta#phi(gen e, sct. e)",4000,-2,2));

  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elmult","gen N e/#gamma",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elpt","gen e p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_eleta","gen e #eta",600,-3,3));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_elphi","gen e #phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"sctmchgenel_dielM","M(e,e)",1000,-10,990));

  all1dhists.push_back(new TH1F(selection+"genmchsct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"genmchsct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"genmchsct_eleta","#eta",600,-3,3));
  all1dhists.push_back(new TH1F(selection+"genmchsct_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"genmchsct_dielM","all M(e,e)",1000,-10,990));
}

// Function to add a set of histograms for gen electrons
void data_robustanalyzer::addgenhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"gen_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"gen_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"gen_eleta","#eta",600,-3,3));
  all1dhists.push_back(new TH1F(selection+"gen_elphi","#phi",66,-3.3,3.3));
  all1dhists.push_back(new TH1F(selection+"gen_dielM","all M(e,e)",1000,-10,990));

}

// Function to add a set of histograms for scouting electrons
void data_robustanalyzer::addhist(TString selection) {

  all1dhists.push_back(new TH1F(selection+"sct_elmult","N e",50,-5,45));
  all1dhists.push_back(new TH1F(selection+"sct_elpt","p_{T} / GeV",1000,-10,990));
  all1dhists.push_back(new TH1F(selection+"sct_eleta","#eta",600,-3,3));
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

// Function to do gen matching
// Currently doing genmatching for lead pT only in mueg type events
// Modify this to return two indices for egeg type events
pair< pair<int,int>,pair<int,int> > data_robustanalyzer::diElecGenMatching(vector<int> genidx, vector<int> elidx) {

  if(genidx.size()>2) { // Less than two is allowed to take care of gen level selections
    cout<<"Warning. Code can do gen matching for 2e at best. Skipping gen matching."<<endl;
    return make_pair(make_pair(-1,-1),make_pair(-1,-1));
  }
  
  if(genidx.size()==0 || elidx.size()==0) {
    return make_pair(make_pair(-1,-1),make_pair(-1,-1));
  }

  // Counter for the total no.of gen matches
  // and 1 history variable for the first gen matched position
  // Logic does not work with more than 2 gen el.
  int genmchcnt=0, genpos1=-1;
  pair< pair<int,int>, pair<int,int> > genelmch = make_pair(make_pair(-1,-1),make_pair(-1,-1));

  // Find the elidx with the best angular match
  for(unsigned int el=0; el<elidx.size(); el++) {

    // Variable to store the position of gen electron if there is a match
    int genposhistory = -1;
    // Variable to store the angle if there is a match
    double anglediffhistory = 0;
    
    // Loop over the gen particles
    for(unsigned int gen=0; gen<genidx.size(); gen++) {
      double diffeta = abs((*ele_eta)->at(elidx[el])-(*gen_eta)->at(genidx[gen]));
      if(gen==genpos1) continue; // Continue to the next gen electron if this electron has already been matched

      if(abs((*ele_eta)->at(elidx[el]))<1.479) {
	if(diffeta<0.06) {
	  if(genposhistory==-1) {
	    genposhistory = gen;
	    anglediffhistory = diffeta;
	    continue;
	  }
	  if(genposhistory!=-1 && diffeta<anglediffhistory) {
	    genposhistory = gen;
	    anglediffhistory = diffeta;
	    continue;
	  }
	}
      }
      else {
	if(diffeta<0.04) {
	  if(genposhistory==-1) {
	    genposhistory = gen;
	    anglediffhistory = diffeta;
	    continue;
	  }
	  if(genposhistory!=-1 && diffeta<anglediffhistory) {
	    genposhistory = gen;
	    anglediffhistory = diffeta;
	    continue;
	  }
	}
      }
    } // End of gen loop

    if(genposhistory!=-1) { // Found a gen match
      genmchcnt++;
      if(genmchcnt==1) {
	genelmch.first = make_pair(elidx[el],genidx[genposhistory]);
	genpos1 = genposhistory;
      }
      if(genmchcnt==2) {
	genelmch.second = make_pair(elidx[el],genidx[genposhistory]);
      }
      if(genmchcnt>2) {
	cout<<"Error! Shouldn't be possible to have more than two gen matches"<<endl;
      }
    }

    if(genmchcnt==2) break; // Break if two gen matches are found
    
  } // End of scouting electron object loop

  return genelmch;
}
