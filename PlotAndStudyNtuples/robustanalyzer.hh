#ifndef ROBUSTANALYZER_H
#define ROBUSTANALYZER_H

#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TLorentzVector.h"

using namespace std;

class robustanalyzer {
  
public:
  robustanalyzer(TString, TString, int, bool);
  ~robustanalyzer();
  
  void sort(int*, TTreeReaderValue<std::vector<float>> *, int);

  void analyzersinglefile(int);

  void addhist(TString);
  void fillhistinevent(TString, vector<int>);
  
  void addGenCollections();
  void addgenhist(TString);
  void fillgenhistinevent(TString, vector<int>, vector<int>);

  
private:
  
  int nC;
  
  TTreeReader* tree;
  TTreeReaderValue<unsigned int> *hlt_isomu27;
  TTreeReaderValue<unsigned int> *hlt_mu50;
  TTreeReaderValue<unsigned int> *hlt_scouting;
  TTreeReaderValue<unsigned int> *hlt_sct_dimu3;
  TTreeReaderValue<unsigned int> *hlt_sct_dieg;
  TTreeReaderValue<unsigned int> *hlt_sct_eg30;
  TTreeReaderValue<unsigned int> *hlt_sct_jet;
  TTreeReaderValue<unsigned int> *hltmu_sct;
  TTreeReaderValue<unsigned int> *hltoth_sctmnt;
  TTreeReaderValue<unsigned int> *n_ele;
  TTreeReaderValue<vector<float>> *ele_pt;
  TTreeReaderValue<vector<float>> *ele_eta;
  TTreeReaderValue<vector<float>> *ele_phi;
  TTreeReaderValue< vector<vector<float>> > *ele_trkpt;
  TTreeReaderValue< vector<vector<float>> > *ele_trketa;
  TTreeReaderValue< vector<vector<float>> > *ele_trkphi;
  TTreeReaderValue< vector<vector<float>> > *ele_trkd0;
  TTreeReaderValue< vector<vector<float>> > *ele_trkdz;
  TTreeReaderValue< vector<vector<int>> > *ele_trkcharge;
  TTreeReaderValue< vector<vector<float>> > *ele_trkchi2;
  TTreeReaderValue<vector<float>> *ele_m;
  TTreeReaderValue<vector<float>> *ele_detain;
  TTreeReaderValue<vector<float>> *ele_dphiin;
  TTreeReaderValue<vector<float>> *ele_sigmaietaieta;
  TTreeReaderValue<vector<float>> *ele_hoe;
  TTreeReaderValue<vector<float>> *ele_ooemoop;
  TTreeReaderValue<vector<int>> *ele_mhits;
  TTreeReaderValue<vector<float>> *ele_ecaliso;
  TTreeReaderValue<vector<float>> *ele_hcaliso;
  TTreeReaderValue<vector<float>> *ele_tkiso;
  TTreeReaderValue<vector<float>> *ele_r9;
  TTreeReaderValue<vector<float>> *ele_smin;
  TTreeReaderValue<vector<float>> *ele_smaj;

  TTreeReaderValue<unsigned int> *n_rho;
  TTreeReaderValue<vector<float>> *rho;
  
  bool isMC;
  TTreeReaderValue<unsigned int> *n_genlep;
  TTreeReaderValue<vector<int>> *genlep_pdg;
  TTreeReaderValue<vector<float>> *genlep_pt;
  TTreeReaderValue<vector<float>> *genlep_eta;
  TTreeReaderValue<vector<float>> *genlep_phi;
  TTreeReaderValue<vector<float>> *genlep_m;
  TTreeReaderValue<vector<float>> *genlep_vx;
  TTreeReaderValue<vector<float>> *genlep_vy;
  TTreeReaderValue<vector<float>> *genlep_vz;
  TTreeReaderValue<unsigned int> *n_genjpsi;
  TTreeReaderValue<vector<int>> *genjpsi_pdg;
  TTreeReaderValue<vector<float>> *genjpsi_pt;
  TTreeReaderValue<vector<float>> *genjpsi_eta;
  TTreeReaderValue<vector<float>> *genjpsi_phi;
  TTreeReaderValue<vector<float>> *genjpsi_m;
  TTreeReaderValue<vector<float>> *genjpsi_vx;
  TTreeReaderValue<vector<float>> *genjpsi_vy;
  TTreeReaderValue<vector<float>> *genjpsi_vz;

  TFile* outfile;

  std::vector<TH1F*> all1dhists;
};

#endif
