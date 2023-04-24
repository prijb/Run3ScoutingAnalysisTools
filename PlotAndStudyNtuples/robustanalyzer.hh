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
  robustanalyzer(TString, TString, int);
  ~robustanalyzer();

  void analyzersinglefile(int);
  void addhist(TString);
  void fillhistinevent(TString, vector<int>);
  void sort(int*, TTreeReaderValue<std::vector<float>> *, int);
    
  private:

  int nC;
  
  TTreeReader* tree;
  TTreeReaderValue<unsigned int> *hlt_isomu27;
  TTreeReaderValue<unsigned int> *hlt_mu50;
  TTreeReaderValue<unsigned int> *hlt_scouting;
  TTreeReaderValue<unsigned int> *n_ele;
  TTreeReaderValue<vector<float>> *ele_pt;
  TTreeReaderValue<vector<float>> *ele_eta;
  TTreeReaderValue<vector<float>> *ele_phi;
  TTreeReaderValue<vector<float>> *ele_m;
  TTreeReaderValue<vector<float>> *ele_d0;
  TTreeReaderValue<vector<float>> *ele_dz;
  TTreeReaderValue<vector<float>> *ele_detain;
  TTreeReaderValue<vector<float>> *ele_dphiin;
  TTreeReaderValue<vector<float>> *ele_sigmaietaieta;
  TTreeReaderValue<vector<float>> *ele_hoe;
  TTreeReaderValue<vector<float>> *ele_ooemoop;
  TTreeReaderValue<vector<int>> *ele_mhits;
  TTreeReaderValue<vector<int>> *ele_charge;
  TTreeReaderValue<vector<float>> *ele_ecaliso;
  TTreeReaderValue<vector<float>> *ele_hcaliso;
  TTreeReaderValue<vector<float>> *ele_tkiso;
  TTreeReaderValue<vector<float>> *ele_r9;
  TTreeReaderValue<vector<float>> *ele_smin;
  TTreeReaderValue<vector<float>> *ele_smaj;

  TTreeReaderValue<unsigned int> *n_rho;
  TTreeReaderValue<vector<float>> *rho;
  
  TFile* outfile;

  std::vector<TH1F*> all1dhists;
};

#endif
