#ifndef ROBUSTANALYZER_H
#define ROBUSTANALYZER_H

#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;

class robustanalyzer {
  
public:
  robustanalyzer(TString, TString, int, bool, bool);
  ~robustanalyzer();

  void analyzersinglefile(int);
  void addhist(TString);
  void addgenhist(TString);
  void addgenmchhist(TString);
  void fillhistinevent(TString, vector<int>);
  void fillgenhistinevent(TString, vector<int>);
  void fillgenmchhistinevent(TString, vector<int>, vector<int>);
  void sort(int*, TTreeReaderValue<std::vector<float>> *, int);
  pair<int,int> inZwindow(vector<int>);
  vector< pair<int,int> > diElecGenMatching(vector<int>, vector<int>);
    
  private:

  int nC;
  bool isDiEl;
  bool isMC;
    
  TTreeReader* tree;
  TTreeReaderValue<unsigned int> *n_gen;
  TTreeReaderValue<vector<int>> *gen_pdg;
  TTreeReaderValue<vector<float>> *gen_pt;
  TTreeReaderValue<vector<float>> *gen_eta;
  TTreeReaderValue<vector<float>> *gen_phi;
  TTreeReaderValue<vector<float>> *gen_vx;
  TTreeReaderValue<vector<float>> *gen_vy;
  TTreeReaderValue<vector<float>> *gen_vz;
  TTreeReaderValue<vector<int>> *gen_nmoms;
  TTreeReaderValue<vector<int>> *gen_mompdg;
  TTreeReaderValue<UInt_t> *n_ele;
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
  TTreeReaderValue<vector<unsigned int>> *ele_seedid;
  TTreeReaderValue<vector<vector<float>>> *ele_enemat;
  TTreeReaderValue<vector<vector<float>>> *ele_timmat;
  TTreeReaderValue<UInt_t> *n_rho;
  TTreeReaderValue<vector<float>> *rho;
  
  TFile* outfile;

  std::vector<TH1F*> all1dhists;
};

#endif
