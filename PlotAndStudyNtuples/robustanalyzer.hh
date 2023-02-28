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
  robustanalyzer(TString, TString, int, bool, bool, bool);
  ~robustanalyzer();

  void analyzersinglefile(int);
  void addhist(TString);
  void addgenhist(TString);
  void addgenmchhist(TString);
  void fillhistinevent(TString, vector<int>);
  void fillgenhistinevent(TString, vector<int>, vector<int>);
  void fillgenmchhistinevent(TString, vector<int>, vector<int>, double);
  void sort(int*, TTreeReaderValue<std::vector<float>> *, int);
  bool inZwind(TLorentzVector, TLorentzVector);
  bool inSideBand1(TLorentzVector, TLorentzVector);
  bool inSideBand2(TLorentzVector, TLorentzVector);
  bool inSideBand3(TLorentzVector, TLorentzVector);
  vector< pair<int,int> > diElecGenMatching(vector<int>, vector<int>);
    
  private:

  int nC;
  bool isDY;
  bool isMC;
  bool isJP;
    
  TTreeReader* tree;
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
  TTreeReaderValue<UInt_t> *n_ele;
  TTreeReaderValue<vector<float>> *ele_pt;
  TTreeReaderValue<vector<float>> *ele_eta;
  TTreeReaderValue<vector<float>> *ele_phi;
  TTreeReaderValue<vector<float>> *ele_m;
  TTreeReaderValue<vector<float>> *ele_trkpt;
  TTreeReaderValue<vector<float>> *ele_trketa;
  TTreeReaderValue<vector<float>> *ele_trkphi;
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
  
  TTreeReaderValue<vector<float>> *rho;

  TFile* outfile;

  std::vector<TH1F*> all1dhists;
};

#endif
