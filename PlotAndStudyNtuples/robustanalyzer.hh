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
  
  void analyzerelectrons();
  void addhist_electron(TString);
  void fillhistinevent_electron(TString, vector<int>);
  void addhist_oflsct_electron(TString);
  void fillhistinevent_oflsct_electron(TString, vector<int>, vector<int>);
  bool offline_scouting_angmch(double, double);
  double effectivearea(double);
  bool vetosel_offline_electron(int);
  bool vetoecalsel_offline_electron(int);

  void analyzerphotons();
  void addhist_photon(TString);
  void fillhistinevent_photon(TString, vector<int>);

  void sort(int*, TTreeReaderValue<std::vector<float>> *, int);
  
  private:

  int nC;
  
  TTreeReader* tree;
  TTreeReaderValue<unsigned int> *l1_2mu;
  TTreeReaderValue<unsigned int> *l1_2musqos;
  TTreeReaderValue<unsigned int> *l1_met;
  TTreeReaderValue<unsigned int> *l1_1jet;
  TTreeReaderValue<unsigned int> *l1_2jet;
  TTreeReaderValue<unsigned int> *l1_1eg;
  TTreeReaderValue<unsigned int> *l1_2eg;
  TTreeReaderValue<unsigned int> *hlt_isomu27;
  TTreeReaderValue<unsigned int> *hlt_mu50;
  TTreeReaderValue<unsigned int> *hlt_scouting;
  
  TTreeReaderValue<unsigned int> *n_ofle;
  TTreeReaderValue<vector<float>> *ofle_pt;
  TTreeReaderValue<vector<float>> *ofle_eta;
  TTreeReaderValue<vector<float>> *ofle_phi;
  TTreeReaderValue<vector<float>> *ofle_d0;
  TTreeReaderValue<vector<float>> *ofle_dz;
  TTreeReaderValue<vector<float>> *ofle_detain;
  TTreeReaderValue<vector<float>> *ofle_dphiin;
  TTreeReaderValue<vector<float>> *ofle_sigmaietaieta;
  TTreeReaderValue<vector<float>> *ofle_hoe;
  TTreeReaderValue<vector<float>> *ofle_ooemoop;
  TTreeReaderValue<vector<int>> *ofle_mhits;
  TTreeReaderValue<vector<int>> *ofle_charge;
  TTreeReaderValue<vector<float>> *ofle_phoiso;
  TTreeReaderValue<vector<float>> *ofle_nthadiso;
  TTreeReaderValue<vector<float>> *ofle_cghadiso;
  TTreeReaderValue<vector<bool>> *ofle_conveto;

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

  TTreeReaderValue<unsigned int> *n_pho;
  TTreeReaderValue<vector<float>> *pho_pt;
  TTreeReaderValue<vector<float>> *pho_eta;
  TTreeReaderValue<vector<float>> *pho_phi;
  TTreeReaderValue<vector<float>> *pho_m;
  TTreeReaderValue<vector<float>> *pho_sigmaietaieta;
  TTreeReaderValue<vector<float>> *pho_hoe;
  TTreeReaderValue<vector<float>> *pho_ecaliso;
  TTreeReaderValue<vector<float>> *pho_hcaliso;
  TTreeReaderValue<vector<float>> *pho_r9;
  TTreeReaderValue<vector<float>> *pho_smin;
  TTreeReaderValue<vector<float>> *pho_smaj;

  TTreeReaderValue<unsigned int> *n_rho;
  TTreeReaderValue<vector<float>> *rho;
  TTreeReaderValue<unsigned int> *n_oflrho;
  TTreeReaderValue<vector<float>> *oflrho;
  
  TFile* outfile;

  std::vector<TH1F*> all1dhists;
};

#endif
