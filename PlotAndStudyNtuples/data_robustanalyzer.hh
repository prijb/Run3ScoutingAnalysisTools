#ifndef DATA_ROBUSTANALYZER_H
#define DATA_ROBUSTANALYZER_H

#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace std;

class data_robustanalyzer {
  
public:
  data_robustanalyzer(TString, TString, bool, bool, bool);
  ~data_robustanalyzer();

  double getL1decision(vector<unsigned int>, vector<unsigned int>);
  void analyzersinglefile(int);
  void addhist(TString);
  void addgenhist(TString);
  void addgenmchhist(TString);
  void fillhistinevent(TString, vector<int>, double);
  void fillgenhistinevent(TString, vector<int>);
  void fillgenmchhistinevent(TString, vector<int>, vector<int>, double);
  void sort(int*, TTreeReaderValue<std::vector<float>> *, int);
  pair<int,int> inZwindow(vector<int>);
  vector< pair<int,int> > diElecGenMatching(vector<int>, vector<int>);
    
  private:

  // Name list of the L1 triggers
  vector<string> l1name{"L1_DoubleMu_12_5", "L1_DoubleMu_15_7", "L1_HTT200er", "L1_HTT255er", "L1_HTT280er", "L1_HTT320er", "L1_HTT360er", "L1_ETT2000", "L1_HTT400er", "L1_HTT450er", "L1_SingleJet180", "L1_SingleJet200", "L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5", "L1_DoubleJet30er2p5_Mass_Min360_dEta_Max1p5", "L1_DoubleMu4p5er2p0_SQ_OS_Mass_Min7", "L1_DoubleMu4_SQ_OS_dR_Max1p2", "L1_SingleEG36er2p5", "L1_SingleLooseIsoEG28er2p1", "L1_SingleEG8er2p5", "L1_SingleEG10er2p5", "L1_SingleEG15er2p5", "L1_SingleEG26er2p5", "L1_SingleEG28_FWD2p5", "L1_DoubleEG4_er1p2_dR_Max0p9", "L1_DoubleEG4p5_er1p2_dR_Max0p9", "L1_DoubleEG5_er1p2_dR_Max0p9", "L1_DoubleEG5p5_er1p2_dR_Max0p8", "L1_DoubleEG7_er1p2_dR_Max0p8", "L1_DoubleEG7p5_er1p2_dR_Max0p7", "L1_DoubleEG_15_10_er2p5", "L1_DoubleEG_20_10_er2p5", "L1_DoubleEG_22_10_er2p5", "L1_DoubleEG_25_12_er2p5", "L1_DoubleEG_25_14_er2p5", "L1_DoubleEG_27_14_er2p5", "L1_DoubleEG_LooseIso22_12_er2p5", "L1_DoubleEG_LooseIso25_12_er2p5", "L1_TripleEG_18_17_8_er2p5", "L1_TripleEG_18_18_12_er2p5", "L1_TripleEG16er2p5", "L1_DoubleEG8er2p5_HTT300er"};
  vector<int> l1prescale{0, 1, 1600, 500, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 44000, 7600, 1350, 780, 45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  
  bool isDiEl;
  bool isDY;
  bool isMC;
    
  TTreeReader* tree;
  TTreeReaderValue<vector<bool>> *l1Result;
  TTreeReaderArray<float> *bsx;
  TTreeReaderArray<float> *bsy;
  TTreeReaderArray<float> *bsz;
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
  TTreeReaderValue<vector<bool>> *gen_isPFS;
  TTreeReaderValue<vector<bool>> *gen_isLastC;
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
  
  TFile* outfile;

  std::vector<TH1F*> all1dhists;
};

#endif
