#include "data_robustanalyzer.hh"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {

  stringstream ss;
  ss << argv[1];
  int cnt;
  ss >> cnt;

  try {
    stringstream ssjpsi;
    ssjpsi<<"hists_JPsi_"<<cnt<<".root";
    data_robustanalyzer drana_JPsi("./data/Scouting_JPsiToEE_nTuple_220206_100k.root",ssjpsi.str(), true, false, true);
    drana_JPsi.analyzersinglefile(cnt);
    /*    
    stringstream ss1;
    ss1<<"hists_DoubleElectronGunPt1To300_"<<cnt<<".root";
    data_robustanalyzer drana_DoubleElectronGunPt1To300("./data/DoubleElectronGun_ScoutingSkimSmall220717.root",ss1.str(), true, false, true);
    drana_DoubleElectronGunPt1To300.analyzersinglefile(cnt);

    stringstream ss1l1v1;
    ss1l1v1<<"hists_DoubleElectronGunPt1To300_1240V74L1v1_"<<cnt<<".root";
    data_robustanalyzer drana_DoubleElectronGunPt1To300L1v1("./data/DoubleElectronGun_ScoutingSkimSmall220717_1240V74L1v1.root",ss1l1v1.str(), true, false, true);
    drana_DoubleElectronGunPt1To300L1v1.analyzersinglefile(cnt);

    stringstream ss1l1v2;
    ss1l1v2<<"hists_DoubleElectronGunPt1To300_1240V74L1v2_"<<cnt<<".root";
    data_robustanalyzer drana_DoubleElectronGunPt1To300L1v2("./data/DoubleElectronGun_ScoutingSkimSmall220717_1240V74L1v2.root",ss1l1v2.str(), true, false, true);
    drana_DoubleElectronGunPt1To300L1v2.analyzersinglefile(cnt);

    stringstream ss2;
    ss2<<"hists_DYToLLM4To50_"<<cnt<<".root";
    data_robustanalyzer drana_DYToLLM4To50("./data/DYToLLM4To50_ScoutingSkimSmall220717.root",ss2.str(), false, true, true);
    drana_DYToLLM4To50.analyzersinglefile(cnt);

    stringstream ss2l1v1;
    ss2l1v1<<"hists_DYToLLM4To50_1240V74L1v1_"<<cnt<<".root";
    data_robustanalyzer drana_DYToLLM4To50l1v1("./data/DYToLLM4To50_ScoutingSkimSmall220717_1240V74L1v1.root",ss2l1v1.str(), false, true, true);
    drana_DYToLLM4To50l1v1.analyzersinglefile(cnt);

    stringstream ss2l1v2;
    ss2l1v2<<"hists_DYToLLM4To50_1240V74L1v2_"<<cnt<<".root";
    data_robustanalyzer drana_DYToLLM4To50l1v2("./data/DYToLLM4To50_ScoutingSkimSmall220717_1240V74L1v2.root",ss2l1v2.str(), false, true, true);
    drana_DYToLLM4To50l1v2.analyzersinglefile(cnt);
    */
    //stringstream ssdata;
    //ssdata<<"hists_ZeroBias2018D_"<<cnt<<".root";
    //data_robustanalyzer drana_ZeroBias2018D("./data/ZeroBias2018D_ScoutingSkimSmall220518.root",ssdata.str(), false, false, false);
    //drana_ZeroBias2018D.analyzersinglefile(cnt);
 
    //stringstream ssephdata;
    //ssephdata<<"hists_Ephemeral1HLTPhysics2018D_"<<cnt<<".root";
    //data_robustanalyzer drana_Ephemeral1HLTPhysics2018D("./data/Ephemeral1HLTPhysics2018D_ScoutingSkimSmall220520_0.root",ssephdata.str(), false, false, false);
    //drana_Ephemeral1HLTPhysics2018D.analyzersinglefile(cnt);
  }
  
  catch (char const* exc) {
    cerr<<exc<<endl;
  }
  
  return -1;
}
