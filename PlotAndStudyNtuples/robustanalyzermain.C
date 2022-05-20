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
    
    stringstream ss1;
    ss1<<"hists_DoubleElectronGunPt1To300_"<<cnt<<".root";
    data_robustanalyzer drana_DoubleElectronGunPt1To300("./data/DoubleElectronGunPt1To300_ScoutingSkimSmall220518.root",ss1.str(), true);
    drana_DoubleElectronGunPt1To300.analyzersinglefile(cnt);
    /*
      stringstream ss2;
      ss2<<"hists_DoubleElectronGunPt1To300Old_"<<cnt<<".root";
      data_robustanalyzer drana_DoubleElectronGunPt1To300Old("./data/DoubleElectronGunPt1To300_ScoutingSkim220411.root",ss2.str(), true);
      drana_DoubleElectronGunPt1To300Old.analyzersinglefile(cnt);
      
      stringstream ss3;
      ss3<<"hists_DoubleElectronGunPt1To300220508Old_"<<cnt<<".root";
      data_robustanalyzer drana_DoubleElectronGunPt1To300220508Old("./data/DoubleElectronGunPt1To300_ScoutingSkim220508Old.root",ss3.str(), true);
      drana_DoubleElectronGunPt1To300220508Old.analyzersinglefile(cnt);
      
      stringstream ss2;
      ss2<<"hists_QCDPt20To30EmEnriched_"<<cnt<<".root";
      data_robustanalyzer drana_QCDPt20To30EmEnriched("./data/QCDPt20To30EmEnriched_ScoutingSkim220127.root",ss2.str(), true);
      drana_QCDPt20To30EmEnriched.analyzersinglefile(cnt);
      
      stringstream ss3;
      ss3<<"hists_QCDPt30To50EmEnriched_"<<cnt<<".root";
      data_robustanalyzer drana_QCDPt30To50EmEnriched("./data/QCDPt30To50EmEnriched_ScoutingSkim220127.root",ss3.str(), true);
      drana_QCDPt30To50EmEnriched.analyzersinglefile(cnt);
    */
  }
  catch (char const* exc) {
    cerr<<exc<<endl;
  }
  
  return -1;
}
