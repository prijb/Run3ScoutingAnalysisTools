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
  
  stringstream ss1;
  ss1<<"hists_DYToLLM50_"<<cnt<<".root";
  data_robustanalyzer drana_DYToLLM50("./data/DYToLLM50_ScoutingSkim220127.root",ss1.str(), true);
  drana_DYToLLM50.analyzersinglefile(cnt);

  stringstream ss2;
  ss2<<"hists_QCDPt20To30EmEnriched_"<<cnt<<".root";
  data_robustanalyzer drana_QCDPt20To30EmEnriched("./data/QCDPt20To30EmEnriched_ScoutingSkim220127.root",ss2.str(), true);
  drana_QCDPt20To30EmEnriched.analyzersinglefile(cnt);

  stringstream ss3;
  ss3<<"hists_QCDPt30To50EmEnriched_"<<cnt<<".root";
  data_robustanalyzer drana_QCDPt30To50EmEnriched("./data/QCDPt30To50EmEnriched_ScoutingSkim220127.root",ss3.str(), true);
  drana_QCDPt30To50EmEnriched.analyzersinglefile(cnt);

  return -1;
}
