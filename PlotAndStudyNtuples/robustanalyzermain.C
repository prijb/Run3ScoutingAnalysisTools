#include "robustanalyzer.hh"
#include <iostream>
#include <string>
#include <sstream>

using namespace std;

int main(int argc, char* argv[]) {
  
  stringstream ss;
  ss << argv[1];
  int cnt;
  ss >> cnt;
  
  stringstream ssargv2;
  ssargv2 << argv[2];
  int numCores;
  ssargv2 >> numCores;

  try {
    /*
    stringstream ssjpsi;
    ssjpsi<<"hists_Eph02022G_"<<cnt<<".root";
    robustanalyzer rana_JPsi("./data/ScoutingNTuple_Eph02022G_230224trial_1509269.root", ssjpsi.str(), numCores, false, false, false);
    rana_JPsi.analyzersinglefile(cnt);
    */
    stringstream ssjpsi;
    ssjpsi<<"hists_JPsi_"<<cnt<<".root";
    robustanalyzer rana_JPsi("./data/ScoutingNTuple_JPsiToEE_1300pre4_Skim230301_1550584.root", ssjpsi.str(), numCores, false, true, true);
    rana_JPsi.analyzersinglefile(cnt);
    
    stringstream ssjpsi2;
    ssjpsi2<<"hists_JPsi_mod_"<<cnt<<".root";
    robustanalyzer rana_JPsi2("./data/ScoutingNTuple_JPsiToEE_1300pre4_Skim230301modScout_1550592.root", ssjpsi2.str(), numCores, false, true, true);
    rana_JPsi2.analyzersinglefile(cnt);
    
}
  
  catch (char const* exc) {
    cerr<<exc<<endl;
  }
  
  return -1;
}
