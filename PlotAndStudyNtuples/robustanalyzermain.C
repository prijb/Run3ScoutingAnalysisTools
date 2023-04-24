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
    
    stringstream ss1;
    ss1<<"hists_data_"<<cnt<<".root";
    robustanalyzer rana_data("./data/CrabSkim230414_ScoutingNTuple_ForScoutingEleEff_fatdummy.root", ss1.str(), numCores);
    rana_data.analyzersinglefile(cnt);
  }
  catch (char const* exc) {
    cerr<<exc<<endl;
  }
  
  return -1;
}
