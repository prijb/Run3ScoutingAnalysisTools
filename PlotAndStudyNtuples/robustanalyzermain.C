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
  
  try {    
    stringstream ss1;
    ss1<<"hists_data_"<<cnt<<".root";
    robustanalyzer rana_data("./data/Run3Scouting_singlefile_Skim220717.root",ss1.str(), false, false);
    rana_data.analyzersinglefile(cnt);
  }
  catch (char const* exc) {
    cerr<<exc<<endl;
  }
  
  return -1;
}
