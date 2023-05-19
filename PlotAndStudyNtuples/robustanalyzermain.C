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
    ss1<<"hists_jpsipt28_"<<cnt<<".root";
    robustanalyzer rana_data("./data/JPsiPt28_EPt15To100_13p6TeV_ntuple_1887590.root", ss1.str(), numCores, true/*isMC*/);
    rana_data.analyzersinglefile(cnt);
  }
  catch (char const* exc) {
    cerr<<exc<<endl;
  }
  
  return -1;
}
