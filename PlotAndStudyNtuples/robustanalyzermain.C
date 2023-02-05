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
    robustanalyzer rana_data("./data/scoutingNTuple.root", ss1.str(), numCores, false, false);
    rana_data.analyzersinglefile(cnt);
    /*
    stringstream ss2;
    ss2<<"hists_DYToLLM50_"<<cnt<<".root";
    robustanalyzer rana_DYM50("./data/DYToLLM50_ntuple_117598.root", ss2.str(), numCores, true, true);
    rana_DYM50.analyzersinglefile(cnt);
    
    stringstream ss3;
    ss3<<"hists_QCDPt30To50_"<<cnt<<".root";
    robustanalyzer rana_QCDPt30To50("./data/QCDPt30To50_ntuple_119416.root", ss3.str(), numCores, false, true);
    rana_QCDPt30To50.analyzersinglefile(cnt);
    stringstream ss4;
    ss4<<"hists_QCDPt50To80_"<<cnt<<".root";
    robustanalyzer rana_QCDPt50To80("./data/QCDPt50To80_ntuple_122416.root", ss4.str(), numCores, false, true);
    rana_QCDPt50To80.analyzersinglefile(cnt);
    stringstream ss5;
    ss5<<"hists_QCDPt80To120_"<<cnt<<".root";
    robustanalyzer rana_QCDPt80To120("./data/QCDPt80To120_ntuple_122417.root", ss5.str(), numCores, false, true);
    rana_QCDPt80To120.analyzersinglefile(cnt);
    */
  }
  catch (char const* exc) {
    cerr<<exc<<endl;
  }
  
  return -1;
}
