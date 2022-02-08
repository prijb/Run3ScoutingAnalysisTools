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

  return -1;
}
