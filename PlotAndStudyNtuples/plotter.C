#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

double datasf = 1.09/19; // Data rate - 1.1 Hz, 22 events from parent trigger selection
double sig3cmsf = 1.0/15;
double sig30cmsf = 1.0/28;
double sig1msf = 1.0/19;
double sig3msf = 1.0/7;
TString cutdeets = "Cut details";
TFile* datahistfile = TFile::Open("hists_data.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};
std::vector<TString> histtype{"p e1", "hist same"};
std::vector<int> markerstyle{20, 24};
std::vector<int> markersize{10, 10};
std::vector<TString> legendmarkerstyle{"lep", "l"};
//std::vector<double> scale{1, 1};

/*
int autoplotter(std::vector<TFile*> file, std::vector<TString> cutname, ) {
}
*/

int invmee_specialplot(TString selection) {

  auto invmeehist = (TH1F*) datahistfile->Get(selection);

  vector<double> xbins;
  xbins.push_back(0.5);
  while(xbins[xbins.size()-1]<100) {
    double val = ceil(110.0*xbins[xbins.size()-1])/100.0;
    xbins.push_back(val);
  }

  invmeehist = (TH1F*) invmeehist->Rebin(xbins.size()-1, selection, &xbins[0]);
  
  TCanvas* c1;
  c1 = new TCanvas();
  gStyle->SetOptStat(0);
  c1->SetLogx(true);
  invmeehist->Draw();
  c1->SaveAs(selection+".png");
  
  return -1;
}

int comparesamevariable(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool underflow=false, bool overflow=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, bool normalize=true, TString xaxistitle="xaxis") {

  // Check if the size of file vector is same as the cutnames
  if(file.size()!=cutname.size()) {
    cout<<"Error! Mismatching size of vectors"<<endl;
    return -1;
  }

  TString foldername = "";
  foldername += ((TString)file[0]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","")+"_"+cutname[0]+"_";

  std::vector<TH1F*> allhists;
  for(unsigned int histctr=0; histctr<cutname.size(); histctr++) {
    TString cutwithvar = cutname[histctr]+"_"+var;
    allhists.push_back((TH1F*)file[histctr]->Get(cutwithvar));
    cutname[histctr] += "_"+((TString)file[histctr]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","");
  }
  
  // Get the title from histogram title
  TString dummytitle = "xaxis";
  TString histtitle = allhists[0]->GetTitle();
  TString xtitle = xaxistitle.CompareTo(dummytitle)?xaxistitle:histtitle;

  if(rebin==-1) rebin = 1;
  xbinlow = xbinlow==-1?(underflow?0:1):xbinlow;
  xbinhigh = xbinhigh==-1?(overflow?allhists[0]->GetNbinsX()+1:allhists[0]->GetNbinsX()):xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(rebin!=1) allhists[histctr]->Rebin(rebin);

    // Make changes to sig and bkg to enable good basic plotting
    double err = 0.0;
    if(underflow) {
      allhists[histctr]->SetBinContent(xbinlow,allhists[histctr]->IntegralAndError(0,xbinlow,err));
      allhists[histctr]->SetBinError(xbinlow,err);
    }
    err = 0.0;
    if(overflow) {
      allhists[histctr]->SetBinContent(xbinhigh,allhists[histctr]->IntegralAndError(xbinhigh,allhists[histctr]->GetNbinsX()+1,err));
      allhists[histctr]->SetBinError(xbinhigh,err);
    }

    allhists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allhists[histctr]->GetXaxis()->SetTitle(xtitle);
    if(normalize) allhists[histctr]->GetYaxis()->SetTitle("normalized number of events (a.u.)");
    else allhists[histctr]->GetYaxis()->SetTitle("number of events");

    allhists[histctr]->SetTitle("");

    allhists[histctr]->SetLineWidth(2);
    allhists[histctr]->SetLineColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalize,histtype,markerstyle,markersize,legendmarkerstyle);
  c1->SaveAs("./dirplots/"+foldername+"/"+var+".png");

  return -1;
}

int efficiency(std::vector<TFile*> file, std::vector<TString> cutnames, int nbins=0, double *rebin=0, TString xaxistitle="p_{T} / GeV") {

  if(cutnames.size()/2!=file.size()) {
    cout<<"Inconsistent entries for file size and cutnames: Suggested cutnameSize = 2* fileSize"<<endl;
    return -1;
  }
  
  std::vector<TEfficiency*> pEff;
  TH1F* demo;

  for(unsigned int filenum=0; filenum<file.size(); filenum++) {
    
    TH1F* nosel = (TH1F*) file[filenum]->Get(cutnames[filenum*2]);
    TH1F* sel = (TH1F*) file[filenum]->Get(cutnames[filenum*2+1]);
    
    nosel = (TH1F*) nosel->Rebin(nbins,"newx",rebin);
    sel = (TH1F*) sel->Rebin(nbins,"newx",rebin);
    if(filenum==0) demo = (TH1F*) sel->Clone();
    
    pEff.push_back(0);
    if(TEfficiency::CheckConsistency((*sel),(*nosel))) {
      pEff[filenum] = new TEfficiency((*sel),(*nosel));
    }
    
    sel->SetLineColor(kRed);
    gStyle->SetOptStat(0);
    pEff[filenum]->SetLineWidth(3);
    pEff[filenum]->SetLineColor(coloropt[filenum]);
  }
  
  std::vector<TH1F*> allhists;
  demo->SetTitle("");
  demo->GetXaxis()->SetTitle(xaxistitle);
  demo->GetYaxis()->SetTitle("RECO eff.");
  demo->SetLineColorAlpha(kWhite,1);
  demo->SetFillColorAlpha(kWhite,1);
  allhists.push_back(demo);
  
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(), (float []){0.7,0.65,0.9,0.95},false,(float []){0.0,1.1},false);
  auto pad = c1->GetPad(3);
  for(unsigned int filenum=0; filenum<file.size(); filenum++) {
    pEff[filenum]->Draw("same");
  }
  c1->SaveAs("./dirplots/"+((TString)file[0]->GetName()).ReplaceAll(".root","")+"/"+cutnames[0]+"_eff.png");
  
  return -1;
}

int plotter() {

  std::vector<TFile*> file;
  std::vector<TString> cutname;
  std::vector<TString> legend;  

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("noselsct");
  coloropt.push_back(kBlack);
  legend.push_back("2021 Scouting Data");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 5, 15, 1, true, true, true, (float []){8e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", 5, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "eleta", 200, 800, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #phi");
  //comparesamevariable(file, cutname, "dielM", -1, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("noselsct");
  coloropt.push_back(kBlack);
  legend.push_back("2021 Scouting Data");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "dielM", -1, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("vetoselsct");
  coloropt.push_back(kBlack);
  legend.push_back("2021 Scouting Data");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 5, 15, 1, true, true, true, (float []){8e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", 5, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "eleta", 200, 800, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #phi");
  //comparesamevariable(file, cutname, "dielM", -1, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("vetoselsct");
  coloropt.push_back(kBlack);
  legend.push_back("2021 Scouting Data");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "dielM", -1, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  invmee_specialplot("noselsct_dielM");
  invmee_specialplot("vetoselsct_dielM");
  invmee_specialplot("looseselsct_dielM");
  invmee_specialplot("mediumselsct_dielM");
  invmee_specialplot("tightselsct_dielM");

  return -1;
}
