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
//TFile* sigDYToLLM50file = TFile::Open("hists_DYToLLM50.root","READ");
TFile* sigDoubleElectronPt1To300file = TFile::Open("hists_DoubleElectronGunPt1To300.root","READ");
//TFile* sigDoubleElectronPt1To300Oldfile = TFile::Open("hists_DoubleElectronGunPt1To300Old.root","READ");
//TFile* sigDoubleElectronPt1To300220508Oldfile = TFile::Open("hists_DoubleElectronGunPt1To300220508Old.root","READ");
//TFile* sigQCDPt20To30file = TFile::Open("hists_QCDPt20To30EmEnriched.root","READ");
//TFile* sigQCDPt30To50file = TFile::Open("hists_QCDPt30To50EmEnriched.root","READ");
//TFile* sigQCDPt50To80file = TFile::Open("hists_QCD_Pt50To80_EMEnriched.root","READ");
//TFile* sigQCDPt15To20bcToEfile = TFile::Open("hists_QCD_Pt15To20bcToE_EMEnriched.root","READ");
//TFile* sigQCDPt20To30bcToEfile = TFile::Open("hists_QCD_Pt20To30bcToE_EMEnriched.root","READ");
//TFile* sigQCDPt30To80bcToEfile = TFile::Open("hists_QCD_Pt30To80bcToE_EMEnriched.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};

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
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalize);
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
  std::vector<int> color;  
  std::vector<TString> legend;  

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen");
  color.push_back(kBlue);
  legend.push_back("gen electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 6, 11, 1, true, true, true, (float []){0.7,9e5}, (float []){0.7,0.65,0.95,0.95}, false, "gen electron multiplicity");
  //comparesamevariable(file, cutname, "diel_M", 50, 150, 2, true, true, true, (float []){50,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen #DeltaM(e_{1}, e_{2}) / GeV");
  //comparesamevariable(file, cutname, "diel_deta", -1, -1, 1, true, true, true, (float []){8e-1,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen #Delta#eta(e_{1}, e_{2}) / GeV");
  //comparesamevariable(file, cutname, "diel_dphi", -1, -1, 1, true, true, true, (float []){8e-1,2e5}, (float []){0.4,0.65,0.65,0.95}, false, "gen #Delta#phi(e_{1}, e_{2}) / GeV");
  //comparesamevariable(file, cutname, "diel_dR", -1, -1, 1, true, true, true, (float []){8e-1,2e5}, (float []){0.4,0.65,0.65,0.95}, false, "gen #Delta#R(e_{1}, e_{2}) / GeV");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_lead");
  color.push_back(kRed);
  legend.push_back("lead gen");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_sublead");
  color.push_back(kGreen);
  legend.push_back("sub-lead gen");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elpt", 10, 320, 1, true, true, true, (float []){5e-1,2e4}, (float []){0.2,0.2,0.45,0.45}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 10, true, true, true, (float []){0.7,1e5}, (float []){0.4,0.2,0.65,0.45}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.4,0.2,0.65,0.45}, false, "gen electron #phi");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselsct");
  color.push_back(kRed+2);
  legend.push_back("sct electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", -1, -1, 1, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "e multiplicity");
  //comparesamevariable(file, cutname, "elpt", -1, 350, 1, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "e p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 1, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "e #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.74,0.85,0.99}, false, "e #phi");
  
  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgenAnoselgenelsctbar");
  color.push_back(kRed+2);
  legend.push_back("sct electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dEta", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){1e2,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.4,0.1,0.65,0.35}, false, "charge#times#Delta#phi (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){1e2,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgenAnoselgenelsctee");
  color.push_back(kRed+2);
  legend.push_back("sct electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dEta", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){5e1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.4,0.1,0.65,0.35}, false, "charge#times#Delta#phi (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){5e1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgenAnoselgenelsctmchbar");
  color.push_back(kRed+2);
  legend.push_back("genmch sct e");
  coloropt = color;
  legendEntries = legend;
  comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgenAnoselgenelsctmchee");
  color.push_back(kRed+2);
  legend.push_back("genmch sct e");
  coloropt = color;
  legendEntries = legend;
  comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  return -1;
}
