#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"

TString cutdeets = "Cut details";
TFile* sigDoubleElectronPt1To300file = TFile::Open("hists_DoubleElectronGunPt1To300.root","READ");
TFile* sigDoubleElectronPt1To300l1v1file = TFile::Open("hists_DoubleElectronGunPt1To300_1240V74L1v1.root","READ");
TFile* sigDoubleElectronPt1To300l1v2file = TFile::Open("hists_DoubleElectronGunPt1To300_1240V74L1v2.root","READ");
TFile* sigDYToLLM4To50file = TFile::Open("hists_DYToLLM4To50.root","READ");
TFile* sigDYToLLM4To50l1v1file = TFile::Open("hists_DYToLLM4To50_1240V74L1v1.root","READ");
TFile* sigDYToLLM4To50l1v2file = TFile::Open("hists_DYToLLM4To50_1240V74L1v2.root","READ");
TFile* datZeroBiasfile = TFile::Open("hists_ZeroBias2018D.root","READ");
TFile* datEphemeralfile = TFile::Open("hists_Ephemeral1HLTPhysics2018D.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};
std::vector<double> scalehist{1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

int makeratehist(std::vector<TFile*> file, std::vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, float yrange[]=(float []){0.1,100}, float legPos[]=(float []){0.7,0.75,0.95,1}, float seltextpos[]=(float []){0.1,1}, TString xaxistitle="xaxis", float drawsignalline=-1.0) {

  // Check if the size of file vector is same as the cutnames
  if(file.size()!=cutname.size()) {
    cout<<"Error! Mismatching size of vectors"<<endl;
    return -1;
  }

  TString foldername = "";
  foldername += ((TString)file[0]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","")+"_"+cutname[0]+"_";

  std::vector<TH1F*> allhists;
  std::vector<TH1F*> allratehists;
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
  xbinlow = xbinlow==-1?1:xbinlow;
  xbinhigh = xbinhigh==-1?allhists[0]->GetNbinsX():xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    allhists[histctr]->Scale(1.0/scalehist[histctr]);
    if(rebin!=1) allhists[histctr]->Rebin(rebin);

    // Make the rate hist
    allratehists.push_back((TH1F*) allhists[histctr]->Clone());
    for(unsigned int bincnt=0; bincnt<allhists[histctr]->GetNbinsX(); bincnt++) {
      double err = 0;
      allratehists[histctr]->SetBinContent(bincnt, allhists[histctr]->IntegralAndError(bincnt,allhists[histctr]->GetNbinsX()+1,err));
      allratehists[histctr]->SetBinError(bincnt, err);
    }
    
    allratehists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allratehists[histctr]->GetXaxis()->SetTitle(xtitle);
    //allratehists[histctr]->GetYaxis()->SetTitle("L1 rate / prev. L1 rate");
    allratehists[histctr]->GetYaxis()->SetTitle("HLT rate / prev. HLT rate");

    allratehists[histctr]->SetTitle("");

    allratehists[histctr]->SetLineWidth(2);
    allratehists[histctr]->SetLineColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allratehists, legendEntries, allratehists[0]->GetXaxis()->GetTitle(),allratehists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,false);
  c1->SaveAs("./dirplots/"+foldername+"/"+var+"_ratehist.png");

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
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalize);
  c1->SaveAs("./dirplots/"+foldername+"/"+var+".png");

  return -1;
}

int efficiency(std::vector<TFile*> file, std::vector<TString> cutnames, float legPos[]=(float []){0.7,0.75,0.95,1}, int nbins=0, double *rebin=0, TString xaxistitle="p_{T} / GeV") {

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
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(), (float []){-1,0.65,0.9,0.95},false,(float []){0.0,1.1},false);
  auto pad = c1->GetPad(3);
  auto leg = new TLegend(legPos[0], legPos[1], legPos[2], legPos[3]);
  for(unsigned int filenum=0; filenum<file.size(); filenum++) {
    pEff[filenum]->Draw("same");
    leg->AddEntry(pEff[filenum],legendEntries[filenum],"l");
    leg->SetTextFont(132);
    leg->SetTextSize(0.065);
    leg->SetBorderSize(0);
  }
  leg->Draw();
  //TLatex latex;
  //latex.DrawLatex((*legpos),(*(legpos+1)),seltext[0]);
  //latex.DrawLatex((*legpos),(*(legpos+1))-0.1,seltext[1]);
  c1->SaveAs("./dirplots/"+((TString)file[0]->GetName()).ReplaceAll(".root","")+"/"+cutnames[0]+"_eff.png");
  
  return -1;
}

int plotter() {

  std::vector<TFile*> file;
  std::vector<TString> cutname;
  std::vector<int> color;  
  std::vector<TString> legend;  

  /////////////////////////////////////////////////////////////////////
  //////////////////////// DoubleElectron /////////////////////////////
  /////////////////////////////////////////////////////////////////////
  
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
  //comparesamevariable(file, cutname, "diel_deta", -1, -1, 1, true, true, true, (float []){8e-1,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen #Delta#eta(e_{1}, e_{2})");
  //comparesamevariable(file, cutname, "diel_dphi", -1, -1, 1, true, true, true, (float []){8e-1,2e5}, (float []){0.4,0.65,0.65,0.95}, false, "gen #Delta#phi(e_{1}, e_{2})");
  //comparesamevariable(file, cutname, "diel_dR", -1, -1, 1, true, true, true, (float []){8e-1,2e5}, (float []){0.4,0.65,0.65,0.95}, false, "gen #Delta#R(e_{1}, e_{2})");
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
  //comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

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
  //comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgenAnoselsctmchgenel");
  color.push_back(kBlue);
  legend.push_back("genmch sct e");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen");
  color.push_back(kBlue-7);
  legend.push_back("gen e");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elpt", 10, 320, 5, false, true, true, (float []){5e-1,2e3}, (float []){0.2,0.2,0.45,0.45}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 5, false, true, true, (float []){5e-1,1200}, (float []){0.4,0.1,0.65,0.35}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, false, true, true, (float []){5e-1,2e3}, (float []){0.2,0.2,0.45,0.45}, false, "gen electron #phi");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();  
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_sublead_elpt");
  cutname.push_back("noselgenAnoselsctmchgenel_sublead_elpt");
  color.push_back(kBlue);
  legend.push_back("new L1");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_sublead_elpt");
  cutname.push_back("noselgenAoldL1selsctmchgenel_sublead_elpt");
  color.push_back(kRed);
  legend.push_back("old L1");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAEG_25_12_L1selsctmchgenel_elpt");
  color.push_back(kYellow+2);
  legend.push_back("old+L1_EG_XX_YY_er2p5");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAEGLIso_22_12_L1selsctmchgenel_elpt");
  color.push_back(kCyan+2);
  legend.push_back("old+L1_EGLooseIso_XX_YY_er2p5");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAEG_18_17_8_L1selsctmchgenel_elpt");
  color.push_back(kGreen+2);
  legend.push_back("old+L1_EG_XX_YY_ZZ_er2p5");
  /*
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAEG_25_14_L1selsctmchgenel_elpt");
  color.push_back(kGreen+2);
  legend.push_back("old+L1_DoubleEG_25_14_er2p5");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAunpreEGL1selsctmchgenel_elpt");
  color.push_back(kGreen+2);
  legend.push_back("un-prescaled L1 only");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenApreEGL1selsctmchgenel_elpt");
  color.push_back(kYellow+2);
  legend.push_back("prescaled L1 only");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAzeroEGL1selsctmchgenel_elpt");
  color.push_back(kCyan+2);
  legend.push_back("zero L1 only");
  */
  coloropt = color;
  legendEntries = legend;
  vector<double> binspt{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60};
  vector<float> legpos = {0.11, 0.75, 0.36, 1};
  //efficiency(file, cutname, &legpos[0], 60, &binspt[0], "e_{2} p_{T} / GeV");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();  
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAnoselsctmchgenel_elpt");
  color.push_back(kBlack);
  legend.push_back("new L1");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAoldL1selsctmchgenel_elpt");
  color.push_back(kRed);
  legend.push_back("old L1");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenAoldscoutselsctmchgenel_elpt");
  color.push_back(kBlue);
  legend.push_back("old scouting sel");
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselgengen_elpt");
  cutname.push_back("noselgenArandselsctmchgenel_elpt");
  color.push_back(kGreen+2);
  legend.push_back("new L1 + e veto ID");
  coloropt = color;
  legendEntries = legend;
  //vector<double> binspt{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60};
  //vector<double> legpos = {40, 0.5};
  //efficiency(file, cutname, &legpos[0], 60, &binspt[0], "p_{T} / GeV");

  /////////////////////////////////////////////////////////////////////
  //////////////////////// DYToLLM4To50 ///////////////////////////////
  /////////////////////////////////////////////////////////////////////
  
  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgengen");
  color.push_back(kBlue);
  legend.push_back("gen electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 6, 15, 1, true, true, true, (float []){0.7,9e5}, (float []){0.7,0.65,0.95,0.95}, false, "gen electron multiplicity");
  //comparesamevariable(file, cutname, "diel_M", 0, 70, 1, true, true, true, (float []){8e-1,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen #DeltaM(e_{1}, e_{2}) / GeV");
  //comparesamevariable(file, cutname, "diel_deta", -1, -1, 10, true, true, true, (float []){8e-1,2e5}, (float []){0.65,0.65,0.9,0.95}, false, "gen #Delta#eta(e_{1}, e_{2}) / GeV");
  //comparesamevariable(file, cutname, "diel_dphi", -1, -1, 10, true, true, true, (float []){8e-1,2e5}, (float []){0.4,0.65,0.65,0.95}, false, "gen #Delta#phi(e_{1}, e_{2}) / GeV");
  //comparesamevariable(file, cutname, "diel_dR", -1, -1, 10, true, true, true, (float []){8e-1,2e5}, (float []){0.4,0.65,0.65,0.95}, false, "gen #Delta#R(e_{1}, e_{2}) / GeV");
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgengen_lead");
  color.push_back(kRed);
  legend.push_back("lead gen");
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgengen_sublead");
  color.push_back(kGreen);
  legend.push_back("sub-lead gen");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elpt", 0, 100, 1, true, true, true, (float []){5e-1,2e4}, (float []){0.6,0.6,0.85,0.85}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 10, true, true, true, (float []){0.7,1e5}, (float []){0.4,0.2,0.65,0.45}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){0.7,1e5}, (float []){0.4,0.2,0.65,0.45}, false, "gen electron #phi");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM4To50file);
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
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgenAnoselgenelsctbar");
  color.push_back(kRed+2);
  legend.push_back("sct electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dEta", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){1e2,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.6,0.85,0.85}, false, "charge#times#Delta#phi (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){1e2,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgenAnoselgenelsctee");
  color.push_back(kRed+2);
  legend.push_back("sct electrons");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dEta", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){5e1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.6,0.85,0.85}, false, "charge#times#Delta#phi (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){5e1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgenAnoselgenelsctmchbar");
  color.push_back(kRed+2);
  legend.push_back("genmch sct e");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgenAnoselgenelsctmchee");
  color.push_back(kRed+2);
  legend.push_back("genmch sct e");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dEta", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta (gen, sct)");
  //comparesamevariable(file, cutname, "qdPhi", 4000, 6000, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "charge#times#Delta#phi (gen, sct)");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgenAnoselsctmchgenel");
  color.push_back(kBlue);
  legend.push_back("genmch sct e");
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgengen");
  color.push_back(kBlue-7);
  legend.push_back("gen e");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elpt", 10, 320, 5, false, true, true, (float []){5e-1,2e3}, (float []){0.2,0.2,0.45,0.45}, false, "gen electron p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 5, false, true, true, (float []){5e-1,1200}, (float []){0.4,0.1,0.65,0.35}, false, "gen electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, false, true, true, (float []){5e-1,2e3}, (float []){0.2,0.2,0.45,0.45}, false, "gen electron #phi");

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgengen_sublead_elpt");
  cutname.push_back("noselgenAnoselsctmchgenel_sublead_elpt");
  color.push_back(kBlue);
  legend.push_back("new L1");
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgengen_sublead_elpt");
  cutname.push_back("noselgenAoldL1selsctmchgenel_sublead_elpt");
  color.push_back(kRed);
  legend.push_back("old L1");
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselgengen_sublead_elpt");
  cutname.push_back("noselgenApreEGL1selsctmchgenel_sublead_elpt");
  color.push_back(kGreen+2);
  legend.push_back("old L1 + L1_EG8_HTT300_er2p5");
  /*
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("basicselgengen_elpt");
  cutname.push_back("basicselgenAnoselsctmchgenel_elpt");
  color.push_back(kBlue);
  legend.push_back("new L1");
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("basicselgengen_elpt");
  cutname.push_back("basicselgenAoldL1selsctmchgenel_elpt");
  color.push_back(kRed);
  legend.push_back("old L1");
  */
  coloropt = color;
  legendEntries = legend;
  legpos = {0.4, 0.2, 0.65, 0.45};
  //efficiency(file, cutname, &legpos[0], 60, &binspt[0], "e_{2} p_{T} / GeV");

  /////////////////////////////////////////////////////////////////////
  //////////////////////////// ZeroBias ///////////////////////////////
  /////////////////////////////////////////////////////////////////////

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(datZeroBiasfile);
  cutname.push_back("noselsctbar");
  color.push_back(kBlack);
  legend.push_back("new L1");
  file.push_back(datZeroBiasfile);
  cutname.push_back("oldL1selsctbar");
  color.push_back(kRed);
  legend.push_back("old L1");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 5, 20, 1, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "e multiplicity");
  //comparesamevariable(file, cutname, "elminpt", 5, 150, 1, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "min. e p_{T} / GeV");
  //comparesamevariable(file, cutname, "elpt", 5, 150, 1, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "e p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 5, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "e #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.74,0.85,0.99}, false, "e #phi");  
  //comparesamevariable(file, cutname, "elhoe", -1, 200, 4, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.74,0.85,0.99}, false, "e #phi");  
  
  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  scalehist.clear();
  file.push_back(datZeroBiasfile);
  cutname.push_back("noselsct");
  color.push_back(kBlack);
  legend.push_back("new L1");
  scalehist.push_back(1992.985);
  file.push_back(datZeroBiasfile);
  cutname.push_back("oldL1selsct");
  color.push_back(kRed);
  legend.push_back("old L1");
  scalehist.push_back(1992.985);
  file.push_back(datZeroBiasfile);
  cutname.push_back("EG_25_12_L1selsct");
  color.push_back(kYellow+2);
  legend.push_back("old+L1_EG_XX_YY_er2p5");
  scalehist.push_back(1992.985);
  file.push_back(datZeroBiasfile);
  cutname.push_back("EGLIso_22_12_L1selsct");
  color.push_back(kCyan+2);
  legend.push_back("old+L1_EGLooseIso_XX_YY_er2p5");
  scalehist.push_back(1992.985);
  file.push_back(datZeroBiasfile);
  cutname.push_back("EG_18_17_8_L1selsct");
  color.push_back(kGreen+2);
  legend.push_back("old+L1_EG_XX_YY_ZZ_er2p5");
  scalehist.push_back(1992.985);
  coloropt = color;
  legendEntries = legend;
  //makeratehist(file, cutname, "elminpt", 5, 150, 1, false, (float []){0,5}, (float []){0.4, 0.6, 0.65, 0.85});
  //makeratehist(file, cutname, "elmaxpt", 5, 150, 1, false, (float []){0,5}, (float []){0.4, 0.6, 0.65, 0.85});

  /////////////////////////////////////////////////////////////////////
  /////////////////////// EphemeralHLTPhysics /////////////////////////
  /////////////////////////////////////////////////////////////////////

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(datEphemeralfile);
  cutname.push_back("noselsct");
  color.push_back(kBlack);
  legend.push_back("new L1");
  file.push_back(datEphemeralfile);
  cutname.push_back("oldL1selsct");
  color.push_back(kRed);
  legend.push_back("old L1");
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 5, 20, 1, true, true, true, (float []){5e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "e multiplicity");
  //comparesamevariable(file, cutname, "elminpt", 5, 150, 1, true, true, true, (float []){5e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "min. e p_{T} / GeV");
  //comparesamevariable(file, cutname, "elpt", 5, 150, 1, true, true, true, (float []){5e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "e p_{T} / GeV");
  //comparesamevariable(file, cutname, "eleta", -1, -1, 5, true, true, true, (float []){5e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "e #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){5e-1,1e7}, (float []){0.6,0.74,0.85,0.99}, false, "e #phi");  
  //comparesamevariable(file, cutname, "dielM", -1, 100, 1, true, true, true, (float []){5e-1,1e7}, (float []){0.6,0.74,0.85,0.99}, false, "invariant mass / GeV");  
  //comparesamevariable(file, cutname, "elhoe", -1, -1, 4, true, true, true, (float []){5e-1,1e5}, (float []){0.6,0.74,0.85,0.99}, false, "e #phi");  
  
  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  scalehist.clear();
  file.push_back(datEphemeralfile);
  cutname.push_back("noselsct");
  color.push_back(kBlack);
  legend.push_back("new L1");
  scalehist.push_back(95855.2);
  file.push_back(datEphemeralfile);
  cutname.push_back("oldL1selsct");
  color.push_back(kRed);
  legend.push_back("old L1");
  scalehist.push_back(95855.2);
  file.push_back(datEphemeralfile);
  cutname.push_back("oldscoutselsct");
  color.push_back(kBlue);
  legend.push_back("old L1 + old sel");
  scalehist.push_back(95855.2);
  file.push_back(datEphemeralfile);
  cutname.push_back("oldscoutvetoselsct");
  color.push_back(kCyan+2);
  legend.push_back("old L1 + e veto ID");
  scalehist.push_back(95855.2);
  file.push_back(datEphemeralfile);
  cutname.push_back("randselsct");
  color.push_back(kGreen+2);
  legend.push_back("new L1 + e veto ID");
  scalehist.push_back(95855.2);
  coloropt = color;
  legendEntries = legend;
  //makeratehist(file, cutname, "elminpt", 5, 150, 1, false, (float []){0,3}, (float []){0.6,0.6,0.85,0.95});
  //makeratehist(file, cutname, "elmaxpt", 5, 150, 1, false, (float []){0,3}, (float []){0.6,0.6,0.85,0.95});

  // For TSG approval
  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDYToLLM4To50file);
  cutname.push_back("noselsct");
  color.push_back(kBlue);
  legend.push_back("old L1");
  scalehist.push_back(1);
  file.push_back(sigDYToLLM4To50l1v1file);
  cutname.push_back("noselsct");
  color.push_back(kGreen+2);
  legend.push_back("new L1v1");
  scalehist.push_back(1);
  file.push_back(sigDYToLLM4To50l1v2file);
  cutname.push_back("noselsct");
  color.push_back(kRed);
  legend.push_back("new L1v2");
  scalehist.push_back(1);
  coloropt = color;
  legendEntries = legend;
  comparesamevariable(file, cutname, "elmult", 5, 10, 1, false, true, true, (float []){0,100}, (float []){0.6,0.7,0.85,0.95}, false, "e multiplicity");  
  comparesamevariable(file, cutname, "elpt", 5, 75, 1, false, true, true, (float []){0,40}, (float []){0.6,0.7,0.85,0.95}, false, "e p_{T} [GeV]");  
  comparesamevariable(file, cutname, "eleta", -1, -1, 20, false, true, true, (float []){0,30}, (float []){0.6,0.7,0.85,0.95}, false, "e #eta");  
  comparesamevariable(file, cutname, "elphi", -1, -1, 3, false, true, true, (float []){0,20}, (float []){0.6,0.7,0.85,0.95}, false, "e #phi");  

  file.clear();
  cutname.clear();
  color.clear();
  legend.clear();
  file.push_back(sigDoubleElectronPt1To300file);
  cutname.push_back("noselsct");
  color.push_back(kBlue);
  legend.push_back("1240_V74 L1");
  scalehist.push_back(1);
  file.push_back(sigDoubleElectronPt1To300l1v1file);
  cutname.push_back("noselsct");
  color.push_back(kGreen+2);
  legend.push_back("add L1 v1");
  scalehist.push_back(1);
  file.push_back(sigDoubleElectronPt1To300l1v2file);
  cutname.push_back("noselsct");
  color.push_back(kRed);
  legend.push_back("add L1 v2");
  scalehist.push_back(1);
  coloropt = color;
  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 5, 15, 1, true, true, true, (float []){1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "e multiplicity");  
  //comparesamevariable(file, cutname, "elpt", 5, 575, 10, false, true, false, (float []){1e-5,8e3}, (float []){0.6,0.7,0.85,0.95}, false, "e p_{T} [GeV]");  
  //comparesamevariable(file, cutname, "eleta", 200, 800, 10, false, true, true, (float []){1e-3,5e3}, (float []){0.6,0.7,0.85,0.95}, false, "e #eta");  
  //comparesamevariable(file, cutname, "elphi", -1, -1, 3, false, true, true, (float []){1e-3,1.5e4}, (float []){0.6,0.7,0.85,0.95}, false, "e #phi");  
  //comparesamevariable(file, cutname, "elminpt", 5, 575, 10, false, true, false, (float []){1e-5,4e3}, (float []){0.6,0.7,0.85,0.95}, false, "e min. p_{T} [GeV]");  
  //comparesamevariable(file, cutname, "elmaxpt", 5, 575, 10, false, true, false, (float []){1e-5,4e3}, (float []){0.6,0.7,0.85,0.95}, false, "e max. p_{T} [GeV]");  
  //comparesamevariable(file, cutname, "elleadeta", 200, 800, 10, false, true, true, (float []){1e-3,3e3}, (float []){0.6,0.7,0.85,0.95}, false, "leading e #eta");  

  return -1;
}
