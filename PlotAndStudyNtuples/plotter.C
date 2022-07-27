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
//TFile* datahistfile = TFile::Open("hists_data.root","READ");
//TFile* datahistfile = TFile::Open("hists_data_12504363.root","READ");
//TFile* datahistfile = TFile::Open("hists_data_12517349.root","READ");
TFile* datahistfile = TFile::Open("hists_data_12521272.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};
std::vector<TString> histtype{"p e1", "hist same"};
std::vector<int> markerstyle{20, 24};
std::vector<int> markersize{10, 10};
std::vector<TString> legendmarkerstyle{"lep", "l"};
//std::vector<double> scale{1, 1};

int autoplotter(std::vector<TFile*> file, std::vector<TString> cutname) {

  TString foldername = "";
  foldername += ((TString)file[0]->GetName()).ReplaceAll(".root","").ReplaceAll("hists_","")+"_"+cutname[0]+"_autoplotted";

  TList* listOfKeys = file[0]->GetListOfKeys(); // Get the list of histograms
  for(unsigned int ctr=0; ctr<listOfKeys->GetEntries(); ctr++) {
    TString keyName = listOfKeys->At(ctr)->GetName(); // Get the histogram name
    bool isWithCut = keyName.BeginsWith(cutname[0]); // Check if the histogram begins with a certain cut

    if(!isWithCut) continue;

    int maxBins = 10000;
    TH1F* histVar = (TH1F*) file[0]->Get(keyName);
    histVar->Sumw2();
    //if(histVar->GetNbinsX()>maxBins) histVar->Rebin(maxBins);
    /*
    // Go over the histogram and rebin it based on the contents
    std::vector<double> binlowedge;
    binlowedge.push_back(histVar->GetBinLowEdge(1));
    double maxvalue = histVar->GetMaximum();
    double binval = 0.0;
    double cutFracWMax = 0.05;
    for(unsigned int binctr=1; binctr<histVar->GetNbinsX(); binctr++) {
      binval += histVar->GetBinContent(binctr);
      if(binval>cutFracWMax*maxvalue) {
	binlowedge.push_back(histVar->GetBinLowEdge(binctr+1));
	binval = 0.0;
      }
    }
    if(binlowedge[binlowedge.size()-1]!=histVar->GetBinLowEdge(histVar->GetNbinsX()+1)) {
      binlowedge.push_back(histVar->GetBinLowEdge(histVar->GetNbinsX()+1));
    }

    histVar = histVar->Rebin(binlowedge.size()-1,histVar->GetTitle(),binlowedge)
    */
    TCanvas* c1;
    c1 = new TCanvas();
    gStyle->SetOptStat(0);
    //c1->SetLogx(true);
    histVar->Draw();
    c1->SaveAs("./dirplots/"+foldername+"/"+keyName+".png");
    
  }

  return -1;
}


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

int fitinvmee(TString selection) {

  auto invmeehist = (TH1F*) datahistfile->Get(selection);

  invmeehist->SetTitle("");
  invmeehist->Rebin(100);
  invmeehist->GetXaxis()->SetRange(30,150);

  TF1* fitbkg = new TF1("fitbkg","exp([0]+[1]*x)",40,65);
  fitbkg->SetParameters(5,-0.002);
  invmeehist->Fit("fitbkg","R");
  fitbkg->SetLineColor(kBlue);
  double b0 = fitbkg->GetParameter(0);
  double be0 = fitbkg->GetParError(0);
  double b1 = fitbkg->GetParameter(1);
  double be1 = fitbkg->GetParError(1);

  TF1* fitpeak = new TF1("fitpeak","[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",80,95);
  fitpeak->SetParameters(702,91,5);
  invmeehist->Fit("fitpeak","R");
  fitpeak->SetLineColor(kGreen+2);
  double p0 = fitpeak->GetParameter(0);
  double pe0 = fitpeak->GetParError(0);
  double p1 = fitpeak->GetParameter(1);
  double pe1 = fitpeak->GetParError(1);
  double p2 = fitpeak->GetParameter(2);
  double pe2 = fitpeak->GetParError(2);

  TF1* fitfunc = new TF1("fitZmass","exp([0]+[1]*x)+[2]*exp(-0.5*((x-[3])/[4])*((x-[3])/[4]))",40,120);
  fitfunc->SetParameters(b0, b1, p0, p1, p2);
  fitfunc->SetParLimits(0, b0-10*be0, b0+10*be0);
  fitfunc->SetParLimits(1, b1-10*be1, b1+10*be1);
  fitfunc->SetParLimits(2, p0-10*pe0, p0+10*pe0);
  fitfunc->SetParLimits(3, p1-10*pe1, p1+10*pe1);
  fitfunc->SetParLimits(4, p2-10*pe2, p2+10*pe2);
  
  TFitResultPtr fitres = invmeehist->Fit("fitZmass","RS");

  cout<<"Chi2/Ndf = "<<fitres->Chi2()<<"/"<<fitres->Ndf()<<endl;

  TF1* Zpeak = new TF1("Zpeak","[0]*exp(-0.5*((x-[1])/[2])*((x-[1])/[2]))",40,120);
  Zpeak->SetParameters(fitfunc->GetParameter(2), fitfunc->GetParameter(3), fitfunc->GetParameter(4));
  cout<<"Signal integral in the (83.73,92.53): "<<Zpeak->Integral(83.73,92.53)<<endl;
  cout<<"Bkg integral in the (83.73,92.53): "<<fitfunc->Integral(83.73,92.53)-Zpeak->Integral(83.73,92.53)<<endl;
  cout<<"Signal integral in the (80,100): "<<Zpeak->Integral(80,100)<<endl;
  cout<<"Bkg integral in the (80,100): "<<fitfunc->Integral(80,100)-Zpeak->Integral(80,100)<<endl;
  cout<<"Signal integral in the (70,100): "<<Zpeak->Integral(70,100)<<endl;
  cout<<"Bkg integral in the (70,100): "<<fitfunc->Integral(70,100)-Zpeak->Integral(70,100)<<endl;
  cout<<"Signal integral in the (50,110): "<<Zpeak->Integral(50,110)<<endl;
  cout<<"Bkg integral in the (50,110): "<<fitfunc->Integral(50,110)-Zpeak->Integral(50,110)<<endl;
  cout<<"Signal integral in the (45,115): "<<Zpeak->Integral(45,115)<<endl;
  cout<<"Bkg integral in the (45,115): "<<fitfunc->Integral(45,115)-Zpeak->Integral(45,115)<<endl;
  cout<<"Signal integral in the (40,120): "<<Zpeak->Integral(40,120)<<endl;
  cout<<"Bkg integral in the (40,120): "<<fitfunc->Integral(40,120)-Zpeak->Integral(40,120)<<endl;

  fitfunc->SetLineWidth(3);
  
  TCanvas* c1;
  //invmeehist->Draw();
  //fitbkg->Draw("SAME");
  //fitpeak->Draw("SAME");
  c1 = enhance_plotter({invmeehist}, {"Scouting"}, "M(e,e)", "number of events", (float []){0.15,0.7,0.4,0.9}, false, (float []){0,1500}, false, {"hist"}, {20}, {2}, {"lep"});
  fitfunc->Draw("same");
  c1->SaveAs(selection+"_invmeefit.png");

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

void group_plotter(std::vector<TFile*> file, std::vector<TString> cutname, bool isBar) {

  //all1dhists.push_back(new TH1F(selection+"sctbar_elm","m / GeV",1000,-1e-5,1e-5));
  //all1dhists.push_back(new TH1F(selection+"sctbar_ellog10d0","log_{10}d_{0} / log_{10}cm",1000,-5,5));
  //all1dhists.push_back(new TH1F(selection+"sctbar_eldz","d_{z} / cm",4000,-20,20));
  //all1dhists.push_back(new TH1F(selection+"sctbar_elcharge","charge",5,-2,3));
  comparesamevariable(file, cutname, "elpt", 5, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron p_{T} [GeV]");
  comparesamevariable(file, cutname, "eld0", 8500, 11500, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron d_{0} [cm]");
  comparesamevariable(file, cutname, "eldetain", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta(SC seed, track)");
  if(isBar) {
    comparesamevariable(file, cutname, "eldphiin", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#phi(SC, track)");
    comparesamevariable(file, cutname, "elsigmaietaieta", -1, 400, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #sigmai#etai#eta");
    comparesamevariable(file, cutname, "elhoe", -1, -1, 50, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron H/E");
    comparesamevariable(file, cutname, "elooemoop", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron E^{-1}-p^{-1} [GeV^{-1}]");
    comparesamevariable(file, cutname, "elmhits", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron missing hits");
    comparesamevariable(file, cutname, "elecaliso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron ECAL iso. [GeV]");
    comparesamevariable(file, cutname, "elhcaliso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron HCAL iso. [GeV]");
    comparesamevariable(file, cutname, "eltkiso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron track iso. [GeV]");
    comparesamevariable(file, cutname, "elr9", -1, 150, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron r9");
    comparesamevariable(file, cutname, "elsmin", -1, 1100, 5, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron smin");
    comparesamevariable(file, cutname, "elsmaj", -1, 700, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron smaj");
  }
  else {
    comparesamevariable(file, cutname, "eldetain", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta(SC seed, track)");
    comparesamevariable(file, cutname, "eldphiin", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#phi(SC, track)");
    comparesamevariable(file, cutname, "elsigmaietaieta", -1, 800, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #sigmai#etai#eta");
    comparesamevariable(file, cutname, "elhoe", -1, -1, 50, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron H/E");
    comparesamevariable(file, cutname, "elooemoop", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron E^{-1}-p^{-1} [GeV^{-1}]");
    comparesamevariable(file, cutname, "elmhits", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron missing hits");
    comparesamevariable(file, cutname, "elecaliso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron ECAL iso. [GeV]");
    comparesamevariable(file, cutname, "elhcaliso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron HCAL iso. [GeV]");
    comparesamevariable(file, cutname, "eltkiso", -1, -1, 10, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron track iso. [GeV]");
    comparesamevariable(file, cutname, "elr9", -1, 300, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron r9");
    comparesamevariable(file, cutname, "elsmin", -1, 2200, 5, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron smin");
    comparesamevariable(file, cutname, "elsmaj", -1, 1400, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron smaj");
  }
    
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
  legend.push_back("2022 Scouting Data");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "elmult", 5, 15, 1, true, true, true, (float []){8e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "electron multiplicity");
  //comparesamevariable(file, cutname, "elpt", 5, 200, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "eleta", 200, 800, 4, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #eta");
  //comparesamevariable(file, cutname, "elphi", -1, -1, 1, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "electron #phi");
  //comparesamevariable(file, cutname, "dielM", -1, 20000, 100, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("noselsctbar");
  coloropt.push_back(kBlack);
  legend.push_back("2022 Scouting Data");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //group_plotter(file, cutname, true);

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("mediumselsctec");
  coloropt.push_back(kBlack);
  legend.push_back("2022 Scouting Data");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //group_plotter(file, cutname, false);

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
  cutname.push_back("looseselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "dielM", -1, 20000, 100, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("mediumselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "dielM", -1, 20000, 100, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("tightselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "dielM", -1, 20000, 100, true, true, true, (float []){8e-1,1e5}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  //invmee_specialplot("noselsct_dielM");
  //invmee_specialplot("vetoselsct_dielM");
  //invmee_specialplot("looseselsct_dielM");
  //invmee_specialplot("mediumselsct_dielM");
  //invmee_specialplot("tightselsct_dielM");

  file.clear();
  cutname.clear();

  file.push_back(datahistfile);
  cutname.push_back("vetoselsct");

  //autoplotter(file, cutname);
  
  //fitinvmee("tightselsct_dielM");

  // Today
  
  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("tightnoneselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "dielM", -1, 20000, 100, false, true, true, (float []){0,50}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("tightlooseselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "dielM", -1, 20000, 100, false, true, true, (float []){8e-1,50}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(datahistfile);
  cutname.push_back("tightmediumselsct");
  coloropt.push_back(kBlack);
  legend.push_back("Scouting");
  histtype.push_back("p e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("lep");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "dielM", -1, 20000, 100, false, true, true, (float []){8e-1,50}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");

  fitinvmee("noselsct_leadsublead_dielM");
  fitinvmee("vetoselsct_leadsublead_dielM");
  fitinvmee("looseselsct_leadsublead_dielM");
  fitinvmee("mediumselsct_leadsublead_dielM");
  fitinvmee("tightselsct_leadsublead_dielM");
  //fitinvmee("tightnoneselsct_leadsublead_dielM");
  //fitinvmee("nonetightselsct_leadsublead_dielM");
  //fitinvmee("tightlooseselsct_leadsublead_dielM");
  //fitinvmee("loosetightselsct_leadsublead_dielM");
  //fitinvmee("tightmediumselsct_leadsublead_dielM");
  //fitinvmee("mediumtightselsct_leadsublead_dielM");

  return -1;
}
