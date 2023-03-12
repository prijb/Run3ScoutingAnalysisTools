#include <string>
#include <sstream>
#include <map>
#include <vector>
#include "enhance_plotter.C"
#include "RooMsgService.h"

TString cutdeets = "Cut details";
TFile* jpsimodhistfile = TFile::Open("hists_JPsi_mod.root","READ");
TFile* jpsihistfile = TFile::Open("hists_JPsi.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};
std::vector<TString> histtype{"p e1", "hist same"};
std::vector<int> markerstyle{20, 24};
std::vector<int> markersize{10, 10};
std::vector<TString> legendmarkerstyle{"lep", "l"};
std::vector<double> scale{-1, 1};

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
    TH1F* temphistholder = (TH1F*) file[histctr]->Get(cutwithvar);
    allhists.push_back((TH1F*)temphistholder->Clone(temphistholder->GetName()));
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

    if(scale[0]!=-1) allhists[histctr]->Scale(scale[histctr]);
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
    allhists[histctr]->SetMarkerColor(coloropt[histctr]);
    if(markerstyle[histctr]==0) allhists[histctr]->SetFillColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos, false, logY,yrange,normalize,histtype,markerstyle,markersize,legendmarkerstyle);
  if(!normalize) c1->SaveAs("./dirplots_jpsi/"+foldername+"/"+var+".png");
  else  c1->SaveAs("./dirplots_jpsi/"+foldername+"/"+var+"_normed.png");

  return -1;
}

int linearfit(TString selection, TString xaxistitle, TString filesavename) {

  auto pthistorig = (TH1F*) jpsihistfile->Get(selection);
  auto pthist = (TH1F*) pthistorig->Clone(pthistorig->GetName());

  pthist->SetTitle("");
  pthist->GetXaxis()->SetRange(0,70);

  TF1* fitbkg1 = new TF1("fitbkg1","exp([0]+[1]*x)",4,10);
  fitbkg1->SetParameters(2,-0.1);
  pthist->Fit("fitbkg1","R");
  double b10 = fitbkg1->GetParameter(0);
  double be10 = fitbkg1->GetParError(0);
  double b11 = fitbkg1->GetParameter(1);
  double be11 = fitbkg1->GetParError(1);
  
  TF1* fitbkg2 = new TF1("fitbkg2","exp([0]+[1]*x)",20,30);
  fitbkg2->SetParameters(2,-0.1);
  pthist->Fit("fitbkg2","R");
  double b20 = fitbkg2->GetParameter(0);
  double be20 = fitbkg2->GetParError(0);
  double b21 = fitbkg2->GetParameter(1);
  double be21 = fitbkg2->GetParError(1);

  TF1* fitbkg = new TF1("fitbkg","exp([0]+[1]*x)+exp([2]+[3]*x)",4,60);
  fitbkg->SetParameters(b10, b11, b20, b21);
  fitbkg->SetParLimits(0, b10-5*be10, b10+5*be10);
  fitbkg->SetParLimits(1, b11-5*be11, b11+5*be11);
  fitbkg->SetParLimits(2, b20-5*be20, b20+5*be20);
  fitbkg->SetParLimits(3, b21-5*be21, b21+5*be21);
  pthist->Fit("fitbkg","R");
  fitbkg->SetLineColor(kRed);

  cout<<"Intgrate between 2 and 60: "<<fitbkg->Integral(2,60)<<endl;
  cout<<"Intgrate between 12 and 60: "<<fitbkg->Integral(12,60)<<endl;
  cout<<"Intgrate between 23 and 60: "<<fitbkg->Integral(23,60)<<endl;
  
  TCanvas* c1;
  c1 = enhance_plotter({pthist}, {"J/#psi#rightarrowee"}, xaxistitle, "number of events", (float []){0.55,0.7,0.75,0.95}, false, true, (float []){8e-1,2e8}, false, {"hist e1"}, {20}, {2}, {"le"});
  fitbkg->Draw("same");
  c1->SaveAs(filesavename);

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
  demo->GetYaxis()->SetTitle("efficiency");
  demo->SetLineColorAlpha(kWhite,1);
  demo->SetFillColorAlpha(kWhite,1);
  allhists.push_back(demo);
  
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legendEntries, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(), (float []){-1,0.65,0.9,0.95},false,false,(float []){0.0,1.1},false,histtype,markerstyle,markersize,legendmarkerstyle);
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
  c1->SaveAs("./dirplots_jpsi/"+((TString)file[0]->GetName()).ReplaceAll(".root","")+"/"+cutnames[0]+"_eff.png");
  
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

  file.push_back(jpsihistfile);
  cutname.push_back("noselgen_gen");
  coloropt.push_back(kBlue);
  legend.push_back("J/#psi#rightarrowee");
  histtype.push_back("hist e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("le");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "jpsi_mult", 5, 10, 1, true, true, true, (float []){8e-1,2e8}, (float []){0.55,0.7,0.75,0.95}, false, "J/#psi multiplicity");
  //comparesamevariable(file, cutname, "jpsi_m", 400, 450, 1, true, true, true, (float []){8e-1,2e8}, (float []){0.55,0.7,0.75,0.95}, false, "J/#psi mass [GeV]");
  //comparesamevariable(file, cutname, "jpsi_pt", -1, 100, 1, true, true, true, (float []){8e-1,2e6}, (float []){0.55,0.7,0.75,0.95}, false, "J/#psi p_{T} [GeV]");
  //comparesamevariable(file, cutname, "jpsi_eta", -1, -1, 10, false, true, true, (float []){8e-1,30000}, (float []){0.55,0.7,0.75,0.95}, false, "J/#psi #eta");
  //comparesamevariable(file, cutname, "jpsi_phi", -1, -1, 1, false, true, true, (float []){8e-1,30000}, (float []){0.55,0.7,0.75,0.95}, false, "J/#psi #phi");
  //comparesamevariable(file, cutname, "el_mult", 5, 15, 1, true, true, true, (float []){8e-1,2e9}, (float []){0.55,0.7,0.75,0.95}, false, "electron multiplicity");
  //comparesamevariable(file, cutname, "el_pt", -1, 70, 1, true, true, true, (float []){8e-1,2e7}, (float []){0.55,0.7,0.75,0.95}, false, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "el_eta", -1, -1, 10, false, true, true, (float []){8e-1,60000}, (float []){0.55,0.7,0.75,0.95}, false, "electron #eta");
  //comparesamevariable(file, cutname, "el_pt", -1, 70, 1, true, true, true, (float []){8e-1,2e7}, (float []){0.55,0.7,0.75,0.95}, false, "electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "el_eta", -1, -1, 10, false, true, true, (float []){8e-1,60000}, (float []){0.55,0.7,0.75,0.95}, false, "electron #eta");
  //comparesamevariable(file, cutname, "el_phi", -1, -1, 1, false, true, true, (float []){0,50000}, (float []){0.55,0.7,0.75,0.95}, false, "electron #phi");
  //comparesamevariable(file, cutname, "leadel_pt", -1, 70, 1, true, true, true, (float []){8e-1,2e7}, (float []){0.55,0.7,0.75,0.95}, false, "electron_{1} p_{T} [GeV]");
  //comparesamevariable(file, cutname, "leadel_eta", -1, -1, 10, false, true, true, (float []){8e-1,60000}, (float []){0.55,0.7,0.75,0.95}, false, "electron_{1} #eta");
  //comparesamevariable(file, cutname, "leadel_phi", -1, -1, 1, false, true, true, (float []){0,50000}, (float []){0.55,0.7,0.75,0.95}, false, "electron_{1} #phi");
  //comparesamevariable(file, cutname, "subleadel_pt", -1, 70, 1, true, true, true, (float []){8e-1,2e7}, (float []){0.55,0.7,0.75,0.95}, false, "electron_{2} p_{T} [GeV]");
  //comparesamevariable(file, cutname, "subleadel_eta", -1, -1, 10, false, true, true, (float []){8e-1,60000}, (float []){0.55,0.7,0.75,0.95}, false, "electron_{2} #eta");
  //comparesamevariable(file, cutname, "subleadel_phi", -1, -1, 1, false, true, true, (float []){0,50000}, (float []){0.55,0.7,0.75,0.95}, false, "electron_{2} #phi");
  //comparesamevariable(file, cutname, "diel_M", 80, 200, -1, true, true, true, (float []){8e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "M(e,e) [GeV]");
  //comparesamevariable(file, cutname, "diel_deta", 650, 750, 1, true, true, true, (float []){8e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#eta(e,e) [GeV]");
  //comparesamevariable(file, cutname, "diel_dphi", 400, 1000, 10, true, true, true, (float []){8e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "#Delta#phi(e,e) [GeV]");
  //comparesamevariable(file, cutname, "diel_dR", -1, 600, 10, true, true, true, (float []){8e-1,1e7}, (float []){0.6,0.7,0.85,0.95}, false, "#DeltaR(e,e) [GeV]");

  //linearfit("noselgen_gen_leadel_pt", "electron_{1} p_{T} [GeV]", "JPsiMC_noselgen_gen_leadel_pt.png");
  //linearfit("noselgen_gen_subleadel_pt", "electron_{2} p_{T} [GeV]", "JPsiMC_noselgen_gen_subleadel_pt.png");
  
  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(jpsihistfile);
  cutname.push_back("nosel_sct");
  coloropt.push_back(kBlue);
  legend.push_back("J/#psi#rightarrowee MC");
  histtype.push_back("hist e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("le");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "el_mult", 5, 15, 1, true, true, true, (float []){8e-1,2e6}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron multiplicity");
  //comparesamevariable(file, cutname, "el_pt", -1, 100, 1, true, true, true, (float []){8e-1,4e6}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron p_{T} [GeV]");
  //comparesamevariable(file, cutname, "el_eta", -1, -1, 10, false, true, true, (float []){8e-1,20000}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron #eta");
  //comparesamevariable(file, cutname, "el_phi", -1, -1, 1, false, true, true, (float []){8e-1,15000}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron #phi");
  //comparesamevariable(file, cutname, "el_trkpt", -1, 100, 1, true, true, true, (float []){8e-1,4e6}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron track p_{T} [GeV]");
  //comparesamevariable(file, cutname, "el_trketa", -1, -1, 10, false, true, true, (float []){8e-1,20000}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron track #eta");
  //comparesamevariable(file, cutname, "el_trkphi", -1, -1, 1, false, true, true, (float []){8e-1,15000}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron track #phi");
  //comparesamevariable(file, cutname, "bar_el_charge", -1, -1, 1, true, true, true, (float []){8e-1,4e7}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron charge");
  //comparesamevariable(file, cutname, "diel_M", 900, 3000, 20, true, true, true, (float []){8e-1,4e6}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron M(e,e) [GeV]");
  //comparesamevariable(file, cutname, "diel_trkM", 900, 3000, 10, true, true, true, (float []){8e-1,4e6}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron track M(e,e) [GeV]");
  //comparesamevariable(file, cutname, "leadsubleaddiel_M", 900, 3000, 10, true, true, true, (float []){8e-1,4e4}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron M(e_{1},e_{2}) [GeV]");
  //comparesamevariable(file, cutname, "leadsubleaddiel_trkM", 900, 3000, 5, true, true, true, (float []){8e-1,4e4}, (float []){0.55,0.7,0.75,0.95}, false, "scouting electron track M(e_{1},e_{2}) [GeV]");
  //comparesamevariable(file, cutname, "leadsubleaddiel_chargeprod", -1, -1, 1, true, true, true, (float []){8e-1,4e7}, (float []){0.55,0.7,0.75,0.95}, false, "q(e_{1})#timesq(e_{2})");
  //comparesamevariable(file, cutname, "leadsubleaddiel_trkchargeprod", -1, -1, 1, true, true, true, (float []){8e-1,4e7}, (float []){0.55,0.7,0.75,0.95}, false, "tracker p_{T} sorted q(e_{1})#timesq(e_{2})");
  //comparesamevariable(file, cutname, "leadbarsubleadbardiel_M", 900, 3000, 10, true, true, true, (float []){8e-1,4e4}, (float []){0.55,0.7,0.75,0.95}, false, "EB scouting electron M(e_{1},e_{2}) [GeV]");
  //comparesamevariable(file, cutname, "leadbarsubleadbardiel_trkM", 900, 3000, 5, true, true, true, (float []){8e-1,4e4}, (float []){0.55,0.7,0.75,0.95}, false, "EE scouting electron track M(e_{1},e_{2}) [GeV]");
  //comparesamevariable(file, cutname, "leadecsubleadecdiel_M", 900, 3000, 10, true, true, true, (float []){8e-1,4e4}, (float []){0.55,0.7,0.75,0.95}, false, "EB scouting electron M(e_{1},e_{2}) [GeV]");
  //comparesamevariable(file, cutname, "leadecsubleadecdiel_trkM", 900, 3000, 5, true, true, true, (float []){8e-1,4e4}, (float []){0.55,0.7,0.75,0.95}, false, "EE scouting electron track M(e_{1},e_{2}) [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(jpsihistfile);
  cutname.push_back("nosel_sct_bar_el");
  coloropt.push_back(kBlue);
  legend.push_back("dEtaIn");
  histtype.push_back("hist e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("le");
  file.push_back(jpsihistfile);
  cutname.push_back("nosel_sct_bar_el_trksc");
  coloropt.push_back(kRed);
  legend.push_back("#Delta#eta(trk, SC)");
  histtype.push_back("hist same e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("le");

  legendEntries = legend;
  //comparesamevariable(file, cutname, "detain", -1, -1, 10, true, true, true, (float []){8e-1,4e6}, (float []){0.65,0.7,0.85,0.95}, false, "value");

  legend.clear();
  legend.push_back("dPhiIn");
  legend.push_back("#Delta#phi(trk, SC)");
  legendEntries = legend;
  //comparesamevariable(file, cutname, "dphiin", -1, -1, 10, true, true, true, (float []){8e-1,4e6}, (float []){0.55,0.7,0.75,0.95}, false, "value");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(jpsihistfile);
  cutname.push_back("noselgen_gen_el");
  coloropt.push_back(kBlue-2);
  legend.push_back("2022 gen p_{T}");
  histtype.push_back("hist e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("le");

  file.push_back(jpsihistfile);
  cutname.push_back("noselgenAnosel_sctmchgenel_el");
  coloropt.push_back(kBlue);
  legend.push_back("2022 scouting");
  histtype.push_back("hist same e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("le");

  file.push_back(jpsimodhistfile);
  cutname.push_back("noselgen_gen_el");
  coloropt.push_back(kRed-2);
  legend.push_back("2023 gen p_{T}");
  histtype.push_back("hist same e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("le");

  file.push_back(jpsimodhistfile);
  cutname.push_back("noselgenAnosel_sctmchgenel_el");
  coloropt.push_back(kRed);
  legend.push_back("2023 scouting");
  histtype.push_back("hist same e1");
  markerstyle.push_back(20);
  markersize.push_back(2);
  legendmarkerstyle.push_back("le");

  legendEntries = legend;
  comparesamevariable(file, cutname, "pt", -1, 50, 1, true, true, true, (float []){8e-1,4e7}, (float []){0.55,0.7,0.75,0.95}, false, "p_{T} [GeV]");

  file.clear();
  cutname.clear();
  coloropt.clear();
  legend.clear();  
  histtype.clear();
  markerstyle.clear();
  markersize.clear();
  legendmarkerstyle.clear();

  file.push_back(jpsihistfile);
  cutname.push_back("noselgen_gen_el_pt");
  cutname.push_back("noselgenAnosel_sctmchgenel_el_pt");
  coloropt.push_back(kBlue);
  legend.push_back("2022 scouting (seeded)");
  histtype.push_back("e1");
  markerstyle.push_back(104);
  markersize.push_back(2);
  legendmarkerstyle.push_back("pe");
  legendmarkerstyle.clear();

  file.push_back(jpsimodhistfile);
  cutname.push_back("noselgen_gen_el_pt");
  cutname.push_back("noselgenAnosel_sctmchgenel_el_pt");
  coloropt.push_back(kRed);
  legend.push_back("2023 Scouting (unseeded)");
  histtype.push_back("same e1");
  markerstyle.push_back(104);
  markersize.push_back(2);
  legendmarkerstyle.push_back("pe");

  legendEntries = legend;
  vector<double> binspt{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40};
  vector<float> legpos = {0.35, 0.15, 0.6, 0.35};
  efficiency(file, cutname, &legpos[0], binspt.size()-1, &binspt[0], "p_{T} / GeV");

  return -1;
}
