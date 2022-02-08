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
TFile* datafile = TFile::Open("hists_Efmrl.root","READ");
TFile* dyfile = TFile::Open("hists_DY.root","READ");
TFile* sig3cmfile = TFile::Open("hists_M200dM20ctau3cm.root","READ");
TFile* sig30cmfile = TFile::Open("hists_M200dM20ctau30cm.root","READ");
TFile* sig1mfile = TFile::Open("hists_M200dM20ctau1m.root","READ");
TFile* sig3mfile = TFile::Open("hists_M200dM20ctau3m.root","READ");

TString seltext[2] = {"line1", "line2"};

std::vector<int> coloropt{1, kGreen-9, kGreen-7, kGreen-3, kGreen+2};
std::vector<TString> legendEntries{"l1", "l2", "l3", "l4", "l5", "l6"};

int comparemultihist(TString cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool isMConly=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float yrange[]=(float []){0.1,100}, bool normalized=false) {

  std::vector<TString> legNam;
  if(!isMConly) legNam.push_back("2018 data");
  legNam.push_back("c#tau = 3 cm");
  //legNam.push_back("c#tau = 30 cm");
  //legNam.push_back("c#tau = 1 m");
  //legNam.push_back("c#tau = 3 m");

  TString histname = cutname+"_"+var;
  TH1F* datahist;
  if(!isMConly) datahist = (TH1F*) datafile->Get(histname);
  TH1F* sig3cmhist = (TH1F*) sig3cmfile->Get(histname);
  //TH1F* sig30cmhist = (TH1F*) sig30cmfile->Get(histname);
  //TH1F* sig1mhist = (TH1F*) sig1mfile->Get(histname);
  //TH1F* sig3mhist = (TH1F*) sig3mfile->Get(histname);

  // Get the title from histogram title
  TString xtitle = sig3cmhist->GetTitle();

  std::vector<TH1F*> allhists;
  if(!isMConly) allhists.push_back(datahist);
  allhists.push_back(sig3cmhist);
  //allhists.push_back(sig30cmhist);
  //allhists.push_back(sig1mhist);
  //allhists.push_back(sig3mhist);

  if(rebin==-1) rebin = 1;
  xbinlow = xbinlow==-1?1:xbinlow;
  xbinhigh = xbinhigh==-1?(overflow?allhists[0]->GetNbinsX()+1:allhists[0]->GetNbinsX()):xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(rebin!=1) allhists[histctr]->Rebin(rebin);

    // Make changes to sig and bkg to enable good basic plotting
    double err = 0.0;
    //allhists[histctr]->SetBinContent(xbinlow,allhists[histctr]->IntegralAndError(0,xbinlow,err));
    //allhists[histctr]->SetBinError(xbinlow,err);
    err = 0.0;
    if(overflow) {
      allhists[histctr]->SetBinContent(xbinhigh,allhists[histctr]->IntegralAndError(xbinhigh,allhists[histctr]->GetNbinsX()+1,err));
      allhists[histctr]->SetBinError(xbinhigh,err);
    }

    allhists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allhists[histctr]->GetXaxis()->SetTitle(xtitle);
    allhists[histctr]->GetYaxis()->SetTitle("normalized number of events (a.u.)");

    allhists[histctr]->SetTitle("");

    allhists[histctr]->SetLineWidth(2);
    if(isMConly) allhists[histctr]->SetLineColor(coloropt[histctr+1]);
    else allhists[histctr]->SetLineColor(coloropt[histctr]);
  }

  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legNam, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,normalized);
  c1->SaveAs("./dirplots/"+cutname+"/"+cutname+"_"+var+".png");
  
  return -1;
}

// Comparison of the type cross-check between two histogram - filled v hollow
int crossChecktwohist(TFile* file, vector<TString> cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float yrange[]=(float []){0.1,100}) {

  std::vector<TString> legNam;
  legNam.push_back("L3#mu filter - DoubleMu33Displaced");
  legNam.push_back("L3#mu filter, p_{T}>16, |#eta|<2.5, d_{0}>0.1mm ");

  TH1F* histfilled = (TH1F*) file->Get(cutname[0]+"_"+var);
  TH1F* histhollow = (TH1F*) file->Get(cutname[1]+"_"+var);

  std::vector<TH1F*> allhists;
  allhists.push_back(histfilled);
  allhists.push_back(histhollow);

  // Get the title from histogram title
  TString xtitle = histfilled->GetTitle();

  if(rebin==-1) rebin = 1;
  xbinlow = xbinlow==-1?0:xbinlow;
  xbinhigh = xbinhigh==-1?allhists[0]->GetNbinsX()+1:xbinhigh;
  xbinlow = xbinlow/rebin;
  xbinhigh = xbinhigh/rebin;

  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {
    if(rebin!=1) allhists[histctr]->Rebin(rebin);

    // Make changes to sig and bkg to enable good basic plotting
    double err = 0.0;
    allhists[histctr]->SetBinContent(xbinlow,allhists[histctr]->IntegralAndError(0,xbinlow,err));
    allhists[histctr]->SetBinError(xbinlow,err);
    err = 0.0;
    allhists[histctr]->SetBinContent(xbinhigh,allhists[histctr]->IntegralAndError(xbinhigh,allhists[histctr]->GetNbinsX()+1,err));
    allhists[histctr]->SetBinError(xbinhigh,err);

    allhists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allhists[histctr]->GetXaxis()->SetTitle(xtitle);
    allhists[histctr]->GetYaxis()->SetTitle("number of events");
    
    allhists[histctr]->SetTitle("");
  }

  allhists[0]->SetFillColor(kBlue);
  allhists[1]->SetLineColor(kRed);
  
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter(allhists, legNam, allhists[0]->GetXaxis()->GetTitle(),allhists[0]->GetYaxis()->GetTitle(),legPos,logY,yrange,false);
  c1->SaveAs("./dirplots/"+cutname[0]+"_"+cutname[1]+"/"+cutname[0]+"_"+cutname[1]+"_"+var+".png");

  return -1;
}

int makeratehist(TString cutname, TString var, int xbinlow, int xbinhigh, int rebin=-1, bool logY=false, bool overflow=false, float legPos[]=(float []){0.7,0.75,0.95,1}, float seltextpos[]=(float []){0.1,1}, float yrange[]=(float []){0.1,1}, float drawsignalline=-1.0) {

  std::vector<TString> legNam;
  legNam.push_back("2018 data");
  legNam.push_back("c#tau = 3 cm");
  //legNam.push_back("c#tau = 30 cm");
  //legNam.push_back("c#tau = 1 m");
  //legNam.push_back("c#tau = 3 m");

  TString histname = cutname+"_"+var;
  TH1F* datahist;
  datahist = (TH1F*) datafile->Get(histname);
  datahist->Scale(datasf);
  TH1F* sig3cmhist = (TH1F*) sig3cmfile->Get(histname);
  sig3cmhist->Scale(sig3cmsf);
  //TH1F* sig30cmhist = (TH1F*) sig30cmfile->Get(histname);
  //sig30cmhist->Scale(sig30cmsf);
  //TH1F* sig1mhist = (TH1F*) sig1mfile->Get(histname);
  //sig1mhist->Scale(sig1msf);
  //TH1F* sig3mhist = (TH1F*) sig3mfile->Get(histname);
  //sig3mhist->Scale(sig3msf);

  // Get the title from histogram title
  TString xtitle = sig3cmhist->GetTitle();

  std::vector<TH1F*> allhists;
  std::vector<TH1F*> allratehists;
  allhists.push_back(datahist);
  allhists.push_back(sig3cmhist);
  //allhists.push_back(sig30cmhist);
  //allhists.push_back(sig1mhist);
  //allhists.push_back(sig3mhist);

  double axisscale = 0.0;
  for(unsigned int histctr=0; histctr<allhists.size(); histctr++) {

    allratehists.push_back((TH1F*) allhists[histctr]->Clone());
    for(unsigned int bincnt=0; bincnt<allhists[histctr]->GetNbinsX(); bincnt++) {
      double err = 0;
      allratehists[histctr]->SetBinContent(bincnt, allhists[histctr]->IntegralAndError(bincnt,allhists[histctr]->GetNbinsX()+1,err));
      allratehists[histctr]->SetBinError(bincnt, err);
    }
    allratehists[histctr]->GetXaxis()->SetRange(xbinlow, xbinhigh);
    allratehists[histctr]->GetXaxis()->SetTitle(xtitle);
    allratehists[histctr]->GetYaxis()->SetTitle("rate / Hz");

    allratehists[histctr]->SetTitle("");

    allratehists[histctr]->SetLineWidth(2);
    allratehists[histctr]->SetLineColor(coloropt[histctr]);
    if(histctr==1) {
      axisscale = allratehists[0]->GetBinContent(0)*0.8/allratehists[1]->GetBinContent(0);
      allratehists[histctr]->Scale(axisscale);
    }
    if(histctr>1) {
      allratehists[histctr]->Scale(axisscale);
    }
  }

  cout<<allratehists[1]->GetMaximum()<<endl;
  TCanvas* c1;
  c1 = new TCanvas();
  c1 = enhance_plotter_rate(allratehists, legNam, allratehists[0]->GetXaxis()->GetTitle(),allratehists[0]->GetYaxis()->GetTitle(),legPos,yrange,logY,false);

  TPad* pad = (TPad*) c1->FindObject("pad3");
  pad->cd();
  TGaxis *axis;
  axis = new TGaxis(allratehists[0]->GetBinLowEdge(xbinhigh+1),yrange[0],allratehists[0]->GetBinLowEdge(xbinhigh+1),yrange[1],yrange[0]/axisscale,yrange[1]/axisscale,510,"-L");
  axis->SetLineColor(coloropt[1]);
  axis->SetLabelColor(coloropt[1]);
  axis->SetLabelFont(132);
  axis->SetLabelSize(0.06);
  axis->SetLabelOffset(-0.035);
  axis->Draw();

  if(drawsignalline!=-1.0) {
    drawsignalline = allratehists[1]->GetMaximum();
    TLine *signalline = new TLine(allratehists[0]->GetBinLowEdge(xbinlow),drawsignalline,allratehists[0]->GetBinLowEdge(xbinhigh+1),drawsignalline);
    signalline->SetLineWidth(2);
    signalline->SetLineColor(coloropt[1]);
    signalline->SetLineStyle(9);
    signalline->Draw();
  }

  TLatex sel;
  sel.SetTextFont(132);
  sel.SetTextSize(0.05);
  sel.DrawLatex(seltextpos[0], seltextpos[1]+0.1*(yrange[1]-yrange[0]), seltext[0]);
  sel.DrawLatex(seltextpos[0], seltextpos[1], seltext[1]);
  
  c1->SaveAs("./dirplots/"+cutname+"/"+cutname+"_"+var+"_ratehist.png");
  
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

int newplotter() {

  std::vector<int> coloroptrate{1, kRed+2, kRed-3, kRed-7, kRed-9};
  coloropt = coloroptrate;
  seltext[0] = "N#mu#geq1, p_{T}#geq38 GeV, |#eta|<2.5, d_{0}>0.01 cm";
  seltext[1] = "Ne/#gamma#geq1, p_{T}>38 GeV, |#eta|<2.65";
  //makeratehist("sel3recoeg", "egpt", 66, 150, 1, false, false, (float []){0.55,0.775,0.75,0.975}, (float []){60,0.8}, (float []){0.1,1.39}, 0.875);

  std::vector<int> coloroptschemeselelevetoid{1, 11, 4, 38};
  coloropt = coloroptschemeselelevetoid;
  std::vector<TFile*> fileselelevetoid;
  std::vector<TString> nameselelevetoid;
  std::vector<TString> legselelevetoid;
  fileselelevetoid.push_back(datafile);
  nameselelevetoid.push_back("selelevetoidusrecoebus");
  legselelevetoid.push_back("data, vetoid");
  fileselelevetoid.push_back(datafile);
  nameselelevetoid.push_back("basicselusrecoebus");
  legselelevetoid.push_back("data, noid");
  fileselelevetoid.push_back(dyfile);
  nameselelevetoid.push_back("selelevetoidusrecoebus");
  legselelevetoid.push_back("DY, egvetoid");
  fileselelevetoid.push_back(dyfile);
  nameselelevetoid.push_back("basicselusrecoebus");
  legselelevetoid.push_back("DY, noid");
  legendEntries = legselelevetoid;  
  //comparesamevariable(fileselelevetoid, nameselelevetoid, "leadsubleadM", -1, 120, 1, true, true, true, (float []){1e-3,3e-1}, (float []){0.2,0.65,0.45,0.95}, true, "M(e/#gamma_{1},e/#gamma_{2}) / GeV");

  std::vector<int> coloroptschemebasicsel{14, 1, 16, 8, 46};
  coloropt = coloroptschemebasicsel;
  std::vector<TFile*> filebasicsel;
  std::vector<TString> namebasicsel;
  std::vector<TString> legbasicsel;
  filebasicsel.push_back(datafile);
  namebasicsel.push_back("basicselrecoeb");
  legbasicsel.push_back("data");
  filebasicsel.push_back(datafile);
  namebasicsel.push_back("selelevetozwindidusrecoebus");
  legbasicsel.push_back("data, egveto, 75<M<95");
  filebasicsel.push_back(datafile);
  namebasicsel.push_back("selelevetozoppoidusrecoebus");
  legbasicsel.push_back("data, egveto, Z veto");
  filebasicsel.push_back(dyfile);
  namebasicsel.push_back("selelevetozwindidusrecoebus");
  legbasicsel.push_back("DY, egvetom 75<M<95");
  filebasicsel.push_back(sig3cmfile);
  namebasicsel.push_back("basicselrecoeb");
  legbasicsel.push_back("signal 3 cm");
  legendEntries = legbasicsel;  
  //comparesamevariable(filebasicsel, namebasicsel, "leadegin5x5noiseclnd", -1, 400, 10, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{1} #sigmai#etai#eta5x5 (noise cleaned)");
  //comparesamevariable(filebasicsel, namebasicsel, "subleadegin5x5noiseclnd", -1, 400, 10, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{2} #sigmai#etai#eta5x5 (noise cleaned)");
  //comparesamevariable(filebasicsel, namebasicsel, "leadeghovereoversupcluse", -1, 150, 2, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{1} H/E");
  //comparesamevariable(filebasicsel, namebasicsel, "subleadeghovereoversupcluse", -1, 150, 2, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{2} H/E");
  //comparesamevariable(filebasicsel, namebasicsel, "leadegecalpfclustisoovere", 50, 150, 4, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{1} ecal iso./E");
  //comparesamevariable(filebasicsel, namebasicsel, "subleadegecalpfclustisoovere", 50, 150, 4, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{2} ecal iso./E");
  //comparesamevariable(filebasicsel, namebasicsel, "leadeghcalpfclustisoovere", 50, 150, 4, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{1} hcal iso./E");
  //comparesamevariable(filebasicsel, namebasicsel, "subleadeghcalpfclustisoovere", 50, 150, 1, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{2} hcal iso./E");
  //comparesamevariable(filebasicsel, namebasicsel, "leadegpixelmchvar_s2", 40, 200, 2, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "pixel match");

  filebasicsel.clear();
  namebasicsel.clear();
  legbasicsel.clear();
  filebasicsel.push_back(datafile);
  namebasicsel.push_back("basicselrecoee");
  legbasicsel.push_back("data");
  filebasicsel.push_back(datafile);
  namebasicsel.push_back("selelevetozwindidusrecoeeus");
  legbasicsel.push_back("data, egveto, 75<M<95");
  filebasicsel.push_back(datafile);
  namebasicsel.push_back("selelevetozoppoidusrecoeeus");
  legbasicsel.push_back("data, egveto, Z veto");
  filebasicsel.push_back(dyfile);
  namebasicsel.push_back("selelevetozwindidusrecoeeus");
  legbasicsel.push_back("DY, egvetom 75<M<95");
  filebasicsel.push_back(sig3cmfile);
  namebasicsel.push_back("basicselrecoee");
  legbasicsel.push_back("signal 3 cm");
  legendEntries = legbasicsel;  
  //comparesamevariable(filebasicsel, namebasicsel, "leadegin5x5noiseclnd", -1, 400, 10, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{1} #sigmai#etai#eta5x5 (noise cleaned)");
  //comparesamevariable(filebasicsel, namebasicsel, "subleadegin5x5noiseclnd", 100, 800, 10, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{2} #sigmai#etai#eta5x5 (noise cleaned)");
  //comparesamevariable(filebasicsel, namebasicsel, "leadeghovereoversupcluse", -1, 150, 2, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{1} H/E");
  //comparesamevariable(filebasicsel, namebasicsel, "subleadeghovereoversupcluse", -1, 150, 2, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{2} H/E");
  //comparesamevariable(filebasicsel, namebasicsel, "leadegecalpfclustisoovere", 50, 150, 4, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{1} ecal iso./E");
  //comparesamevariable(filebasicsel, namebasicsel, "subleadegecalpfclustisoovere", 50, 150, 4, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{2} ecal iso./E");
  //comparesamevariable(filebasicsel, namebasicsel, "leadeghcalpfclustisoovere", 50, 150, 4, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{1} hcal iso./E");
  //comparesamevariable(filebasicsel, namebasicsel, "subleadeghcalpfclustisoovere", 50, 150, 1, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma_{2} hcal iso./E");
  //comparesamevariable(filebasicsel, namebasicsel, "leadegpixelmchvar_s2", 40, 200, 2, true, true, true, (float []){1e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "pixel match");

  std::vector<int> coloroptschemegennosel{4, 2, 46};
  coloropt = coloroptschemegennosel;
  std::vector<TFile*> filegennosel;
  std::vector<TString> namegennosel;
  std::vector<TString> leggennosel;
  filegennosel.push_back(dyfile);
  namegennosel.push_back("gennoselgeneg");
  leggennosel.push_back("DY, mother Z");
  filegennosel.push_back(sig3cmfile);
  namegennosel.push_back("gennoselgeneg");
  leggennosel.push_back("signal 3 cm");
  filegennosel.push_back(sig1mfile);
  namegennosel.push_back("gennoselgeneg");
  leggennosel.push_back("signal 1 m");
  legendEntries = leggennosel;  
  //comparesamevariable(filegennosel, namegennosel, "t1", 1000, 2000, 5, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma ecal time / ns");
  //comparesamevariable(filegennosel, namegennosel, "log10d0", -1, -1, 1, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "e/#gamma log_{10}d_{0} / log_{10}cm");

  std::vector<TFile*> filegenbarsel;
  std::vector<TString> namegenbarsel;
  std::vector<TString> leggenbarsel;
  filegenbarsel.push_back(dyfile);
  namegenbarsel.push_back("genbarselgeneg");
  leggenbarsel.push_back("DY, mother Z");
  filegenbarsel.push_back(sig3cmfile);
  namegenbarsel.push_back("genbarselgeneg");
  leggenbarsel.push_back("signal 3 cm");
  filegenbarsel.push_back(sig1mfile);
  namegenbarsel.push_back("genbarselgeneg");
  leggenbarsel.push_back("signal 1 m");
  legendEntries = leggenbarsel;  
  //comparesamevariable(filegenbarsel, namegenbarsel, "t1", 1000, 1200, 2, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen e/#gamma ecal time / ns");
  //comparesamevariable(filegenbarsel, namegenbarsel, "t0", 1000, 1200, 2, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen e/#gamma ecal time light / ns");
  //comparesamevariable(filegenbarsel, namegenbarsel, "t1mt0", -1, -1, 1, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen e/#gamma ecal time diff / ns");

  std::vector<int> coloroptgenbarselptgt10{kBlue, kRed+2, kRed-3, kRed-7, kRed-9};
  coloropt = coloroptgenbarselptgt10;
  std::vector<TFile*> filegenbarselptgt10;
  std::vector<TString> namegenbarselptgt10;
  std::vector<TString> leggenbarselptgt10;
  filegenbarselptgt10.push_back(dyfile);
  namegenbarselptgt10.push_back("genbarselptgt10geneg");
  leggenbarselptgt10.push_back("DY, mother Z");
  filegenbarselptgt10.push_back(sig3cmfile);
  namegenbarselptgt10.push_back("genbarselptgt10geneg");
  leggenbarselptgt10.push_back("signal 3 cm");
  filegenbarselptgt10.push_back(sig30cmfile);
  namegenbarselptgt10.push_back("genbarselptgt10geneg");
  leggenbarselptgt10.push_back("signal 30 cm");
  filegenbarselptgt10.push_back(sig1mfile);
  namegenbarselptgt10.push_back("genbarselptgt10geneg");
  leggenbarselptgt10.push_back("signal 1 m");
  filegenbarselptgt10.push_back(sig3mfile);
  namegenbarselptgt10.push_back("genbarselptgt10geneg");
  leggenbarselptgt10.push_back("signal 3 m");
  legendEntries = leggenbarselptgt10;  
  //comparesamevariable(filegenbarselptgt10, namegenbarselptgt10, "t1", 1000, 1200, 2, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen e/#gamma ecal t_{e}+t_{LLP} / ns");
  //comparesamevariable(filegenbarselptgt10, namegenbarselptgt10, "t0", 1000, 1200, 2, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen e/#gamma ecal t_{prompt} / ns");
  //comparesamevariable(filegenbarselptgt10, namegenbarselptgt10, "t1mt0", 900, 1600, 10, true, true, true, (float []){4e-3,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen e/#gamma ecal time delay / ns");

  /* The code for analyzing different 
     tracking seeds to displaced electrons 
     starts from here 
  */

  std::vector<int> coloroptgenbarnewtracks{kBlue, kRed+2, kRed-3, kRed-7, kRed-9};
  coloropt = coloroptgenbarnewtracks;
  std::vector<TFile*> filegennewtr;
  std::vector<TString> namegennewtr;
  std::vector<TString> leggennewtr;
  filegennewtr.push_back(sig3cmfile);
  namegennewtr.push_back("gennoselgeneg");
  leggennewtr.push_back("signal 3 cm");
  legendEntries = leggennewtr;  
  //comparesamevariable(filegennewtr, namegennewtr, "egmult", 5, 10, 1, true, true, true, (float []){5e-1,2}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron multiplicity");
  //comparesamevariable(filegennewtr, namegennewtr, "pt", 50, 95, 1, true, true, true, (float []){5e-4,2e-1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron p_{T} / GeV");
  //comparesamevariable(filegennewtr, namegennewtr, "eta", -1, -1, 1, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron #eta");
  //comparesamevariable(filegennewtr, namegennewtr, "phi", -1, -1, 1, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron #phi");
  //comparesamevariable(filegennewtr, namegennewtr, "log10d0", 250, -1, 10, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron log_{10} d_{0} / log_{10} cm");
  //comparesamevariable(filegennewtr, namegennewtr, "leadpt", 50, 95, 1, true, true, true, (float []){5e-4,2e-1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron_{1} p_{T} / GeV");
  //comparesamevariable(filegennewtr, namegennewtr, "leadeta", -1, -1, 1, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron_{1} #eta");
  //comparesamevariable(filegennewtr, namegennewtr, "leadphi", -1, -1, 1, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron_{1} #phi");
  //comparesamevariable(filegennewtr, namegennewtr, "leadlog10d0", 250, -1, 10, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron_{1} log_{10} d_{0} / log_{10} cm");
  //comparesamevariable(filegennewtr, namegennewtr, "subleadpt", 50, 95, 1, true, true, true, (float []){5e-4,2e-1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron_{2} p_{T} / GeV");
  //comparesamevariable(filegennewtr, namegennewtr, "subleadeta", -1, -1, 1, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron_{2} #eta");
  //comparesamevariable(filegennewtr, namegennewtr, "subleadphi", -1, -1, 1, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron_{2} #phi");
  //comparesamevariable(filegennewtr, namegennewtr, "subleadlog10d0", 250, -1, 10, true, true, true, (float []){5e-4,1}, (float []){0.5,0.65,0.75,0.95}, true, "gen electron_{2} log_{10} d_{0} / log_{10} cm");
  
  
  return -1;
}
