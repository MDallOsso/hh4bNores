// M.Dall'Osso
// last change....
// basic 2 distributions comparison
//....................................

#include <TROOT.h>
#include <TStyle.h>
#include <TH1F.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TGaxis.h>
#include <iostream>
#include <TLatex.h>
#include <THStack.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip> 
//kBlue+2, kCyan+2, kBlack, kViolet+2,  debug
int colors[] = {kBlue, kOrange+2, kCyan+2, kBlack, kViolet+2, kGreen+2, kPink,kMagenta,
                kSpring+8,   kYellow+2, kMagenta-6, kGreen-6, kOrange-6, kBlue-6, kCyan-6, kGreen-6, kPink-6, kRed-6,
                kMagenta-8, kOrange-8, kBlue-8, kCyan-8, kGreen-8, kPink-8, kRed-8,
                kMagenta-4, kViolet-8, kAzure-8, kSpring-8, kYellow-8, kGreen-8,
                kViolet-4, kAzure-4, kSpring-4, kYellow-4, kGreen-4};
const int maxCol = 36;

// input files - parameters
//******************************************

// 13Tev - 3 parameters values - 636 files

string L [] = {"10p0","12p5","12p0","15p0","1p0","2p4","3p5","5p0","7p5","m10p0","m12p5","m15p0","m1p0","m2p4","m3p5","m5p0","m5","m7p5"};
string Y[]  = {"0p5","0p75","1p0","1p25","1p5","1p75","2p0","2p25","2p5"};
string C2[] = {"c0p0","c0p5","c1p0","c2p0","c3p0","cm0p5","cm1p0","cm3p0","cm2p0"};
string Totsamples = "636";
int CMenergy = 13;   //tev
int pars = 3;        //space parameters dimension
bool Privat = false; //true
int Maxtotclu = 20;
string Inputroot  = "Distros_3par_20000ev_13Tev.root";
string Inputfolder = "Results/FCL_";
TString Outfolder = "plots_3par_13Tev/"; //debug  _2var
string testoption = ""; //to add 'debug'  _2var

//******************************************

TString filename; 
TString mainDIR;
TString suffix;
bool initialized = false;
bool lgRIGHT = true;
bool lgTOP = true;
int totNclu;
string Outname;
std::vector< std::vector<string> > clu;

bool init() {

  //dafault file name format
  std::stringstream inputclusters;
  std::vector<string> samples;

  filename = Inputroot;

  //read cluster result
  inputclusters << Inputfolder << pars << "par_" << CMenergy << "Tev_Nclu" << totNclu << testoption << ".dat"; //<< testoption
  string infname = inputclusters.str();
  ifstream inresfile;
  inresfile.open(infname.c_str());

  return name.str();
}

void drawCMS() {
  TLatex* text=new TLatex(0.2, 0.95, "CMS Data, 2012, #sqrt{s}=8 TeV, Preliminary");
  text->SetNDC();
  text->SetTextSize(0.04);
  text->Draw();
}

void drawCMSsimulation() {
  TLatex* text=new TLatex(0.18, 0.95, "CMS Simulation, 2014, #sqrt{s}=8 TeV, Preliminary");
  text->SetNDC();
  text->SetTextSize(0.04);
  text->Draw();
}

void drawCMSprivate() {
  TLatex* text=new TLatex(0.18, 0.95, "CMS Simulation, 2014, #sqrt{s}=8 TeV, Private");
  text->SetNDC();
  text->SetTextSize(0.04);
  text->Draw();
}

void drawPrivate(float textsize = 0.035) {
  stringstream header;
  if(Privat) header << "Private simulation 2015, #sqrt{s}=" << CMenergy << " TeV, " << Totsamples << " samples, " << totNclu << " tot clusters, " << testoption ;
  else     header << "Simulation 2015, #sqrt{s}=" << CMenergy << " TeV, " << Totsamples << " samples, " << totNclu << " tot clusters, " << testoption ;
  TLatex* text=new TLatex(0.18, 0.95, header.str().c_str());
  text->SetNDC();
  text->SetTextSize(textsize);
  text->Draw(); 
}

void drawPU(int pu = 40, bool right = true) {
  double xmin = 0.77;
  if (!right) xmin = 0.18;
  TLatex* text=new TLatex(xmin, 0.87, Form("<PU> ~ %d",pu));
  text->SetNDC();
  text->SetTextSize(0.04);
  text->Draw();
}

void draw(std::vector<TH1F*> h,
	  TString name, TString xTitle,
	  double xmin, double xmax, double ymin, double ymax,
	  TString legHeader = "", bool legRIGHT = true, bool legTOP = true,
	  bool logX = false, bool logY = false, bool stat = false,
	  double scale = -9., int rebin = -1, int orbin = -1,
	  TString option = "", int nclus = 99) {  //double ymin_ratio, double ymax_ratio,

  TCanvas* can = new TCanvas(name,name,900,900);
  can->cd();
  if (logX) gPad->SetLogx();
  if (logY) gPad->SetLogy();

  double legymax, legymin, legxmin, legxmax;
  if(legHeader == "BENCHMARKS COMPARISON") {
    legxmin = (legRIGHT ? 0.58 : 0.18);
    legxmax = legxmin+0.20;
    legymax = (legTOP ? 0.90 : 0.15);
    legymin = legymax-0.35;    
  }
  TLegend* leg = new TLegend(legxmin,legymin,legxmax,legymax);
  leg->SetTextSize(0.024);
  leg->SetFillColor(0);
  leg->SetLineColor(0);

  TString options = (option=="" ? "pe" : option);
  
  //normalize and set y range
  ymax=0.;
  h[0]->Sumw2();  
  double norm = 10000.;     
  for (size_t i=0; i<h.size(); i++) {  
    if(h[i]->GetNbinsX() != orbin) cout << "WARNING: orbin for " << h[i]->GetName() << " are " << h[i]->GetNbinsX() << endl; //debug - shift of h[][] wrt clu[][]
    if (rebin>0) h[i]->Rebin(rebin);
    scale = norm/(h[i]->Integral());
    h[i]->Scale(scale);
    if (h[i]->GetMaximum() > ymax) ymax = h[i]->GetMaximum();
  }
  ymax = ymax*1.2;

  for (size_t i=0; i<h.size(); i++) {
    h[i]->GetXaxis()->SetRangeUser(xmin,xmax);
    h[i]->GetXaxis()->SetTitle(xTitle);
    h[i]->SetMinimum(ymin);
    h[i]->SetMaximum(ymax);
    h[i]->GetYaxis()->SetLabelSize(0.03);
    h[i]->GetYaxis()->SetTitleOffset(1);
    h[i]->GetYaxis()->SetTitleSize(0.04);    
    string nam = "";
    if(legHeader == "BENCHMARKS COMPARISON"){
      if(i==(h.size()-1)){
        nam = translate(clu[nclus][0].c_str());
        leg->AddEntry(h[i],nam.c_str(),"l"); //to print all for bench comp
      }
      else {
        nam = translate(clu[nclus][i+1].c_str());
        leg->AddEntry(h[i],nam.c_str(),"l"); //to print all for bench comp    
      }
    }
    else if(i==(h.size()-1)) {
      h[i]->SetLineColor(kRed); //to print benchmark always red.
      nam = translate(clu[nclus][0].c_str());
      leg->AddEntry(h[i],nam.c_str(),"l"); //to print only benchmark
    }
    if (i==1) options = options + (stat ? "sames" : "same"); //once is enought
    h[i]->Draw(options);
  }  
  leg->Draw("same");
  drawPrivate(0.025); //Ncl..
  can->Update();
  can->SaveAs(Outfolder+name+".png");
}

void draw_ratio(std::vector<TH1F*> h,
	  TString name, TString xTitle,
	  double xmin, double xmax,
	  TString legHeader = "", bool legRIGHT = true, bool legTOP = true,
	  bool logX = false, bool stat = false, int rebin = -1, int orbin = -1,
	  TString option = "", int nclus = 99) {  //double ymin_ratio, double ymax_ratio,
  TCanvas* can = new TCanvas(name+"_ratio",name+"_ratio",900,450);
  can->cd();

  double legxmin = (legRIGHT ? 0.55 : 0.18);
  double legxmax = legxmin+0.25;
  double legymin = (legTOP ? 0.70 : 0.15);
  double legymax = legymin+0.15;
  TLegend* leg = new TLegend(legxmin,legymin,legxmax,legymax);
  if (legHeader!="") leg->SetHeader(legHeader);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);
  leg->SetLineColor(0);

  TString options = (option=="" ? "pe" : option);

  if (rebin>0) h[h.size()-1]->Rebin(rebin); //to rebin benchmark before divide
  for (size_t i=0; i<h.size(); i++) {
    //if(h[i]->GetNbinsX() != orbin) cout << "WARNING: orbin for " << h[i]->GetName() << " are " << h[i]->GetNbinsX() << endl; //debug - shift of h[][] wrt clu[][]
    if (rebin>0 && i<(h.size()-1)) h[i]->Rebin(rebin);
    TH1F* ratio = (TH1F*)h[i]->Clone("ratio_"+name);
    ratio->Sumw2();
    ratio->Divide(h[h.size()-1]); //benchmark is at the end.

    if (logX) gPad->SetLogx();
    ratio->SetMarkerStyle(20+i);
    ratio->SetMarkerSize(1.0); //1.2
    ratio->GetXaxis()->SetRangeUser(xmin,xmax);
    ratio->SetMinimum(-0.1);
    ratio->SetMaximum(4);
    if (i==0){    //just for the first one
     ratio->GetXaxis()->SetLabelSize(0.05);
     ratio->GetXaxis()->SetTitle(xTitle);
     ratio->GetXaxis()->SetTitleOffset(1);
     ratio->GetXaxis()->SetTitleSize(0.06);
     ratio->GetYaxis()->SetTitle("ratio");
     ratio->GetYaxis()->SetTitleSize(0.06);
     ratio->GetYaxis()->SetTitleOffset(0.7);
     ratio->GetYaxis()->SetLabelSize(0.05);
    }
    string nam = "";
    nam = translate(clu[nclus][0].c_str());
    if(i==(h.size()-1)) leg->AddEntry(h[i],nam.c_str(),"l"); //to print only benchmark (first in the list)
    if (i==1) options = options + (stat ? "sames" : "same"); //once is enought
    ratio->Draw(options);
  }
  leg->Draw("same");
  drawPrivate(0.04);
  can->Update();
  can->SaveAs(Outfolder+name+"_ratio.png");
}

void performancePlot1D(bool ratio, int nclust, TString hName,
		       double xmin, double xmax, double ymin, double ymax,
		       TString xaxis, TString yaxis,
		       bool logX = false, bool logY = false, bool stat = false,
		       double scale = -9., int rebin = -1, int orbin = -1, TString option = "")
{

  gROOT ->Reset();
  //=========  settings ====================
  gROOT->SetStyle("Plain");
 
  //if (!initialized) init();
  //=============================================
  if(nclust < 0  || nclust >= (int)clu.size()){
    std::cout << "ERROR: bad ncluster value(" << nclust+1 << "), ncluster between 1 " << "and " << clu.size() << std::endl;
    return;
  }
  std::vector<TH1F*>  h;   
  TFile* f = new TFile(filename);
  if (f == NULL) return;
  
  //debug  
  //cout << (sizeof(L)/sizeof(std::string)) << " " << (sizeof(Y)/sizeof(std::string))<< " "  << (sizeof(C2)/sizeof(std::string)) << endl;

  //create samples name --
  vector<string> samplename;
  for (unsigned int l = 0; l<(sizeof(L)/sizeof(std::string)); ++l){
    for (unsigned int y = 0; y<(sizeof(Y)/sizeof(std::string)); ++y){
      if(sizeof(C2)>0){  //for 3par       
        for (unsigned int c = 0; c<(sizeof(C2)/sizeof(std::string)); ++c){
          string ss = "L" + L[l] + "y" + Y[y] + C2[c];
          samplename.push_back(ss);
        }
      }
      else { //for 2par     
        string ss = "L" + L[l] + "y" + Y[y];
        samplename.push_back(ss);
      }
    }
  }
  //-----------

  int nc = nclust; //ehi!
  int size = clu[nc].size();
  std::cout << "# Start reading distros from file" << std::endl;
  for(int nsam=1; nsam<size; nsam++) { //on samples - skip benchmark.   

    TH1F* histo = NULL;
    string sample = clu[nc][nsam];
    sample = sample + "_" + hName;
    TString fname = sample;
    //std::cout << " Getting " << fname << std::endl;
    histo = (TH1F*)f->Get(fname); 
    if(nsam<maxCol){ histo->SetLineColor(colors[nsam]);
                     histo->SetMarkerColor(colors[nsam]);
    }
    else { histo->SetLineColor(kGray+1);
           histo->SetMarkerColor(kGray+1);
    }
    histo->SetMarkerSize(1.0);
    histo->SetMarkerStyle(20);
    histo->SetLineWidth(2);
    histo->GetXaxis()->SetTitle(xaxis);
    histo->GetYaxis()->SetTitle(yaxis);
    h.push_back(histo);
  }
  //to append benchmark as last histo. 
  TH1F* histo = NULL;
  string sample = clu[nc][0];
  sample = sample + "_" + hName;
  TString fname = sample;
  std::cout << " Getting the benchmark: " << fname << std::endl;
  histo = (TH1F*)f->Get(fname);        
  histo->SetMarkerSize(1.0);
  histo->SetMarkerStyle(20);
  if((unsigned int)(nc+1) == clu.size()) {
    histo->SetLineWidth(2);
    histo->SetLineColor(colors[0]); //to print benchmark always red.
  }
  else {
    histo->SetLineWidth(3);
    histo->SetLineColor(kRed); //to print benchmark always red.  
  }
  histo->GetXaxis()->SetTitle(xaxis);
  histo->GetYaxis()->SetTitle(yaxis);
  h.push_back(histo);

  cout << " Cluster size: " << size << " sample" << endl << endl;
  stringstream sst, sc;
  sc << "CLUSTER " << (nc+1) << " - " << "#samples: " << size;
  TString clusterLabel = sc.str();
  if ((unsigned int)(nc+1) == clu.size()) {
    clusterLabel = "BENCHMARKS COMPARISON";
    sst << "BC";   //benchmark comp
  }
  else sst << nc+1;   

  if(h.size()>0){
      if(!ratio){ draw(h, Outname+testoption+"_Clus"+sst.str()+"_"+hName,
           xaxis,xmin,xmax,ymin,ymax,
  	   clusterLabel,lgRIGHT,lgTOP, logX, logY, stat,
 	   scale,rebin,orbin,option,nc);  //,ymin_ratio,ymax_ratio
      }
      else {
         if((unsigned int)(nc+1) != clu.size()){   //to skip ratio of benchmark comparison
          draw_ratio(h, Outname+testoption+"_Clus"+sst.str()+"_"+hName,
             xaxis,xmin,xmax,clusterLabel,lgRIGHT,lgTOP, logX, stat, rebin,orbin,option,nc);
        }
      }
  }  
  else cout << "WARNING: empty cluster!" << endl;

}
 
//plot call for different variables:
void plot_pt(bool rat, int ncluster = 99, int rebin = 1, TString opt="") {
  performancePlot1D(rat,ncluster-1,"pt",0.,450.,0.,800.,"higgs p_{T} [GeV]","a.u.", false, false, false,-9.,rebin,100,opt);
}

void plot_pzh(bool rat, int ncluster = 99, int rebin = 4, TString opt="") {
  performancePlot1D(rat,ncluster-1,"pzh",0.,1000.,0.,800.,"higgs p_{z higher} [GeV]","a.u.", false, false, false,-9.,rebin,500,opt);
}

void plot_pzl(bool rat, int ncluster = 99, int rebin = 4, TString opt="") {
  performancePlot1D(rat,ncluster-1,"pzl",-300.,500.,0.,800.,"higgs p_{z lower} [GeV]","a.u.", false, false, false,-9.,rebin,500,opt);
}

void plot_mhh(bool rat, int ncluster = 99, int rebin = 2, TString opt="") {
  performancePlot1D(rat,ncluster-1,"mhh",240.,900.,0.,800.,"di-higgs mass [GeV]","a.u.", false, false, false,-9.,rebin,200,opt);
}

void plot_hth(bool rat, int ncluster = 99, int rebin = 2, TString opt="") {
  performancePlot1D(rat,ncluster-1,"hth",0.,3.2,0.,800.,"higgs #theta","a.u.", false, false, false,-9.,rebin,200,opt);
}

void plot_hcth(bool rat, int ncluster = 99, int rebin = 2, TString opt="") {
  performancePlot1D(rat,ncluster-1,"hcth",-1.,1.,0.,800.,"higgs cos#theta","a.u.", false, false, false,-9.,rebin,200,opt);
}

void plot_hths(bool rat, int ncluster = 99, int rebin = 2, TString opt="") {
  performancePlot1D(rat,ncluster-1,"hths",0.,3.2,0.,800.,"higgs #theta*","a.u.", false, false, false,-9.,rebin,200,opt);
}

void plot_hcths(bool rat, int ncluster = 99, int rebin = 2, TString opt="") {
  performancePlot1D(rat,ncluster-1,"hcths",0.,1.,0.,800.,"higgs #cbar cos#theta*#cbar","a.u.", false, false, false,-9.,rebin,100,opt);
}

//genaral function to plot all or single variable, single cluster.
void plot(int totclu = 20, bool r = false, int var = 0, int nclu = 0, int reb = 99, TString opt="hist") { 
// r: 0-distros, 1-ratio ; var: 1-pt, 2-pzh, 3-pzl, 4-mhh ;  nclu = 0 to do all the clusters

  bool doall = false;
  if(totclu == 0) doall = true;
  string app;

  for(int i=2; i<=Maxtotclu; i++) { //mintotclu .. continued!

    if(!doall) totNclu = totclu;
    else   totNclu = i;

    std::stringstream ss;
    ss << pars << "par_" << CMenergy << "Tev_Nclu" << totNclu;
    Outname = ss.str();
    
    //default rebin
    int reb1 = (reb == 99) ? 1 : reb;
    int reb2 = (reb == 99) ? 4 : reb;
    int reb3 = (reb == 99) ? 4 : reb;
    int reb4 = (reb == 99) ? 2 : reb;
    int reb5 = (reb == 99) ? 2 : reb;
    int reb6 = (reb == 99) ? 2 : reb;
    int reb7 = (reb == 99) ? 2 : reb;
    int reb8 = (reb == 99) ? 2 : reb;
     
    if(!init()) return;
     
    cout << clu.size() << endl; //DEBUG
    if(nclu == 0){
      int size = clu.size();
      for(int nc=1; nc<(size+1); nc++) { //debug  -- on clusters
        if(var == 0) { plot_pt(r,nc, reb1, opt); 
                     plot_mhh(r,nc, reb4, opt); plot_hcths(r,nc, reb8, opt);} // plot_pzh(r,nc, reb2, opt); plot_pzl(r,nc, reb3, opt); plot_hths(r,nc, reb7, opt); plot_hcths(r,nc, reb8, opt);
        else if(var == 1) plot_pt(r,nc, reb1, opt);
        else if(var == 2) plot_pzh(r,nc, reb2, opt);
        else if(var == 3) plot_pzl(r,nc, reb3, opt);
        else if(var == 4) plot_mhh(r,nc, reb4, opt);
        else if(var == 5) plot_hth(r,nc, reb5, opt);
        else if(var == 6) plot_hcth(r,nc, reb6, opt);
        else if(var == 7) plot_hths(r,nc, reb7, opt);
        else if(var == 8) plot_hcths(r,nc, reb8, opt);
        else cout << "wrong var value (0-4)" << endl;
      }
    }
    else {
        if(var == 0) { plot_pt(r,nclu, reb1, opt); 
		     plot_mhh(r,nclu, reb4, opt); plot_hcths(r,nclu, reb8, opt); } // plot_pzh(r,nclu, reb2, opt); plot_pzl(r,nclu, reb3, opt); plot_hths(r,nclu, reb7, opt);
        else if(var == 1) plot_pt(r,nclu, reb1, opt);
        else if(var == 2) plot_pzh(r,nclu, reb2, opt);
        else if(var == 3) plot_pzl(r,nclu, reb3, opt);
        else if(var == 4) plot_mhh(r,nclu, reb4, opt);
        else if(var == 5) plot_hth(r,nclu, reb5, opt);
        else if(var == 6) plot_hcth(r,nclu, reb6, opt);
        else if(var == 7) plot_hths(r,nclu, reb7, opt);
        else if(var == 8) plot_hcths(r,nclu, reb8, opt);
        else cout << "wrong var value (0-4)" << endl;
    }

    clu.clear(); //needed!
    if(!doall) return;
  }//end for  
}

