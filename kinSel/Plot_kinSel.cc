// Code for non resonant HH->bbbb in CMS Run2 
// Make plots from Step1 -> KINEMATIC SELECTION
//  Author: Martino Dall'Osso
//    Oct 14th, 2015
//----------------------------------------------------------------

#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TKey.h>

// matrix parameters  --> needed only for the header..
static const int binPt = 10;
static const int binCSV = 8;
static const int binEta = 4;

static const int colors[] = {kBlack, kOrange+2, kBlue, kCyan+2, kGreen-1, kRed, kPink,kMagenta,kViolet+2,
                kSpring+8,   kYellow+2, kMagenta-6, kGreen-6, kOrange-6, kBlue-6, kCyan-6, kGreen-6, kPink-6, kRed-6,
                kMagenta-8, kOrange-8, kBlue-8, kCyan-8, kGreen-8, kPink-8, kRed-8,
                kMagenta-4, kViolet-8, kAzure-8, kSpring-8, kYellow-8, kGreen-8,
                kViolet-4, kAzure-4, kSpring-4, kYellow-4, kGreen-4};
const int maxCol = 36;

static const std::string plotsFld="plots/";
static const TString OutFld="plots/"; //debug!!
TString plotFld = "";
#include "../utils/hh4bStructs.h"
#include "hh4b_13TeV_kinSel.h" //Histos are here

TH1F *h1 = new TH1F("h1", "", 100, 0., 1000.);


//just to have order
//--------------------
class Plot_kinSel{

  private:
    std::vector<Jet4Plot>  h1v;

  public:
    Plot_kinSel(std::string ,bool , TString,  std::string , int =0);
    ~Plot_kinSel();
    
//    void setBranches(bool, TTree*);   
//    bool readMatrix( std::string, std::string );
//    double selectBestDiJets (int );
//    void fillHistos(bool );   
    void drawCMSprel(TString);
    //void getHistos(TFile * ,bool );
    //void getH1(TFile *, TH1F* );
    void drawH1(TString, std::vector<Jet4Plot> , TString , TString ="", int = -99, double = -99., 
          bool = true, bool = true,
          TString ="", double =-99., double =-99., double =-99., double =-99., TString = "", bool = false,
	  int = -99, TString = "", bool = false, bool = false);
    void saveHistos(std::string, bool );   
//    string jet4SelectionMethod(int );   
};


Plot_kinSel::Plot_kinSel(std::string sample, bool isData, TString dataset, std::string opt, int whatdisplay)
{
  //drawing style  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);

  //read File with Histos
  std::string inputfilename;  
  inputfilename = plotsFld+"Histograms_"+sample+"_"+opt+".root"; 
  TFile *f=new TFile(inputfilename.c_str());

  plotFld = sample+"_"+opt+"/";

  Jet4Plot jn;
  if (whatdisplay == 0 || whatdisplay == 1){ //Jets variable plots

    //4jets pt comparison
//    h1=(TH1F*)f->Get("h_naJets");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
//    h1=(TH1F*)f->Get("h_nfJets");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
//    drawH1(dataset, h1v, "nJets_comp", "# jets", 1, 1000);      h1v.clear();

    h1=(TH1F*)f->Get("h_nJetsAll");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_nJetsAcc");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_nfJets");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_naJets");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "nJets_comp", "# jets", 0, 10000);      h1v.clear();

    h1=(TH1F*)f->Get("h_JetsAll_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_JetsAcc_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJets_mass");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_aJets_mass");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "Jets_mass_comp", "jet m (GeV/c^{2})", 2, 10000);      h1v.clear();

    h1=(TH1F*)f->Get("h_JetsAll_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_JetsAcc_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJets_pT");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_aJets_pT");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "Jets_pT_comp", "jet p_{T} (GeV/c)", 2, 10000);      h1v.clear();

    h1=(TH1F*)f->Get("h_JetsAll_eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_JetsAcc_eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJets_eta");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_aJets_eta");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "Jets_eta_comp", "jet #eta", 2, 10000);      h1v.clear();

    h1=(TH1F*)f->Get("h_JetsAll_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_JetsAcc_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJets_CSV");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_aJets_CSV");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "Jets_CSV_comp", "jet CSV", 4, 10000, 0);      h1v.clear();
  }
  else if (whatdisplay == 0 || whatdisplay == 2){ //Jets comparison

    //4 final jets
    h1=(TH1F*)f->Get("h_fJet1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet3_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet4_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_mass_comp", "jet m (GeV/c^{2})", 2);      h1v.clear();

    h1=(TH1F*)f->Get("h_fJet1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet3_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet4_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_pT_comp", "jet p_{T} (GeV/c)", 2);      h1v.clear();

    h1=(TH1F*)f->Get("h_fJet1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet3_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet4_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_Eta_comp", "jet #eta", 2);      h1v.clear();

    h1=(TH1F*)f->Get("h_fJet1_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet2_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet3_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJet4_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_CSV_comp", "jet CSV", 4, 0, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_fJet3avg_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets3avg_CSV", "avg 3jets CSV", 4, 0, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_fJet3min_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets3min_CSV", "min 3jet CSV", 4, 0, 0);      h1v.clear();

   // h1=(TH1F*)f->Get("h_fJets_Centr");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
//    drawH1(dataset, h1v, "fJets_Centr", "jet C", 4, 0, 0);      h1v.clear();

  }
  else if (whatdisplay == 0 || whatdisplay == 3){  //diJets comparison

    //di-Jets
    h1=(TH1F*)f->Get("h_H1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diJets_mass_comp", "di-jet m (GeV/c^{2})", 4 ,10000);      h1v.clear();

    h1=(TH1F*)f->Get("h_H1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diJets_pT_comp", "di-jet p_{T} (GeV/c)", 2, 10000);      h1v.clear();

    h1=(TH1F*)f->Get("h_H1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diJets_eta_comp", "di-jet #eta", 2, 10000);      h1v.clear();

    h1=(TH1F*)f->Get("h_H1_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H2_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diJets_phi_comp", "di-jet #phi", 2, 10000, 1, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_H1_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H2_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diJets_CosThSt_comp", "di-jet |cos#theta*|", 2, 10000);      h1v.clear();

    h1=(TH1F*)f->Get("h_H1_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H2_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diJets_DR_comp", "di-jet #Delta R", 2, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_H1_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H2_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diJets_DPhi_comp", "di-jet #Delta #phi", 2, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_H1_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_H2_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diJets_DEta_comp", "di-jet #Delta #eta", 2, 0);      h1v.clear();
  }
  else if (whatdisplay == 4){ //diHiggs system

    h1=(TH1F*)f->Get("h_HH_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diHiggs_mass", "", 2, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_HH_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diHiggs_pT", "", 2, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_HH_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diHiggs_Eta", "", 2, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_HH_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diHiggs_Phi", "", 2, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_HH_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diHiggs_DR", "", 2, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_HH_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diHiggs_DPhi", "", 2, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_HH_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "diHiggs_DEta", "", 2, 0);      h1v.clear();

//    h1=(TH1F*)f->Get("h_HH_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
//h_H1_H2_mass
//    drawH1(dataset, h1v, "diHiggs_mass", "", 2, 10000);      h1v.clear();

  }
  else if (!isData && (whatdisplay == 0 || whatdisplay == 5)){ 



  }
}

Plot_kinSel::~Plot_kinSel(){
}
//---------------

void Plot_kinSel::drawCMSprel(TString dataset) {
  TLatex* text=new TLatex(0.18, 0.95, "CMS preliminary 2015, #sqrt{s}=13 TeV, "+dataset+", non-res hh->4b");
  text->SetNDC();
  text->SetTextSize(0.02);
  text->Draw();
}
//---------------

//from clustering...
void Plot_kinSel::drawH1(TString dataset, std::vector<Jet4Plot> hn, TString name, TString xTitle, int rebin, double norm,
          bool legRIGHT , bool legTOP , TString yTitle, double xmin, double xmax, double ymin, double ymax,
	  TString legHeader ,  bool stat , int orbin , TString option , 
          bool logX , bool logY ) {

  TCanvas* c = new TCanvas(name,name,900,900);
  c->cd();
  if (logX) gPad->SetLogx();
  if (logY) gPad->SetLogy();

  double legymax, legymin, legxmin, legxmax;
  legxmin = (legRIGHT ? 0.70 : 0.25);
  legxmax = legxmin+0.15;
  legymax = (legTOP ? 0.85 : 0.55);
  legymin = legymax-0.15;     
  TLegend* leg = new TLegend(legxmin,legymin,legxmax,legymax);
  if (legHeader!="") leg->SetHeader(legHeader);
  leg->SetTextSize(0.025);

  TString options = (option=="" ? "HIST" : option);
  
  //normalize and set y range
  ymax=0.;
  hn[0].h->Sumw2();  
  //cout << h.size() << endl;
  for (size_t i=0; i<hn.size(); i++) {  
  //cout << i << endl;
    //if(h[i]->GetNbinsX() != orbin) cout << "WARNING: orbin for " << h[i]->GetName() << " are " << h[i]->GetNbinsX() << endl; //debug - shift of h[][] wrt clu[][]
    if(rebin >0) hn[i].h->Rebin(rebin);
    if(norm > 0){ 
      double scale = norm/(hn[i].h->Integral());
      hn[i].h->Scale(scale);
    }
    else if(hn[i].norm >0) {  //scale histos if required..
      double scale = hn[i-1].h->Integral()/hn[1].h->Integral();
      hn[i].h->Scale(scale);     
    } 
    if (hn[i].h->GetMaximum() > ymax) ymax = hn[i].h->GetMaximum();
  }
  ymax = ymax*1.2;

  for (size_t i=0; i<hn.size(); i++) {
    hn[i].h->SetLineColor(colors[i]);
    hn[i].h->SetMarkerColor(colors[i]);
    if(xmin>0 && xmax>0)hn[i].h->GetXaxis()->SetRangeUser(xmin,xmax);
    if(xTitle!="")hn[i].h->GetXaxis()->SetTitle(xTitle);
    if(ymin>0)hn[i].h->SetMinimum(ymin);
    if(ymax>0)hn[i].h->SetMaximum(ymax);
    if(yTitle!="")hn[i].h->GetYaxis()->SetTitle(yTitle);
    else hn[i].h->GetYaxis()->SetTitle("a.u.");
    hn[i].h->GetXaxis()->SetLabelSize(0.03);
    hn[i].h->GetXaxis()->SetTitleOffset(1);
    hn[i].h->GetXaxis()->SetTitleSize(0.04);    
    hn[i].h->GetYaxis()->SetLabelSize(0.03);
    hn[i].h->GetYaxis()->SetTitleOffset(1);
    hn[i].h->GetYaxis()->SetTitleSize(0.04);    
    string nam = "";
    leg->AddEntry(hn[i].h); //to print all for bench comp      
    if (i==1) options = options + (stat ? "sames" : "same"); //once is enought
    hn[i].h->Draw(options);
  }  
  leg->Draw("same");
  drawCMSprel(dataset);
  c->Update();
  c->SaveAs(OutFld+plotFld+name+".png");
  c->SaveAs(OutFld+plotFld+name+".pdf");
}
//---------------

