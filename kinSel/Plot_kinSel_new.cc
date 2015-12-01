// Code for non resonant HH->bbbb in CMS Run2 
// Make plots from Step1 -> KINEMATIC SELECTION
//  Author: Martino Dall'Osso
//    Oct 14th, 2015
//     To get plot of set of variables from Histograms produced with Step1 ntuples.
//      .L Plot_kinSel_new.cc++
//       Plot_kinSel("BTagCSVRun2015C-D","BTagCSVRun2015C-D_4bTag", "pT20CSVL", "", 2)
//        Option: Sample1, Sample2, optionSample1, optionSample2, VariablesToBePlotted.
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
#include <TPaletteAxis.h>
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

static const std::string frameworkVersionFld="V15/";  //WARNING: CHANGE ME!! -- equal to data version
bool saveCanvas = true;
//debug -- arranged for b-Tag ratio...
//std::string subFld = "BTagCSVRun2015C-D_Trg2b_pT20CSVL_2-4CSV/"; //debug
std::string subFld = "BTagCSVRun2015C-D_Trg2b_pT30CSVL_2-4CSV/"; //debug
//std::string subFld = "BTagCSVRun2015C-D_Trg2b_pT20CSVM_2-4CSV/"; //debug
//std::string subFld = "BTagCSVRun2015C-D_Trg3b_pT20CSVL_3-4CSV/"; //debug
//std::string subFld = "BTagCSVRun2015C-D_Trg3b_pT30CSVL_3-4CSV/"; //debug
//std::string subFld = "BTagCSVRun2015C-D_Trg3b_pT20CSVM_3-4CSV/"; //debug

static const std::string plotsFld="plots/"+frameworkVersionFld;
static const TString OutFld="plots/"+frameworkVersionFld; //debug!!
TString plotSubFld = "";
#include "../utils/hh4bStructs.h"
#include "hh4b_13TeV_kinSel.h" //Histos are here

TH1F *h1 = NULL; //dummy for the get.
TH2F *h2a = NULL;
TH2F *h2b = NULL;
TH2F *h2ra = NULL;
TH2F *h2rb = NULL;

//just to order
//--------------------
class Plot_kinSel{

  private:
    std::vector<Jet4Plot>  h1v;
    bool twoSamples = false;
    TString Legend[2];

  public:
    Plot_kinSel(std::string ,std::string ,std::string ,  std::string , int =0);
    ~Plot_kinSel();
    
    void drawCMSprel(TString, float);
    void drawH2(TString , TH2F* , TString , TString ="", 
          TString ="", double= -99 , double= -99 , double = -99, double = -99,
	  bool  = false, int = -99 , TString ="" , bool = false , bool = false );
    void drawH2ratio(TString , TH2F* , TH2F* , TString , TString ="", 
          TString ="", double= -99 , double= -99 , double = -99, double = -99,
	  bool  = false, int = -99 , TString ="" , bool = false , bool = false );
    void drawH1(TString, std::vector<Jet4Plot> , TString , TString ="", int = -99, double = -99., 
          bool = true, bool = true,
          TString ="", double =-99., double =-99., double =-99., double =-99., TString = "", bool = false,
	  int = -99, TString = "", bool = false, bool = false);
    void saveHistos(std::string, bool );   
};


Plot_kinSel::Plot_kinSel(std::string sample1, std::string sample2, std::string opt1, std::string opt2, int whatdisplay)
{
  TString dataset = sample1; //debug

  //drawing style  
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  std::string inputfilename, sample;  
  TFile *f, *f1, *f2;

  //read File with Histos
  if(sample1 == "" || sample2 == ""){
    sample = (sample1 == "") ? sample2 : sample1;
    inputfilename = plotsFld+"Histograms_"+sample+"_"+opt1+".root"; 
    f = TFile::Open(inputfilename.c_str());
    if(!f) return;
    //plotSubFld = sample+"_"+opt1+"/";
  }
  else if(sample1 == sample2){  //1 sample - different option
    inputfilename = plotsFld+"Histograms_"+sample1+"_"+opt1+".root"; 
    f1 = TFile::Open(inputfilename.c_str());
    if(!f1) return;
    Legend[0]=opt1;
    inputfilename = plotsFld+"Histograms_"+sample1+"_"+opt2+".root"; 
    f2 = TFile::Open(inputfilename.c_str());
    if(!f2) return;
    Legend[1]=opt2;
    plotSubFld = subFld;//"comp_"+sample1+"_"+opt1+"-"+opt2"/" //debug
    twoSamples = true;
  }
  else {
    inputfilename = plotsFld+"Histograms_"+sample1+"_"+opt1+".root"; 
    f1 = TFile::Open(inputfilename.c_str());
    if(!f1) return;
    Legend[0]=sample1;
    if(opt2=="") opt2 = opt1;
    inputfilename = plotsFld+"Histograms_"+sample2+"_"+opt2+".root"; 
    f2 = TFile::Open(inputfilename.c_str());
    if(!f2) return;
    Legend[1]=sample2;
//    plotSubFld = "comp_"+sample1+"_"+opt1+"-"+opt2"/";
    plotSubFld = "";//"comp"+sample2+"_"+opt+"/" - debug
    twoSamples = true;
  }

  Jet4Plot jn;
  if (1 == whatdisplay || whatdisplay == 0){ //Jets variable plots
   if(!twoSamples){
    h1=(TH1F*)f->Get("h_nJetsAll");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_nJetsAcc");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_nfJets");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_naJets");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "nJets_comp", "# jets", 0, 1.);      h1v.clear();

    h1=(TH1F*)f->Get("h_JetsAll_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_JetsAcc_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJets_mass");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_aJets_mass");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "Jets_mass_comp", "jet m (GeV/c^{2})", 2, 1.);      h1v.clear();

    h1=(TH1F*)f->Get("h_JetsAll_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_JetsAcc_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJets_pT");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_aJets_pT");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "Jets_pT_comp", "jet p_{T} (GeV/c)", 2, 1.);      h1v.clear();

    h1=(TH1F*)f->Get("h_JetsAll_eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_JetsAcc_eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJets_eta");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_aJets_eta");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "Jets_eta_comp", "jet #eta", 2, 1.);      h1v.clear();

    h1=(TH1F*)f->Get("h_JetsAll_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_JetsAcc_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_fJets_CSV");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f->Get("h_aJets_CSV");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "Jets_CSV_comp", "jet CSV", 4, 1., 0);      h1v.clear();

    //h1=(TH1F*)f->Get("h_fJets_Centr");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    //drawH1(dataset, h1v, "fJets_Centr", "jet C", 4, 0, 0);      h1v.clear();

    //h1=(TH1F*)f->Get("h_MET");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    //drawH1(dataset, h1v, "met", "met (GeV)", 4, 0 );      h1v.clear();
   }
   else {
    h1=(TH1F*)f1->Get("h_nJetsAcc");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_nJetsAcc");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "nJetsAcc_comp", "# jets in Acc", 0, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_nfJets");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_nfJets");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "nfJets_comp", "# final jets", 0, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJets_mass");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJets_mass");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_mass_comp", "final jets m (GeV/c^{2})", 2, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJets_pT");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJets_pT");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_pT_comp", "final jets p_{T} (GeV/c)", 2, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJets_eta");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJets_eta");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_eta_comp", "final jets #eta", 2, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJets_CSV");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJets_CSV");      jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_CSV_comp", "final jets CSV", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJets_Centr");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJets_Centr");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets_Centr", "final jets C", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_MET");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_MET");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "met", "met (GeV)", 4, 1.);      h1v.clear();
   }
  }
  else if (2 == whatdisplay || whatdisplay == 0){ //Jets comparison
   if(!twoSamples){
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

    h1=(TH1F*)f->Get("h_fJets3avg_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets3avg_CSV", "avg 3jets CSV", 4, 0, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_fJets3min_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets3min_CSV", "min 3jet CSV", 4, 0, 0);      h1v.clear();

    h1=(TH1F*)f->Get("h_fJets4avg_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets4avg_CSV", "avg 4jets CSV", 4, 0, 0);      h1v.clear();
   }
   else{
    h1=(TH1F*)f1->Get("h_fJet1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_mass_comp", "final jet1 m (GeV/c^{2})", 4, 1.);      h1v.clear();  //re-scale..
    h1=(TH1F*)f1->Get("h_fJet2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_mass_comp", "final jet2 m (GeV/c^{2})", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_mass_comp", "final jet3 m (GeV/c^{2})", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_mass_comp", "final jet4 m (GeV/c^{2})", 4, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_pT_comp", "final jet1 p_{T} (GeV/c)", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_pT_comp", "final jet2 p_{T} (GeV/c)", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_pT_comp", "final jet3 p_{T} (GeV/c)", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_pT_comp", "final jet4 p_{T} (GeV/c)", 4, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_Eta_comp", "final jet1 #eta", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_Eta_comp", "final jet2 #eta", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_Eta_comp", "final jet3 #eta", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_Eta_comp", "final jet4 #eta", 4, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_CSV_comp", "final jet1 CSV", 4, 1., 0, 1, "", 0.5, 1.0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_CSV_comp", "final jet2 CSV", 4, 1., 0, 1, "", 0.5, 1.0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_CSV_comp", "final jet3 CSV", 4, 1., 0, 1, "", 0.5, 1.0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_CSV_comp", "final jet4 CSV", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJets3avg_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJets3avg_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets3avg_CSV", "final avg 3jets CSV", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJets3min_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJets3min_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets3min_CSV", "final min 3jet CSV", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJets4avg_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJets4avg_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJets4avg_CSV", "avg 4jets CSV", 4, 1., 0);      h1v.clear();
   }
  }
  else if (3 == whatdisplay || whatdisplay == 0){  //diJets comparison

    //di-Jets
    if(!twoSamples){
      h1=(TH1F*)f->Get("h_H1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_mass_comp", "di-jet 1 m (GeV/c^{2})", 4 ,1.);      h1v.clear();

      h1=(TH1F*)f->Get("h_H1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_pT_comp", "di-jet p_{T} (GeV/c)", 2, 1.);      h1v.clear();

      h1=(TH1F*)f->Get("h_H1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_eta_comp", "di-jet #eta", 2, 1.);      h1v.clear();

      h1=(TH1F*)f->Get("h_H1_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H2_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_phi_comp", "di-jet #phi", 2, 1., 1, 0);      h1v.clear();

      h1=(TH1F*)f->Get("h_H1_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H2_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f->Get("h_H_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_CosThSt_comp", "di-jet |cos#theta*|", 2, 1.);      h1v.clear();

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
    else { //debug
      h1=(TH1F*)f1->Get("h_H1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "H1_mass_comp", "H1 m (GeV/c^{2})", 4 ,1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_pT_comp", "H1 p_{T} (GeV/c)", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_eta_comp", "H1 #eta", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H1_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H1_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_phi_comp", "H1 #phi", 2, 1., 1, 0);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H1_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H1_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_CosThSt_comp", "H1 |cos#theta*|", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H1_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H1_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_DR_comp", "H1 #Delta R", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H1_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H1_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_DPhi_comp", "H1 #Delta #phi", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H1_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H1_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_DEta_comp", "H1 #Delta #eta", 2, 1.);      h1v.clear();
    }
  }
  else if (4 == whatdisplay || whatdisplay == 0){  //diJets comparison

    if(twoSamples){
      h1=(TH1F*)f1->Get("h_H2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "H2_mass_comp", "H2 m (GeV/c^{2})", 4 ,1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_pT_comp", "H2 p_{T} (GeV/c)", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_eta_comp", "H2 #eta", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H2_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H2_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_phi_comp", "H2 #phi", 2, 1., 1, 0);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H2_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H2_CosTheta*");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_CosThSt_comp", "H2 |cos#theta*|", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H2_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H2_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_DR_comp", "H2 #Delta R", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H2_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H2_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_DPhi_comp", "H2 #Delta #phi", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_H2_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_H2_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diJets_DEta_comp", "H2 #Delta #eta", 2, 1.);      h1v.clear();
    }
  }
  else if (5 == whatdisplay || whatdisplay == 0 ){ //diHiggs system

    if(!twoSamples){
      h1=(TH1F*)f->Get("h_HH_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diHiggs_mass", "", 2, 0);      h1v.clear();

      h1=(TH1F*)f->Get("h_HH_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diHiggs_pT", "", 2, 0);      h1v.clear();

      h1=(TH1F*)f->Get("h_HH_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diHiggs_Eta", "", 2, 0);      h1v.clear();

      h1=(TH1F*)f->Get("h_HH_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diHiggs_Phi", "", 2, 0, 1, 0);      h1v.clear();

      h1=(TH1F*)f->Get("h_HH_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diHiggs_DR", "", 2, 0);      h1v.clear();

      h1=(TH1F*)f->Get("h_HH_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diHiggs_DPhi", "", 2, 0);      h1v.clear();

      h1=(TH1F*)f->Get("h_HH_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "diHiggs_DEta", "", 2, 0);      h1v.clear();

      h2a=(TH2F*)f->Get("h_H1_H2_mass");
      drawH2(dataset, h2a, "diHiggs_H1H2_mass", "");
    }
    else {
      h1=(TH1F*)f1->Get("h_HH_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_HH_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "HH_mass", "", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_HH_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_HH_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "HH_pT", "", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_HH_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_HH_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "HH_Eta", "", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_HH_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_HH_Phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "HH_Phi", "", 2, 1., 1, 0);      h1v.clear();

      h1=(TH1F*)f1->Get("h_HH_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_HH_DR");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "HH_DR", "", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_HH_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_HH_DPhi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "HH_DPhi", "", 2, 1.);      h1v.clear();

      h1=(TH1F*)f1->Get("h_HH_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      h1=(TH1F*)f2->Get("h_HH_DEta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
      drawH1(dataset, h1v, "HH_DEta", "", 2, 1.);      h1v.clear();
    }
  }
  else if (6 == whatdisplay || whatdisplay == 0 ){ //efficiencies

    if(!twoSamples){
      h2a=(TH2F*)f->Get("h_H1_H2_mass");
      drawH2(dataset, h2a, "diHiggs_H1H2_mass_1", "");   
    }
    else{
      h2a=(TH2F*)f1->Get("h_H1_H2_mass");
      drawH2(dataset, h2a, "diHiggs_H1H2_mass_1", "");

      h2b=(TH2F*)f2->Get("h_H1_H2_mass");
      drawH2(dataset, h2b, "diHiggs_H1H2_mass_2", "");

      h2ra=(TH2F*)f1->Get("h_H1_H2_mass");
      h2rb=(TH2F*)f2->Get("h_H1_H2_mass");
      drawH2ratio(dataset, h2ra, h2rb, "diHiggs_H1H2_mass_ratio", "");
    }
  }
  else if (7 == whatdisplay || whatdisplay == 0 ){ //efficiencies
    if(twoSamples){
      h1=(TH1F*)f1->Get("h_Cuts");
      //h1=(TH1F*)f2->Get("h_Cuts");
      //drawEff(dataset, h2, "diHiggs_H1H2_mass_MC", "");
    }
  }
  else if (8 == whatdisplay || whatdisplay == 0){ //Jets comparison
   if(twoSamples){
    /*h1=(TH1F*)f1->Get("h_fJet1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_mass_comp", "final jet1 m (GeV/c^{2})", 4, 1.);      h1v.clear();  //re-scale..
    h1=(TH1F*)f1->Get("h_fJet2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_mass_comp", "final jet2 m (GeV/c^{2})", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_mass_comp", "final jet3 m (GeV/c^{2})", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_mass");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_mass_comp", "final jet4 m (GeV/c^{2})", 4, 1.);      h1v.clear();*/

    h1=(TH1F*)f1->Get("h_fJet1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_pT_comp", "final jet1 p_{T} (GeV/c)", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_pT_comp", "final jet2 p_{T} (GeV/c)", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_pT_comp", "final jet3 p_{T} (GeV/c)", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_pT");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_pT_comp", "final jet4 p_{T} (GeV/c)", 4, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_Eta_comp", "final jet1 #eta", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_Eta_comp", "final jet2 #eta", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_Eta_comp", "final jet3 #eta", 4, 1.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_Eta");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_Eta_comp", "final jet4 #eta", 4, 1.);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_CSV_comp", "final jet1 CSV", 4, 1., 0, 1, "", 0.5, 1.0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_CSV_comp", "final jet2 CSV", 4, 1., 0, 1, "", 0.5, 1.0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_CSV_comp", "final jet3 CSV", 4, 1., 0, 1, "", 0.5, 1.0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_CSV");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_CSV_comp", "final jet4 CSV", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_Et");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_Et");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_Et_comp", "final jet1 Et", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_Et");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_Et");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_Et_comp", "final jet2 Et", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_Et");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_Et");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_Et_comp", "final jet3 Et", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_Et");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_Et");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_Et_comp", "final jet4 Et", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_chMult");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_chMult");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_chMult_comp", "final jet1 chMult", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_chMult");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_chMult");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_chMult_comp", "final jet2 chMult", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_chMult");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_chMult");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_chMult_comp", "final jet3 chMult", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_chMult");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_chMult");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_chMult_comp", "final jet4 chMult", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_phi_comp", "final jet1 phi", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_phi_comp", "final jet2 phi", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_phi_comp", "final jet3 phi", 4, 1., 0);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_phi");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_phi_comp", "final jet4 phi", 4, 1., 0);      h1v.clear();

    h1=(TH1F*)f1->Get("h_fJet1_leadTrackPt");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet1_leadTrackPt");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet1_leadTrackPt_comp", "final jet1 leadTrackPt", 0, 1., 0, 1, "", 0., 200.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet2_leadTrackPt");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet2_leadTrackPt");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet2_leadTrackPt_comp", "final jet2 leadTrackPt", 0, 1., 0, 1, "", 0., 200.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet3_leadTrackPt");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet3_leadTrackPt");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet3_leadTrackPt_comp", "final jet3 leadTrackPt", 0, 1., 0, 1, "", 0., 200.);      h1v.clear();
    h1=(TH1F*)f1->Get("h_fJet4_leadTrackPt");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    h1=(TH1F*)f2->Get("h_fJet4_leadTrackPt");    jn.h=h1; jn.norm=0.; h1v.push_back(jn);
    drawH1(dataset, h1v, "fJet4_leadTrackPt_comp", "final jet4 leadTrackPt", 0, 1., 0, 1, "", 0., 200.);      h1v.clear();
  
   }
  }
  //else if (MCsample != "" && !twoSamples && (whatdisplay == 0 || whatdisplay == 5)){ 
  //}
}

Plot_kinSel::~Plot_kinSel(){
}
//---------------

void Plot_kinSel::drawCMSprel(TString dataset, float textsize) {
  TLatex* text=new TLatex(0.18, 0.93, "CMS preliminary 2015, #sqrt{s}=13 TeV, "+dataset+", non-res hh->4b");
  text->SetNDC();
  text->SetTextSize(textsize);
  text->Draw();
}
//---------------

void Plot_kinSel::drawH1(TString dataset, std::vector<Jet4Plot> hn, TString name, TString xTitle, int rebin, double norm,
          bool legRIGHT , bool legTOP , TString yTitle, double xmin, double xmax, double ymin, double ymax,
	  TString legHeader ,  bool stat , int orbin , TString option , 
          bool logX , bool logY ) {

  TCanvas* c = new TCanvas(name,name,600,600);
  c->cd();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current p

  if (logX) gPad->SetLogx();
  if (logY) gPad->SetLogy();

  double legymax, legymin, legxmin, legxmax;
  legxmin = (legRIGHT ? 0.67 : 0.25);
  legxmax = legxmin+0.15;
  legymax = (legTOP ? 0.85 : 0.55);
  legymin = legymax-0.15;     
  TLegend* leg = new TLegend(legxmin,legymin,legxmax,legymax);
  if (legHeader!="") leg->SetHeader(legHeader);
  leg->SetTextSize(0.035);

  TString options = (option=="" ? "HIST" : option);
 
  //normalize and set y range
  ymax=0.;
  //cout << h.size() << endl;
  for (size_t i=0; i<hn.size(); i++) {  
  //cout << i << endl;
    //if(h[i]->GetNbinsX() != orbin) cout << "WARNING: orbin for " << h[i]->GetName() << " are " << h[i]->GetNbinsX() << endl; //debug - shift of h[][] wrt clu[][]
    if(rebin >0) hn[i].h->Rebin(rebin);//
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

  //enough for the first histo
  if(xTitle!="")hn[0].h->GetXaxis()->SetTitle(xTitle);
  if(ymin>0)hn[0].h->SetMinimum(ymin);
  if(ymax>0)hn[0].h->SetMaximum(ymax);
  if(yTitle!="")hn[0].h->GetYaxis()->SetTitle(yTitle);
  else hn[0].h->GetYaxis()->SetTitle("a.u.");
  hn[0].h->GetXaxis()->SetLabelSize(0.03);
  hn[0].h->GetXaxis()->SetTitleOffset(1);
  hn[0].h->GetXaxis()->SetTitleSize(0.04);    
  hn[0].h->GetYaxis()->SetLabelSize(0.05);
  hn[0].h->GetYaxis()->SetTitleOffset(0.8);
  hn[0].h->GetYaxis()->SetTitleSize(0.06);    

  for (size_t i=0; i<hn.size(); i++) {
    hn[i].h->SetLineColor(colors[i]);
    hn[i].h->SetMarkerColor(colors[i]);
    if(xmax>0)hn[i].h->GetXaxis()->SetRangeUser(xmin,xmax);
    string nam = "";
    if(!twoSamples) leg->AddEntry(hn[i].h); //to print all for bench comp      
    else leg->AddEntry(hn[i].h, Legend[i]);
    if (i==1) options = options + (stat ? "sames" : "same"); //once is enought
    hn[i].h->Draw(options);
  }  
  leg->Draw("same");
  drawCMSprel(dataset, 0.035);
  c->Update();

//debug -- add veto if more than 2 histo..
  c->cd();
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd();
  TH1F* hr = (TH1F*)hn[0].h->Clone("hr");
  hr->Divide(hn[1].h);
  hr->SetMaximum(2.);
  hr->SetMinimum(0.);

  /*TH1F* Err=((TH1F*)gROOT->FindObject("BkgErr")); //debug
  delete Err;
  Err=(TH1F*)hn[1].h->Clone("BkgErr");
  for(int i=1; i<=Err->GetNbinsX(); i++) Err->SetBinError(i, Err->GetBinError(i)/Err->GetBinContent(i));
  for(int i=1; i<=Err->GetNbinsX(); i++) Err->SetBinContent(i, 1); // Bin content */

// Y axis ratio plot settings
  hr->GetYaxis()->SetTitle("R");
  hr->GetYaxis()->SetNdivisions(505);
  hr->GetYaxis()->SetTitleSize(20);
  hr->GetYaxis()->SetTitleFont(43);
  hr->GetYaxis()->SetTitleOffset(1.35);
  hr->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hr->GetYaxis()->SetLabelSize(15);
  // X axis ratio plot settings
  //hr->GetXaxis()->SetTitleSize(20);
  //hr->GetXaxis()->SetTitleFont(43);
  //hr->GetXaxis()->SetTitleOffset(4.);
  hr->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  hr->GetXaxis()->SetLabelSize(15);
  hr->GetXaxis()->SetTitleOffset(0.6);
  hr->GetXaxis()->SetTitleSize(0.15); 

  hr->SetMarkerStyle(20);
  hr->SetMarkerSize(0.8);

  //Err->Draw("E3");   //understand utility...
  hr->Draw("PE");

  c->Update();

  if(saveCanvas){
    c->SaveAs(OutFld+plotSubFld+name+".png");
    c->SaveAs(OutFld+plotSubFld+name+".pdf");
  }
}
//---------------

void Plot_kinSel::drawH2(TString dataset, TH2F* h2, TString name, TString xTitle, 
          TString yTitle, double xmin, double xmax, double ymin, double ymax,
           bool stat , int orbin , TString option , bool logX , bool logY ) {

  gStyle->SetPadRightMargin(0.09);

  TCanvas* c = new TCanvas(name,name,600,600);
  c->cd();
  if (logX) gPad->SetLogx();
  if (logY) gPad->SetLogy();

  TString options = (option=="" ? "COLZ" : option);
    
  //normalize and set y range
  //if(h[i]->GetNbinsX() != orbin) cout << "WARNING: orbin for " << h[i]->GetName() << " are " << h[i]->GetNbinsX() << endl; //debug - shift of h[][] wrt clu[][]

  if(xmin>0 && xmax>0)h2->GetXaxis()->SetRangeUser(xmin,xmax);
  if(xTitle!="")h2->GetXaxis()->SetTitle(xTitle);
  if(ymin>0)h2->SetMinimum(ymin);
  if(ymax>0)h2->SetMaximum(ymax);
  if(yTitle!="")h2->GetYaxis()->SetTitle(yTitle);
  h2->GetXaxis()->SetLabelSize(0.03);
  h2->GetXaxis()->SetTitleOffset(1);
  h2->GetXaxis()->SetTitleSize(0.04);    
  h2->GetYaxis()->SetLabelSize(0.03);
  h2->GetYaxis()->SetTitleOffset(1);
  h2->GetYaxis()->SetTitleSize(0.04);    
  string nam = "";
  h2->Draw(options);
  c->Update(); 
  TPaletteAxis *palette = (TPaletteAxis*)h2->GetListOfFunctions()->FindObject("palette");
  palette->GetAxis()->SetLabelSize(0.005);
  palette->GetAxis()->SetTitleSize(0.005);

  drawCMSprel(dataset, 0.02);
  c->Update();

  if(saveCanvas){
    c->SaveAs(OutFld+plotSubFld+name+".png");
    c->SaveAs(OutFld+plotSubFld+name+".pdf");
  }

}
//---------------

void Plot_kinSel::drawH2ratio(TString dataset, TH2F* ha, TH2F* hb, TString name, TString xTitle, 
          TString yTitle, double xmin, double xmax, double ymin, double ymax,
           bool stat , int orbin , TString option , bool logX , bool logY ) {

  gStyle->SetPadRightMargin(0.09);

  TCanvas* c1 = new TCanvas(name,name,600,600);
  c1->cd();
  if (logX) gPad->SetLogx();
  if (logY) gPad->SetLogy();

  TString options = (option=="" ? "COLZ" : option);
  if(xmin>0 && xmax>0)ha->GetXaxis()->SetRangeUser(xmin,xmax);
  if(xTitle!="")ha->GetXaxis()->SetTitle(xTitle);
  if(ymin>0)ha->SetMinimum(ymin);
  if(ymax>0)ha->SetMaximum(ymax);
  if(yTitle!="")ha->GetYaxis()->SetTitle(yTitle);
  ha->GetXaxis()->SetLabelSize(0.03);
  ha->GetXaxis()->SetTitleOffset(1);
  ha->GetXaxis()->SetTitleSize(0.04);    
  ha->GetYaxis()->SetLabelSize(0.03);
  ha->GetYaxis()->SetTitleOffset(1);
  ha->GetYaxis()->SetTitleSize(0.04);    
  string nam = "";
  
  ha->Divide(hb);
  ha->Draw(options);
  TPaletteAxis *palette = (TPaletteAxis*)ha->GetListOfFunctions()->FindObject("palette");
  palette->GetAxis()->SetLabelSize(0.005);
  palette->GetAxis()->SetTitleSize(0.005);

  drawCMSprel(dataset, 0.02);
  c1->Update(); 

  if(saveCanvas){
    c1->SaveAs(OutFld+plotSubFld+name+".png");
    c1->SaveAs(OutFld+plotSubFld+name+".pdf");
  }

}
//---------------

