// Code for non resonant HH->bbbb in CMS Run2 
// Step3 - 2D shape fit --> to create 2D plots: BDT output vs CSV
//  Author: Martino Dall'Osso
//   from Carlo Alberto Gottardo
//    Nov 03rd, 2015
//     g++ `root-config --libs  --cflags` hh4b_13TeV_make2D.cpp -o make2D
//      ./make2D BTagCSVRun2015C-D M-260 pT20CSVL_def
//----------------------------------------------------------------

#include <iostream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TAxis.h"
#include <vector>
#include <cmath>
#include <string>
#include <fstream>

//WARNING: CHANGE ME!!
static const std::string frameworkVersionFld = "V15/";   //equal to the Data version   
static const std::string outoption="";
float CSVcut = 0.605;  //Run2: low 0.605; medium 0.890; high 0.970. 
int massbins = 23; //was 23
int bdtbins = 20;  //was 20

static const std::string dataFld="../data/"+frameworkVersionFld;;
static const std::string resultsFld="results/"+frameworkVersionFld;
static const std::string outFile = resultsFld+"scatters"+outoption+".root";
//INPUT FILE NAMEs DEFINED IN THE MAIN (bottom)


void domake2D(std::string inputFile, std::string inputFile_CR, std::string inputMC){
  TH2F *mass_scatt_ALL = new TH2F("BDT_vs_mass", "BDT vs mass", bdtbins, 0.25, 0.75, massbins, 40, 500); 
  TH2F *mass_scatt_3CSV = new TH2F("BDT_vs_mass_3CSV", "BDT vs mass (3btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
  TH2F *mass_scatt_4CSV = new TH2F("BDT_vs_mass_4CSV", "BDT vs mass (4btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
  TH2F *mass_scatt_2CSV = new TH2F("BDT_vs_mass_2CSV", "BDT vs mass (2btag)", bdtbins, 0.25, 0.75, massbins, 40, 500);
  TH2F *mass_scatt_TEST = new TH2F("BDT_vs_mass_TEST", "2-3 btag syst. error / 4 btag stat. unc. ", bdtbins, 0.25, 0.75, massbins, 40, 500);  // not used...
  TH2F *MC_mass_scatt_4CSV = new TH2F("MC_BDT_vs_mass_4CSV", "MC BDT vs mass (4btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
  TH2F *MC_mass_scatt_3CSV = new TH2F("MC_BDT_vs_mass_3CSV", "MC BDT vs mass (3btag)", bdtbins, 0.25, 0.75, massbins, 40, 500);
  TH2F *MC_mass_scatt_ALL = new TH2F("MC_BDT_vs_mass_TEST", "MC BDT vs mass (ALL)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
  TH2F *diff = new TH2F("twothree_btag_Residuals", "2-3 btag Residuals", bdtbins, 0.25, 0.75, massbins, 40, 500);
  gStyle->SetOptStat(0);
  mass_scatt_ALL->GetXaxis()->SetTitle("BDT response"); mass_scatt_ALL->GetYaxis()->SetTitle("h mass");
  mass_scatt_3CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_3CSV->GetYaxis()->SetTitle("h mass");
  mass_scatt_4CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_4CSV->GetYaxis()->SetTitle("h mass");
  mass_scatt_2CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_2CSV->GetYaxis()->SetTitle("h mass");
  mass_scatt_TEST->GetXaxis()->SetTitle("BDT response"); mass_scatt_TEST->GetYaxis()->SetTitle("h mass");
  MC_mass_scatt_4CSV->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_4CSV->GetYaxis()->SetTitle("h mass");
  MC_mass_scatt_3CSV->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_3CSV->GetYaxis()->SetTitle("h mass");
  MC_mass_scatt_ALL->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_ALL->GetYaxis()->SetTitle("h mass");

  // Control Region -- debug - 2btag
  //------------------------------------
  TChain *upperTree = new TChain("tree");
  upperTree->Add(inputFile_CR.c_str());
  //upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-24.root");
  //upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-Pr.root");
  //upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012D-Pr.root");
  float u_tr, u_eta, u_pt, u_BDT, u_csv, u_m1, u_m2;
  //upperTree->SetBranchAddress("QPt_3", &u_pt);
  //upperTree->SetBranchAddress("QEta_3", &u_eta);
  //upperTree->SetBranchAddress("ThirdJetTracks", &u_tr);
  upperTree->SetBranchAddress("fJet4_CSV", &u_csv);
  upperTree->SetBranchAddress("BDT", &u_BDT);
  upperTree->SetBranchAddress("H1_mass", &u_m1);
  upperTree->SetBranchAddress("H2_mass", & u_m2);
  std::cout << "BKG templ Entries: " << upperTree->GetEntries() << std::endl;
  if(upperTree->GetEntries() == 0) return;
  for(int i=0; i<upperTree->GetEntries(); i++) {
    upperTree->GetEntry(i);
    float mass = (u_m1+u_m2)/2.0;
    mass_scatt_ALL->Fill(u_BDT, mass);
    mass_scatt_2CSV->Fill(u_BDT, mass);
  }

  //3btag - 4btag DATA
  //------------------------------------
  TChain *lowerTree = new TChain("tree");
  lowerTree->Add(inputFile.c_str());
  //lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug.root");
  //lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-Promp.root");
  //lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012D-Promp.root");	
  float d_tr, d_eta, d_pt, d_BDT, d_csv, d_m1, d_m2;
  //lowerTree->SetBranchAddress("QPt_3", &d_pt);
  //lowerTree->SetBranchAddress("QEta_3", &d_eta);
  //lowerTree->SetBranchAddress("ThirdJetTracks", &d_tr);
  lowerTree->SetBranchAddress("fJet4_CSV", &d_csv);
  lowerTree->SetBranchAddress("BDT", &d_BDT);
  lowerTree->SetBranchAddress("H1_mass", & d_m1);
  lowerTree->SetBranchAddress("H2_mass", & d_m2);
  std::cout << "SIG templ Entries: " << lowerTree->GetEntries() << std::endl;
  if(lowerTree->GetEntries() == 0) return;
  for(int i=0; i<lowerTree->GetEntries(); i++) {
    lowerTree->GetEntry(i);
    float mass = (d_m1+d_m2)/2.0;
    mass_scatt_ALL->Fill(d_BDT, mass);
    if(d_csv < CSVcut) mass_scatt_3CSV->Fill(d_BDT, mass);
    if(d_csv >= CSVcut) mass_scatt_4CSV->Fill(d_BDT, mass);
  }

  //check 2btag vs 3 btag - fill residuals
  //------------------------------------
  int TOTbins = mass_scatt_2CSV->GetXaxis()->GetNbins()*mass_scatt_2CSV->GetYaxis()->GetNbins();
  float delta, unc;
  float N2 = mass_scatt_2CSV->GetEntries();
  float N3 = mass_scatt_3CSV->GetEntries();
  float bin2, bin3, bin4;
  float test;
  for(int i=1; i<=TOTbins; i++) {
    bin2 = mass_scatt_2CSV->GetBinContent(i);
    bin3 = mass_scatt_3CSV->GetBinContent(i);
    bin4 = mass_scatt_4CSV->GetBinContent(i);
    if (bin2 == 0 || bin3 == 0) {
      delta = 0.0; 
      unc = 1.1;
    }
    else { 
      delta = bin2 - (bin3*(N2/N3));  //WARNING!!
      unc = sqrt(bin2 + (N2/N3)*(N2/N3)*bin3 + (bin3/N3)*(bin3/N3)*N2 + ((bin3*N2)/(N3*N3))*((bin3*N2)/(N3*N3))*N3 );
      if (bin4 != 0) {test = ((((N3/N2)*bin2)-bin3)/bin3)*sqrt(bin4);}
      if (bin4 == 0) {test = 0.0;}
    }
    diff->SetBinContent(i,delta/unc);
    mass_scatt_TEST->SetBinContent(i,test);
  }

  //3btag - 4btag MC
  //------------------------------------
  TFile * MC = new TFile(inputMC.c_str(), "READ");
  //BDT_NOCUT_hh4b-nores-8Tev300k_step2_L1y1.root", "READ");
  TTree * leggo = (TTree*)MC->Get("tree");
  float mc_BDT, mc_csv, mc_m1, mc_m2;
  //lowerTree->SetBranchAddress("QPt_3", &d_pt);
  //lowerTree->SetBranchAddress("QEta_3", &d_eta);
  //lowerTree->SetBranchAddress("ThirdJetTracks", &d_tr);
  leggo->SetBranchAddress("fJet4_CSV", &mc_csv);  //4th Jet CSV
  leggo->SetBranchAddress("BDT", &mc_BDT);
  leggo->SetBranchAddress("H1_mass", & mc_m1); //diJet1 mass
  leggo->SetBranchAddress("H2_mass", & mc_m2); //diJet2 mass
  std::cout << "SIG MC Entries: " << leggo->GetEntries() << std::endl;
  if(leggo->GetEntries() == 0) return;
  for(int i=0; i<leggo->GetEntries(); i++) {
    leggo->GetEntry(i);
    float mass = (mc_m1+mc_m2)/2.0;
    MC_mass_scatt_ALL->Fill(mc_BDT, mass);
    if(mc_csv < CSVcut) MC_mass_scatt_3CSV->Fill(mc_BDT, mass);
    if(mc_csv >= CSVcut) MC_mass_scatt_4CSV->Fill(mc_BDT, mass);
  }
  MC->Close();

  TFile *outhisto = new TFile(outFile.c_str(), "RECREATE");
  mass_scatt_ALL->Write();
  mass_scatt_2CSV->Write();
  mass_scatt_3CSV->Write();
  mass_scatt_4CSV->Write();
  mass_scatt_TEST->Write();
  MC_mass_scatt_3CSV->Write();
  MC_mass_scatt_4CSV->Write();
  MC_mass_scatt_ALL->Write();
  diff->Write();
  outhisto->Close();

  return;
}

//-----------------------------------
int main (int argc, char **argv) {  //DEBUG!!!!

  std::string dataSample, MCSample, opt;  
  for(int i=1; i<argc; i++) {
    if(i==1) dataSample = argv[i];
    if(i==2) MCSample = argv[i];
    if(i==3) opt = argv[i];
  }  

  std::string inputFile_CR_= dataFld+"tree_Step2_"+dataSample+"_"+opt+"_CR.root"; // control region (2-btag only?)
  std::string inputFile_= dataFld+"tree_Step2_"+dataSample+"_"+opt+".root";
  std::string inputMC_= dataFld+"tree_Step2_"+MCSample+"_"+opt+".root";
  
  std::cout << "Make 2D plots reading: " << std::endl;
  std::cout << inputFile_CR_ << std::endl;
  std::cout << inputFile_ << std::endl;
  std::cout << inputMC_ << std::endl;

  domake2D(inputFile_ , inputFile_CR_, inputMC_);
  return 0;
}


