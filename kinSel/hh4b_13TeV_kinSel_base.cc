// Code for non resonant HH->bbbb in CMS Run2 
// Step1 - KINEMATIC SELECTION
//  Author: Martino Dall'Osso
//   thanks Souvik Das (Univ. of Florida) & Caterina Vernieri (FNAL) 
//    Sept 21th, 2015
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

double pi=3.14159265;

// selection parameters:
bool yesTrgA = true;
bool yesTrgB = true;
double Jet_pt_cut_low = 20.; //30.;
double deltaRCut = 0.5; //0.5
int nJets_cut = 4; //3
double Jet_eta_cut = 10; //2.5 - 10=noCut
double CSV_cut = 0.1; //Run2: low 0.605; medium 0.890; high 0.970. Run1: low 0.679 (used in the first presentation).
double CMVA_cut_med = 0.1; //0.71; //0.71;
double H_mass = 115.0; //for withinRegion
double dijetM_cut_low = 100.;
double dijetM_cut_high = 150.;

// matrix parameters
int binPt = 10;
int binCSV = 8;
int binEta = 4;

// environment
//------------------
int finalIndex = 2; //0 MCTruth, 1 matrix, 2 CSV, 3 minMass
//std::string rootfolder= "/lustre/cmswork/dallosso/hh2bbbb/non-resonant/event_generation/CMSSW_7_2_2_patch2/src/nores_ntuples13TeV/Ntuples_L1y1C0/Loop_2/";
std::string utilsFld="../utils/";
std::string dataFld="../data/";
std::string plotsFld="plots/";
std::string matrixFld="jetsMatch_matrix/";
std::string subdataFld=dataFld+"vhbb_heppy_V13/";

typedef struct{
  double CMVA;
  double CSV;
  double mass;
  double pT;
  double eta;
  double phi;
} Jet;

//get TLorentzVector from Jet
TLorentzVector get_jetVector(Jet* jj){ 
  TLorentzVector jetv;
  jetv.SetPtEtaPhiM(jj->pT,jj->eta,jj->phi,jj->mass);
  return jetv;
}

//Jets sorting by CSV
bool cmp_CSV(Jet jet1, Jet jet2){
  return jet1.CSV > jet2.CSV;
}

//Jets sorting by CSV
bool cmp_CMVA(Jet jet1, Jet jet2){
  return jet1.CMVA > jet2.CMVA;
}

//star angle computation    
double computeCosThetaStar(TLorentzVector dijet, TLorentzVector dihiggs){     
  dijet.Boost(-dihiggs.BoostVector()); //equal to (-P1.px()/P1.e(),-P1.py()/P1.e(),-P1.pz()/P1.e())
  double costhetast = dijet.CosTheta();
  return costhetast;
}

//check...
int withinElipse(double mH1, double mH2, double a=30., double b=60., double mH1_c=H_mass, double mH2_c=H_mass)
{
 //45deg antirotation
  double fact = (sqrt(2)/2); //1; //
  double mH1_ = (mH1-mH1_c)*fact - (mH2-mH2_c)*fact;
  double mH2_ = (mH1-mH1_c)*fact + (mH2-mH2_c)*fact;
  double a_ = a*fact;
  double b_ = b*fact ;
  
  double elips = (mH1_/a_)*(mH1_/a_) + (mH2_/b_)*(mH2_/b_);

  int ret=-1;
  if(elips == 1) ret = 1;
  else if(elips > 1) ret = 0;
  return ret;
}

//check...
int withinRegion(double mH1, double mH2, double r1=15., double r2=30., double mH1_c=H_mass, double mH2_c=H_mass)
{
  double r=pow(pow(mH1-mH1_c, 2)+pow(mH2-mH2_c, 2), 0.5);
  double angle=atan2(mH2-mH2_c, mH1-mH1_c);
  int ret=-1;
  if (r<r1) ret=0;
  else if (r>r1 && r<r2){
    if (angle>=0 && angle<pi/2.) ret=1;
    else if (angle>=pi/2. && angle<pi) ret=4;
    else if (angle<0 && angle>=-pi/2.) ret=2;
    else if (angle<pi/2.&& angle>=-pi) ret=3;
    else std::cout<<"This is within annulus but not within any CR!"<<std::endl;
  }
  else ret=5;
  return ret;
}

//-------------
//book Histos
//-------------
TH1F *h_nJets_InAcc=new TH1F("h_nJets_InAcc", "# Central Jets with p_{T}> 30 GeV; # jets", 18, 0., 18.);
TH1F *h_nJets=new TH1F("h_nJets", "# Jets; # jets", 18, 0., 18.);
TH1F *h_nJets_90=new TH1F("h_nJets_90", "# Central Jets with p_{T}> 90 GeV;  n", 10, 0., 10.); 
//  TH1F *h_nCand=new TH1F("h_nCand", "# HH candidates", 20, 0., 20.);
//  TH1F *h_nCand_true=new TH1F("h_nCand_true", "# HH candidates if true", 20, 0., 20.);	
TH1F *h_nPV=new TH1F("h_nPV", "# of Primary Vertices; nPV", 10, 0., 10.);
TH1F *h_nPV_weighted=new TH1F("h_nPV_weighted", "# of Primary Vertices after Reweighting; nPV", 50, 0., 50.);
TH1F *h_HLT_HH4b=new TH1F("h_HLT_HH4b", "h_HLT_HH4b; Quad - Double", 2, 0, 2);
TH1F *h_HLT_HH4ball=new TH1F("h_HLT_HH4ball", "h_HLT_HH4ball; Quad - Double", 1, 0, 1);

TH1F *h_Jet1_pT=new TH1F("h_Jet1_pT", "Jet1_pT; p_{T} (GeV/c)", 100, 0., 900.);
TH1F *h_Jet2_pT=new TH1F("h_Jet2_pT", "Jet2_pT; p_{T} (GeV/c)", 100, 0., 500.);
TH1F *h_Jet3_pT=new TH1F("h_Jet3_pT", "Jet3_pT; p_{T} (GeV/c)", 100, 0., 350.);
TH1F *h_Jet4_pT=new TH1F("h_Jet4_pT", "Jet4_pT; p_{T} (GeV/c)", 100, 0., 350.);
TH1F *h_Jet1_Eta=new TH1F("h_Jet1_Eta", "Jet1_Eta; #eta", 100, -4., 4.);
TH1F *h_Jet2_Eta=new TH1F("h_Jet2_Eta", "Jet2_Eta; #eta", 100, -4., 4.);
TH1F *h_Jet3_Eta=new TH1F("h_Jet3_Eta", "Jet3_Eta; #eta", 100, -4., 4.);
TH1F *h_Jet4_Eta=new TH1F("h_Jet4_Eta", "Jet4_Eta; #eta", 100, -4., 4.);

TH1F *h_Jet1_CSV=new TH1F("h_Jet1_CSV", "Jet1_CSV; CSV", 50, 0.6, 1.);
TH1F *h_Jet2_CSV=new TH1F("h_Jet2_CSV", "Jet2_CSV; CSV", 50, 0.6, 1.);
TH1F *h_Jet3_CSV=new TH1F("h_Jet3_CSV", "Jet3_CSV; CSV", 50, 0.6, 1.);
TH1F *h_Jet4_CSV=new TH1F("h_Jet4_CSV", "Jet4_CSV; CSV", 100, 0., 1.);
TH1F *h_Jet1_CMVA=new TH1F("h_Jet1_CMVA", "h_Jet1_CMVA; CMVA", 100, 0., 1.);
TH1F *h_Jet2_CMVA=new TH1F("h_Jet2_CMVA", "h_Jet2_CMVA; CMVA", 100, 0., 1.);
TH1F *h_Jet3_CMVA=new TH1F("h_Jet3_CMVA", "h_Jet3_CMVA; CMVA", 100, 0., 1.);
TH1F *h_Jet4_CMVA=new TH1F("h_Jet4_CMVA", "h_Jet4_CMVA; CMVA", 100, 0., 1.);

TH1F *h_MET=new TH1F("h_MET", "MET; MET (GeV)", 100, 0, 200.);
//  TH1F *h_MET_sig=new TH1F("h_MET_sig", "MET Significance; Sig", 20, 0., 20.);

TH1F *h_3Jets_avgCSV=new TH1F("h_3Jets_avgCSV", "3Jets_avgCSV; CSV", 50, 0.6, 1.);
TH1F *h_3Jets_minCSV=new TH1F("h_3Jets_minCSV", "3Jets_minCSV; CSV", 50, 0.6, 1.);

TH1F *h_Jets_mass=new TH1F("h_Jets_mass", "Jets_mass; m (GeV/c^{2})", 100, 0., 150.);
TH1F *h_Jets_pT=new TH1F("h_Jets_pT", "Jets_pT; p_{T} (GeV/c)", 250, 0., 500.);
TH1F *h_Jets_eta = new TH1F("h_Jets_eta","h_Jets_eta; #eta", 100, -4., 4.);
TH1F *h_Jets_Phi=new TH1F("h_Jets_Phi", "h_Jets_phi; #phi", 100, -3.5, 3.5);
TH1F *h_Jets_CSV=new TH1F("h_Jets_CSV", "Jets_CSV; CSV", 200, 0., 1.);
TH1F *h_Jets_CMVA=new TH1F("h_Jets_CMVA", "Jets_CMVA; CMVA", 50, 0., 1.);
TH1F *h_Jets_Centr=new TH1F("h_Jets_Centrality", "Jets_Centrality; C", 100, 0., 1.);
TH1F *h_Jets_HT=new TH1F("h_Jets_HT", "h_Jets_HT; HT (GeV/c)", 50, 0., 900.);

TH1F *h_H_mass=new TH1F("h_H_mass", "H mass; m (GeV/c^{2})", 100, 50., 180.);
TH1F *h_H_pT=new TH1F("h_H_pT", "H p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_H_Phi=new TH1F("h_H_Phi", "H phi; #phi", 100, -3.5, 3.5);
TH1F *h_H_Eta=new TH1F("h_H_Eta", "H eta; #eta", 100, -6., 6.);
TH1F *h_H1_mass=new TH1F("h_H1_mass", "H1 mass; m (GeV/c^{2})", 100, 50., 250.);
TH1F *h_H1_pT=new TH1F("h_H1_pT", "H1 p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_H1_Phi=new TH1F("h_H1_Phi", "H1 phi; #phi", 100, -3.5, 3.5);
TH1F *h_H1_Eta=new TH1F("h_H1_Eta", "H1 eta; #eta", 100, -6., 6.);
TH1F *h_H1_CosThSt=new TH1F("h_H1_CosTheta*", "h_H1_CosTheta*; |cos#theta*|", 50, 0., 1.);
TH1F *h_H1_deltaR = new TH1F("h_H1_DR","h_H1_DR; #DeltaR", 100, 0., 7.);
TH1F *h_H1_deltaPhi = new TH1F("h_H1_DPhi","h_H1_DPhi; #Delta#phi", 100, -3.5, 3.5);
TH1F *h_H1_deltaEta = new TH1F("h_H1_DEta","h_H1_DEta; #Delta#eta", 100, -6., 6.);
TH2F *h_H1_deltaPhiVSpT = new TH2F("h_H1_deltaPhiVSpT", "h_H1_deltaPhiVSpT", 200, 0., 3.5, 200, 0., 600.);

TH1F *h_H2_mass=new TH1F("h_H2_mass", "H2 mass; m (GeV/c^{2})", 100, 50., 250.);
TH1F *h_H2_pT=new TH1F("h_H2_pT", "H2 p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_H2_CosThSt=new TH1F("h_H2_CosTheta*", "h_H2_CosTheta*; |cos#theta*|", 50, 0., 1.);
TH1F *h_H2_deltaR = new TH1F("h_H2_DR","h_H2_DR; #DeltaR", 100, 0., 7.);
TH1F *h_H2_deltaPhi = new TH1F("h_H2_DPhi","h_H2_DPhi; #Delta#phi", 100, -3.5, 3.5);
TH1F *h_H2_deltaEta = new TH1F("h_H2_DEta","h_H2_DEta; #Delta#eta", 100, -6., 6.);
TH2F *h_H2_deltaPhiVSpT = new TH2F("h_H2_deltaPhiVSpT", "h_H2_deltaPhiVSpT; #Delta#phi; p_{T} (GeV/c)", 200, 0., 3.5, 200, 0., 600.);

TH1F *h_HH_mass = new TH1F("h_HH_mass"," HH mass; m (GeV/c^{2})" , 100, 0., 1500.); 
TH1F *h_HH_pT=new TH1F("h_HH_pT", "HH p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_HH_Eta=new TH1F("h_HH_Eta", "HH eta; #eta", 100, -6., 6.);
TH1F *h_HH_Phi=new TH1F("h_HH_Phi", "HH phi; #phi", 100, -3.5, 3.5);
TH1F *h_HH_deltaR = new TH1F("h_HH_DR","h_HH_DR; #DeltaR", 100, 0., 7.);
TH1F *h_HH_deltaPhi = new TH1F("h_HH_DPhi","h_HH_DPhi; #Delta#phi", 100, -3.5, 3.5);
TH1F *h_HH_deltaEta = new TH1F("h_HH_DEta","h_HH_DEta; #Delta#eta", 100, -6., 6.);
TH2F *h_H1_H2_mass = new TH2F("h_H1_H2_mass", " mh mh; m_{H_{lead}} (GeV/c^{2}); m_{H_{sublead}} (GeV/c^{2})", 300, 0., 600., 300, 0., 600.);

TH1F *h_HH_massInReg = new TH1F("h_HH_massInReg","HH_mass SR; m (GeV/c^{2}) SR" , 200, 0., 1500.);
TH1F *h_HH_pTInReg=new TH1F("h_HH_pTInReg", "HH p_{T} SR; p_{T} (GeV/c) SR", 50, 0., 900.);
TH2F *h_H1_H2_massInReg = new TH2F("h_H1_H2_massInReg", " all comb if min(#sigma (#delta m)^{2}) SR; m (GeV/c^{2}) SR; m (GeV/c^{2}) SR ", 300, 0., 600., 300, 0., 600.);
TH2F *h_H1_H2_massInReg2 = new TH2F("h_H1_H2_massInReg2", " all comb if min(#sigma (#delta m)^{2}) SR2; m (GeV/c^{2}) SR; m (GeV/c^{2}) SR ", 300, 0., 600., 300, 0., 600.);
TH2F *h_H1_H2_massInReg3 = new TH2F("h_H1_H2_massInReg3", " all comb if min(#sigma (#delta m)^{2}) SR3; m (GeV/c^{2}) SR; m (GeV/c^{2}) SR ", 300, 0., 600., 300, 0., 600.);

TH1F *h_genHH_mass=new TH1F("h_genHH_mass", "Generator HH mass; m_{genHH} (GeV/c^{2})", 100, 0., 1000.);
TH1F *h_genH1_mass=new TH1F("h_genH1_mass", "Generator H1 mass; m_{genH1} (GeV/c^{2})", 100, 0., 1000.);
TH1F *h_genH2_mass=new TH1F("h_genH2_mass", "Generator H2 mass; m_{genH2} (GeV/c^{2})", 100, 0., 1000.);
TH1F *h_nGenH=new TH1F("h_nGenH", "# generated H per event; # genH", 8, 0., 8.);
TH1F *h_genH_pT=new TH1F("h_genH_pT", "Generate H p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_genH_mass=new TH1F("h_genH_mass", "Generate H mass; mass (GeV/c^{2})", 100, 50., 180.);

TH1F *h_genB1Jets_DR = new TH1F("h_genB1Jets_DeltaR","h_genB1Jets_DeltaR; #DeltaR", 100, 0., 7.);
TH1F *h_genBfH1_DR=new TH1F("h_genBfH1_DR", "DR between genBfromH1 and jet; #DeltaR", 100, 0., 7.);
TH1F *h_genBfH2_DR=new TH1F("h_genBfH2_DR", "DR between genBfromH2 and jet; #DeltaR", 100, 0., 7.);

//  TH1F *h_genBfromH_pT=new TH1F("h_genH_pT", "Generate H p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
//  TH1F *h_genBfromH_pT=new TH1F("h_genH_pT", "Generate H p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
//  TH1F *h_genBfromH_pT=new TH1F("h_genH_pT", "Generate H p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
//  TH1F *h_genBfromH_pT=new TH1F("h_genH_pT", "Generate H p_{T}; p_{T} (GeV/c)", 50, 0., 900.);
//  TH1F *h_genH_mass=new TH1F("h_genH_mass", "Generate H mass; mass (GeV)", 100, 50., 180.);

TH1F *h_HH_mass_diff_cand=new TH1F("h_HH_mass_diff_cand", "|#Deltam| between Higgs masses - all candidates", 50, 0., 200.);
TH1F *h_HH_massNorm_diff=new TH1F("h_HH_massNorm_diff", "|#Deltam| between Higgs masses", 50, 0., 2.);
TH1F *h_HH_CSV = new TH1F("h_HH_CSV"," Sum CSV | between the two higgs ", 70, -4., 4.); 
//  TH1F *h_HHmass_right = new TH1F("h_HHmass_right"," h_HHmass right" , 100, 0., 1000.); 
//  TH1F *h_HHmass_wrong = new TH1F("h_HHmass_wrong"," h_HHmass wrong" , 100, 0., 1000.);
//  TH1F *h_HHmass_first = new TH1F("h_HHmass_first"," h_HHmass first" , 100, 0., 1000.);
//  TH1F *h_HHmass_other = new TH1F("h_HHmass_other"," h_HHmass other" , 100, 0., 1000.);		

TH1F *h_Cuts=new TH1F("h_Cuts", "Cut flow", 20, 0, 20);

//new histos
TH1F *  h_jet1b_dr = new TH1F("h_jet1b_dr","h_jet1b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet2b_dr = new TH1F("h_jet2b_dr","h_jet2b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet3b_dr = new TH1F("h_jet3b_dr","h_jet3b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet4b_dr = new TH1F("h_jet4b_dr","h_jet4b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet4b_drMatrix = new TH1F("h_jet4b_drMatrix","h_jet4b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet4b_drAll = new TH1F("h_jet4b_drAll","h_jet4b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet4b_drNotMatched = new TH1F("h_jet4b_drNotMatched","h_jet4b_dr; jet-parton #DeltaR", 100, 0., 1.);

TH1F *h_Jet4match_pT=new TH1F("h_Jet4match_pT", "Jet4match_pT; 4th jet p_{T} (GeV/c)", 100, 0., 900.);
TH1F *h_Jet4all_pT=new TH1F("h_Jet4all_pT", "Jet4All_pT; 4th jet p_{T} (GeV/c)", 100, 0., 900.);
TH1F *h_Jet4match_eta=new TH1F("h_Jet4match_eta", "Jet4match_Eta; 4th jet #eta", 100, -4., 4.);
TH1F *h_Jet4all_eta=new TH1F("h_Jet4all_eta", "Jet4all_Eta; 4th jet #eta", 100, -4., 4.);
TH1F *h_Jet4match_CSV=new TH1F("h_Jet4match_CSV", "Jet4match_CSV; 4th jet CSV", 50, 0., 1.);
TH1F *h_Jet4all_CSV=new TH1F("h_Jet4all_CSV", "Jet4all_CSV; 4th jet CSV", 50, 0., 1.);

TH1F *h_Jet4match_DpT3=new TH1F("h_Jet4match_DpT3", "Jet4match_DpT3; 3rd-4th jets #Deltap_{T} (GeV/c)", 100, -300., 300.);
TH1F *h_Jet4all_DpT3=new TH1F("h_Jet4all_DpT3", "Jet4All_DpT3; 3rd-4th jets #Delta p_{T} (GeV/c)", 100, -300., 300.);
TH1F *h_Jet4match_DCSV3=new TH1F("h_Jet4match_DCSV3", "Jet4match_DCSV3; 3rd-4th jets #Delta CSV", 50, 0., 1.);
TH1F *h_Jet4all_DCSV3=new TH1F("h_Jet4all_DCSV3", "Jet4all_DCSV3; 3rd-4th jets #Delta CSV", 50, 0., 1.);
TH1F *h_Jet4match_Deta3=new TH1F("h_Jet4match_Deta3", "h_Jet4match_Deta3; 3rd-4th jets #Delta#eta", 100, 0., 4.);
TH1F *h_Jet4all_Deta3=new TH1F("h_Jet4all_Deta3", "h_Jet4all_Deta3; 3rd-4th jets #Delta#eta", 100, 0., 4.);

TH1F *h_M_true=new TH1F("h_M_true", "<mass> true; <m_{H}> (GeV/c^{2})", 150, 0., 300.);
TH1F *h_M_matrix=new TH1F("h_M_matrix", "<mass> matrix; <m_{H}> (GeV/c^{2})", 150, 0., 300.);
TH1F *h_M_csv=new TH1F("h_M_csv", "<mass> csv; <m_{H}> (GeV/c^{2})", 150, 0., 300.);
TH1F *h_M_highpt=new TH1F("h_M_highpt", "<mass> highpt; <m_{H}> (GeV/c^{2})", 150, 0., 300.);
TH1F *h_M_matrixCSV=new TH1F("h_M_matrixCSV", "<mass> matrixCSV; <m_{H}> (GeV/c^{2})", 150, 0., 300.);

TH1F *h_HH_mass_tr = new TH1F("h_HH_mass_tr"," HH mass; m_{HH} (GeV/c^{2})" , 100, 0., 1500.); 
TH1F *h_HH_mass_matr = new TH1F("h_HH_mass_matr"," HH mass; m_{HH} (GeV/c^{2})" , 100, 0., 1500.); 
TH1F *h_HH_mass_csv = new TH1F("h_HH_mass_csv"," HH mass; m_{HH} (GeV/c^{2})" , 100, 0., 1500.); 

TH2F *h_CSVmass_no4thCSV = new TH2F("h_CSVmass_no4thCSV", " no4thCSV; <m_{H}> (GeV/c^{2}); CSV", 300, 0., 600., 200, 0., 1.);
TH2F *h_CSVmass_4thCSV = new TH2F("h_CSVmass_4thCSV", " 4thCSV; <m_{H}> (GeV/c^{2}); CSV", 300, 0., 600., 200, 0., 1.);
TH2F *h_CSVmass_4thCSV_b = new TH2F("h_CSVmass_4thCSV_b", " 4thCSV; <m_{H}> (GeV/c^{2}); CSV", 300, 0., 600., 200, 0., 1.);

TH2F *h_HCSV_HRL_mass = new TH2F("h_HCSV_HRL_mass", " mh mh; <m_{H_{RL}}> (GeV/c^{2}); <m_{H_{CSV}}> (GeV/c^{2})", 200, 0., 600., 200, 0., 600.);
TH2F *h2_HCSV_HRL_mass = new TH2F("h2_HCSV_HRL_mass", " mh mh; <m_{H_{RL}}> (GeV/c^{2}); <m_{H_{CSV}}> (GeV/c^{2})", 200, 0., 600., 200, 0., 600.);

//  TH2F *h_R_pTeta = new TH2F("h_R_pTeta", " R; pT (GeV/c); eta)", 100, 0., 500., 100, -4., 4.);
//  TH2F *h_R_pTCSV = new TH2F("h_R_pTCSV", " R; pT (GeV/c); CSV)", 100, 0., 500., 50, 0., 1.);
//  TH2F *h_R_CSVeta = new TH2F("h_R_CSVeta", " R; CSV; eta)", 50, 0., 1., 100, -4., 4.);

//-------------
//   MAIN
//-------------
void hh4b_kinSel(std::string sample, std::string opt, int maxEvents =0, bool isData = false, bool isSignalMC = true, bool draw = false)
{

  std::string inputfilename;  
  //std::ifstream inf(filesList);
  //getline(inf,inputfilename);
  inputfilename = subdataFld+"tree_"+sample+".root"; 

  TFile * f = new TFile(inputfilename.c_str());
  f->cd();	
  TTree *tree=(TTree*)f->Get("tree");
  std::cout<<"Opened input file "<<inputfilename<<std::endl;
   // Book variables
  ULong64_t evt;
  unsigned int run;
  int nJet, nGenH, nPV, nGenBfHafterISR;    
  double Jet_id[20],Jet_puId[20],Jet_btagCSV[20],Jet_btagCMVA[20],Jet_rawPt[20],Jet_mcPt[20],Jet_mcFlavour[20];
  double Jet_corr_JECUp[20], Jet_corr_JECDown[20], Jet_corr[20];
  double Jet_mcMatchId[20],Jet_pt[20],Jet_eta[20],Jet_phi[20],Jet_mass[20],Jet_hadronFlavour[20],Jet_btagProb[20];
  double Jet_btagBProb[20], Jet_btagnew[20],Jet_btagCSVV0[20], Jet_chHEF[20],Jet_neHEF[20],Jet_chEmEF[20],Jet_neEmEF[20];
  double Jet_chMult[20],Jet_leadTrackPt[20],Jet_mcEta[20],Jet_mcPhi[20],Jet_mcM[20], ht;
  double GenH_pt[2], GenH_mass[2], GenH_eta[2], GenH_phi[2], GenH_status[2];
  double nGenBfH[4], GenBfH_pdgId[4],GenBfH_pt[4],GenBfH_eta[4],GenBfH_phi[4],GenBfH_mass[4], GenBfH_status[4], GenBfH_charge[4];
  double GenBfHafterISR_pdgId[4],GenBfHafterISR_pt[4],GenBfHafterISR_eta[4],GenBfHafterISR_phi[4],GenBfHafterISR_mass[4];
  float weightPU, Vtype_, met_pt,met_eta,met_phi,met_mass;
  float HLT_HH4bAll, HLT_BIT_QuadTriple, HLT_BIT_QuadDouble, HLT_BIT_DoubleTriple, HLT_BIT_DoubleDouble; 


  //Retrieve variables
  //--------------------
  tree->SetBranchAddress("nprimaryVertices", &(nPV));
  tree->SetBranchAddress("Vtype", &(Vtype_));
  tree->SetBranchAddress("evt",&evt); 
  tree->SetBranchAddress("run",&run);

  if(isData){
    tree->SetBranchAddress("HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v",&HLT_BIT_QuadTriple);          //new
    tree->SetBranchAddress("HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v",&HLT_BIT_QuadDouble);          //new
    tree->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v",&HLT_BIT_DoubleTriple);          //new
    tree->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v",&HLT_BIT_DoubleDouble);          //new
    tree->SetBranchAddress("HLT_HH4bAll",&HLT_HH4bAll);          //new
  }
  else {
    tree->SetBranchAddress("puWeight", &(weightPU)); //added
    tree->SetBranchAddress("HLT_BIT_HLT_QuadJet45_TripleCSV0p5_v",&HLT_BIT_QuadTriple);          //new
    tree->SetBranchAddress("HLT_BIT_HLT_QuadJet45_DoubleCSV0p5_v",&HLT_BIT_QuadDouble);          //new
    tree->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_TripleCSV0p5_v",&HLT_BIT_DoubleTriple);          //new
    tree->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_DoubleCSV0p5_v",&HLT_BIT_DoubleDouble);          //new
    tree->SetBranchAddress("HLT_HH4bAll",&HLT_HH4bAll);          //new

    tree->SetBranchAddress("Jet_hadronFlavour",&Jet_hadronFlavour); 
    tree->SetBranchAddress("Jet_mcPt",&Jet_mcPt);          
    tree->SetBranchAddress("Jet_mcFlavour",&Jet_mcFlavour);     
    tree->SetBranchAddress("Jet_mcMatchId",&Jet_mcMatchId);     
    tree->SetBranchAddress("Jet_mcEta",&Jet_mcEta);         
    tree->SetBranchAddress("Jet_mcPhi",&Jet_mcPhi);         
    tree->SetBranchAddress("Jet_mcM",&Jet_mcM);           
    tree->SetBranchAddress("Jet_corr_JECUp",&Jet_corr_JECUp);     
    tree->SetBranchAddress("Jet_corr_JECDown",&Jet_corr_JECDown);     
    tree->SetBranchAddress("Jet_corr",&Jet_corr); 

    tree->SetBranchAddress("nGenHiggsBoson",&nGenH);
    tree->SetBranchAddress("GenHiggsBoson_pt",&GenH_pt);
    tree->SetBranchAddress("GenHiggsBoson_mass",&GenH_mass);
    tree->SetBranchAddress("GenHiggsBoson_eta",&GenH_eta);
    tree->SetBranchAddress("GenHiggsBoson_phi",&GenH_phi);
    tree->SetBranchAddress("GenHiggsBoson_status",&GenH_status);
    tree->SetBranchAddress("nGenBQuarkFromH",&nGenBfH);
    tree->SetBranchAddress("GenBQuarkFromH_pdgId",&GenBfH_pdgId);
    tree->SetBranchAddress("GenBQuarkFromH_pt",&GenBfH_pt);
    tree->SetBranchAddress("GenBQuarkFromH_eta",&GenBfH_eta);
    tree->SetBranchAddress("GenBQuarkFromH_phi",&GenBfH_phi);
    tree->SetBranchAddress("GenBQuarkFromH_mass",&GenBfH_mass);
    tree->SetBranchAddress("GenBQuarkFromH_charge",&GenBfH_charge);
    tree->SetBranchAddress("GenBQuarkFromH_status",&GenBfH_status);
    tree->SetBranchAddress("nGenBQuarkFromHafterISR",&nGenBfHafterISR);
    tree->SetBranchAddress("GenBQuarkFromHafterISR_pdgId",&GenBfHafterISR_pdgId);
    tree->SetBranchAddress("GenBQuarkFromHafterISR_pt",&GenBfHafterISR_pt);
    tree->SetBranchAddress("GenBQuarkFromHafterISR_eta",&GenBfHafterISR_eta);
    tree->SetBranchAddress("GenBQuarkFromHafterISR_phi",&GenBfHafterISR_phi);
    tree->SetBranchAddress("GenBQuarkFromHafterISR_mass",&GenBfHafterISR_mass);  
  }

  tree->SetBranchAddress("nJet",&nJet);              
  tree->SetBranchAddress("Jet_id",&Jet_id);            
  tree->SetBranchAddress("Jet_puId",&Jet_puId);          
  tree->SetBranchAddress("Jet_btagCSV",&Jet_btagCSV);       
  tree->SetBranchAddress("Jet_btagCMVA",&Jet_btagCMVA);      
  tree->SetBranchAddress("Jet_rawPt",&Jet_rawPt);            
  tree->SetBranchAddress("Jet_pt",&Jet_pt);            
  tree->SetBranchAddress("Jet_eta",&Jet_eta);           
  tree->SetBranchAddress("Jet_phi",&Jet_phi);           
  tree->SetBranchAddress("Jet_mass",&Jet_mass);          
  tree->SetBranchAddress("Jet_btagProb",&Jet_btagProb);      
  tree->SetBranchAddress("Jet_btagBProb",&Jet_btagBProb);     
  tree->SetBranchAddress("Jet_btagnew",&Jet_btagnew);       
  tree->SetBranchAddress("Jet_btagCSVV0",&Jet_btagCSVV0);     
  tree->SetBranchAddress("Jet_chHEF",&Jet_chHEF);         
  tree->SetBranchAddress("Jet_neHEF",&Jet_neHEF);         
  tree->SetBranchAddress("Jet_chEmEF",&Jet_chEmEF);        
  tree->SetBranchAddress("Jet_neEmEF",&Jet_neEmEF);        
  tree->SetBranchAddress("Jet_chMult",&Jet_chMult);        
  tree->SetBranchAddress("Jet_leadTrackPt",&Jet_leadTrackPt);   
  tree->SetBranchAddress("met_pt",&met_pt);            
  tree->SetBranchAddress("met_eta",&met_eta);           
  tree->SetBranchAddress("met_phi",&met_phi);           
  tree->SetBranchAddress("met_mass",&met_mass);       
  //std::cout<<Jet_pt[0]<<"  "<<Jet_pt[1]<<"  "<<Jet_pt[2]<< "   "<<std::endl;

  //Histos options
  h_nPV->Sumw2();
  h_nPV_weighted->Sumw2();

  //output file
  std::string eve = "";
  if(maxEvents > 0) eve = "_" + to_string(maxEvents);
  std::string outfilename=dataFld+"tree_S1_"+sample+eve+".root";
  TFile *outfile=new TFile(outfilename.c_str(), "recreate");
  TTree *outtree=tree->CloneTree(0);
  int H1jet1_i, H1jet2_i;
  int H2jet1_i, H2jet2_i;
  //new branches... not jet implemented
  outtree->Branch("H1jet1_i", &H1jet1_i, "H1jet1_i/I");
  outtree->Branch("H1jet2_i", &H1jet2_i, "H1jet2_i/I");
  outtree->Branch("H2jet1_i", &H2jet1_i, "H2jet1_i/I");
  outtree->Branch("H2jet2_i", &H2jet2_i, "H2jet2_i/I");
  outtree->Branch("Ht", &ht, "ht/F");

  int nEvents;
  if(maxEvents>0) nEvents=maxEvents;
  else nEvents=tree->GetEntries();

  int nCut0=0, nCut1=0, nCut2=0, nCut3=0, nCut4=0, nCut4a=0, nCut5=0, nCut5a=0, nCut5b=0, HHf=0;
  int nJets_InAcc=0;
  double deltaR_HH=-99.,   deltaR_H1=-99.,   deltaR_H2=-99.;
  double deltaPhi_HH=-99., deltaPhi_H1=-99., deltaPhi_H2=-99.; 
  double deltaEta_HH=-99., deltaEta_H1=-99., deltaEta_H2=-99.; 

  int effB1=0, effB2=0;
  int effBH1_wind=0, effBH2_wind=0;

  double N_all [8]={}; //ok
  double Nt_all [8]={};
  double N_all_b [8]={};
  double Nt_all_b [8]={};
  double Nmatched [8]={};
  double Ntr [8]={};

  double nM[binPt][binCSV][binEta];
  double nA[binPt][binCSV][binEta];
  double R[binPt][binCSV][binEta];
  int nRmaxOk = 0, nMCTruth =0, nCSVOk =0, nMatrCSVOk=0;
  double Rave = 9999;;

  for(int ipT=0; ipT<binPt; ipT++){
    for(int iCSV=0; iCSV<binCSV; iCSV++){
      for(int ieta=0; ieta<binEta; ieta++){
        nM[ipT][iCSV][ieta] = 0;
        nA[ipT][iCSV][ieta] = 0;
      }
    }
  }
  //read matrix and calculate R
    ifstream inmatrix;
    std::string fn = utilsFld+matrixFld+"nMnA_"+opt+".asc";
    inmatrix.open(fn);
    if (inmatrix) {
      int i1, i2, i3;
      double sumnA=0, sumnM=0; 
      for(int ipT=0; ipT<binPt; ipT++){
        for(int iCSV=0; iCSV<binCSV; iCSV++){
          for(int ieta=0; ieta<binEta; ieta++){
            inmatrix >> i1 >> i2 >> i3 >> nM[ipT][iCSV][ieta] >> nA[ipT][iCSV][ieta];
            if(i1!=ipT ||i2!=iCSV ||i3!=ieta) cout << "WARNING: matrix mismatch" << endl;
            if(nA[ipT][iCSV][ieta]>0) R[ipT][iCSV][ieta] = nM[ipT][iCSV][ieta]/nA[ipT][iCSV][ieta]; 
            else R[ipT][iCSV][ieta] = 1;
            //cout << nM[ipT][iCSV][ieta] << "  " << nA[ipT][iCSV][ieta] <<endl; 
            sumnM+=nM[ipT][iCSV][ieta];
            sumnA+=nA[ipT][iCSV][ieta];
            if(nM[ipT][iCSV][ieta]!=0 || nM[ipT][iCSV][ieta]!=0)  {
              cout << ipT << " " << iCSV << " " << ieta << " " << R[ipT][iCSV][ieta] << " +- " << R[ipT][iCSV][ieta]*sqrt(1/nM[ipT][iCSV][ieta]+ 1/nA[ipT][iCSV][ieta]) <<endl; 
            }
            else 
              cout << "nM or nA null" <<endl;
	  }
        }
      }
      Rave=sumnM/sumnA;
      inmatrix.close();
    }
    else cout << "Unable to open matrix file " << fn << endl;      

  //------------------  
  // Loop over events
  //------------------
  cout << "nEvents: " << nEvents << endl << endl;
  for (int i=0; i<nEvents; ++i){
    bool jet_inAcc [20];
    vector<Jet> jet;
    bool foundHH=false; 
    ++nCut0;
    tree->GetEvent(i); //READ EVENT
    //std::cout<<Vtype_<<"  " <<Jet_pt[0]<<"  "<<Jet_pt[1]<<"  "<<Jet_pt[2]<< "   "<<std::endl;
    nJets_InAcc = 0;
    h_nJets->Fill(nJet);

    if(isSignalMC && i<=0.5*nEvents) continue; // to use only second part of the sample if Signal MC //to not run on events used to fill the matrix
    else {    
      vector<TLorentzVector> genB_P, genBISR_P;
      for(int p=0; p<4; p++){     
        TLorentzVector b;
        TLorentzVector bISR;
        bISR.SetPtEtaPhiM(GenBfH_pt[p],GenBfH_eta[p],GenBfH_phi[p],GenBfH_mass[p]);
        b.SetPtEtaPhiM(GenBfH_pt[p],GenBfH_eta[p],GenBfH_phi[p],GenBfH_mass[p]);
        genB_P.push_back(b);     
        genBISR_P.push_back(b);     
      }

      //-----------------------
      // 0.Acceptance cuts - pT & |eta| - id puid
      //-----------------------
      for (int j=0; j<nJet; ++j){
        // h_Jets_HT->Fill(Jet_pt[j]); //debug - check if it better after kin cuts..
        // cout << Jet_id[j] << endl;
        if ((fabs(Jet_eta[j])<Jet_eta_cut) && (Jet_pt[j]>Jet_pt_cut_low) ) { //debug //&& Jet_puId[j]>0 && && (Jet_id[j]>0)
          //if(yesTrgA || yesTrgB) {        
            //if (fabs(Jet_eta[j])<2.5) { //only central jets for trigger
            //  ++nJets_InAcc;
            //  jet_inAcc[j] = true;
            //}
            //else jet_inAcc[j] = false;
          //}
          //else {
            ++nJets_InAcc;
            jet_inAcc[j] = true;
          //}
        }
        else jet_inAcc[j] = false;
      }

      //-----------------------
      // 1.cut on trigger - new trigger path
      //-----------------------
      if(yesTrgA || (yesTrgB)){ 
        h_HLT_HH4ball->Fill(HLT_HH4bAll);
//cout << HLT_HH4bAll << endl;
        if(HLT_HH4bAll) ++nCut1;
        else continue;  //break loop if not right trigger path 
      }
 
      h_nJets_InAcc->Fill(nJets_InAcc); //weightPU? - debug
      //----------------------------------------------
      // 3.cut on number of Jets in acceptance region
      //----------------------------------------------
      if(nJets_InAcc>=nJets_cut){ 
        ++nCut3;
        //----------------------------------------
        // fill Jet vector with jet sorted in pT
        //----------------------------------------
	for (int j=0; j<nJet; ++j){    //loop on all the jets (only inAcc are saved)
            //additional info...
           // TLorentzVector jj;
           // jj.SetPtEtaPhiM(Jet_pt[j],Jet_eta[j],Jet_phi[j],Jet_mass[j]); 
           // h_genB1Jets_DR->Fill(jj.DeltaR(gen_B1)); //DR between all jets and B quark at jen level

          if(jet_inAcc[j]){
            Jet je;         
            je.CSV = Jet_btagCSV[j];
            je.CMVA = Jet_btagCMVA[j];
            je.pT  = Jet_pt[j];
            je.mass = Jet_mass[j];
            je.eta = Jet_eta[j];
            je.phi = Jet_phi[j];
            jet.push_back(je);  
                h_Jet4all_pT->Fill(jet[j].pT); //debug
                h_Jet4all_eta->Fill(jet[j].eta); //debug
                h_Jet4all_CSV->Fill(jet[j].CSV);//debug
          }              

        }
        //----------------------------------------------
        // 4. jet sort in CSV + CSV cut on first 3 jets
        //----------------------------------------------
        std::sort (jet.begin(), jet.end(), cmp_CSV );      
        if(jet[0].CSV>0.1){ //cut on jets sorted in cmva //CSV_cut 1
  	  ++nCut4;
           //debug - check sorting --> OK!!
            std::cout<< "size " << jet.size() << " ";                
            for (int j=0; j<jet.size(); ++j){  
              std::cout<<jet[j].CSV << " ";
              //std::cout<<jet[j].pT << " ";
            }
            cout << endl; 	
           //vector<double> c_jet1_CSV, c_jet2_CSV, c_jet3_CSV, c_jet4_CSV;
           //  vector<double> c_jet1_CMVA, c_jet2_CMVA, c_jet3_CMVA, c_jet4_CMVA;
           //  vector<TLorentzVector> c_diJet1, c_diJet2;                  
           //  double Mdiff[25];
           //  double Mdiff_min = -1;
           //  int Ind = 99;

          //all jets vector in acc sorted by csv
          vector<TLorentzVector> jets_P;            
          int NJetInAcc = jet.size();
          for(int l=0; l<NJetInAcc;l++){jets_P.push_back(get_jetVector(&jet[l]));}
          //.............

          // Match first 3 jets with b and fill deltaR
          //------------
          vector<TLorentzVector> genB_P_match;            
          double dR=0;
          bool noJetMatch = false;
          std::vector<TLorentzVector>::iterator index;
          for(int l=0; l<3;l++){
            double dRFin = deltaRCut;
            for(std::vector<TLorentzVector>::iterator it = genB_P.begin() ; it != genB_P.end(); ++it){
              dR = jets_P[l].DeltaR(*it);
              if(dR<dRFin) {
                dRFin=dR;
                index = it;            
              }
            }
            if(dRFin!=deltaRCut){
              genB_P_match.push_back(*index);
              genB_P.erase(index);
              //debug - improve!!
              if(l==0)h_jet1b_dr->Fill(dRFin);
              else if(l==1)h_jet2b_dr->Fill(dRFin);
              else if(l==2)h_jet3b_dr->Fill(dRFin);
            }
            else noJetMatch = true;
          }

          // Fill deltaR for all the fourth jet and find the 'MC Truth' (jet closest to remaining genB)
          //------------
          double dRFinal = deltaRCut;
          int InTrue = 0;
          for(int j=3; j<NJetInAcc; ++j){
            dR = jets_P[j].DeltaR(genB_P[0]);
            h_jet4b_drAll->Fill(dR,1./(NJetInAcc-3));
            if(dR<dRFinal) {
                dRFinal=dR;
                InTrue = j; //indice del jet piu' vicino al partone rimasto           
            }
          }
          if(dRFinal==deltaRCut) noJetMatch = true;
          if(noJetMatch) continue; //skip event if there are no match for all the 4 jets

          h_jet4b_dr->Fill(dRFinal);
          //fill deltaR for all 4th jets not matched
          for(int j=3; j<NJetInAcc; ++j){
            dR = jets_P[j].DeltaR(genB_P[0]);
            if(j!=InTrue) h_jet4b_drNotMatched->Fill(dR,1./(NJetInAcc-4));
          }

          // fill discriminating variables
           //------------
          for(int j=3; j<NJetInAcc; ++j){
              double DpT = jet[2].pT - jet[j].pT;  //jet has same ordering of jets_P ?..
              double DCSV = jet[2].CSV - jet[j].CSV;
              double Deta = fabs(jet[2].eta) - fabs(jet[j].eta);
              if(j==InTrue){
                h_Jet4match_pT->Fill(jet[j].pT);
                h_Jet4match_eta->Fill(jet[j].eta);
                h_Jet4match_CSV->Fill(jet[j].CSV);
                h_Jet4match_DpT3->Fill(DpT);
                h_Jet4match_DCSV3->Fill(DCSV);
                h_Jet4match_Deta3->Fill(Deta);
              }
              else{
             //   h_Jet4all_pT->Fill(jet[j].pT);  debug
//                h_Jet4all_eta->Fill(jet[j].eta);
//                h_Jet4all_CSV->Fill(jet[j].CSV);
                h_Jet4all_DpT3->Fill(DpT);
                h_Jet4all_DCSV3->Fill(DCSV);
                h_Jet4all_Deta3->Fill(Deta);
              }
          }

          // read matrix and assign max R
          //------------
            double Rmax = 0; 
            int In =0;
            for(int j=3; j<NJetInAcc; ++j){
              int ipT = jet[j].pT*(binPt-1)/300;
              if(ipT>binPt-1) ipT = binPt-1;
              int iCSV = jet[j].CSV*binCSV;
              if(iCSV<0) iCSV =0;
              if(iCSV>(binCSV-1)) iCSV = binCSV-1;
              int ieta = fabs(jet[j].eta*2); // //(binEta-1)/2 jet[j].eta*(binEta-1)/2;
              //  if(ieta>binEta-1) ieta = binEta-1;
              //    if(ieta>binEta-1) ieta = binEta-1;
              if(ieta<2) ieta = 0;
              else if(ieta<3) ieta = 1;
              else if(ieta<4) ieta = 2;
              else ieta = 3;
              double thisR = R[ipT][iCSV][ieta];
              if(i<0.5*nEvents+20){ 
//                cout << j << "  " << thisR << "  " << jets_P[j].DeltaR(genB_P[0]) << endl;
//                cout << j << "  " << jet[j].pT << "  " << ipT << "  " << iCSV << "  " << ieta << endl; 
              }
              if(thisR>Rmax) {
                Rmax = thisR;
                In=j;
              }
            }                  
            dR = jets_P[In].DeltaR(genB_P[0]); //matched?!
            h_jet4b_drMatrix->Fill(dR);
            if(In==InTrue)nRmaxOk++;

            // compute di-jets mass and plot it
            // matrix method
            // ----------------------------------
            int k;
	    k = In;
            double M12=0,M13=0,M14=0, M23=0, M24=0, M34=0;
            M12 = (jets_P[0] + jets_P[1]).M();
            M13 = (jets_P[0] + jets_P[2]).M();
            M14 = (jets_P[0] + jets_P[k]).M();
            M23 = (jets_P[1] + jets_P[2]).M();
            M24 = (jets_P[1] + jets_P[k]).M();
            M34 = (jets_P[2] + jets_P[k]).M();
            double DM1, DM2, DM3;
            double mAve_matrix = 9999;
            DM1= fabs(M12 - M34);
            DM2= fabs(M13 - M24);
            DM3= fabs(M14 - M23);
            if(DM1 < DM2 && DM1<DM3) mAve_matrix=(M12+M34)/2.;
            if(DM2 < DM1 && DM2<DM3) mAve_matrix=(M13+M24)/2.;
            if(DM3 < DM1 && DM3<DM2) mAve_matrix=(M14+M23)/2.;
            h_M_matrix->Fill(mAve_matrix);

            // MC truth
            // ----------------------------------
	    k=InTrue;
            nMCTruth++;
            M12 = (jets_P[0] + jets_P[1]).M();
            M13 = (jets_P[0] + jets_P[2]).M();
            M14 = (jets_P[0] + jets_P[k]).M();
            M23 = (jets_P[1] + jets_P[2]).M();
            M24 = (jets_P[1] + jets_P[k]).M();
            M34 = (jets_P[2] + jets_P[k]).M();
            DM1= fabs(M12 - M34);
            DM2= fabs(M13 - M24);
            DM3= fabs(M14 - M23);
            double mAve = 9999;;
            if(DM1 < DM2 && DM1<DM3) mAve=(M12+M34)/2.;
            if(DM2 < DM1 && DM2<DM3) mAve=(M13+M24)/2.;
            if(DM3 < DM1 && DM3<DM2) mAve=(M14+M23)/2.;
            h_M_true->Fill(mAve);

            // 4th CSV
            // ----------------------------------
	    k=3;
            if(3==InTrue)nCSVOk++;
            M12 = (jets_P[0] + jets_P[1]).M();
            M13 = (jets_P[0] + jets_P[2]).M();
            M14 = (jets_P[0] + jets_P[k]).M();
            M23 = (jets_P[1] + jets_P[2]).M();
            M24 = (jets_P[1] + jets_P[k]).M();
            M34 = (jets_P[2] + jets_P[k]).M();
            DM1= fabs(M12 - M34);
            DM2= fabs(M13 - M24);
            DM3= fabs(M14 - M23);
            double mAve_csv = 9999;;
            if(DM1 < DM2 && DM1<DM3) mAve_csv=(M12+M34)/2.;
            if(DM2 < DM1 && DM2<DM3) mAve_csv=(M13+M24)/2.;
            if(DM3 < DM1 && DM3<DM2) mAve_csv=(M14+M23)/2.;
            h_M_csv->Fill(mAve_csv);

           // compute variables and fill all histos
           // ----------------------------------
            TLorentzVector H1, H2, HH;
            H1 = jets_P[0]+jets_P[1];
            //MC truth
            H2 = jets_P[2]+jets_P[InTrue];
            HH = H1+H2;
  	    h_HH_mass_tr->Fill(HH.M());
            //matrix
            H2 = jets_P[2]+jets_P[In];
            HH = H1+H2;
  	    h_HH_mass_matr->Fill(HH.M());
            //CSV
            H2 = jets_P[2]+jets_P[3];
            HH = H1+H2;
  	    h_HH_mass_csv->Fill(HH.M());
       
            foundHH = true;
            HHf++;       
            vector<TLorentzVector> fJets;
            H1 = jets_P[0]+jets_P[1];  //leading in CSV
            H2 = jets_P[2]+jets_P[finalIndex];
            HH = H1+H2;
            fJets.push_back(jets_P[0]);
            fJets.push_back(jets_P[1]);
            fJets.push_back(jets_P[2]);
            fJets.push_back(jets_P[finalIndex]);
            //double H1_CosThSt = computeCosThetaStar(H1,HH);
            //double H2_CosThSt = computeCosThetaStar(H2,HH);

            //delta angles computation:     
            deltaR_HH = H1.DeltaR(H2);  
            deltaR_H1 = fJets[0].DeltaR(fJets[1]);
            deltaR_H2 = fJets[2].DeltaR(fJets[3]);
            deltaPhi_HH = H1.DeltaPhi(H2);  
            deltaPhi_H1 = fJets[0].DeltaPhi(fJets[1]);
            deltaPhi_H2 = fJets[2].DeltaPhi(fJets[3]);
            deltaEta_HH = H1.Eta() - H2.Eta();  
            deltaEta_H1 = fJets[0].Eta() - fJets[1].Eta();
            deltaEta_H2 = fJets[2].Eta() - fJets[3].Eta();

          //-------------------------  
          // Fill histos with final variables
          //-------------------------
	    if (foundHH){	             
  	      double C = 0.;
	      //improve...
	      /*h_Jets_CSV->Fill(c_jet1_CSV[Ind]);
	      h_Jets_CSV->Fill(c_jet2_CSV[Ind]);
	      h_Jets_CSV->Fill(c_jet3_CSV[Ind]);
	      h_Jets_CSV->Fill(c_jet4_CSV[Ind]);
	      h_Jets_CMVA->Fill(c_jet1_CMVA[Ind]);
	      h_Jets_CMVA->Fill(c_jet2_CMVA[Ind]);
	      h_Jets_CMVA->Fill(c_jet3_CMVA[Ind]);
	      h_Jets_CMVA->Fill(c_jet4_CMVA[Ind]);*/
	      //
	      //bool efB1 = false;
	      //bool efB2 = false;
              //----------------
              // Loop over 4 selected jets
              //----------------
	      for(std::vector<TLorentzVector>::iterator it = fJets.begin() ; it != fJets.end(); ++it){
	        h_Jets_mass->Fill((*it).M());
 	        h_Jets_pT->Fill((*it).Pt());
	        h_Jets_eta->Fill((*it).Eta()); 
	        h_Jets_Phi->Fill((*it).Phi()); 
	        C += ((*it).Pt()/(*it).E());
    	      //count on deltaR to get efficiency.. 
	        /*if(fabs(gen_B1.Eta())<2.5 && gen_B1.Pt()>20 && fabs(gen_B2.Eta())<2.5 && gen_B2.Pt()>20){

	          if((*it).DeltaR(gen_B1)<0.5 && !efB1){
	            efB1 = true;
	            effB1++;
	          }
	          if((*it).DeltaR(gen_B2)<0.5 && !efB2){
	            efB2 = true;
	            effB2++;
	          }
	          h_genBfH1_DR->Fill((*it).DeltaR(gen_B1));
	          h_genBfH2_DR->Fill((*it).DeltaR(gen_B2));      
	        }*/
	      }
	      h_Jets_Centr->Fill(C/4); //sum over 4 jets

	      h_Jet1_pT->Fill(fJets[0].Pt());
	      h_Jet2_pT->Fill(fJets[1].Pt());
	      h_Jet3_pT->Fill(fJets[2].Pt());
	      h_Jet4_pT->Fill(fJets[3].Pt());
	      h_Jet1_Eta->Fill(fJets[0].Eta());
	      h_Jet2_Eta->Fill(fJets[1].Eta());
	      h_Jet3_Eta->Fill(fJets[2].Eta());
	      h_Jet4_Eta->Fill(fJets[3].Eta());
	     /* h_Jet1_CSV->Fill(c_jet1_CSV[Ind]);
	      h_Jet2_CSV->Fill(c_jet2_CSV[Ind]);
	      h_Jet3_CSV->Fill(c_jet3_CSV[Ind]);
	      h_Jet4_CSV->Fill(c_jet4_CSV[Ind]);
	      h_Jet1_CMVA->Fill(c_jet1_CMVA[Ind]);
	      h_Jet2_CMVA->Fill(c_jet2_CMVA[Ind]);
	      h_Jet3_CMVA->Fill(c_jet3_CMVA[Ind]);
	      h_Jet4_CMVA->Fill(c_jet4_CMVA[Ind]);
	
	      double Jets_avgCSV = (c_jet1_CSV[Ind]+c_jet2_CSV[Ind]+c_jet3_CSV[Ind])/3;    
	      double Jets_minCSV = std::min({c_jet1_CSV[Ind],c_jet2_CSV[Ind],c_jet3_CSV[Ind]});
	      h_3Jets_avgCSV->Fill(Jets_avgCSV);
	      h_3Jets_minCSV->Fill(Jets_minCSV);
	*/
	      h_H1_mass->Fill(H1.M());
	      h_H1_pT->Fill(H1.Pt());
	      h_H1_Eta->Fill(H1.Eta());
	      h_H1_Phi->Fill(H1.Phi());
	      //h_H1_CosThSt->Fill(fabs(H1_CosThSt));
	      h_H1_deltaR->Fill(deltaR_H1);
	      h_H1_deltaPhi->Fill(deltaPhi_H1);
	      h_H1_deltaEta->Fill(deltaEta_H1);
	      h_H1_deltaPhiVSpT->Fill(deltaPhi_H1,H1.Pt());
	      h_H2_mass->Fill(H2.M());
	      h_H2_pT->Fill(H2.Pt());
	      //h_H2_CosThSt->Fill(fabs(H2_CosThSt));
	      h_H2_deltaR->Fill(deltaR_H2);
	      h_H2_deltaPhi->Fill(deltaPhi_H2);
	      h_H2_deltaEta->Fill(deltaEta_H2);
	      h_H2_deltaPhiVSpT->Fill(deltaPhi_H2,H2.Pt());
	      h_H_mass->Fill(H1.M());
	      h_H_mass->Fill(H2.M());
	      h_H_pT->Fill(H1.Pt());
	      h_H_pT->Fill(H2.M());
	      h_H_Eta->Fill(H1.Eta());
	      h_H_Eta->Fill(H2.Eta());
	      h_H_Phi->Fill(H1.Phi());
	      h_H_Phi->Fill(H2.Phi());
	      h_HH_mass->Fill((H1+H2).M());
	      h_HH_pT->Fill((H1+H2).Pt());
	      h_HH_Eta->Fill((H1+H2).Eta());
	      h_HH_Phi->Fill((H1+H2).Phi());
              h_H1_H2_mass->Fill(H1.M(), H2.M());  //leading H first
	      h_HH_deltaR->Fill(deltaR_HH);
	      h_HH_deltaPhi->Fill(deltaPhi_HH);
	      h_HH_deltaEta->Fill(deltaEta_HH);
	      //int region=withinRegion(H1_P.M(), H2_P.M(), 17.5, 37.5, 125, 125);
	
	      //h_MET->Fill(met.E());
         
	    //----------------------------------- 
	    // WINDOWS : cut on dijet masses (rectangle)
	    //-----------------------------------
	     /* if(dijetM_cut_low<c_diJet1[Ind].M() && c_diJet1[Ind].M()<dijetM_cut_high){
	        nCut5a++;
	        if(dijetM_cut_low<c_diJet2[Ind].M() && c_diJet2[Ind].M()<dijetM_cut_high){ 
	          nCut5++;
	          //  cout << "inWindows " << nCut5 << " " << c_diJet1[Ind].M() << " " << c_diJet2[Ind].M() << endl;
	          h_HH_massInReg->Fill((H1+H2).M());
	          h_HH_pTInReg->Fill((H1+H2).Pt());
	          h_H1_H2_massInReg->Fill(H1.M(), H2.M());

  	        // count on deltaR to get efficiency.. 
        	  if(fabs(gen_B1.Eta())<2.5 && gen_B1.Pt()>20 && fabs(gen_B2.Eta())<2.5 && gen_B2.Pt()>20){
	            if(fJets[0].DeltaR(gen_B1)<0.5 || fJets[0].DeltaR(gen_B2)<0.5) effBH1_wind++;
                    else if(fJets[1].DeltaR(gen_B1)<0.5 || fJets[1].DeltaR(gen_B2)<0.5) effBH1_wind++;
         	    if(fJets[2].DeltaR(gen_B1)<0.5 || fJets[2].DeltaR(gen_B2)<0.5) effBH2_wind++;
	            else if(fJets[3].DeltaR(gen_B1)<0.5 || fJets[3].DeltaR(gen_B2)<0.5) effBH2_wind++;
	          }
	        }
	      } // windows
              */
  	      //ELIPTICAL WINDOWS -- TO BE UNDERSTOOD
              if(withinElipse(H1.M(), H2.M())){
	        nCut5b++;
                h_H1_H2_massInReg2->Fill(H1.M(), H2.M());
              }
              else h_H1_H2_massInReg3->Fill(H1.M(), H2.M());
	      // windows

	      outtree->Fill();	

            }//if FOUND
        } // CSV cut
      } // nJets_InAcc
    jet.clear();
    } // sample splitting
  } // Event loop

  cout << endl;
  cout << "'ACCEPTANCE' (Ntrue/Nevents): " << (double)nMCTruth/nEvents*100 << endl;
  cout << "N true jets " << Ntr[0] << endl;
  cout << "nMCTruth " << nMCTruth << endl;
  cout << endl;

 if(draw){
//DRAW...
  TCanvas * c001 = new TCanvas("c001","scatter",600,500);
  c001->cd();
  h_H1_H2_mass->Draw("colz");

  TCanvas * c000 = new TCanvas("c000","#DeltaR",600,500);
  c000->cd();
cout << h_jet4b_drMatrix->Integral() << " " << h_jet4b_dr->Integral() << endl;
  h_jet4b_dr->Draw();
  h_jet4b_drMatrix->SetLineColor(kRed);
  h_jet4b_drMatrix->Draw("same");
//  h_jet4b_dr->Scale(h_jet4b_drMatrix->Integral()/h_jet4b_dr->Integral());
//  h_jet4b_drMatrix->Scale(1/h_jet4b_drMatrix->Integral());
//  h_jet4b_drAll->Scale(h_jet4b_drMatrix->Integral()/h_jet4b_drAll->Integral());
  h_jet4b_drNotMatched ->SetLineColor(kBlue);
  h_jet4b_drNotMatched ->Draw("same");
//  h_jet4b_drAll->SetLineColor(kBlue);
//  h_jet4b_drAll->Draw("same");


  TCanvas * c00 = new TCanvas("c00","nJets",500,500);  //not normalized . Okay
  c00->cd();
  h_nJets_InAcc->Scale(1/h_nJets_InAcc->Integral());
  h_nJets_InAcc->GetYaxis()->SetTitle("1/Nevents");
  h_nJets_InAcc->SetLineColor(kRed);
  h_nJets->Scale(h_nJets_InAcc->Integral()/h_nJets->Integral());
  h_nJets->GetYaxis()->SetTitle("1/Nevents");
  h_nJets_InAcc->Draw();
  h_nJets->Draw("same");


  TCanvas * c0 = new TCanvas("c0","<m>",600,500);
  c0->cd();
  h_M_true->Rebin(3);
  h_M_true->Draw();
  h_M_matrix->Rebin(3);
  h_M_matrix->SetLineColor(kRed);
  h_M_matrix->Draw("same");
  h_M_csv->Rebin(3);
  h_M_csv->SetLineColor(kBlue);
  h_M_csv->Draw("same");
//  h_M_matrixCSV->SetLineColor(kGray);
//  h_M_matrixCSV->Draw("same");

  TCanvas * c = new TCanvas("c","nome",600,500);
  c->Divide(2,3);
  c->cd(1);
  h_Jet4all_pT->Scale(h_Jet4match_pT->Integral()/h_Jet4all_pT->Integral());
  h_Jet4all_pT->Draw();
  h_Jet4match_pT->SetLineColor(kRed);
  h_Jet4match_pT->Draw("same");
  c->cd(2);
  h_Jet4match_eta->SetLineColor(kRed);
  h_Jet4match_eta->Draw();
  h_Jet4all_eta->Scale(h_Jet4match_eta->Integral()/h_Jet4all_eta->Integral());
  h_Jet4all_eta->Draw("same");
  c->cd(3);
  h_Jet4all_CSV->Scale(h_Jet4match_CSV->Integral()/h_Jet4all_CSV->Integral());
  h_Jet4all_CSV->Draw();
  h_Jet4match_CSV->SetLineColor(kRed);
  h_Jet4match_CSV->Draw("same");
  c->cd(4);
  h_Jet4all_DpT3->Scale(h_Jet4match_DpT3->Integral()/h_Jet4all_DpT3->Integral());
  h_Jet4all_DpT3->Draw();
  h_Jet4match_DpT3->SetLineColor(kRed);
  h_Jet4match_DpT3->Draw("same");
  c->cd(5);
  h_Jet4match_DCSV3->SetLineColor(kRed);
  h_Jet4match_DCSV3->Draw();
  h_Jet4all_DCSV3->Scale(h_Jet4match_DCSV3->Integral()/h_Jet4all_DCSV3->Integral());
  h_Jet4all_DCSV3->Draw("same");

  c->cd(6);
  h_Jet4match_Deta3->SetLineColor(kRed);
  h_Jet4match_Deta3->Draw();
  h_Jet4all_Deta3->Scale(h_Jet4match_Deta3->Integral()/h_Jet4all_Deta3->Integral());
  h_Jet4all_Deta3->Draw("same");

 TCanvas * c1 = new TCanvas("c1","discr",1000,250);
  c1->Divide(3,1);
  c1->cd(1);
  h_Jet4all_pT->Scale(h_Jet4match_pT->Integral()/h_Jet4all_pT->Integral());
  h_Jet4all_pT->Draw();
  h_Jet4match_pT->SetLineColor(kRed);
  h_Jet4match_pT->Draw("same");
  c1->cd(2);
  h_Jet4match_eta->SetLineColor(kRed);
  h_Jet4match_eta->Draw();
  h_Jet4all_eta->Scale(h_Jet4match_eta->Integral()/h_Jet4all_eta->Integral());
  h_Jet4all_eta->Draw("same");
  c1->cd(3);
  h_Jet4all_CSV->Scale(h_Jet4match_CSV->Integral()/h_Jet4all_CSV->Integral());
  h_Jet4all_CSV->Draw();
  h_Jet4match_CSV->SetLineColor(kRed);
  h_Jet4match_CSV->Draw("same");

 TCanvas * c2 = new TCanvas("c2","Mhh",600,500);
  c2->cd();
  h_HH_mass_tr->Draw();
  h_HH_mass_matr->SetLineColor(kRed);
  h_HH_mass_matr->Draw("same");
  h_HH_mass_csv->SetLineColor(kBlue);
  h_HH_mass_csv->Draw("same");
 }
 
  h_Cuts->Fill(1,nCut0);
  h_Cuts->Fill(3,nCut1);
  h_Cuts->Fill(5,nCut2);
  h_Cuts->Fill(7,nCut3);
  h_Cuts->Fill(9,nCut4);  
  h_Cuts->Fill(11,nCut4a);  
  h_Cuts->Fill(13,HHf);
  h_Cuts->Fill(15,nCut5a);
  h_Cuts->Fill(17,nCut5);

  outtree->Write();
  outfile->Close();

  std::string histfilename= plotsFld+"Histograms_"+sample+"_"+opt+".root";
  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");

//new histos
  h_jet1b_dr->Write();
  h_jet2b_dr->Write();
  h_jet3b_dr->Write();
  h_jet4b_dr->Write();
////

  h_nJets->Write();
  h_nJets_InAcc->Write();

  h_nPV->Write();
  h_nPV_weighted->Write();
  h_HLT_HH4b->Write();
  h_nGenH->Write();
  h_genH_pT->Write();
  h_genH_mass->Write();
  h_genHH_mass->Write();
  h_genH1_mass->Write();
  h_genH2_mass->Write();
  h_Jets_HT->Write();
  h_nJets_InAcc->Write();
  h_Jet1_pT->Write();
  h_Jet2_pT->Write();
  h_Jet3_pT->Write();
  h_Jet4_pT->Write();
  h_Jet1_Eta->Write();
  h_Jet2_Eta->Write();
  h_Jet3_Eta->Write();
  h_Jet4_Eta->Write();
  h_Jet1_CSV->Write();
  h_Jet2_CSV->Write();
  h_Jet3_CSV->Write();
  h_Jet4_CSV->Write();
  h_Jet1_CMVA->Write();
  h_Jet2_CMVA->Write();
  h_Jet3_CMVA->Write();
  h_Jet4_CMVA->Write();
  h_3Jets_avgCSV->Write();
  h_3Jets_minCSV->Write();
  h_nJets_90->Write();
  h_Jets_mass->Write();
  h_Jets_pT->Write();
  h_Jets_eta->Write();
  h_Jets_Phi->Write();
  h_Jets_Centr->Write();
  h_Jets_CSV->Write();
  h_Jets_CMVA->Write();
  h_H_mass->Write();
  h_H_pT->Write();
  h_H_Eta->Write();
  h_H_Phi->Write();
  h_H1_mass->Write();
  h_H1_pT->Write();
  h_H1_Eta->Write();
  h_H1_Phi->Write();
  h_H1_CosThSt->Write();
  h_H1_deltaR->Write();
  h_H1_deltaPhi->Write();
  h_H1_deltaEta->Write();
  h_H1_deltaPhiVSpT->Write();
  h_H2_mass->Write();
  h_H2_pT->Write();
  h_H2_CosThSt->Write();
  h_H2_deltaR->Write();
  h_H2_deltaPhi->Write();
  h_H2_deltaEta->Write();
  h_H2_deltaPhiVSpT->Write();
  h_HH_mass->Write();
  h_HH_pT->Write();
  h_HH_Eta->Write();
  h_HH_Phi->Write();
  h_HH_deltaR->Write();
  h_HH_deltaPhi->Write();
  h_HH_deltaEta->Write();
  h_H1_H2_mass->Write();
  h_HH_massInReg->Write();
  h_HH_pTInReg->Write();
  h_H1_H2_massInReg->Write();
  h_H1_H2_massInReg2->Write();
  h_H1_H2_massInReg3->Write();
  h_genB1Jets_DR->Write();
  h_genBfH1_DR->Write();
  h_genBfH2_DR->Write();
  h_MET->Write();

  h_Cuts->Write();
  tFile->Write();
  tFile->Close();

//-----------------
// Some output
//-----------------
  std::cout<<std::setprecision(2);
  std::cout<<std::fixed;  
  std::cout<<"Wrote output file "<<histfilename<<std::endl;
  std::cout<<"0.Number of events from input file = "<< nCut0 << std::endl;
  std::cout<<"1.Number of events after trigger   = "<< nCut1 << " ,  " <<(float)nCut1/nCut0*100<<"%"<<std::endl;
  std::cout<<"2.Number of events after Vtype     = "<< nCut2 << " ,  " <<(float)nCut2/nCut0*100<<"%"<<std::endl;
  std::cout<<"3.Number of events after nJets_InAcc  = "<< nCut3 << " ,  " <<(float)nCut3/nCut0*100<<"%"<<std::endl;
  std::cout<<"4.Number of events after bTag cut  = "<< nCut4 << " ,  " <<(float)nCut4/nCut0*100<<"%"<<std::endl;
  std::cout<<"Number of events after HHfound     = "<< HHf << " ,  " <<(float)HHf/nCut0*100<<"%"<<std::endl;
  std::cout<<"Number of events with only diJet1 in mass windows  = "<< (nCut5a-nCut5) << " ,  " <<(float)(nCut5a-nCut5)/nCut0*100<<"%"<<std::endl;
  std::cout<<"Number of events in mass windows   = "<< nCut5 << " ,  " <<(float)nCut5/nCut0*100<<"%"<<std::endl;
  std::cout<<"Number of events in mass windows 2 = "<< nCut5b << " ,  " <<(float)nCut5b/nCut0*100<<"%"<<std::endl;
  std::cout<<endl;
  std::cout<<"eff B1 for HHfound = "<<effB1<< "  " << (float)effB1/HHf*100 <<"%"<<std::endl;
  std::cout<<"eff B2 for HHfound = "<<effB2<< "  " << (float)effB2/HHf*100 <<"%"<<std::endl;
  std::cout<<"eff Bh1 in windows  = "<<effBH1_wind<< "  " << (float)effBH1_wind/nCut5*100 <<"%"<<std::endl;
  std::cout<<"eff Bh2 in windows  = "<<effBH2_wind<< "  " << (float)effBH2_wind/nCut5*100 <<"%"<<std::endl;

  return;
}

