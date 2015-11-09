// Code for non resonant HH->bbbb in CMS Run2 
// Step2 - TMVA for S/B discrimination
//  Author: Martino Dall'Osso
//   thanks Carlo Alberto Gottardo
//    Oct 19th, 2015
//     g++ `root-config --libs --cflags` -lTMVA hh4b_13TeV_tmva.cc -o doTmva
//      ./doTmva BTagCSVRun2015C-D M-260 pT20CSVL_def 0 1
//----------------------------------------------------------------

#define hh4b_tmva_C

#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include <cmath>
#include <string>

//using namespace std;
// NOTES
//***********************************************
//  weightfile 	= names for .C and .xml weights file namely "window", "all"
//  outTMVAfileName = names for .root training and test files "TMVAinWindowTraining.root", "TMVAallTraining.root"
//  inputSignal = training() signal sample
//  inputBkg = training() bkg sample
//  training() trains both methods and produces two C and xml weights files
//	applications: one TH1F per method and per sample choice
//  merits: one TH1F per method and per application choice
//************************************************

static const std::string frameworkVersionFld = "V15/";   //equal to the Data version
static const int numVars = 9;   //9

static const std::string weightfile = "all";
static const std::string dataFld="../data/"+frameworkVersionFld;;
static const std::string resultsFld="results/"+frameworkVersionFld;;

//#include "../utils/hh4bStructs.h"
#include "hh4b_13TeV_tmva.h"

//-------------
int main (int argc, char **argv) {  //DEBUG!!!!
 
  std::string sample_, opt_, MC_;
  int appOn_ = 1;   // 0 == MC, 1 == Data , 2 == Data_CR
  int wtd_ = 1;
  
  for(int i=1; i<argc; i++) {
    if(i==1) sample_ = argv[i];
    if(i==2) MC_ = argv[i];
    if(i==3) opt_ = argv[i];
    if(i==4) appOn_ = std::atoi(argv[i]);
    if(i==5) wtd_ = std::atoi(argv[i]);
  }

  hh4b_tmva hhtmva (sample_,MC_, opt_, appOn_, wtd_);
  hhtmva.doTmva();
  return 0;
}


hh4b_tmva::hh4b_tmva(std::string sam, std::string MC, std::string op, int appO_, int whatToD){
  datasample=sam;
  MCsample=MC;
  opt=op;
  appOn=appO_; 
  whatToDo=whatToD;
}

//--------------------
//  CORE FUNCTION
//--------------------
void hh4b_tmva::doTmva(){
  //init filenames:
  std::string inputSignal = dataFld+"tree_Step1_"+MCsample+"_"+opt+".root";  //debug - same sample for train and appl
  std::string inputBkg = dataFld+"tree_Step1_"+datasample+"_"+opt+".root";  //debug - same sample for train and appl
  std::string outfileTrain = resultsFld+"TMVATrain_"+datasample+"_"+opt+".root";
  std::string outfileApp = resultsFld+"TMVATest_"+datasample+"_"+opt+".root";
  std::string infileApp;

  if(whatToDo == 1 || whatToDo == 0){
    //_______TRAINING____________ _______________________________________________________________
    std::cout << "Training all" << std::endl; 
    training(weightfile, outfileTrain, inputSignal, inputBkg);
  }
  if(whatToDo == 2 || whatToDo == 0){   //sistemare!!!!
    //_______FIGURES OF MERIT____________________________________________________________________
    float dummy1, dummy2;
    std::cout << "calculating merits " << std::endl;
    TFile *f = new TFile(outfileApp.c_str(), "RECREATE");
    if(CutScanL(outfileTrain) && CutScanBDT(outfileTrain)){  //execute functions.
      dummy1 = merit_L->GetMaximum(); dummy2 = merit_L->GetMaximumBin();
      std::cout << "All Likelihood \n max Q: " <<dummy1<< " at " << dummy2 << std::endl;
      dummy1 = merit_BDT->GetMaximum(); dummy2 = merit_BDT->GetMaximumBin();
      std::cout << "All window Likelihood" << "\n" << "max Q: " << dummy1 << " at " << dummy2 << std::endl;
      merit_L->SetName("TrainA_merit_L");
      merit_BDT->SetName("TrainA_merit_BDT");
      f->cd();
      merit_L->Write();
      merit_BDT->Write();
    }
    f->Close(); //debug
    return;
  }
  if(whatToDo == 3 || whatToDo == 0){
    //_______APPLICATION__________________________________________________________________________
    if(appOn == 1){
      infileApp = dataFld+"tree_Step1_"+datasample+"_"+opt+".root";   //debug - same sample for train and appl
      outNtuples = dataFld+"tree_Step2_"+datasample+"_"+opt+".root";
      std::cout << "APPLYING ON DATA " << infileApp << std::endl;
    }
    else if(appOn == 2){
      infileApp = dataFld+"tree_Step1_"+datasample+"_"+opt+"_CR.root";   //debug - same sample for train and appl
      outNtuples = dataFld+"tree_Step2_"+datasample+"_"+opt+"_CR.root";
      std::cout << "APPLYING ON DATA_CR " << infileApp << std::endl;
    }
    else if(appOn == 0){
      infileApp = dataFld+"tree_Step1_"+MCsample+"_"+opt+".root";   //debug - same sample for train and appl
      outNtuples = dataFld+"tree_Step2_"+MCsample+"_"+opt+".root";
      std::cout << "APPLYING ON MC " << infileApp << std::endl;
    }
    else { std::cout << "NO SAMPLE TO APPLY ON" << std::endl; return;}

    TFile *f = new TFile(outfileApp.c_str(), "UPDATE"); //debug
    if(!Lapplication(infileApp, weightfile)) return;
    else {
      std::cout << std::endl << " #######################  Lapplication done ################### " << std::endl << std::endl;
      if(!BDTapplication(infileApp, weightfile)) return;
      else {  //execute functions.
        std::cout << std::endl << " #######################  BDTapplication done ################### " << std::endl << std::endl;
        like_out->SetName("CR_MC_Likel");
        bdt_out->SetName("CR_MC_BDT");
        f->cd();
        like_out->Write();
        bdt_out->Write();
      } 
    }
    f->Close(); //debug
    return;
/*
infileApp = "../Samples/CR_NOCUT_hh4b-nores-8Tev300k_step2_L1y1.root";
TH1D * S2_MCLikelihood = Lapplication(infileApp, weightfile);
TH1D * S2_MCBDT = BDTapplication(infileApp, weightfile);
S2_MCLikelihood->SetName("MC_Likel");
S2_MCBDT->SetName("MC_BDT");
f->cd();
S2_MCLikelihood->Write();
S2_MCBDT->Write();

infileApp = "../Samples/NOCUT_DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1.root";
TH1D * B_Likelihood = Lapplication(infileApp, weightfile);
TH1D * B_BDT = BDTapplication(infileApp, weightfile);
B_Likelihood->SetName("BKG_Likel");
B_BDT->SetName("B_BKG_BDT");
f->cd();
B_Likelihood->Write();
B_BDT->Write();

infileApp = "../Samples/NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1.root";
TH1D * C_Likelihood = Lapplication(infileApp, weightfile);
TH1D * C_BDT = BDTapplication(infileApp, weightfile);
C_Likelihood->SetName("BKG_Likel");
C_BDT->SetName("C_BKG_BDT");
f->cd();
C_Likelihood->Write();
C_BDT->Write();

infileApp = "../Samples/NOCUT_DiJetPt_BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1.root";
TH1D * CP_Likelihood = Lapplication(infileApp, weightfile);
TH1D * CP_BDT = BDTapplication(infileApp, weightfile);
CP_Likelihood->SetName("CP_BKG_Likel");
CP_BDT->SetName("CP_BKG_BDT");
f->cd();
CP_Likelihood->Write();
CP_BDT->Write();

infileApp = "../Samples/NOCUT_DiJetPt_BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1.root";
TH1D * D_Likelihood = Lapplication(infileApp, weightfile);
TH1D * D_BDT = BDTapplication(infileApp, weightfile);
D_Likelihood->SetName("D_BKG_Likel");
D_BDT->SetName("D_BKG_BDT");
f->cd();
D_Likelihood->Write();
D_BDT->Write();

infileApp = "../Samples/CR_NOCUT_DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1.root";
TH1D * CR_B_Likelihood = Lapplication(infileApp, weightfile);
TH1D * CR_B_BDT = BDTapplication(infileApp, weightfile);
CR_B_Likelihood->SetName("CR_B_Likel");
CR_B_BDT->SetName("CR_B_BDT");
f->cd();
CR_B_Likelihood->Write();
CR_B_BDT->Write();

infileApp = "../Samples/CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1.root";
TH1D * CR_C_Likelihood = Lapplication(infileApp, weightfile);
TH1D * CR_C_BDT = BDTapplication(infileApp, weightfile);
CR_C_Likelihood->SetName("CR_C_Likel");
CR_C_BDT->SetName("CR_C_BDT");
f->cd();
CR_C_Likelihood->Write();
CR_C_BDT->Write();

infileApp = "../Samples/CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1.root";
TH1D * CR_CP_Likelihood = Lapplication(infileApp, weightfile);
TH1D * CR_CP_BDT = BDTapplication(infileApp, weightfile);
CR_CP_Likelihood->SetName("CR_CP_Likel");
CR_CP_BDT->SetName("CR_CP_BDT");
f->cd();
CR_CP_Likelihood->Write();
CR_CP_BDT->Write();

infileApp = "../Samples/CR_NOCUT_DiJetPt_BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1.root";
TH1D * CR_D_Likelihood = Lapplication(infileApp, weightfile);
TH1D * CR_D_BDT = BDTapplication(infileApp, weightfile);
CR_D_Likelihood->SetName("CR_D_Likel");
CR_D_BDT->SetName("CR_D_BDT");
f->cd();
CR_D_Likelihood->Write();
CR_D_BDT->Write();*/
    }
  }
  hh4b_tmva::~hh4b_tmva(){
}
//---------------

//=====================================TRAINING=FUNCTION=====================================
void hh4b_tmva::training(std::string weightfile, std::string outTMVAfileName, std::string inputSignal, std::string inputBkg) {

  TFile* outputFile = new TFile( outTMVAfileName.c_str(), "RECREATE" );
  TMVA::Factory *factory = new TMVA::Factory(weightfile.c_str(), outputFile, "!V:AnalysisType=Classification");
  TFile *sign = TFile::Open(inputSignal.c_str());
  if(!sign) return;
  TFile *bkg = TFile::Open(inputBkg.c_str());
  if(!bkg) return;

  factory->AddSignalTree((TTree*)sign->Get("tree"),1.0);
  factory->AddBackgroundTree((TTree*)bkg->Get("tree"),1.0);

  std::string varible_name;
  for(int i=0; i<numVars; ++i){ //debug
    //	cos*	dj_1 R   dj_2 R   dj_1 pt dj_2_pt	 Qcent    HHM   QPt_3   TDJ deltaR
    //if(i == 41| i== 27 | i== 32 | i==24 |  i==29 | i== 10 | i==37 | i==13 | i==36)
    //{
    factory->AddVariable(setTitle(i).c_str(), 'F');
     //}
  } 	
  //debug -- check -- CHANGE ME
  factory->PrepareTrainingAndTestTree("", "!V:nTrain_Signal=10000:nTrain_Background=10000:nTest_Signal=5000:nTest_Background=5000:SplitMode=Random");
  factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood"); //debug -- reduce correlation.. , "VarTransform=N+P(H1_pT,H1_dR)+P(H2_pT,H2_dR)"   , 
  factory->BookMethod(TMVA::Types::kBDT, "BDT_GiniIndex", "NTrees=1200:MaxDepth=3:SeparationType=GiniIndex:AdaBoostR2Loss=Quadratic"); //:VarTransform=N+P(H1_pT,H1_dR)+P(H2_pT,H2_dR)
  //--------------
  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->Close();
  delete factory;
  sign->Close();
  bkg->Close();
  return;
}

//=====================================BDT=APPLICATION=====================================
bool hh4b_tmva::BDTapplication(std::string appfile, std::string weightfile) {

  std::cout << "APPLICOOOO A " << appfile << std::endl;
  // Ouput of a new TFile which is a clone of the previous but has the new BDT branch
  //copy TTree and add new branch

  float Bout_orig = 0.0;
  float Bout = 0.0;
  TFile *old = TFile::Open(appfile.c_str(), "READ");
  if(!old) return 0;
  TTree *oldtree = (TTree*)old->Get("tree");
  oldtree->SetBranchStatus("*",1);
  //Create a new file + a clone of old tree in new file
  TFile *newtfile = new TFile(outNtuples.c_str(),"recreate");
  TTree *newtree = oldtree->CloneTree();
  TBranch *br = newtree->Branch("BDT", &Bout, "Bout/F");
	
  TMVA::Reader *reader = new TMVA::Reader("!Color");
  float var[46];
  TFile *input = TFile::Open(appfile.c_str());
  if(!input) return 0;
  TTree *theTree = (TTree*)input->Get("tree");

  std::string varible_name;
  for(int i=0; i<numVars; ++i){
    //cos*	dj_1 R   dj_2 R   dj_1 pt dj_2_pt	 Qcent    HHM   QPt_3   TDJ deltaR
    //if(i == 41| i== 27 | i== 32 | i==24 |  i==29 | i== 10 | i==37 | i==13 | i==36)
    //{
    reader->AddVariable(setTitle(i).c_str(), &var[i]);
    //}
  }
  std::string readweighfile = "weights/"+weightfile+"_BDT_GiniIndex.weights.xml";
  reader->BookMVA("BDT_GiniIndex", readweighfile.c_str() );
  for(int i=0; i<numVars; ++i){
    theTree->SetBranchAddress(setTitle(i).c_str(), &var[i]);
  }	
  for(long i=0; i<theTree->GetEntries(); i++) {
    theTree->GetEntry(i);
    Bout_orig = reader->EvaluateMVA( "BDT_GiniIndex" );
    Bout = (Bout_orig+1.0)/2.0;
    bdt_out->Fill(Bout);
    br->Fill();
  }
  //newtree->SetEntries(br->GetEntries());
  newtree->Print();
  newtfile->Write();
  return 1;
}

//=====================================LIKELIHOOD=APPLICATION=====================================
bool hh4b_tmva::Lapplication(std::string appfile, std::string weightfile) {

  TMVA::Reader *reader = new TMVA::Reader("!Color");
  float var[46];
  TFile *input = TFile::Open(appfile.c_str());
  if(!input) return 0;
  TTree *theTree = (TTree*)input->Get("tree");
  std::string varible_name;
  for(int i=0; i<numVars; ++i){
    //cos*	dj_1 R   dj_2 R   dj_1 pt dj_2_pt	 Qcent    HHM   QPt_3   TDJ deltaR
    //if(i == 41| i== 27 | i== 32 | i==24 |  i==29 | i== 10 | i==37 | i==13 | i==36)
    //{
    reader->AddVariable(setTitle(i).c_str(), &var[i]);
    //}
  }
  std::string readweighfile = "weights/"+weightfile+"_Likelihood.weights.xml";
  reader->BookMVA("Likelihood", readweighfile.c_str() );
  for(int i=0; i<numVars; ++i){
   //if(i==21 | i==22) continue;
   theTree->SetBranchAddress(setTitle(i).c_str(), &var[i]);
  }	
  float Lout;
  for(long i=0; i<theTree->GetEntries(); i++) {
    theTree->GetEntry(i);
    Lout = reader->EvaluateMVA( "Likelihood" );
    like_out->Fill(Lout);
  }
  return 1;
}

//=====================================CUT=SCAN=L============================================
bool hh4b_tmva::CutScanL(std::string outTMVAfile) {

  TFile *MVAtrain = new TFile(outTMVAfile.c_str(),"READ");
  if(!MVAtrain) return 0;
  TH1F *Likelihood_S_original = (TH1F*)MVAtrain->Get("Method_Likelihood/Likelihood/MVA_Likelihood_S");
  TH1F *Likelihood_B_original = (TH1F*)MVAtrain->Get("Method_Likelihood/Likelihood/MVA_Likelihood_B");
  TH1F *Likelihood_S = new TH1F("L_S", "Likelihood Signal", 40, 0, 1);
  TH1F *Likelihood_B = new TH1F("L_B", "Likelihood Background", 40, 0, 1);

  float tempS = 0;float tempB = 0;
  for(int i=0; i<=Likelihood_B_original->GetXaxis()->GetNbins(); i++){
    tempB = Likelihood_B_original->GetBinContent(i);
    Likelihood_B->SetBinContent(i, tempB);
  }
  float integral = Likelihood_B->Integral();
  Likelihood_B->Scale(100./integral);
  for(int i=0; i<=Likelihood_S_original->GetXaxis()->GetNbins(); i++){
    tempS = Likelihood_S_original->GetBinContent(i);
    Likelihood_S->SetBinContent(i, tempS);
  }
  integral = Likelihood_S->Integral();
  Likelihood_S->Scale(100./integral);
  float cut = 0.0;
  float QL, sl, bl;
  while(cut <= 1) {
    sl = Likelihood_S->Integral(Likelihood_S->FindBin(cut),Likelihood_S->FindBin(1.0));
    bl = Likelihood_B->Integral(Likelihood_B->FindBin(cut),Likelihood_B->FindBin(1.0));
    QL = 2*(sqrt(sl+bl)-sqrt(bl));
    merit_L->SetBinContent(merit_L->FindBin(cut), QL);
    cut+=0.025;
  }
  gStyle->SetOptStat(0);
  Likelihood_S->SetLineColor(kBlue);
  Likelihood_S->SetLineWidth(2);
  Likelihood_B->SetLineColor(kRed);
  Likelihood_B->SetLineWidth(2);
  merit_L->SetXTitle("cut");
  merit_L->SetYTitle("2(sqrt(s+b)-sqrt(b))");
  merit_L->SetTitle("Figures of merit");
  merit_L->SetMinimum(0);
  merit_L->SetLineWidth(2);
  merit_L->SetLineColor(kRed);

  return 1;
}

//=====================================CUT=SCAN=BDT============================================
bool hh4b_tmva::CutScanBDT(std::string outTMVAfile) {

  TFile *MVAtrain = new TFile(outTMVAfile.c_str(),"READ");
  if(!MVAtrain) return 0;
  TH1F *BDT_S_original = (TH1F*)MVAtrain->Get("Method_BDT/BDT_GiniIndex/MVA_BDT_GiniIndex_S");
  TH1F *BDT_B_original = (TH1F*)MVAtrain->Get("Method_BDT/BDT_GiniIndex/MVA_BDT_GiniIndex_B");
  TH1F *BDT_S = new TH1F("BDT_S", "BDT signal", 40, 0, 1);
  TH1F *BDT_B = new TH1F("BDT_B", "BDT bkg", 40, 0, 1);
  float tempS = 0;float tempB = 0;
  for(int i=0; i<=BDT_B_original->GetXaxis()->GetNbins(); i++){
    tempB = BDT_B_original->GetBinContent(i);
    BDT_B->SetBinContent(i, tempB);
  }
  std::cout << " Nbins " << BDT_B_original->GetXaxis()->GetNbins() << std::endl;
  //	float integral = BDT_B->Integral();
  BDT_B->Scale(2.83222748815166);     //debug -- scale?
  for(int i=0; i<=BDT_S_original->GetXaxis()->GetNbins(); i++){
    tempS = BDT_S_original->GetBinContent(i);
    BDT_S->SetBinContent(i, tempS);
  }
  //integral = BDT_S->Integral();
  BDT_S->Scale(5.99E-4);     //debug -- scale?

  float cut = 0.0;
  float QBDT, sb, bb, QBDT_err;
  while(cut <= 1) {
    sb = BDT_S->Integral(BDT_S->FindBin(cut), BDT_S->FindBin(1.0));
    bb = BDT_B->Integral(BDT_B->FindBin(cut), BDT_B->FindBin(1.0));
    QBDT_err = sqrt(bb/(bb + sb) + 4*sb*pow(-1/(2.*sqrt(sb)) + 1/(2.*sqrt(bb + sb)),2));
    QBDT = 2*(sqrt(sb+bb)-sqrt(bb));	
    merit_BDT->SetBinContent(merit_BDT->FindBin(cut), QBDT);
    merit_BDT->SetBinError(merit_BDT->FindBin(cut), QBDT_err);
    cut+=0.025;
  }
  gStyle->SetOptStat(0);
  BDT_S->SetLineColor(kGreen);
  BDT_S->SetMaximum(0.25);
  BDT_S->SetLineWidth(2);
  BDT_B->SetLineColor(kMagenta);
  BDT_B->SetLineWidth(2);
  merit_BDT->SetLineWidth(2);
  merit_BDT->SetLineColor(kBlue);
  merit_BDT->SetXTitle("cut");
  merit_BDT->SetYTitle("2(sqrt(s+b)-sqrt(b))");
  merit_BDT->SetTitle("Figures of merit");

  return 1;
}

//=====================================VARIABLES=NAMES============================================
std::string hh4b_tmva::setTitle(int nh) {
  
  std::string title;
  if(nh==0)  {title = "H1_pT";}
  if(nh==1)  {title = "H2_pT";}
  if(nh==2)  {title = "H1_dR";}
  if(nh==3)  {title = "H2_dR";}
  if(nh==4)  {title = "Centr";}
  if(nh==5)  {title = "fJet3_pT";}
  if(nh==6)  {title = "HH_mass";}
  if(nh==7)  {title = "HH_dR";}
  if(nh==8)  {title = "H1_cos";}  //equal to H2
/*
if(nh==0)  {title = "APt_min" ;}
if(nh==1)  {title = "APt_mean";}
if(nh==2)  {title = "APt_max" ;}
if(nh==3)  {title =	"AEta_min" ;}
if(nh==4)  {title =	"AEta_mean" ;}
if(nh==5)  {title =	"AEta_max" ;}
if(nh==6)  {title =	"ACSV_min";}
if(nh==7)  {title =	"ACSV_mean";}
if(nh==8)  {title =	"ACSV_max";}
if(nh==9)  {title = "Acent";}
if(nh==10) {title = "Qcent";}
if(nh==11) {title = "QPt_1";}
if(nh==12) {title = "QPt_2";}
if(nh==13) {title = "QPt_3";}
if(nh==14) {title = "QPt_4";}
if(nh==15) {title = "QEta_1";}
if(nh==16) {title = "QEta_2";}
if(nh==17) {title = "QEta_3";}
if(nh==18) {title = "QEta_4";}
if(nh==19) {title = "QCSV_1";}
if(nh==20) {title = "QCSV_2";}
if(nh==21) {title = "QCSV_3";}
if(nh==22) {title = "QCSV_4";}
if(nh==23) {title = "DJ_1_mass";}
if(nh==24) {title = "DJ_1_pt";}
if(nh==25) {title =	"DJ_1_Phi_aperture";}
if(nh==26) {title = "DJ_1_Eta_aperture";}
if(nh==27) {title = "DJ_1_R_aperture";}
if(nh==28) {title = "DJ_2_mass";}
if(nh==29) {title = "DJ_2_pt";}
if(nh==30) {title = "DJ_2_Phi_aperture";}
if(nh==31) {title = "DJ_2_Eta_aperture";}
if(nh==32) {title = "DJ_2_R_aperture";}
if(nh==33) {title = "TDJ_pt";}
if(nh==34) {title = "TDJ_deltaPhi";}
if(nh==35) {title = "TDJ_deltaEta";}
if(nh==36) {title = "TDJ_deltaR";}
if(nh==37) {title = "HHM";}
if(nh==38) {title = "met";}
if(nh==39) {title = "min_3csv";}
if(nh==40) {title = "avg_3csv";}
if(nh==41) {title = "costhetast";}
if(nh==42) {title = "costhetaCS";}
if(nh==43) {title = "tau_1";}
if(nh==44) {title = "tau_2";}
if(nh==45) {title = "JetsN";}*/
  return title;
}

