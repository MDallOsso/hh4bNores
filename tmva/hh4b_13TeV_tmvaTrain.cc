// Code for non resonant HH->bbbb in CMS Run2 
// Step2 - TMVA for S/B discrimination
//  Author: Martino Dall'Osso
//   thanks Carlo Alberto Gottardo
//    Octt 19th, 2015
//----------------------------------------------------------------

#define hh4b_tmvaTrain_C

///lustre/cmswork/dallosso/hh2bbbb/non-resonant/analysis/13TeV/hh4bNores/tmva/TMVA-v4.2.0

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
// NAMES
//***********************************************
//  weightfile 	= names for .C and .xml weights file namely "window", "all"
//  outTMVAfileName = names for .root training and test files "TMVAinWindowTraining.root", "TMVAallTraining.root"
//  inputSignal = training() signal sample
//  inputBkg = training() bkg sample
//************************************************
// NOTES
//  training() trains both methods and produces two C and xml weights files
// 
//	applications: one TH1F per method and per sample choice
//  merits: one TH1F per method and per application choice
// 	

static const std::string dataFld="../data/";

#include "hh4b_13TeV_tmvaTrain.h"

hh4b_tmvaTrain::hh4b_tmvaTrain(){

//_______TRAINING____________ _______________________________________________________________
cout << "Training all" << endl; 
training(weightfile, outTMVAfileName, inputSignal, inputBkg);


//_______FIGURES OF MERIT____________________________________________________________________
//float dummy1, dummy2;
//cout << "calculating merits " << endl;
//TFile *f = new TFile("TMVA_testing.root", "RECREATE");
// all 

/*
TH1F * Amerit_L = CutScanL("TMVAallTraining.root");
TH1F * Amerit_BDT = CutScanBDT("TMVAallTraining.root");
dummy1 = Amerit_L->GetMaximum(); dummy2 = Amerit_L->GetMaximumBin();
cout << "All Likelihood \n max Q: " <<dummy1<< " at " << dummy2 << endl;
dummy1 = Amerit_BDT->GetMaximum(); dummy2 = Amerit_BDT->GetMaximumBin();
cout << "All window Likelihood" << "\n" << "max Q: " << dummy1 << " at " << dummy2 << endl;
Amerit_L->SetName("TrainA_merit_L");
Amerit_BDT->SetName("TrainA_merit_BDT");
f->cd();
Amerit_L->Write();
Amerit_BDT->Write();
*/

//_______APPLICATION__________________________________________________________________________
//cout << "Application in all trained in all" << endl;
//weightfile = "all";

//appfile = dataFld+"tree_Step1_BTagCSVRun2015C.root"; //debug----
//TH1D * S_MCLikelihood = Lapplication(appfile, weightfile);
//TH1D * S_MCBDT = BDTapplication(appfile, weightfile);
//S_MCLikelihood->SetName("CR_MC_Likel");
//S_MCBDT->SetName("CR_MC_BDT");
//f->cd();
//S_MCLikelihood->Write();
//S_MCBDT->Write();

/*
appfile = "../Samples/CR_NOCUT_hh4b-nores-8Tev300k_step2_L1y1.root";
TH1D * S2_MCLikelihood = Lapplication(appfile, weightfile);
TH1D * S2_MCBDT = BDTapplication(appfile, weightfile);
S2_MCLikelihood->SetName("MC_Likel");
S2_MCBDT->SetName("MC_BDT");
f->cd();
S2_MCLikelihood->Write();
S2_MCBDT->Write();

appfile = "../Samples/NOCUT_DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1.root";
TH1D * B_Likelihood = Lapplication(appfile, weightfile);
TH1D * B_BDT = BDTapplication(appfile, weightfile);
B_Likelihood->SetName("BKG_Likel");
B_BDT->SetName("B_BKG_BDT");
f->cd();
B_Likelihood->Write();
B_BDT->Write();

appfile = "../Samples/NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1.root";
TH1D * C_Likelihood = Lapplication(appfile, weightfile);
TH1D * C_BDT = BDTapplication(appfile, weightfile);
C_Likelihood->SetName("BKG_Likel");
C_BDT->SetName("C_BKG_BDT");
f->cd();
C_Likelihood->Write();
C_BDT->Write();

appfile = "../Samples/NOCUT_DiJetPt_BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1.root";
TH1D * CP_Likelihood = Lapplication(appfile, weightfile);
TH1D * CP_BDT = BDTapplication(appfile, weightfile);
CP_Likelihood->SetName("CP_BKG_Likel");
CP_BDT->SetName("CP_BKG_BDT");
f->cd();
CP_Likelihood->Write();
CP_BDT->Write();

appfile = "../Samples/NOCUT_DiJetPt_BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1.root";
TH1D * D_Likelihood = Lapplication(appfile, weightfile);
TH1D * D_BDT = BDTapplication(appfile, weightfile);
D_Likelihood->SetName("D_BKG_Likel");
D_BDT->SetName("D_BKG_BDT");
f->cd();
D_Likelihood->Write();
D_BDT->Write();

appfile = "../Samples/CR_NOCUT_DiJetPt_BJetPlusX_Run2012B-13Jul2012-v1_EDMNtuple_V42_ProcV1.root";
TH1D * CR_B_Likelihood = Lapplication(appfile, weightfile);
TH1D * CR_B_BDT = BDTapplication(appfile, weightfile);
CR_B_Likelihood->SetName("CR_B_Likel");
CR_B_BDT->SetName("CR_B_BDT");
f->cd();
CR_B_Likelihood->Write();
CR_B_BDT->Write();

appfile = "../Samples/CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug2012-v2_EDMNtuple_V42_ProcV1.root";
TH1D * CR_C_Likelihood = Lapplication(appfile, weightfile);
TH1D * CR_C_BDT = BDTapplication(appfile, weightfile);
CR_C_Likelihood->SetName("CR_C_Likel");
CR_C_BDT->SetName("CR_C_BDT");
f->cd();
CR_C_Likelihood->Write();
CR_C_BDT->Write();

appfile = "../Samples/CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-PromptReco-v2_EDMNtuple_V42_ProcV1.root";
TH1D * CR_CP_Likelihood = Lapplication(appfile, weightfile);
TH1D * CR_CP_BDT = BDTapplication(appfile, weightfile);
CR_CP_Likelihood->SetName("CR_CP_Likel");
CR_CP_BDT->SetName("CR_CP_BDT");
f->cd();
CR_CP_Likelihood->Write();
CR_CP_BDT->Write();

appfile = "../Samples/CR_NOCUT_DiJetPt_BJetPlusX_Run2012D-PromptReco-v1_EDMNtuple_V42_ProcV1.root";
TH1D * CR_D_Likelihood = Lapplication(appfile, weightfile);
TH1D * CR_D_BDT = BDTapplication(appfile, weightfile);
CR_D_Likelihood->SetName("CR_D_Likel");
CR_D_BDT->SetName("CR_D_BDT");
f->cd();
CR_D_Likelihood->Write();
CR_D_BDT->Write();*/

cout << "ho finito!" << endl;
//f->Close();

}

hh4b_tmvaTrain::~hh4b_tmvaTrain(){
}
//---------------


//******************************************************************************************************************
//		FUNCTIONS FROM NOW ON
//******************************************************************************************************************

//=====================================TRAINING=FUNCTION=====================================
void hh4b_tmvaTrain::training(string weightfile, string outTMVAfileName, string inputSignal, string inputBkg)
{
	TFile* outputFile = TFile::Open( outTMVAfileName.c_str(), "RECREATE" );
	TMVA::Factory *factory = new TMVA::Factory(weightfile.c_str(), outputFile, "!V:AnalysisType=Classification");

	TFile *sign = TFile::Open(inputSignal.c_str());
	TFile *bkg = TFile::Open(inputBkg.c_str());

	factory->AddSignalTree((TTree*)sign->Get("tree"),1.0);
	factory->AddBackgroundTree((TTree*)bkg->Get("tree"),1.0);

	string varible_name;
	for(int i=0; i<8; ++i){ 
//			cos*	dj_1 R   dj_2 R   dj_1 pt dj_2_pt	 Qcent    HHM   QPt_3   TDJ deltaR
		//if(i == 41| i== 27 | i== 32 | i==24 |  i==29 | i== 10 | i==37 | i==13 | i==36)
		//{
          factory->AddVariable(setTitle(i).c_str(), 'F');
          //}
	} 	
      //debug
	factory->PrepareTrainingAndTestTree("", "!V:nTrain_Signal=15000:nTrain_Background=15000:nTest_Signal=5000:nTest_Background=5000:SplitMode=Random");
	factory->BookMethod(TMVA::Types::kLikelihood, "Likelihood", "VarTransform=N+P(H1.pT,H1.dR)+P(H2.pT,H2.dR)");
	factory->BookMethod(TMVA::Types::kBDT, "BDT_GiniIndex", "NTrees=1200:MaxDepth=3:SeparationType=GiniIndex:AdaBoostR2Loss=Quadratic:VarTransform=N+P(H1.pT,H1.dR)+P(H2.pT,H2.dR)");
      //--------------

	factory->TrainAllMethods();
	factory->TestAllMethods();
	factory->EvaluateAllMethods();

	outputFile->Close();
	delete factory;
	sign->Close();
	bkg->Close();	// return;
}





//=====================================BDT=APPLICATION=====================================
TH1D* hh4b_tmvaTrain::BDTapplication(string appfile, string weightfile)
{
cout << "APPLICOOOO A " << appfile << endl;
	// Ouput of a new TFile which is a clone of the previous but has the new BDT branch
	//copy TTree and add new branch
	float Bout_orig = 0.0;
	float Bout = 0.0;
	TFile *old = new TFile(appfile.c_str(), "READ");
	TTree *oldtree = (TTree*)old->Get("tree");
	oldtree->SetBranchStatus("*",1);
	//Create a new file + a clone of old tree in new file
	string culino = appfile.substr(11,38);
	string Nntuple = "BDT_"+culino+".root";
  	TFile *newtfile = new TFile(Nntuple.c_str(),"recreate");
  	TTree *newtree = oldtree->CloneTree();
   	TBranch *br = newtree->Branch("BDT", &Bout, "Bout/F");

	TH1D *bdt_out = new TH1D("bdt_application","bdt_application",100,0,1);
	TMVA::Reader *reader = new TMVA::Reader("!Color");
	float var[46];
	TFile *input = TFile::Open(appfile.c_str());
	TTree *theTree = (TTree*)input->Get("tree");

	string varible_name;
	for(int i=0; i<8; ++i){
//			cos*	dj_1 R   dj_2 R   dj_1 pt dj_2_pt	 Qcent    HHM   QPt_3   TDJ deltaR
		//if(i == 41| i== 27 | i== 32 | i==24 |  i==29 | i== 10 | i==37 | i==13 | i==36)
		//{
          reader->AddVariable(setTitle(i).c_str(), &var[i]);
        //}
	}
	string readweighfile = "weights/"+weightfile+"_BDT_GiniIndex.weights.xml";
	reader->BookMVA("BDT_GiniIndex", readweighfile.c_str() );
	for(int i=0; i<8; ++i){
	theTree->SetBranchAddress(setTitle(i).c_str(), &var[i]);
	}	
	for(long i=0; i<theTree->GetEntries(); i++)
		{
			theTree->GetEntry(i);
			Bout_orig = reader->EvaluateMVA( "BDT_GiniIndex" );
			Bout = (Bout_orig+1.0)/2.0;
			bdt_out->Fill(Bout);
			br->Fill();
		}
	//newtree->SetEntries(br->GetEntries());
  	newtree->Print();
   	newtfile->Write();
	return bdt_out;
}



//=====================================LIKELIHOOD=APPLICATION=====================================
TH1D* hh4b_tmvaTrain::Lapplication(string appfile, string weightfile)
{

	TH1D *like_out = new TH1D("like_application","like_application",100,0,1);
	TMVA::Reader *reader = new TMVA::Reader("!Color");
	float var[46];
	TFile *input = TFile::Open(appfile.c_str());
	TTree *theTree = (TTree*)input->Get("tree");

	string varible_name;
	for(int i=0; i<8; ++i){
//			cos*	dj_1 R   dj_2 R   dj_1 pt dj_2_pt	 Qcent    HHM   QPt_3   TDJ deltaR
		//if(i == 41| i== 27 | i== 32 | i==24 |  i==29 | i== 10 | i==37 | i==13 | i==36)
		//{
          reader->AddVariable(setTitle(i).c_str(), &var[i]);
         //}
	}
	string readweighfile = "weights/"+weightfile+"_Likelihood.weights.xml";
	reader->BookMVA("Likelihood", readweighfile.c_str() );
	//for(int i=0; i<8; ++i){ //debug -46
	//if(i==21 | i==22) continue;
	  //theTree->SetBranchAddress(setTitle(i).c_str(), &var[i]);
	//}	
	double Lout;
	for(long i=0; i<theTree->GetEntries(); i++)
		{
			theTree->GetEntry(i);
			Lout = reader->EvaluateMVA( "Likelihood" );
			like_out->Fill(Lout);
		}
return like_out;
}




//=====================================CUT=SCAN=L============================================
TH1F* hh4b_tmvaTrain::CutScanL(string outTMVAfile)
{
	TFile *MVAtrain = new TFile(outTMVAfile.c_str(),"READ");
	TH1F *Likelihood_S_original = (TH1F*)MVAtrain->Get("Method_Likelihood/Likelihood/MVA_Likelihood_S");
	TH1F *Likelihood_B_original = (TH1F*)MVAtrain->Get("Method_Likelihood/Likelihood/MVA_Likelihood_B");
	TH1F *Likelihood_S = new TH1F("L_S", "Likelihood Signal", 40, 0, 1);
	TH1F *Likelihood_B = new TH1F("L_B", "Likelihood Background", 40, 0, 1);
	TH1F *merit_L = new TH1F("Q_L", "Likelihood Q", 40, 0, 1);

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
	while(cut <= 1)
	{
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


return merit_L;
}



//=====================================CUT=SCAN=BDT============================================
TH1F* hh4b_tmvaTrain::CutScanBDT(string outTMVAfile)
{
	TFile *MVAtrain = new TFile(outTMVAfile.c_str(),"READ");
	TH1F *BDT_S_original = (TH1F*)MVAtrain->Get("Method_BDT/BDT_GiniIndex/MVA_BDT_GiniIndex_S");
	TH1F *BDT_B_original = (TH1F*)MVAtrain->Get("Method_BDT/BDT_GiniIndex/MVA_BDT_GiniIndex_B");
	TH1F *BDT_S = new TH1F("BDT_S", "BDT signal", 40, 0, 1);
	TH1F *BDT_B = new TH1F("BDT_B", "BDT bkg", 40, 0, 1);
	TH1F *merit_BDT = new TH1F("Q_BDT", "BDT Q", 40, 0, 1);

	float tempS = 0;float tempB = 0;
	for(int i=0; i<=BDT_B_original->GetXaxis()->GetNbins(); i++){
		tempB = BDT_B_original->GetBinContent(i);
		BDT_B->SetBinContent(i, tempB);
		
	}
	cout << " Nbins " << BDT_B_original->GetXaxis()->GetNbins() << endl;
//	float integral = BDT_B->Integral();
	BDT_B->Scale(2.83222748815166);
	for(int i=0; i<=BDT_S_original->GetXaxis()->GetNbins(); i++){
		tempS = BDT_S_original->GetBinContent(i);
		BDT_S->SetBinContent(i, tempS);
	}
	//integral = BDT_S->Integral();
	BDT_S->Scale(5.99E-4);

	float cut = 0.0;
	float QBDT, sb, bb, QBDT_err;
	while(cut <= 1)
	{
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

	return merit_BDT;
}

//=====================================VARIABLES=NAMES============================================
string hh4b_tmvaTrain::setTitle(int nh)
{
string title;
if(nh==0)  {title = "H2.pT";}
if(nh==1)  {title = "H2.dR";}
if(nh==2)  {title = "H1.dR";}
if(nh==3)  {title = "H1.pT";}
if(nh==4)  {title = "H1.eta";}//fJet3.CSV
if(nh==5)  {title = "H2.eta";} //fJet3avgCSV
if(nh==6)  {title = "H1.dPhi";}
if(nh==7)  {title = "H2.dPhi";}  //fJet4.CSV

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
