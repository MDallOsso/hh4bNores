// Code for non resonant HH->bbbb in CMS Run2 
// Step3 - 2D shape fit --> to create 2D plots: BDT output vs CSV
//  Author: Martino Dall'Osso
//   from Carlo Alberto Gottardo
//    Nov 03rd, 2015
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

using namespace std;
//float Validation_cut[11] = {0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52};
//float lower_cut[14] = {0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51};

//WARNING: CHANGE ME!!
static const std::string frameworkVersionFld = "V14/";   //V13
static const std::string dataFld="../data/"+frameworkVersionFld;;
static const std::string inputFile1= dataFld+"tree_Step2_BTagCSVRun2015DOct5_pT20CSVL_2nd.root";  //DEBUG
static const std::string inputFile2= dataFld+"tree_Step2_BTagCSVRun2015DOct5_pT20CSVL.root";
static const std::string inputMC= dataFld+"tree_Step2_M-260_pT20CSVL.root";

static const std::string resultsFld="results/"+frameworkVersionFld;
static const std::string outoption="";
static const std::string outFile = resultsFld+"scatters"+outoption+".root";

float Validation_cut[1] = {0.45};
float Q_cut = 0.58;
float lower_cut = 0.39;
//float Validation_cut = 0.50;
float CSVcut = 0.679;

struct abcd
{
float nTr;
float eta;
float pt;
float bv;
float qcsv;
float mass;
int region; // A=1, A_val=2, B=3, C = 4 , C_val = 5, D = 6; 
abcd()
	{
		nTr=0;
		eta=0;
		pt=0;
		bv=0;
		qcsv=0;
		region = 0; // 0 = unassigned
		mass = 0;
	}
};
vector <double> Dcount(vector<abcd> &, int, int, int, int, float *);



int main()
{

float N[10];
vector <abcd> evU;
vector <abcd> evL;
int massbins = 23;
int bdtbins = 20;
TH2F *mass_scatt_ALL = new TH2F("BDT_vs_mass", "BDT vs mass", bdtbins, 0.25, 0.75, massbins, 40, 500); 
TH2F *mass_scatt_3CSV = new TH2F("BDT_vs_mass_3CSV", "BDT vs mass (3btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
TH2F *mass_scatt_4CSV = new TH2F("BDT_vs_mass_4CSV", "BDT vs mass (4btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
TH2F *mass_scatt_2CSV = new TH2F("BDT_vs_mass_2CSV", "BDT vs mass (2btag)", bdtbins, 0.25, 0.75, massbins, 40, 500);
TH2F *mass_scatt_TEST = new TH2F("BDT_vs_mass_TEST", "2-3 btag syst. error / 4 btag stat. unc. ", bdtbins, 40.25, 0.75, massbins, 40, 500); 
TH2F *MC_mass_scatt_4CSV = new TH2F("MC_BDT_vs_mass_4CSV", "MC BDT vs mass (4btag)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
TH2F *MC_mass_scatt_3CSV = new TH2F("MC_BDT_vs_mass_3CSV", "MC BDT vs mass (3btag)", bdtbins, 0.25, 0.75, massbins, 40, 500);
TH2F *MC_mass_scatt_ALL = new TH2F("MC_BDT_vs_mass_TEST", "MC BDT vs mass (ALL)", bdtbins, 0.25, 0.75, massbins, 40, 500); 
gStyle->SetOptStat(0);
mass_scatt_ALL->GetXaxis()->SetTitle("BDT response"); mass_scatt_ALL->GetYaxis()->SetTitle("h mass");
mass_scatt_3CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_3CSV->GetYaxis()->SetTitle("h mass");
mass_scatt_4CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_4CSV->GetYaxis()->SetTitle("h mass");
mass_scatt_2CSV->GetXaxis()->SetTitle("BDT response"); mass_scatt_2CSV->GetYaxis()->SetTitle("h mass");
mass_scatt_TEST->GetXaxis()->SetTitle("BDT response"); mass_scatt_TEST->GetYaxis()->SetTitle("h mass");
MC_mass_scatt_4CSV->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_4CSV->GetYaxis()->SetTitle("h mass");
MC_mass_scatt_3CSV->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_3CSV->GetYaxis()->SetTitle("h mass");
MC_mass_scatt_ALL->GetXaxis()->SetTitle("BDT response"); MC_mass_scatt_ALL->GetYaxis()->SetTitle("h mass");



//-------------------Reading from files---------

TChain *upperTree = new TChain("tree");
upperTree->Add(inputFile1.c_str());
//upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-24.root");
//upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012C-Pr.root");
//upperTree->Add("../Samples_BDT/BDT_CR_NOCUT_DiJetPt_BJetPlusX_Run2012D-Pr.root");
	float u_tr, u_eta, u_pt, u_BDT, u_csv, u_m1, u_m2;
//	upperTree->SetBranchAddress("QPt_3", &u_pt);
//	upperTree->SetBranchAddress("QEta_3", &u_eta);
//	upperTree->SetBranchAddress("ThirdJetTracks", &u_tr);
	upperTree->SetBranchAddress("fJet4_CSV", &u_csv);
	upperTree->SetBranchAddress("BDT", &u_BDT);
	upperTree->SetBranchAddress("H1_mass", &u_m1);
	upperTree->SetBranchAddress("H2_mass", & u_m2);
cout << upperTree->GetEntries() << endl;
	for(int i=0; i<upperTree->GetEntries(); i++)
	{
		upperTree->GetEntry(i);
		evU.push_back(abcd());
//		evU[i].nTr = u_tr;
//		evU[i].eta = u_eta;
//		evU[i].pt = u_pt;
		evU[i].bv = u_BDT;
		evU[i].qcsv = u_csv;
		evU[i].mass = (u_m1+u_m2)/2.0;
		mass_scatt_ALL->Fill(evU[i].bv, evU[i].mass);
		mass_scatt_2CSV->Fill(evU[i].bv, evU[i].mass);
	}

TChain *lowerTree = new TChain("tree");
lowerTree->Add(inputFile2.c_str());
//lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-24Aug.root");
//lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012C-Promp.root");
//lowerTree->Add("../Samples_BDT/BDT_NOCUT_DiJetPt_BJetPlusX_Run2012D-Promp.root");
	
	float d_tr, d_eta, d_pt, d_BDT, d_csv, d_m1, d_m2;
//	lowerTree->SetBranchAddress("QPt_3", &d_pt);
//	lowerTree->SetBranchAddress("QEta_3", &d_eta);
//	lowerTree->SetBranchAddress("ThirdJetTracks", &d_tr);
	lowerTree->SetBranchAddress("fJet4_CSV", &d_csv);
	lowerTree->SetBranchAddress("BDT", &d_BDT);
	lowerTree->SetBranchAddress("H1_mass", & d_m1);
	lowerTree->SetBranchAddress("H2_mass", & d_m2);
cout << lowerTree->GetEntries() << endl;
	for(int i=0; i<lowerTree->GetEntries(); i++)
	{
		lowerTree->GetEntry(i);
		evL.push_back(abcd());
//		evL[i].nTr = d_tr;
//		evL[i].eta = fabs(d_eta);
//		evL[i].pt = d_pt;
		evL[i].bv = d_BDT;
		evL[i].qcsv = d_csv;
		evL[i].mass = (d_m1+d_m2)/2.0;
		mass_scatt_ALL->Fill(evL[i].bv, evL[i].mass);
		if(evL[i].qcsv < CSVcut) mass_scatt_3CSV->Fill(evL[i].bv, evL[i].mass);
		if(evL[i].qcsv >= CSVcut) mass_scatt_4CSV->Fill(evL[i].bv, evL[i].mass);
	}

TH2F *diff = new TH2F("twothree_btag_Residuals", "2-3 btag Residuals", bdtbins, 0.25, 0.75, massbins, 40, 500);
int TOTbins = mass_scatt_2CSV->GetXaxis()->GetNbins()*mass_scatt_2CSV->GetYaxis()->GetNbins();

float delta, unc;
float N2 = mass_scatt_2CSV->GetEntries();
float N3 = mass_scatt_3CSV->GetEntries();
float bin2, bin3, bin4;
float test;
for(int i=1; i<=TOTbins; i++)
{
	bin2 = mass_scatt_2CSV->GetBinContent(i);
	bin3 = mass_scatt_3CSV->GetBinContent(i);
	bin4 = mass_scatt_4CSV->GetBinContent(i);

	if (bin2 == 0 || bin3 == 0) {delta = 0.0; unc = 1.1;}
	else { 
		delta = bin2 - (bin3*(N2/N3));
		unc = sqrt(bin2 + (N2/N3)*(N2/N3)*bin3 + (bin3/N3)*(bin3/N3)*N2 + ((bin3*N2)/(N3*N3))*((bin3*N2)/(N3*N3))*N3 );
		if (bin4 != 0) {test = ((((N3/N2)*bin2)-bin3)/bin3)*sqrt(bin4);}
		if (bin4 == 0) {test = 0.0;}
		}
		diff->SetBinContent(i,delta/unc);
		mass_scatt_TEST->SetBinContent(i,test);
}


TFile * MC = new TFile(inputMC.c_str(), "READ");
//BDT_NOCUT_hh4b-nores-8Tev300k_step2_L1y1.root", "READ");
TTree * leggo = (TTree*)MC->Get("tree");

vector <abcd> evMC;

	float mc_BDT, mc_csv, mc_m1, mc_m2;
//	lowerTree->SetBranchAddress("QPt_3", &d_pt);
//	lowerTree->SetBranchAddress("QEta_3", &d_eta);
//	lowerTree->SetBranchAddress("ThirdJetTracks", &d_tr);
	leggo->SetBranchAddress("fJet4_CSV", &mc_csv);  //4th Jet CSV
	leggo->SetBranchAddress("BDT", &mc_BDT);
	leggo->SetBranchAddress("H1_mass", & mc_m1); //diJet1 mass
	leggo->SetBranchAddress("H2_mass", & mc_m2); //diJet2 mass

	for(int i=0; i<leggo->GetEntries(); i++)
	{
		leggo->GetEntry(i);
		evMC.push_back(abcd());
		evMC[i].bv = mc_BDT;
		evMC[i].qcsv = mc_csv;
		evMC[i].mass = (mc_m1+mc_m2)/2.0;
		MC_mass_scatt_ALL->Fill(evMC[i].bv, evMC[i].mass);
		if(evMC[i].qcsv < CSVcut) MC_mass_scatt_3CSV->Fill(evMC[i].bv, evMC[i].mass);
		if(evMC[i].qcsv >= CSVcut) MC_mass_scatt_4CSV->Fill(evMC[i].bv, evMC[i].mass);
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

return 0;
}

