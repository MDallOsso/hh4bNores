// Code for non resonant HH->bbbb in CMS Run2 
// Step3 - 2D shape fit --> to split 2D plots (BDT output vs CSV) in 1D slices
//  Author: Martino Dall'Osso
//   from Carlo Alberto Gottardo
//    Nov 03rd, 2015
//----------------------------------------------------------------

#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"

using namespace std;

static const std::string outoption="";
static const std::string frameworkVersionFld = "V14/";   //V13
static const std::string resultsFld="results/"+frameworkVersionFld;
static const std::string inputFile= resultsFld+"scatters"+outoption+".root";  //DEBUG
static const std::string outFile= resultsFld+"toCombine"+outoption+".root";

int main()
{
TFile * apro = new TFile(inputFile.c_str(), "READ");
TH2D * sig_templ = (TH2D*)apro->Get("MC_BDT_vs_mass_4CSV");
sig_templ->Scale(2.01461e-4); //see notebook calculation
TH2D * bkg_templ = (TH2D*)apro->Get("BDT_vs_mass_3CSV");
TH2D * test = (TH2D*)apro->Get("BDT_vs_mass_4CSV");
TH2D * Twobtag = (TH2D*)apro->Get("BDT_vs_mass_2CSV");



//systematic scale error 0.63e-4
int binX_bdt = sig_templ->GetXaxis()->GetNbins();
int binY_mass = sig_templ->GetYaxis()->GetNbins();
int binTOT = binX_bdt*binY_mass;
for(int i=1; i<=binTOT; i++)
{
	if(sig_templ->GetBinContent(i) == 0) sig_templ->SetBinContent(i,0.000000001); 
	if(bkg_templ->GetBinContent(i) == 0) bkg_templ->SetBinContent(i,0.000000001); 
}

TH1D * sig_1D = new TH1D("sig1D_temp", "sig1D_temp", binTOT, 0, binTOT);
TH1D * bkg_1D = new TH1D("bkg1D_temp", "bkg1D_temp", binTOT, 0, binTOT);
TH1D * test_1D = new TH1D("data_obs_temp", "data_obs_temp", binTOT, 0, binTOT);
TH1D * err_1D_up = new TH1D("bkg1D_systUp_temp", "bkg1D_systUp_temp", binTOT, 0, binTOT);
TH1D * err_1D_down = new TH1D("bkg1D_systDown_temp", "bkg1D_systDown_temp", binTOT, 0, binTOT);

double delta = 0;
double SdeltaQ = 0;
double err = 0;
double N2 = Twobtag->Integral();
double N3 = bkg_templ->Integral();
double f = N3/N2;
double b3;
double b2;
double up = 0;
double down = 0;

for(int i=1; i<=binY_mass; i++)
{
	for(int j=1; j<=binX_bdt; j++)
	{
	sig_1D->SetBinContent((i-1)*binX_bdt+j, sig_templ->GetBinContent(i,j)); 
	bkg_1D->SetBinContent((i-1)*binX_bdt+j, bkg_templ->GetBinContent(i,j)); 
	test_1D->SetBinContent((i-1)*binX_bdt+j, test->GetBinContent(i,j)); 
	
	b2=Twobtag->GetBinContent(i,j);
	b3=bkg_templ->GetBinContent(i,j);
	delta = (b3-(f*b2));
	SdeltaQ = (b3 + (f*f)*b2 + (b2/N2)*(b2/N2)*N3 + ((b2*N3)/(N2*N2))*((b2*N3)/(N2*N2))*N2);

	if((delta*delta - SdeltaQ)>0.0) err = sqrt(delta*delta - SdeltaQ);
	if((delta*delta - SdeltaQ)<=0.0) err = 0.0;
	up = bkg_templ->GetBinContent(i,j)+err;
	down = bkg_templ->GetBinContent(i,j)-err;
	err_1D_up->SetBinContent((i-1)*binX_bdt+j, up );
	err_1D_down->SetBinContent((i-1)*binX_bdt+j, down );

	}
}


TH1D * sig_1D_crop = new TH1D("sig1D", "sig1D", 330, 0, 330);
TH1D * bkg_1D_crop = new TH1D("bkg1D", "bkg1D", 330, 0, 330);
TH1D * test_1D_crop = new TH1D("data_obs", "data_obs", 330, 0, 330);
TH1D * err_1D_crop_up = new TH1D("bkg1D_systUp", "bkg1D_systUp", 330, 0, 330);
TH1D * err_1D_crop_down = new TH1D("bkg1D_systDown", "bkg1D_systDown", 330, 0, 330);

for(int i=20; i<=binTOT; i++)
{
	if(i>=351)continue;
	sig_1D_crop  -> SetBinContent(i-19, sig_1D->GetBinContent(i));
	bkg_1D_crop  -> SetBinContent(i-19, bkg_1D->GetBinContent(i));
	test_1D_crop -> SetBinContent(i-19, test_1D->GetBinContent(i));
	err_1D_crop_up -> SetBinContent(i-19, err_1D_up->GetBinContent(i));
	err_1D_crop_down -> SetBinContent(i-19, err_1D_down->GetBinContent(i));
}

float scalefactor = test_1D_crop->Integral()/bkg_1D_crop->Integral();
bkg_1D_crop->Scale(scalefactor );
err_1D_crop_up->Scale(scalefactor );
err_1D_crop_down->Scale(scalefactor );

cout << "Integrals for Combine datacard" << endl;
cout << "data_obs: \t" << test_1D_crop->Integral() << endl;
cout << "signal: \t" << sig_1D_crop->Integral() << endl;
cout << "background: \t" << bkg_1D_crop->Integral() << endl;

TFile * output = new TFile(outFile.c_str(), "RECREATE");
sig_1D_crop->Write();
bkg_1D_crop->Write();
test_1D_crop->Write();
err_1D_crop_up->Write();
err_1D_crop_down->Write();
output->Close();

apro->Close();

return 0;
}

