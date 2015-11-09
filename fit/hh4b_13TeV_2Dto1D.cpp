// Code for non resonant HH->bbbb in CMS Run2 
// Step3 - 2D shape fit --> to split 2D plots (BDT output vs CSV) in 1D slices
//  Author: Martino Dall'Osso
//   from Carlo Alberto Gottardo
//    Nov 03rd, 2015
//     g++ `root-config --libs  --cflags` hh4b_13TeV_2Dto1D.cpp -o toCombine
//----------------------------------------------------------------

#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TAxis.h"

//WARNING: CHANGE ME!!
static const std::string outoption="";
static const std::string frameworkVersionFld = "V15/";   //V13
static const float scaleFactorSig = 1.; //11,54 = 34.3*0.58*0.58 //was 2.01461e-4 -- DEBUG  --what to be ued?
// 13TeV XS:
// X->HH narrow M-260 --> 
// nonres HH SM --> 
float scaleFactorBkg = 0; //computed inside code

static const std::string resultsFld="results/"+frameworkVersionFld;
static const std::string inputFile= resultsFld+"scatters"+outoption+".root";  //DEBUG
static const std::string outFile= resultsFld+"toCombine"+outoption+".root";

int main()
{
  TFile * apro = new TFile(inputFile.c_str(), "READ");
  TH2D * sig_templ = (TH2D*)apro->Get("MC_BDT_vs_mass_4CSV");
  TH2D * bkg_templ = (TH2D*)apro->Get("BDT_vs_mass_3CSV");
  TH2D * test = (TH2D*)apro->Get("BDT_vs_mass_4CSV");
  TH2D * Twobtag = (TH2D*)apro->Get("BDT_vs_mass_2CSV");

  // Apply Signal Scale Factor (f'=Lbkg/Lsig)
 //--------------------------------------------
  sig_templ->Scale(scaleFactorSig);

  // Set Bin Content and compute shape sys uncert.
  //-------------------------------------------------
//systematic scale error 0.63e-4 ?? -- DEBUG
  int binX_bdt = sig_templ->GetXaxis()->GetNbins();
  int binY_mass = sig_templ->GetYaxis()->GetNbins();
  int binTOT = binX_bdt*binY_mass;
  for(int i=1; i<=binTOT; i++) {  // to avoid zeros for the fit? - debug
    if(sig_templ->GetBinContent(i) == 0) sig_templ->SetBinContent(i,0.000000001); 
    if(bkg_templ->GetBinContent(i) == 0) bkg_templ->SetBinContent(i,0.000000001); 
  }

  TH1D * sig_1D = new TH1D("sig1D_temp", "sig1D_temp", binTOT, 0, binTOT); //MC 4b
  TH1D * bkg_1D = new TH1D("bkg1D_temp", "bkg1D_temp", binTOT, 0, binTOT); //data 3b
  TH1D * test_1D = new TH1D("data_obs_temp", "data_obs_temp", binTOT, 0, binTOT);  //data 4b
  TH1D * err_1D_up = new TH1D("bkg1D_systUp_temp", "bkg1D_systUp_temp", binTOT, 0, binTOT);  //shape syst computed below
  TH1D * err_1D_down = new TH1D("bkg1D_systDown_temp", "bkg1D_systDown_temp", binTOT, 0, binTOT);  //shape syst computed below
  double delta = 0, SdeltaQ = 0, err = 0, up = 0, down = 0;
  double N2 = Twobtag->Integral();
  double N3 = bkg_templ->Integral();
  double f = N3/N2;
  double b3, b2;
  for(int i=1; i<=binY_mass; i++) {
    for(int j=1; j<=binX_bdt; j++) {
      sig_1D->SetBinContent((i-1)*binX_bdt+j, sig_templ->GetBinContent(i,j)); 
      bkg_1D->SetBinContent((i-1)*binX_bdt+j, bkg_templ->GetBinContent(i,j)); 
      test_1D->SetBinContent((i-1)*binX_bdt+j, test->GetBinContent(i,j));
	
      b2=Twobtag->GetBinContent(i,j);
      b3=bkg_templ->GetBinContent(i,j);
      delta = (b3-(f*b2));
      SdeltaQ = (b3 + (f*f)*b2 + (b2/N2)*(b2/N2)*N3 + ((b2*N3)/(N2*N2))*((b2*N3)/(N2*N2))*N2);  // debug - to be checked
      if((delta*delta - SdeltaQ)>0.0) err = sqrt(delta*delta - SdeltaQ);
      if((delta*delta - SdeltaQ)<=0.0) err = 0.0;  // DEBUG
      up = bkg_templ->GetBinContent(i,j)+err;
      down = bkg_templ->GetBinContent(i,j)-err;
      err_1D_up->SetBinContent((i-1)*binX_bdt+j, up );
      err_1D_down->SetBinContent((i-1)*binX_bdt+j, down );
    }
  }

  // Final Histos
  //-----------------
  int Xmax = 330, finalBins = 330; //was 330
  TH1D * sig_1D_crop = new TH1D("sig1D", "sig1D", finalBins, 0, Xmax);
  TH1D * bkg_1D_crop = new TH1D("bkg1D", "bkg1D", finalBins, 0, Xmax);
  TH1D * test_1D_crop = new TH1D("data_obs", "data_obs", finalBins, 0, Xmax);
  TH1D * err_1D_crop_up = new TH1D("bkg1D_systUp", "bkg1D_systUp", finalBins, 0, Xmax);
  TH1D * err_1D_crop_down = new TH1D("bkg1D_systDown", "bkg1D_systDown", finalBins, 0, Xmax);
  int bin0 = 19; //was 19
  for(int i=bin0+1; i<=binTOT; i++) {
    if(i>=351)continue; //debug--
    sig_1D_crop  -> SetBinContent(i-bin0, sig_1D->GetBinContent(i));
    bkg_1D_crop  -> SetBinContent(i-bin0, bkg_1D->GetBinContent(i));
    test_1D_crop -> SetBinContent(i-bin0, test_1D->GetBinContent(i));
    err_1D_crop_up -> SetBinContent(i-bin0, err_1D_up->GetBinContent(i));
    err_1D_crop_down -> SetBinContent(i-bin0, err_1D_down->GetBinContent(i));
  }

  // Apply Scale Factor (f=tot4bTag/tot3bTag)
  //--------------------------------------------
  scaleFactorBkg = test_1D_crop->Integral()/bkg_1D_crop->Integral();
  bkg_1D_crop->Scale(scaleFactorBkg );
  err_1D_crop_up->Scale(scaleFactorBkg );
  err_1D_crop_down->Scale(scaleFactorBkg );


  // Output for combine
  //-----------------------
  std::cout << "Integrals for Combine datacard" << std::endl;
  std::cout << "data_obs: \t" << test_1D_crop->Integral() << std::endl;
  std::cout << "signal: \t" << sig_1D_crop->Integral() << std::endl;
  std::cout << "background: \t" << bkg_1D_crop->Integral() << std::endl;  //with f scale factor

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

