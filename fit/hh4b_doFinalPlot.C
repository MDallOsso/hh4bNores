// Code for non resonant HH->bbbb in CMS Run2 
// Step3 - 2D shape fit --> to plot Fit result from Combine
//  Author: Martino Dall'Osso
//   from Carlo Alberto Gottardo
//    Nov 03rd, 2015
//----------------------------------------------------------------

static const std::string outoption="";  //TO BE CHANGED
static const std::string frameworkVersionFld = "V15/";   //equal to the Data version
static const std::string resultsFld="results/"+frameworkVersionFld;;
static const std::string inputFile= "combine/mlfit.root";
static const std::string errFile= resultsFld+"toCombine.root";
static const std::string outFile= resultsFld+"fit2D_"+outoption+".root";

void hh4b_doFinalPlot()
{
  TFile *infile = new TFile(inputFile.c_str(), "READ");

  TH1F* data  = new TH1F("observed","observed",330,0,330);
 //data->SetFillColor(kBlack);

  RooPlot* pre=(RooPlot*)infile->Get("hh4b_prefit");
  RooPlot* post_bsFit=(RooPlot*)infile->Get("hh4b_fit_s");
  RooPlot* post_bFit=(RooPlot*)infile->Get("hh4b_fit_b");

  // Get PreFit Background
  TDirectory* root=gDirectory;
  gDirectory->cd("shapes_prefit/hh4b/");
  TH1F * bkgPre = (TH1F*)gDirectory->Get("bkg1D");
   bkgPre->SetName("bkg pre");
   bkgPre->SetLineColor(kYellow-7);
   //bkgPre->SetLineWidth(5);
   //bkgPre->SetFillColor(kYellow-9);
   //bkgPre->SetFillStyle(3003);

  //DATA
  RooHist* rooData=pre->getHist("h_hh4b"); // this is the rooData
  double* y_arr = rooData->GetY();
  for (int i=0; i<=rooData->GetN(); ++i) {
    if (y_arr[i] == 0) rooData->SetPointError(i,0.,0.,0.,0.);
  }
  //double* y1_arr = rooData->GetY();
  double* x_arr = rooData->GetX();
  //cout << "rooData " << rooData->GetN()<< endl;
  for (int idata=0; idata<rooData->GetN(); ++idata) {
    data->Fill(x_arr[idata],y_arr[idata]);
  }
  data->SetLineColor(kBlack);
  data->SetLineWidth(1);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1.);
  data->GetYaxis()->SetTitleOffset(0.7);

  // Background Post Fit (B only)
  gDirectory->cd("../../shapes_fit_b/hh4b/");
  TH1F * bkgPostB = (TH1F*)gDirectory->Get("bkg1D");
   bkgPostB->SetLineColor(kGreen+1);
   bkgPostB->SetLineWidth(4);

  // Background Post Fit (B+S)
  gDirectory->cd("../../shapes_fit_s/hh4b/");
  TH1F * bkgPostSB = (TH1F*)gDirectory->Get("bkg1D");
   bkgPostSB->SetLineColor(kViolet+7);
   bkgPostSB->SetLineWidth(3);
   //bkgPostSB->SetFillColor(kViolet+1);

  TH1F * sigPostSB = (TH1F*)gDirectory->Get("sig1D");
   sigPostSB->SetLineColor(kRed+1);
   sigPostSB->SetLineWidth(3);

  TH1F * sigbkgPostSB = (TH1F*) sigPostSB->Clone();
   sigbkgPostSB->Add(bkgPostSB, 1.0);
   sigbkgPostSB->SetLineColor(kBlue+1);
   sigbkgPostSB->SetLineWidth(3);
   //sigbkgPostSB->SetFillColor(kBlue-2);
	
  TLegend* l=new TLegend(0.62,0.67,.90,.90);
   l->SetFillColor(kWhite);
   l->SetFillStyle(0);
   l->SetTextSize(0.03);
   l->AddEntry(data,"Data","pe");
   l->AddEntry(bkgPre,"Bkg (pre fit)","l");
   l->AddEntry(bkgPostB,"Bkg (B only fit)","l");
   l->AddEntry(bkgPostSB,"Bkg (S+B fit)","l");
   l->AddEntry(sigPostSB,"Signal (S+B fit)","l");
   l->AddEntry(sigbkgPostSB,"Sig+Bkg (S+B fit)","l");

  TFile * errfile = new TFile(errFile.c_str(), "READ");
  TH1F * up = (TH1F*)errfile->Get("bkg1D_systUp");
  TH1F * down = (TH1F*)errfile->Get("bkg1D_systDown");

  TCanvas * c = new TCanvas("c", "c", 800, 600);
  c->SetLogy();
  gStyle->SetOptStat("0");
  up->SetLineColor(kBlue-9);
  up->SetFillColor(kBlue-9);
  down->SetLineColor(kBlue-9);

  down->SetFillColor(kYellow-7);
  up->SetMinimum(1);
  up->Draw();
  down->Draw("SAME");
  up->GetXaxis()->SetTitle("Bin number");
  up->GetYaxis()->SetTitle("Count");
  up->GetYaxis()->SetTitleOffset(1.4);
  up->SetTitle("BDT vs mass Fit");
  bkgPre->SetLineWidth(4);
  bkgPre->Draw("SAME");
  bkgPostB->Draw("SAME");
  bkgPostSB->Draw("SAME");
  sigbkgPostSB->Draw("SAME");
  sigPostSB->SetFillColor(kRed+1);
  sigPostSB->Draw("SAME");

  data->SetMarkerStyle(8);
  data->SetMarkerSize(0.7);
  data->Draw("SAME, PE");

  l->Draw();
  c->Print(outFile.c_str());
  c->Print(outFile.c_str());

  //infile->Close();
  return;
}

