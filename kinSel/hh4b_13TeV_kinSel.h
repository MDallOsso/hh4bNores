#define hh4b_kinSel_H

//just to have order
//--------------------
class hh4b_kinSel{

  private:
//variables
//-------------------
  ULong64_t evt;
  unsigned int run; //lumi
  int nPV, nJet, nGenH, nGenBfHafterISR;    
  float weightPU, Vtype_, met_pt, met_eta, met_phi, met_mass, ht;
  float HLT_HH4bAll, HLT_BIT_QuadTriple, HLT_BIT_QuadDouble, HLT_BIT_DoubleTriple, HLT_BIT_DoubleDouble;
//float Jet_btagCMVA[25], Jet_rawPt[25],Jet_mcPt[25],Jet_mcFlavour[25]
//float Jet_corr_JECUp[25], Jet_corr_JECDown[25], Jet_corr[25], Jet_mcMatchId[25], Jet_hadronFlavour[25],Jet_btagProb[25];
//float Jet_btagBProb[25], Jet_btagnew[25],Jet_btagCSVV0[25], Jet_chHEF[25],Jet_neHEF[25],Jet_chEmEF[25],Jet_neEmEF[25]
//float Jet_chMult[25],Jet_leadTrackPt[25],Jet_mcEta[25],Jet_mcPhi[25],Jet_mcM[25];
  float Jet_id[25],Jet_puId[25],Jet_btagCSV[25];
  float Jet_pt[25],Jet_eta[25],Jet_phi[25],Jet_mass[25];  
  float GenH_pt[2], GenH_mass[2], GenH_eta[2], GenH_phi[2], GenH_status[2];
  float nGenBfH[4], GenBfH_pdgId[4],GenBfH_pt[4],GenBfH_eta[4],GenBfH_phi[4],GenBfH_mass[4], GenBfH_status[4], GenBfH_charge[4];
  float GenBfHafterISR_pdgId[4],GenBfHafterISR_pt[4],GenBfHafterISR_eta[4],GenBfHafterISR_phi[4],GenBfHafterISR_mass[4];

//debug...
  int nCut0=0, nCut1=0, nCut2=0, nCut3=0, nCut4=0, nCut4a=0, nCut4b=0, nCut5=0, nCut5a=0, nCut5b=0, HHf=0;
  int nJets_InAcc=0;
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
  int nRmaxOk = 0, nMCTruth =0, nCSVOk =0, nMatrCSVOk=0, nMixJets=0, errfJets=0;
  double Rave = 9999;
  int InTrue=99, In=99, jet4index=99;
  vector<int> fJetsIndex;
//----------

  //variables for histos
  diJet H1, H2, HH; //Higgs and di-Higgs candidates per event
  TLorentzVector H1_P, H2_P, HH_P, met; //Higgs and di-Higgs candidates per event
  std::vector<TLorentzVector> jets_inAcc_P;   //jets in acc sorted by csv - per event
  std::vector<TLorentzVector> fJets_P; //final 4 jets per event
  std::vector<Jet>  fJets_CSV; //final 4 jets sorted in CSV
  std::vector<Jet>  fJets; //final 4 jets per event
  std::vector<Jet>  aJets; //additional jets per event
  std::vector<Jet>  jets_inAcc; //jet vector -> baseline in the event loop - per event
  std::vector<Jet>  jets_all; //all jets after trigger - per event
  Jet fJet1_CSV, fJet2_CSV, fJet3_CSV, fJet4_CSV;

//----------------------------

  public:
    hh4b_kinSel(std::string , bool , std::string , int =2, int =0, std::string ="");
// void    doKinSel(std::string , bool , std::string , int =2, int =0, std::string ="");
    ~hh4b_kinSel();

    double fJets3avgCSV, fJets3minCSV, fJets4avgCSV;
    double Centr;

    void anglesComputation();
    void setBranches(bool, TTree*);   
    bool readMatrix( std::string, std::string );
    double selectBestDiJets (int );
    void fillHistos(bool );   
    void writeHistos(std::string, bool );   
    string jet4SelectionMethod(int );  
    void createTree(TTree*);
    void printOutput(std::string, bool );
    //int dokinSel(std::string , bool , std::string , int =2, int =0, std::string ="");
};


// book Histos
//--------------
TH1F *h_nPV = new TH1F("h_nPV", "# of Primary Vertices; nPV", 10, 0., 10.);
TH1F *h_Cuts=new TH1F("h_Cuts", "Cut flow", 20, 0, 20);

TH1F *h_nJets = new TH1F("h_nJets", "# Jets; # jets", 18, 0., 18.);
TH1F *h_nJets_InAcc = new TH1F("h_nJets_InAcc", "# Jets in Acc; # jets in Acc", 18, 0., 18.);

//TH1F *h_nPV_weighted=new TH1F("h_nPV_weighted", "# of Primary Vertices after Reweighting; nPV", 50, 0., 50.);
//TH1F *h_HLT_HH4b=new TH1F("h_HLT_HH4b", "h_HLT_HH4b; Quad - Double", 2, 0, 2);
//TH1F *h_HLT_HH4ball=new TH1F("h_HLT_HH4ball", "h_HLT_HH4ball; Quad - Double", 1, 0, 1);

TH1F *h_MET = new TH1F("h_MET", "MET; MET (GeV)", 100, 0, 200.);

TH1F *h_nJetsAll = new TH1F("h_nJetsAll", "all jets; # All jets", 18, 0., 18.);
TH1F *h_JetsAll_mass = new TH1F("h_JetsAll_mass", "all jets; Jets All m (GeV/c^{2})", 100, 0., 150.);
TH1F *h_JetsAll_pT = new TH1F("h_JetsAll_pT", "all jets; Jets All p_{T} (GeV/c)", 100, 0., 900.);
TH1F *h_JetsAll_eta = new TH1F("h_JetsAll_eta","all jets; Jets All #eta", 100, -4., 4.);
TH1F *h_JetsAll_Phi= new TH1F("h_JetsAll_Phi", "all jets; Jets All #phi", 100, -3.5, 3.5);
TH1F *h_JetsAll_CSV = new TH1F("h_JetsAll_CSV", "all jets; Jets All CSV", 100, 0., 1.);

TH1F *h_nJetsAcc = new TH1F("h_nJetsAcc", "jets in acc.; # jets in Acc", 18, 0., 18.);
TH1F *h_JetsAcc_mass = new TH1F("h_JetsAcc_mass", "jets in acc.; Jets in Acc m (GeV/c^{2})", 100, 0., 150.);
TH1F *h_JetsAcc_pT = new TH1F("h_JetsAcc_pT", "jets in acc.; Jets in Acc p_{T} (GeV/c)", 100, 0., 900.);
TH1F *h_JetsAcc_eta = new TH1F("h_JetsAcc_eta","jets in acc.; Jets in Acc #eta", 100, -4., 4.);
TH1F *h_JetsAcc_Phi= new TH1F("h_JetsAcc_Phi", "jets in acc.; Jets in Acc #phi", 100, -3.5, 3.5);
TH1F *h_JetsAcc_CSV = new TH1F("h_JetsAcc_CSV", "jets in acc.; Jets in Acc CSV", 100, 0., 1.);
//TH1F *h_JetsAcc_Centr = new TH1F("h_Jets_Centrality", "Jets_Centrality; allJets C", 100, 0., 1.);
//TH1F *h_JetsAcc_HT = new TH1F("h_Jets_HT", "h_Jets_HT; allJets HT (GeV/c)", 50, 0., 900.);

TH1F *h_nfJets = new TH1F("h_nfJets", "# Final Jets; # jets", 5, 0., 5.);
TH1F *h_fJets_mass = new TH1F("h_fJets_mass", "final jets; Final Jets m (GeV/c^{2})", 100, 0., 100.);
TH1F *h_fJets_pT = new TH1F("h_fJets_pT", "final jets; Final Jets p_{T} (GeV/c)", 100, 0., 700.);
TH1F *h_fJets_eta = new TH1F("h_fJets_eta","final jets; Final Jets #eta", 100, -4., 4.);
TH1F *h_fJets_CSV = new TH1F("h_fJets_CSV", "final jets; Final Jets CSV", 100, 0., 1.);
TH1F *h_fJets_Centr = new TH1F("h_fJets_Centr", "final jets; Final Jets C", 100, 0., 1.);

//sorted in CSV ----
TH1F *h_fJet1_mass = new TH1F("h_fJet1_mass", "Jet 1; jet1 m (GeV/c^{2})", 100, 0., 100.);
TH1F *h_fJet2_mass = new TH1F("h_fJet2_mass", "Jet 2; jet2 m (GeV/c^{2})", 100, 0., 100.);
TH1F *h_fJet3_mass = new TH1F("h_fJet3_mass", "Jet 3; jet3 m (GeV/c^{2})", 100, 0., 100.);
TH1F *h_fJet4_mass = new TH1F("h_fJet4_mass", "Jet 4; jet4 m (GeV/c^{2})", 100, 0., 100.);
TH1F *h_fJet1_pT = new TH1F("h_fJet1_pT", "Jet 1; jet1 p_{T} (GeV/c)", 100, 0., 700.);
TH1F *h_fJet2_pT = new TH1F("h_fJet2_pT", "Jet 2; jet2 p_{T} (GeV/c)", 100, 0., 700.);
TH1F *h_fJet3_pT = new TH1F("h_fJet3_pT", "Jet 3; jet3 p_{T} (GeV/c)", 100, 0., 700.);
TH1F *h_fJet4_pT = new TH1F("h_fJet4_pT", "Jet 4; jet4 p_{T} (GeV/c)", 100, 0., 700.);
TH1F *h_fJet1_Eta=new TH1F("h_fJet1_Eta", "Jet 1; jet1 #eta", 100, -4., 4.);
TH1F *h_fJet2_Eta=new TH1F("h_fJet2_Eta", "Jet 2; jet2 #eta", 100, -4., 4.);
TH1F *h_fJet3_Eta=new TH1F("h_fJet3_Eta", "Jet 3; jet3 #eta", 100, -4., 4.);
TH1F *h_fJet4_Eta=new TH1F("h_fJet4_Eta", "Jet 4; jet4 #eta", 100, -4., 4.);
TH1F *h_fJet1_CSV=new TH1F("h_fJet1_CSV", "Jet 1; jet1 CSV", 100, 0., 1.);
TH1F *h_fJet2_CSV=new TH1F("h_fJet2_CSV", "Jet 2; jet2 CSV", 100, 0., 1.);
TH1F *h_fJet3_CSV=new TH1F("h_fJet3_CSV", "Jet 3; jet3 CSV", 100, 0., 1.);
TH1F *h_fJet4_CSV=new TH1F("h_fJet4_CSV", "Jet 4; jet4 CSV", 100, 0., 1.);
TH1F *h_fJets3avg_CSV=new TH1F("h_fJets3avg_CSV", "first 3 Jets; avg 3 jets CSV", 100, 0., 1.);
TH1F *h_fJets3min_CSV=new TH1F("h_fJets3min_CSV", "first 3 Jets; min 3 jets CSV", 100, 0., 1.);
TH1F *h_fJets4avg_CSV=new TH1F("h_fJets4avg_CSV", "4 Jets; avg 4 jets CSV", 100, 0., 1.);

TH1F *h_naJets = new TH1F("h_naJets", "additional jets; # additional jets", 18, 0., 18.);
TH1F *h_aJets_mass = new TH1F("h_aJets_mass", "additional jets; additional Jets m (GeV/c^{2})", 100, 0., 150.);
TH1F *h_aJets_pT = new TH1F("h_aJets_pT", "additional jets; additional Jets p_{T} (GeV/c)", 100, 0., 900.);
TH1F *h_aJets_eta = new TH1F("h_aJets_eta","additional jets; additional Jets #eta", 100, -4., 4.);
TH1F *h_aJets_CSV = new TH1F("h_aJets_CSV", "additional jets; additional Jets CSV", 100, 0., 1.);

TH1F *h_H_mass=new TH1F("h_H_mass", "all H candidates; allH m (GeV/c^{2})", 200, 0., 1000.);
TH1F *h_H_pT=new TH1F("h_H_pT", "all H candidates; allH p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_H_Phi=new TH1F("h_H_Phi", "all H candidates; allH #phi", 100, -3.5, 3.5);
TH1F *h_H_Eta=new TH1F("h_H_Eta", "all H candidates; allH #eta", 100, -6., 6.);
TH1F *h_H_CosThSt=new TH1F("h_H_CosTheta*", "all H candidates; all H |cos#theta*|", 50, 0., 1.);

TH1F *h_H1_mass=new TH1F("h_H1_mass", "H1 candidates; H1 m (GeV/c^{2})", 200, 0., 1000.);
TH1F *h_H1_pT=new TH1F("h_H1_pT", "H1 candidates; H1 p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_H1_Phi=new TH1F("h_H1_Phi", "H1 candidates; H1 #phi", 100, -3.5, 3.5);
TH1F *h_H1_Eta=new TH1F("h_H1_Eta", "H1 candidates; H1 #eta", 100, -6., 6.);
TH1F *h_H1_CosThSt=new TH1F("h_H1_CosTheta*", "H1 candidates; H1 |cos#theta*|", 50, 0., 1.);
TH1F *h_H1_deltaR = new TH1F("h_H1_DR","H1 candidates; H1 #DeltaR", 100, 0., 7.);
TH1F *h_H1_deltaPhi = new TH1F("h_H1_DPhi","H1 candidates; H1 #Delta#phi", 100, -3.5, 3.5);
TH1F *h_H1_deltaEta = new TH1F("h_H1_DEta","H1 candidates; H1 #Delta#eta", 100, -6., 6.);
TH1F *h_H1_deltaPhi_abs = new TH1F("h_H1_DPhi_abs","H1 candidates; H1 |#Delta#phi|", 50, 0., 3.5);
TH1F *h_H1_deltaEta_abs = new TH1F("h_H1_DEta_abs","H1 candidates; H1 |#Delta#eta|", 50, 0., 6.);
TH2F *h_H1_deltaPhiVSpT = new TH2F("h_H1_deltaPhiVSpT", "H1 candidates; H1 #Delta#phi; H1 p_{T} (GeV/c)", 200, 0., 3.5, 200, 0., 600.);

TH1F *h_H2_mass=new TH1F("h_H2_mass", "H2 candidates; H2 m (GeV/c^{2})", 200, 0., 1000.);
TH1F *h_H2_pT=new TH1F("h_H2_pT", "H2 candidates; H2 p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_H2_Phi=new TH1F("h_H2_Phi", "H2 candidates; H2 #phi", 100, -3.5, 3.5);
TH1F *h_H2_Eta=new TH1F("h_H2_Eta", "H2 candidates; H2 #eta", 100, -6., 6.);
TH1F *h_H2_CosThSt=new TH1F("h_H2_CosTheta*", "H2 candidates; H2 |cos#theta*|", 50, 0., 1.);
TH1F *h_H2_deltaR = new TH1F("h_H2_DR","H2 candidates; H2 #DeltaR", 100, 0., 7.);
TH1F *h_H2_deltaPhi = new TH1F("h_H2_DPhi","H2 candidates; H2 #Delta#phi", 100, -3.5, 3.5);
TH1F *h_H2_deltaEta = new TH1F("h_H2_DEta","H2 candidates; H2 #Delta#eta", 100, -6., 6.);
TH1F *h_H2_deltaPhi_abs = new TH1F("h_H2_DPhi_abs","H2 candidates; H2 |#Delta#phi|", 50, 0., 3.5);
TH1F *h_H2_deltaEta_abs = new TH1F("h_H2_DEta_abs","H2 candidates; H2 |#Delta#eta|", 50, 0., 6.);
TH2F *h_H2_deltaPhiVSpT = new TH2F("h_H2_deltaPhiVSpT", "H2 candidates; H2 #Delta#phi; H2 p_{T} (GeV/c)", 200, 0., 3.5, 200, 0., 600.);

TH1F *h_HH_mass = new TH1F("h_HH_mass"," di-Higgs candidates; HH m (GeV/c^{2})" , 100, 0., 1500.); 
TH1F *h_HH_pT=new TH1F("h_HH_pT", "di-Higgs candidates; HH p_{T} (GeV/c)", 50, 0., 900.);
TH1F *h_HH_Eta=new TH1F("h_HH_Eta", "di-Higgs candidates; HH #eta", 100, -6., 6.);
TH1F *h_HH_Phi=new TH1F("h_HH_Phi", "di-Higgs candidates; HH #phi", 100, -3.5, 3.5);
TH1F *h_HH_deltaR = new TH1F("h_HH_DR","di-Higgs candidates; HH #DeltaR", 100, 0., 7.);
TH1F *h_HH_deltaPhi = new TH1F("h_HH_DPhi","di-Higgs candidates; HH #Delta#phi", 100, -3.5, 3.5);
TH1F *h_HH_deltaEta = new TH1F("h_HH_DEta","di-Higgs candidates; HH #Delta#eta", 100, -6., 6.);
TH1F *h_HH_deltaPhi_abs = new TH1F("h_HH_DPhi_abs","di-Higgs candidates; HH |#Delta#phi|", 50, 0., 3.5);
TH1F *h_HH_deltaEta_abs = new TH1F("h_HH_DEta_abs","di-Higgs candidates; HH |#Delta#eta|", 50, 0., 6.);
//TH1F *h_HH_CSV = new TH1F("h_HH_CSV","di-Higgs candidates; Sum CSV | between the two higgs ", 70, -4., 4.); 
TH2F *h_H1_H2_mass = new TH2F("h_H1_H2_mass", " mh mh; m_{H_{lead}} (GeV/c^{2}); m_{H_{sublead}} (GeV/c^{2})", 300, 0., 600., 300, 0., 600.);

TH1F *h_HH_massInReg = new TH1F("h_HH_massInReg","di-Higgs candidates in SR; m (GeV/c^{2}) SR" , 200, 0., 1500.);
TH1F *h_HH_pTInReg=new TH1F("h_HH_pTInReg", "di-Higgs candidates in SR; p_{T} (GeV/c) SR", 50, 0., 900.);
TH2F *h_H1_H2_massInReg = new TH2F("h_H1_H2_massInReg", " all comb if min(#sigma (#delta m)^{2}) SR; m (GeV/c^{2}) SR; m (GeV/c^{2}) SR ", 300, 0., 600., 300, 0., 600.);
TH2F *h_H1_H2_massInReg2 = new TH2F("h_H1_H2_massInReg2", " all comb if min(#sigma (#delta m)^{2}) SR2; m (GeV/c^{2}) SR; m (GeV/c^{2}) SR ", 300, 0., 600., 300, 0., 600.);
TH2F *h_H1_H2_massInReg3 = new TH2F("h_H1_H2_massInReg3", " all comb if min(#sigma (#delta m)^{2}) SR3; m (GeV/c^{2}) SR; m (GeV/c^{2}) SR ", 300, 0., 600., 300, 0., 600.);
TH1F *h_HH_mass_diff_cand=new TH1F("h_HH_mass_diff_cand", "|#Deltam| between Higgs masses - all candidates", 50, 0., 200.);
TH1F *h_HH_massNorm_diff=new TH1F("h_HH_massNorm_diff", "|#Deltam| between Higgs masses", 50, 0., 2.);
//  TH1F *h_HHmass_right = new TH1F("h_HHmass_right"," h_HHmass right" , 100, 0., 1000.); 
//  TH1F *h_HHmass_wrong = new TH1F("h_HHmass_wrong"," h_HHmass wrong" , 100, 0., 1000.);
//  TH1F *h_HHmass_first = new TH1F("h_HHmass_first"," h_HHmass first" , 100, 0., 1000.);
//  TH1F *h_HHmass_other = new TH1F("h_HHmass_other"," h_HHmass other" , 100, 0., 1000.);		

//generator level histos
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

//new histos
TH1F *  h_jet1b_dr = new TH1F("h_jet1b_dr","h_jet1b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet2b_dr = new TH1F("h_jet2b_dr","h_jet2b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet3b_dr = new TH1F("h_jet3b_dr","h_jet3b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet4b_dr = new TH1F("h_jet4b_dr","h_jet4b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet4b_drMatrix = new TH1F("h_jet4b_drMatrix","h_jet4b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet4b_drAll = new TH1F("h_jet4b_drAll","h_jet4b_dr; jet-parton #DeltaR", 100, 0., 1.);
TH1F *  h_jet4b_drNotMatched = new TH1F("h_jet4b_drNotMatched","h_jet4b_dr; jet-parton #DeltaR", 100, 0., 1.);

//da eliminare....
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

TH1F *h_Jets_CSV_all=new TH1F("h_Jets_CSV_all", "Jets_CSV_all; all jets CSV", 50, 0., 1.); //debug--

//--------------------------------------
