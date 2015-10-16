// Code for non resonant HH->bbbb in CMS Run2 
// Step1 - KINEMATIC SELECTION
//  Author: Martino Dall'Osso
//   thanks Souvik Das (Univ. of Florida) & Caterina Vernieri (FNAL) 
//    Sept 21th, 2015
//----------------------------------------------------------------

#define hh4b_kinSel_C

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
//-------------------------
bool yesTrg_hh4bAll = true;
bool yesTrg_hh4bQuadTriple = false;
bool yesTrg_hh4bQuadDouble = false; //for trg efficiency
bool yesTrg_hh4bDoubleTriple = false;
bool yesTrg_hh4bDoubleDouble = false; //for trg efficiency
bool yesTrg_hh4bFinal = false; //only (QuadTriple || DoubleTriple)

bool isSignalMC = false; //To take only second part with RL method!!

double Jet_pt_cut_low = 20.; //30.;
double deltaRCut = 0.5; //0.5
int nJets_cut = 4; //3
double Jet_eta_cut = 10; //2.5 - 10=noCut
double CSV_cut = 0.605;  //Run2: low 0.605; medium 0.890; high 0.970. Run1: low 0.679 (used in the first presentation).

//for withinRegion
double H_mass = 115.0;
double dijetM_cut_low = 100.;
double dijetM_cut_high = 150.;
// matrix parameters  --> needed?...
static const int binPt = 10;
static const int binCSV = 8;
static const int binEta = 4;

// environment
//---------------
static const std::string utilsFld="../utils/";
static const std::string dataFld="../data/";
static const std::string plotsFld="plots/";
static const std::string matrixFld="jetsMatch_matrix/";
static const std::string subdataFld=dataFld+"vhbb_heppy_V13/";

// struct and common functions
#include "../utils/hh4bStructs.h"
#include "../utils/HelperFunctions.h"

#include "hh4b_13TeV_kinSel.h" //Histos are here

//just to have order
//--------------------
class hh4b_kinSel{
  public:
    hh4b_kinSel(std::string , bool , std::string , int =2, int =0, std::string ="");
    ~hh4b_kinSel();

    void setBranches(bool, TTree*);   
    bool readMatrix( std::string, std::string );
    double selectBestDiJets (int );
    void fillHistos(bool );   
    void writeHistos(std::string, bool );   
    string jet4SelectionMethod(int );   
};


//-------------
//   MAIN
//-------------
hh4b_kinSel::hh4b_kinSel(std::string sample, bool isData, std::string opt, int finalIndex, int maxEvents , std::string MCsample_RL)
//final index -> to decide which method to use for the 4th jet matching: 0 MCTruth, 1 matrix, 2 CSV, 3 minMass
//MCsample_RL -> MC sample used to fill RL matrix
{
  //Histos
  h_nPV->Sumw2();
  //h_nPV_weighted->Sumw2();

  std::string inputfilename;  
  inputfilename = subdataFld+"tree_"+sample+".root"; 

  TFile * f = new TFile(inputfilename.c_str());
  f->cd();	
  TTree *tree=(TTree*)f->Get("tree");
  std::cout<<"Opened input file "<<inputfilename<<std::endl;

  //retrieve variables
  setBranches(isData, tree);

  //std::cout<<Jet_pt[0]<<"  "<<Jet_pt[1]<<"  "<<Jet_pt[2]<< "   "<<std::endl;

  //output file
  std::string eve = "";
  if(maxEvents > 0) eve = "_" + to_string(maxEvents);
  std::string outfilename=dataFld+"tree_S1_"+sample+eve+".root";
  TFile *outfile=new TFile(outfilename.c_str(), "recreate");

  //output tree and new branches
  TTree *outtree=tree->CloneTree(0);
  /*int H1jet1_i, H1jet2_i;
  int H2jet1_i, H2jet2_i;
  outtree->Branch("H1jet1_i", &H1jet1_i, "H1jet1_i/I");
  outtree->Branch("H1jet2_i", &H1jet2_i, "H1jet2_i/I");
  outtree->Branch("H2jet1_i", &H2jet1_i, "H2jet1_i/I");
  outtree->Branch("H2jet2_i", &H2jet2_i, "H2jet2_i/I");
  outtree->Branch("Ht", &ht, "ht/F");*/

  int nEvents;
  if(maxEvents>0) nEvents=maxEvents;
  else nEvents=tree->GetEntries();

  //read Matrix to select 4th jet (only if RL method is selected)
  if(!readMatrix(MCsample_RL, opt)) {
    if(finalIndex==1)return;
  }

  //------------------  
  // Loop over events
  //------------------
  cout << "nEvents: " << nEvents << endl;
  for (int i=0; i<nEvents; ++i){
    bool jet_inAcc [30];
    bool foundHH=false; 
    ++nCut0;
    tree->GetEvent(i); //READ EVENT
    nJets_InAcc = 0;
    h_nJets->Fill(nJet);

    if(isSignalMC && i<=0.5*nEvents) continue; // to use only second part of the sample if Signal MC //to not run on events used to fill the matrix
    else {    
      vector<TLorentzVector> genB_P, genBISR_P;
//      if(!isData){  //to check jet matching
        for(int p=0; p<4; p++){     
          TLorentzVector b;
          //TLorentzVector bISR;
          //bISR.SetPtEtaPhiM(GenBfH_pt[p],GenBfH_eta[p],GenBfH_phi[p],GenBfH_mass[p]);
          b.SetPtEtaPhiM(GenBfH_pt[p],GenBfH_eta[p],GenBfH_phi[p],GenBfH_mass[p]);
          genB_P.push_back(b);     
          //genBISR_P.push_back(b);     
        }
  //    }

      //-----------------------
      // 1.cut on trigger - new trigger path
      //-----------------------
      if(yesTrg_hh4bAll){ 
        if(HLT_HH4bAll) ++nCut1;
        else continue;  //break loop if not right trigger path 
      }
      else if(yesTrg_hh4bFinal){ 
        if(HLT_BIT_QuadTriple || HLT_BIT_DoubleTriple) ++nCut1;
        else continue;  //break loop if not right trigger path 
      } 
      else if(yesTrg_hh4bQuadTriple){ 
        if(HLT_BIT_QuadTriple) ++nCut1;
        else continue;  //break loop if not right trigger path 
      } 
      else if(yesTrg_hh4bDoubleTriple){ 
        if(HLT_BIT_DoubleTriple) ++nCut1;
        else continue;  //break loop if not right trigger path 
      } 
      h_nJets_InAcc->Fill(nJets_InAcc); //weightPU? - debug

      //-----------------------
      // 2.cut on vtype
      //-----------------------
      if(1) ++nCut2; //debug - to be implemented..
      else continue;

      //----------------------------------------------
      // 3.cut on number of Jets in acceptance region
      //----------------------------------------------
      for (int j=0; j<nJet; ++j){
        if ((fabs(Jet_eta[j])<Jet_eta_cut) && (Jet_pt[j]>Jet_pt_cut_low) ) { //debug //&& Jet_puId[j]>0 && && (Jet_id[j]>0)
            ++nJets_InAcc;
            jet_inAcc[j] = true;
        }
        else jet_inAcc[j] = false;
        Jet je;         
        je.CSV = Jet_btagCSV[j];
        je.pT  = Jet_pt[j];
        je.mass = Jet_mass[j];
        je.eta = Jet_eta[j];
        je.phi = Jet_phi[j];
        jets_all.push_back(je);  //all jets after trigger
      }
      if(nJets_InAcc>=nJets_cut){ 
        ++nCut3;
        //----------------------------------------
        // fill Jet vector with jet sorted in pT
        //----------------------------------------
	for (int j=0; j<nJet; ++j){ //loop on all the jets (only inAcc are saved)
          if(jet_inAcc[j]){
            Jet je;         
            je.CSV = Jet_btagCSV[j];
            //je.CMVA = Jet_btagCMVA[j];
            je.pT  = Jet_pt[j];
            je.mass = Jet_mass[j];
            je.eta = Jet_eta[j];
            je.phi = Jet_phi[j];
            jets_inAcc.push_back(je);  
            //h_Jet4all_pT->Fill(jets_inAcc[j].pT); //debug
            //h_Jet4all_eta->Fill(jets_inAcc[j].eta); //debug
            //h_Jet4all_CSV->Fill(jets_inAcc[j].CSV);//debug
          }              
        }

        //----------------------------------------------
        // 4. jet sort in CSV + CSV cut on first 3 jets
        //----------------------------------------------
        std::sort (jets_inAcc.begin(), jets_inAcc.end(), cmp_CSV );      
        if(jets_inAcc[2].CSV>CSV_cut){ //cut on jets sorted in cmva
  	  ++nCut4;
           /*
           //debug - check sorting --> OK!!
            std::cout<< "size " << jet.size() << " ";                
            for (int j=0; j<jet.size(); ++j){  
             // std::cout<<jet[j].CSV << " ";
              //std::cout<<jet[j].pT << " ";
            }
            cout << endl; */

          //fill jets vector (in acc sorted by csv):
          int NJetInAcc = jets_inAcc.size();
          for(int l=0; l<NJetInAcc;l++){
            jets_inAcc_P.push_back(get_jetVector(&jets_inAcc[l]));  //TLorentz - same order of jet_inAcc!!
            met.SetPtEtaPhiM(met_pt,met_eta,met_phi,met_mass);
          }
          //.............

          double dR=0;
          //check matching if MC
          //-------------------------
          if(!isData){ 
            // Match first 3 jets with b and fill deltaR
            //------------
            vector<TLorentzVector> genB_P_match;            
            bool noJetMatch = false;
            std::vector<TLorentzVector>::iterator index;
            for(int l=0; l<3;l++){
              double dRFin = deltaRCut;
              for(std::vector<TLorentzVector>::iterator it = genB_P.begin() ; it != genB_P.end(); ++it){
                dR = jets_inAcc_P[l].DeltaR(*it);
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
            for(int j=3; j<NJetInAcc; ++j){
              dR = jets_inAcc_P[j].DeltaR(genB_P[0]);
              h_jet4b_drAll->Fill(dR,1./(NJetInAcc-3));
              if(dR<dRFinal) {
                  dRFinal=dR;
                  InTrue = j; //indice del jet piu' vicino al partone rimasto           
              }
            }
            if(dRFinal==deltaRCut) noJetMatch = true;
            if(noJetMatch) continue; //skip event if there are no match for all the 4 jets  -- debug !!!??? still valid?...

            h_jet4b_dr->Fill(dRFinal);
            //fill deltaR for all 4th jets not matched
            for(int j=3; j<NJetInAcc; ++j){
              dR = jets_inAcc_P[j].DeltaR(genB_P[0]);
              if(j!=InTrue) h_jet4b_drNotMatched->Fill(dR,1./(NJetInAcc-4));
            }

            // fill discriminating variables
            //------------
            for(int j=3; j<NJetInAcc; ++j){
              double DpT = jets_inAcc[2].pT - jets_inAcc[j].pT;  //debug --- jet has same ordering of jets_P ?..
              double DCSV = jets_inAcc[2].CSV - jets_inAcc[j].CSV;
              double Deta = fabs(jets_inAcc[2].eta) - fabs(jets_inAcc[j].eta);
              if(j==InTrue){
                h_Jet4match_pT->Fill(jets_inAcc[j].pT);
                h_Jet4match_eta->Fill(jets_inAcc[j].eta);
                h_Jet4match_CSV->Fill(jets_inAcc[j].CSV);
                h_Jet4match_DpT3->Fill(DpT);
                h_Jet4match_DCSV3->Fill(DCSV);
                h_Jet4match_Deta3->Fill(Deta);
              }
              else{
               // h_Jet4all_pT->Fill(jet[j].pT);  debug
               // h_Jet4all_eta->Fill(jet[j].eta);
               // h_Jet4all_CSV->Fill(jet[j].CSV);
                h_Jet4all_DpT3->Fill(DpT);
                h_Jet4all_DCSV3->Fill(DCSV);
                h_Jet4all_Deta3->Fill(Deta);
              }
            }
          }//!isData

          //select 4th jet and match di-jets
          //-----------------------------------------
          double mAve;
          // with RL method
          //----------------
          if(finalIndex == 1){
            //read matrix and assign max R
            double Rmax = 0; 
            for(int j=3; j<NJetInAcc; ++j){
              int ipT = jets_inAcc[j].pT*(binPt-1)/300;
              if(ipT>binPt-1) ipT = binPt-1;
              int iCSV = jets_inAcc[j].CSV*binCSV;
              if(iCSV<0) iCSV =0;
              if(iCSV>(binCSV-1)) iCSV = binCSV-1;
              int ieta = fabs(jets_inAcc[j].eta*2); // //(binEta-1)/2 jet[j].eta*(binEta-1)/2;
              //  if(ieta>binEta-1) ieta = binEta-1;
              //  if(ieta>binEta-1) ieta = binEta-1;
              if(ieta<2) ieta = 0;
              else if(ieta<3) ieta = 1;
              else if(ieta<4) ieta = 2;
              else ieta = 3;
              double thisR = R[ipT][iCSV][ieta];
              //if(i<0.5*nEvents+20){ 
               // cout << j << "  " << thisR << "  " << jets_P[j].DeltaR(genB_P[0]) << endl;
               // cout << j << "  " << jet[j].pT << "  " << ipT << "  " << iCSV << "  " << ieta << endl; 
              //}
              if(thisR>Rmax) {
                Rmax = thisR;
                In=j;
              }
            }                  
            dR = jets_inAcc_P[In].DeltaR(genB_P[0]); //matched?!
            h_jet4b_drMatrix->Fill(dR);
            if(In==InTrue)nRmaxOk++;

	    jet4index = In;
            mAve = selectBestDiJets(jet4index); //choose best di-jets combination
            h_M_matrix->Fill(mAve);
  	    h_HH_mass_matr->Fill(HH.M());
          }

          // with MC truth
          // ----------------------------------       
          if(!isData){
	    jet4index = InTrue;
            nMCTruth++; //1 truth per event..
            mAve = selectBestDiJets(jet4index); //choose best di-jets combination
            h_M_true->Fill(mAve);
  	    h_HH_mass_tr->Fill(HH.M());
          }

          // with CSV
          // ----------------------------------
          jet4index=3;
          if(3==InTrue)nCSVOk++;
          mAve = selectBestDiJets(jet4index); //choose best di-jets combination
          h_M_csv->Fill(mAve);
          h_HH_mass_csv->Fill(HH.M());

          // with minMass
          // ----------------------------------
          jet4index=999;
             //to be implemented....

          //now do it with the selected method (repetition..)        
          foundHH = true;
          HHf++;                  
          jet4SelectionMethod(finalIndex); //to initialize jet4index
          selectBestDiJets(jet4index);     //to chose best di-jets
          //initialize higgs vectors
          for(int i=0;i<4;i++) fJets.push_back(jets_inAcc[fJetsIndex[i]]);
          for(int i=0;i<4;i++) fJets_P.push_back(jets_inAcc_P[fJetsIndex[i]]);
          for(std::vector<Jet>::iterator it = jets_inAcc.begin()+4 ; it != jets_inAcc.end(); ++it){ //additional jets
              aJets.push_back(*it);
          }
//does not work?..they are equal!!
//if((double)fJets_P[0].M() != fJets[0].mass) {
//  errfJets++;
//  cout << (double)fJets_P[0].M() << "  " << fJets[0].mass << endl;
//}
if(jets_inAcc_P[1].M() != fJets_P[1].M()) nMixJets++;

          //angles computation:     
          H1_CosThSt = computeCosThetaStar(H1,HH);
          H2_CosThSt = computeCosThetaStar(H2,HH);
          deltaR_HH = H1.DeltaR(H2);  
          deltaR_H1 = fJets_P[0].DeltaR(fJets_P[1]);
          deltaR_H2 = fJets_P[2].DeltaR(fJets_P[3]);
          deltaPhi_HH = H1.DeltaPhi(H2);  
          deltaPhi_H1 = fJets_P[0].DeltaPhi(fJets_P[1]);
          deltaPhi_H2 = fJets_P[2].DeltaPhi(fJets_P[3]);
          deltaEta_HH = H1.Eta() - H2.Eta();  
          deltaEta_H1 = fJets_P[0].Eta() - fJets_P[1].Eta();
          deltaEta_H2 = fJets_P[2].Eta() - fJets_P[3].Eta();

          //-------------------------  
          // Fill histos with final variables
          //-------------------------
          if (foundHH){	             //debug
             //bool efB1 = false;
	      //bool efB2 = false;
              //----------------
	      //for(std::vector<TLorentzVector>::iterator it = fJets_P.begin() ; it != fJets_P.end(); ++it){
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
	      //}

              //fill the Histos
              fillHistos(isData);
	   
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
	      } // windows*/

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
    } // sample splitting
    jets_all.clear();
    jets_inAcc.clear();
    jets_inAcc_P.clear();
    fJets_P.clear();
    fJets.clear();
  } // Event loop

  if(!isData){
    cout << "'ACCEPTANCE' (Ntrue/Nevents): " << (double)nMCTruth/nEvents*100 << endl;
    cout << "N true jets " << Ntr[0] << endl;
    cout << "nMCTruth " << nMCTruth << endl;
    cout << endl;
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

cout << "# error fJets: " << errfJets << endl;

  //-----------------
  // Some output
  //-----------------
  std::cout<<std::setprecision(2);
  std::cout<<std::fixed;  
  std::cout<< sample << " " << opt << " 4thJet:" << jet4SelectionMethod(finalIndex) << std::endl;
  std::cout<<"00. Number of event with 'mixed' jets = " << nMixJets << std::endl;
  std::cout<<"0. Number of events from input file = "<< nCut0 << std::endl;
  std::cout<<"1. Number of events after trigger   = "<< nCut1 << " ,  " <<(float)nCut2/nCut0*100<<"%  " <<(float)nCut1/nCut0*100<<"%"<<std::endl;
  std::cout<<"2. Number of events after Vtype     = "<< nCut2 << " ,  " <<(float)nCut2/nCut1*100<<"%  " <<(float)nCut2/nCut0*100<<"%"<<std::endl;
  std::cout<<"3. Number of events after nJets_InAcc  = "<< nCut3 << " ,  " <<(float)nCut3/nCut2*100<<"%  " <<(float)nCut3/nCut0*100<<"%"<<std::endl;
  std::cout<<"4. Number of events after bTag cut  = "<< nCut4 << " ,  " <<(float)nCut4/nCut3*100<<"%  " <<(float)nCut4/nCut0*100<<"%"<<std::endl;
  std::cout<<"Number of events after HHfound     = "<< HHf << " ,  " <<(float)HHf/nCut4*100<<"%  " <<(float)HHf/nCut0*100<<"%"<<std::endl;
//  std::cout<<"Number of events with only diJet1 in mass windows  = "<< (nCut5a-nCut5) << " ,  " << (float)(nCut5a-nCut5)/HHf*100<<"%  " << (float)(nCut5a-nCut5)/nCut0*100<<"%"<<std::endl;
//  std::cout<<"Number of events in mass windows   = "<< nCut5 << " ,  " <<(float)nCut5/(nCut5a-nCut5)*100<<"%  " <<(float)nCut5/nCut0*100<<"%"<<std::endl;
  std::cout<<"Number of events in mass windows 2 (elipse) = "<< nCut5b << " ,  " <<(float)nCut5b/nCut5*100<<"%  " <<(float)nCut5b/nCut0*100<<"%"<<std::endl;
  std::cout<<endl;

/*  std::cout<<"eff B1 for HHfound = "<<effB1<< "  " << (float)effB1/HHf*100 <<"%"<<std::endl;
  std::cout<<"eff B2 for HHfound = "<<effB2<< "  " << (float)effB2/HHf*100 <<"%"<<std::endl;
  std::cout<<"eff Bh1 in windows  = "<<effBH1_wind<< "  " << (float)effBH1_wind/nCut5*100 <<"%"<<std::endl;
  std::cout<<"eff Bh2 in windows  = "<<effBH2_wind<< "  " << (float)effBH2_wind/nCut5*100 <<"%"<<std::endl;*/

  //write tree and close
  outtree->Write();
  outfile->Close(); 

  //create Hist file and write histos
  std::string histfilename= plotsFld+"Histograms_"+sample+"_"+opt+".root";
  writeHistos(histfilename, isData);
  std::cout<<"Wrote output file "<<histfilename<<std::endl <<std::endl;

  return;
}
//-----------

hh4b_kinSel::~hh4b_kinSel(){
}
//---------------

//choose best di-jets combination, compute avg Higgs mass and initialize Higgs vector
double hh4b_kinSel::selectBestDiJets(int k = 0){
  fJetsIndex.clear(); //new at every call...
  double M12=0,M13=0,M14=0, M23=0, M24=0, M34=0;
            M12 = (jets_inAcc_P[0] + jets_inAcc_P[1]).M();
            M13 = (jets_inAcc_P[0] + jets_inAcc_P[2]).M();
            M14 = (jets_inAcc_P[0] + jets_inAcc_P[k]).M();
            M23 = (jets_inAcc_P[1] + jets_inAcc_P[2]).M();
            M24 = (jets_inAcc_P[1] + jets_inAcc_P[k]).M();
            M34 = (jets_inAcc_P[2] + jets_inAcc_P[k]).M();
            double DM1, DM2, DM3;
            double Higgs_mAve = 9999;
            DM1= fabs(M12 - M34);
            DM2= fabs(M13 - M24);
            DM3= fabs(M14 - M23);
            if(DM1 < DM2 && DM1<DM3) {
              H1 = jets_inAcc_P[0]+jets_inAcc_P[1];
              H2 = jets_inAcc_P[2]+jets_inAcc_P[k];
              fJetsIndex.push_back(0);
              fJetsIndex.push_back(1);
              fJetsIndex.push_back(2);
              fJetsIndex.push_back(k);
              Higgs_mAve=(M12+M34)/2.;
            }
            if(DM2 < DM1 && DM2<DM3) {
              H1 = jets_inAcc_P[0]+jets_inAcc_P[2];
              H2 = jets_inAcc_P[1]+jets_inAcc_P[k];
              fJetsIndex.push_back(0);
              fJetsIndex.push_back(2);
              fJetsIndex.push_back(1);
              fJetsIndex.push_back(k);
              Higgs_mAve=(M13+M24)/2.;
            }
            if(DM3 < DM1 && DM3<DM2) {
              H1 = jets_inAcc_P[0]+jets_inAcc_P[k];
              H2 = jets_inAcc_P[1]+jets_inAcc_P[2];
              fJetsIndex.push_back(0);
              fJetsIndex.push_back(k);
              fJetsIndex.push_back(1);
              fJetsIndex.push_back(2);
              Higgs_mAve=(M14+M23)/2.;
            }
  HH = H1+H2;        
  return Higgs_mAve;
}
//----------------

bool hh4b_kinSel::readMatrix( std::string MCsample_RL, std::string opt){
  //read matrix and calculate R - for 4th jet selection
  for(int ipT=0; ipT<binPt; ipT++){
    for(int iCSV=0; iCSV<binCSV; iCSV++){
      for(int ieta=0; ieta<binEta; ieta++){
        nM[ipT][iCSV][ieta] = 0;
        nA[ipT][iCSV][ieta] = 0;
      }
    }
  }
  ifstream inmatrix;
  std::string fn = utilsFld+matrixFld+"nMnA_"+MCsample_RL+"_"+opt+".asc";
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
      return true;
  }
  else {
    cout << "WARNING: Unable to open matrix file " << fn << endl;      
    return false;
  }
}
//----------------------

std::string hh4b_kinSel::jet4SelectionMethod(int index){
 std::string method;
 if(index==0) {
   jet4index = InTrue;
   method = "MCTruth";
 }
 else if(index==1) {
   jet4index = In;
   method = "RL";
 }
 else if(index==2) {
   jet4index = 3;
   method = "CSV";
 }
 else if(index==3) {
   jet4index = 99; //debug
   method = "minMass";
 }
 return method;
}

void hh4b_kinSel::setBranches(bool isData, TTree* tree){

  //Retrieve variables
  //--------------------
  tree->SetBranchAddress("nprimaryVertices", &(nPV));
  tree->SetBranchAddress("Vtype", &(Vtype_));
//tree->SetBranchAddress("evt",&evt); 
//tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("nJet",&nJet);              
  tree->SetBranchAddress("Jet_id",&Jet_id);            
  tree->SetBranchAddress("Jet_puId",&Jet_puId);          
  tree->SetBranchAddress("Jet_btagCSV",&Jet_btagCSV);       
  //tree->SetBranchAddress("Jet_btagCMVA",&Jet_btagCMVA);      
  tree->SetBranchAddress("Jet_pt",&Jet_pt);            
  tree->SetBranchAddress("Jet_eta",&Jet_eta);           
  tree->SetBranchAddress("Jet_phi",&Jet_phi);           
  tree->SetBranchAddress("Jet_mass",&Jet_mass);          
  /*tree->SetBranchAddress("Jet_rawPt",&Jet_rawPt);            
  tree->SetBranchAddress("Jet_btagProb",&Jet_btagProb);      
  tree->SetBranchAddress("Jet_btagBProb",&Jet_btagBProb);     
  tree->SetBranchAddress("Jet_btagnew",&Jet_btagnew);       
  tree->SetBranchAddress("Jet_btagCSVV0",&Jet_btagCSVV0);     
  tree->SetBranchAddress("Jet_chHEF",&Jet_chHEF);         
  tree->SetBranchAddress("Jet_neHEF",&Jet_neHEF);         
  tree->SetBranchAddress("Jet_chEmEF",&Jet_chEmEF);        
  tree->SetBranchAddress("Jet_neEmEF",&Jet_neEmEF);        
  tree->SetBranchAddress("Jet_chMult",&Jet_chMult);        
  tree->SetBranchAddress("Jet_leadTrackPt",&Jet_leadTrackPt);*/   
  tree->SetBranchAddress("met_pt",&met_pt);            
  tree->SetBranchAddress("met_eta",&met_eta);           
  tree->SetBranchAddress("met_phi",&met_phi);           
  tree->SetBranchAddress("met_mass",&met_mass);       

  if(isData){
    tree->SetBranchAddress("HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v",&HLT_BIT_QuadTriple);
    tree->SetBranchAddress("HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v",&HLT_BIT_QuadDouble);
    tree->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v",&HLT_BIT_DoubleTriple);
    tree->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v",&HLT_BIT_DoubleDouble);
    tree->SetBranchAddress("HLT_HH4bAll",&HLT_HH4bAll);
  }
  else {
    tree->SetBranchAddress("puWeight", &(weightPU));
    tree->SetBranchAddress("HLT_BIT_HLT_QuadJet45_TripleCSV0p5_v",&HLT_BIT_QuadTriple);
    tree->SetBranchAddress("HLT_BIT_HLT_QuadJet45_DoubleCSV0p5_v",&HLT_BIT_QuadDouble);
    tree->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_TripleCSV0p5_v",&HLT_BIT_DoubleTriple);
    tree->SetBranchAddress("HLT_BIT_HLT_DoubleJet90_Double30_DoubleCSV0p5_v",&HLT_BIT_DoubleDouble);
    tree->SetBranchAddress("HLT_HH4bAll",&HLT_HH4bAll);

    /*tree->SetBranchAddress("Jet_hadronFlavour",&Jet_hadronFlavour); 
    tree->SetBranchAddress("Jet_mcPt",&Jet_mcPt);          
    tree->SetBranchAddress("Jet_mcFlavour",&Jet_mcFlavour);     
    tree->SetBranchAddress("Jet_mcMatchId",&Jet_mcMatchId);     
    tree->SetBranchAddress("Jet_mcEta",&Jet_mcEta);         
    tree->SetBranchAddress("Jet_mcPhi",&Jet_mcPhi);         
    tree->SetBranchAddress("Jet_mcM",&Jet_mcM);           
    tree->SetBranchAddress("Jet_corr_JECUp",&Jet_corr_JECUp);     
    tree->SetBranchAddress("Jet_corr_JECDown",&Jet_corr_JECDown);     
    tree->SetBranchAddress("Jet_corr",&Jet_corr); */

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
}
//---------------

void hh4b_kinSel::fillHistos(bool isData){

   for(std::vector<Jet>::iterator it = jets_all.begin() ; it != jets_all.end(); ++it){ //additional jets
     h_JetsAll_pT->Fill((*it).pT);
     h_JetsAll_eta->Fill((*it).eta);
     h_JetsAll_mass->Fill((*it).mass);
     h_JetsAll_CSV->Fill((*it).CSV);
   }
   h_nJetsAll->Fill(jets_all.size());

   for(std::vector<Jet>::iterator it = jets_inAcc.begin() ; it != jets_inAcc.end(); ++it){ //additional jets
     h_JetsAcc_pT->Fill((*it).pT);
     h_JetsAcc_eta->Fill((*it).eta);
     h_JetsAcc_mass->Fill((*it).mass);
     h_JetsAcc_CSV->Fill((*it).CSV);
   }
   h_nJetsAcc->Fill(jets_inAcc.size());

   double C = 0.;
   for(unsigned int i=0;i<fJets.size();i++) {  //4 final jets
     h_fJets_pT->Fill(fJets[i].pT);
     h_fJets_eta->Fill(fJets[i].eta);
     h_fJets_mass->Fill(fJets[i].mass);
     h_fJets_CSV->Fill(fJets[i].CSV);
     //C += (fJets_P[i].Pt/fJets_P[i].E()); //debug!
   }
   h_fJets_Centr->Fill(C/4); //sum over 4 jets
   h_nfJets->Fill(fJets.size());

   for(std::vector<Jet>::iterator it = aJets.begin() ; it != aJets.end(); ++it){ //additional jets
     h_aJets_pT->Fill((*it).pT);
     h_aJets_eta->Fill((*it).eta);
     h_aJets_mass->Fill((*it).mass);
     h_aJets_CSV->Fill((*it).CSV);
   }
   h_naJets->Fill(aJets.size());

   h_fJet1_pT->Fill(fJets_P[0].Pt());
   h_fJet2_pT->Fill(fJets_P[1].Pt());
   h_fJet3_pT->Fill(fJets_P[2].Pt());
   h_fJet4_pT->Fill(fJets_P[3].Pt());
   h_fJet1_Eta->Fill(fJets_P[0].Eta());
   h_fJet2_Eta->Fill(fJets_P[1].Eta());
   h_fJet3_Eta->Fill(fJets_P[2].Eta());
   h_fJet4_Eta->Fill(fJets_P[3].Eta());
   h_fJet1_CSV->Fill(fJets[0].CSV);
   h_fJet2_CSV->Fill(fJets[1].CSV);
   h_fJet3_CSV->Fill(fJets[2].CSV);
   h_fJet4_CSV->Fill(fJets[3].CSV);	
//	      double Jets_avgCSV = (c_jet1_CSV[Ind]+c_jet2_CSV[Ind]+c_jet3_CSV[Ind])/3;    
	      //double Jets_minCSV = std::min({c_jet1_CSV[Ind],c_jet2_CSV[Ind],c_jet3_CSV[Ind]});
//	      h_3Jets_avgCSV->Fill(Jets_avgCSV);
	      //h_3Jets_minCSV->Fill(Jets_minCSV);	
   h_H1_mass->Fill(H1.M());
   h_H1_pT->Fill(H1.Pt());
   h_H1_Eta->Fill(H1.Eta());
   h_H1_Phi->Fill(H1.Phi());
   h_H1_CosThSt->Fill(fabs(H1_CosThSt));
   h_H1_deltaR->Fill(deltaR_H1);
   h_H1_deltaPhi->Fill(deltaPhi_H1);
   h_H1_deltaEta->Fill(deltaEta_H1);
	      h_H1_deltaPhiVSpT->Fill(deltaPhi_H1,H1.Pt());
	      h_H2_mass->Fill(H2.M());
	      h_H2_pT->Fill(H2.Pt());
	      h_H2_CosThSt->Fill(fabs(H2_CosThSt));
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
	      h_HH_mass->Fill(HH.M());
	      h_HH_pT->Fill(HH.Pt());
	      h_HH_Eta->Fill(HH.Eta());
	      h_HH_Phi->Fill(HH.Phi());
              h_H1_H2_mass->Fill(H1.M(), H2.M());  //leading H first
	      h_HH_deltaR->Fill(deltaR_HH);
	      h_HH_deltaPhi->Fill(deltaPhi_HH);
	      h_HH_deltaEta->Fill(deltaEta_HH);
	      //int region=withinRegion(H1_P.M(), H2_P.M(), 17.5, 37.5, 125, 125);	
	      h_MET->Fill(met.E());     
}
//--------------

void hh4b_kinSel::writeHistos(std::string histfilename, bool isData){

  TFile *tFile=new TFile(histfilename.c_str(), "RECREATE");

  h_nJets->Write();
  h_nJets_InAcc->Write();
  h_nPV->Write();
  //h_nPV_weighted->Write();
  //h_HLT_HH4b->Write();
  h_nJetsAll->Write();
  h_JetsAll_mass->Write();
  h_JetsAll_pT->Write();
  h_JetsAll_eta->Write();
  h_JetsAll_CSV->Write();
  h_nJetsAcc->Write();
  h_JetsAcc_mass->Write();
  h_JetsAcc_pT->Write();
  h_JetsAcc_eta->Write();
  h_JetsAcc_CSV->Write();
  h_nfJets->Write();
  h_fJets_mass->Write();
  h_fJets_pT->Write();
  h_fJets_eta->Write();
  h_fJets_CSV->Write();
  h_naJets->Write();
  h_aJets_mass->Write();
  h_aJets_pT->Write();
  h_aJets_eta->Write();
  h_aJets_CSV->Write();
  h_fJet1_pT->Write();
  h_fJet2_pT->Write();
  h_fJet3_pT->Write();
  h_fJet4_pT->Write();
  h_fJet1_Eta->Write();
  h_fJet2_Eta->Write();
  h_fJet3_Eta->Write();
  h_fJet4_Eta->Write();
  h_fJet1_CSV->Write();
  h_fJet2_CSV->Write();
  h_fJet3_CSV->Write();
  h_fJet4_CSV->Write();
  //h_3Jets_avgCSV->Write();
  //h_3Jets_minCSV->Write();
  //h_nJets_90->Write();
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
  h_MET->Write();
  h_Cuts->Write();

  if(!isData){
    h_jet1b_dr->Write();
    h_jet2b_dr->Write();
    h_jet3b_dr->Write();
    h_jet4b_dr->Write();
    h_nGenH->Write();
    h_genH_pT->Write();
    h_genH_mass->Write();
    h_genHH_mass->Write();
    h_genH1_mass->Write();
    h_genH2_mass->Write();
    h_genB1Jets_DR->Write();
    h_genBfH1_DR->Write();
    h_genBfH2_DR->Write();
  }

  tFile->Write();
  tFile->Close();
}


