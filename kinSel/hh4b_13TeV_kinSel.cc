// Code for non resonant HH->bbbb in CMS Run2 
// Step1 - KINEMATIC SELECTION
//  Author: Martino Dall'Osso
//   thanks Souvik Das (Univ. of Florida) & Caterina Vernieri (FNAL) 
//    Sept 21th, 2015
//     g++  `root-config --libs  --cflags` hh4b_13TeV_kinSel.cc -o kinSel
//      ./kinSel BTagCSVRun2015C-D 1 pT20CSVL_def 2
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
#include <vector>

using namespace std;

// selection parameters:f
//-------------------------
bool yesTrg_hh4bFinal = true; //only (QuadTrsiple || DoubleTriple)

bool yesTrg_hh4bQuadTriple = false;
bool yesTrg_hh4bDoubleTriple = false;
bool yesTrg_hh4bQuadDouble = false; //for trg efficiency
bool yesTrg_hh4bDoubleDouble = false; //for trg efficiency
bool yesTrg_hh4bAll = false; //include all the 4th paths

float Jet_pt_cut_low = 20.;
float deltaRCut = 0.5; //0.5 - for matching
int nJets_cut = 4;
float Jet_eta_cut = 10.;
float CSV_cut = 0.605;  //Run2: low 0.605; medium 0.890; high 0.970. Run1: low 0.679 (used in the first presentation).
float nJets_CSVcut = 3;  //number of Jets on which apply CSV cut

// matrix parameters  --> needed?...
static const int binPt = 10;
static const int binCSV = 8;
static const int binEta = 4;

// environment
//---------------
//WARNING: CHANGE ME!!
static const std::string dataframeworkVersionFld="V15/";
static const std::string MCframeworkVersionFld="V14/";

static const std::string utilsFld="../utils/";
static const std::string dataFld="../data/"+dataframeworkVersionFld;
static const std::string subdataFld=dataFld+"officialNtuples/";
static const std::string MCFld="../data/"+MCframeworkVersionFld;
static const std::string subMCFld=MCFld+"officialNtuples/";

static const std::string plotsFld="plots/"+dataframeworkVersionFld;   //dataVersion is the one that rules for next steps
static const std::string matrixFld="jetsMatch_matrix/";

// struct and common functions
//#include "../utils/HelperFunctions.h"
#include "hh4b_13TeV_kinSel.h" //Histos are here

#include "../utils/hh4bStructs.h"

//-------------
int main (int argc, char **argv) {  //DEBUG!!!!
  
   std::string sample_, opt_, MCsample_RL_ = "";
  bool isData_; 
  int finalIndex_ = 2, maxEvents_ = 0;
  
  for(int i=1; i<argc; i++) {
    if(i==1) sample_ = argv[i];
    if(i==2) isData_ = std::atoi(argv[i]);
    if(i==3) opt_ = argv[i];
    if(i==4) finalIndex_ = std::atoi(argv[i]);
    if(i==5) maxEvents_ = std::atoi(argv[i]);
    if(i==6) MCsample_RL_ = argv[i];
  }

  hh4b_kinSel kinSel (sample_,isData_, opt_, finalIndex_, maxEvents_, MCsample_RL_);
  kinSel.dokinSel();
  return 0;
}

hh4b_kinSel::hh4b_kinSel(std::string sam, bool isD, std::string op, int fi, int mx , std::string MC){
 
  sample=sam;
  opt=op;
  MCsample_RL=MC;
  isData=isD; 
  finalIndex=fi;
  maxEvents=mx;

}

//-------------
//   CORE FUNCTION
//-------------
void hh4b_kinSel::dokinSel()
//final index -> to decide which method to use for the 4th jet matching: 0 MCTruth, 1 matrix, 2 CSV, 3 minMass
//MCsample_RL -> MC sample used to fill RL matrix
{
  //Histos
  h_nPV->Sumw2();
  //h_nPV_weighted->Sumw2();

  //read all Files
  if(!readFiles(fch)) return; //check--

  //TTree *tree=(TTree*)fch->LoadTree("tree");

  //retrieve variables
  setBranches(fch);

  //std::cout<<Jet_pt[0]<<"  "<<Jet_pt[1]<<"  "<<Jet_pt[2]<< "   "<<std::endl;

  //output file
  std::string eve = "";
  if(maxEvents > 0) eve = "_" + std::to_string(maxEvents);
  std::string outfilename=dataFld+"tree_Step1_"+sample+"_"+opt+eve+".root";  //dataVersion is the one that rules for next steps
  TFile *outfile=new TFile(outfilename.c_str(), "RECREATE"); //debug

  //output tree with brand new branches
  TTree *outtree = new TTree("tree", "tree_step1");
  createTree(outtree);

  int nEvents;
  if(maxEvents>0) nEvents=maxEvents;
  else nEvents=fch->GetEntries();

  //read Matrix to select 4th jet (only if RL method is selected)
  if(!readMatrix()) {
    if(finalIndex==1)return;
  }

  //------------------  
  // Loop over events
  //------------------
  std::cout << "nEvents: " << nEvents << std::endl;
  for (int i=0; i<nEvents; ++i){
    ++nCut0;

    jets_all.clear();
    jets_inAcc.clear();
    jets_inAcc_P.clear();
    fJets_P.clear();
    fJets_CSV.clear();
    fJets.clear();
    aJets.clear();
    if((i%500000) == 0) std::cout << i <<std::endl;
    bool jet_inAcc [30];
    bool foundHH=false; 
    nJets_InAcc = 0;

    fch->GetEvent(i); //READ EVENT
    h_nJets->Fill(nJet);

    if(finalIndex==1 && !isData && i<=0.5*nEvents) continue; // to use only second part of the sample if Signal MC //to not run on events used to fill the matrix
    else {    
      std::vector<TLorentzVector> genB_P, genBISR_P;
      //if(!isData){  //to check jet matching
        for(int p=0; p<4; p++){     
          TLorentzVector b;
          //TLorentzVector bISR;
          //bISR.SetPtEtaPhiM(GenBfH_pt[p],GenBfH_eta[p],GenBfH_phi[p],GenBfH_mass[p]);
          b.SetPtEtaPhiM(GenBfH_pt[p],GenBfH_eta[p],GenBfH_phi[p],GenBfH_mass[p]);
          genB_P.push_back(b);     
          //genBISR_P.push_back(b);     
        }
      //}

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

      //----------------------------------------
      // fill Jet vector with jet sorted in pT
      //----------------------------------------
      for (int j=0; j<nJet; ++j){ //loop on all the jets (only inAcc are saved)
        Jet je;         
        je.CSV = Jet_btagCSV[j];
        je.pT  = Jet_pt[j];
        je.mass = Jet_mass[j];
        je.eta = Jet_eta[j];
        je.phi = Jet_phi[j];
        jets_all.push_back(je);  //all jets after trigger
        if ((fabs(Jet_eta[j])<Jet_eta_cut) && (Jet_pt[j]>Jet_pt_cut_low) ) {//debug //&& Jet_puId[j]>0 && && (Jet_id[j]>0)
          ++nJets_InAcc; 
          jets_inAcc.push_back(je);     
        }        
      }
      //----------------------------------------------
      // 3.cut on number of Jets in acceptance region
      //----------------------------------------------
      if(nJets_InAcc>=nJets_cut) ++nCut3;
      else continue;

      //----------------------------------------------
      // 4. jet sort in CSV + CSV cut on first 3 jets
      //----------------------------------------------
      std::sort (jets_inAcc.begin(), jets_inAcc.end() );      // by csv - defined in declaration , cmp_CSV
      if(jets_inAcc[nJets_CSVcut-1].CSV>CSV_cut) ++nCut4; //cut on jets sorted in btag
      else continue;
           /*c
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
//std::cout << jets_inAcc[l].CSV << std::endl;
            jets_inAcc_P.push_back(get_jetVector(&jets_inAcc[l]));  //TLorentz - same order of jet_inAcc!!
            met.SetPtEtaPhiM(met_pt,met_eta,met_phi,met_mass);
          }
          //.............

//debug--- move....
          float dR=0;
          //check matching if MC --- DEBUG!!!!!!!!!!
          //-------------------------
          if(!isData){ 
            // Match first 3 jets with b and fill deltaR
            //------------
            std::vector<TLorentzVector> genB_P_match;            
            bool noJetMatch = false;
            std::vector<TLorentzVector>::iterator index;
            for(int l=0; l<3;l++){
              float dRFin = deltaRCut;
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
            float dRFinal = deltaRCut;
            for(int j=3; j<NJetInAcc; ++j){
              dR = jets_inAcc_P[j].DeltaR(genB_P[0]);
              h_jet4b_drAll->Fill(dR,1./(NJetInAcc-3));
              if(dR<dRFinal) {
                  dRFinal=dR;
                  InTrue = j; //indice del jet piu' vicino al partone rimasto           
              }
            }
            if(dRFinal==deltaRCut) noJetMatch = true;
            // debug - removed....  
            //if(noJetMatch) continue; //skip event if there are no match for all the 4 jets  -- debug !!!??? still valid?...

            h_jet4b_dr->Fill(dRFinal);
            //fill deltaR for all 4th jets not matched
            for(int j=3; j<NJetInAcc; ++j){
              dR = jets_inAcc_P[j].DeltaR(genB_P[0]);
              if(j!=InTrue) h_jet4b_drNotMatched->Fill(dR,1./(NJetInAcc-4));
            }

            // fill discriminating variables
            //------------
            for(int j=3; j<NJetInAcc; ++j){
              float DpT = jets_inAcc[2].pT - jets_inAcc[j].pT;  //debug --- jet has same ordering of jets_P ?..
              float DCSV = jets_inAcc[2].CSV - jets_inAcc[j].CSV;
              float Deta = fabs(jets_inAcc[2].eta) - fabs(jets_inAcc[j].eta);
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
//debug--- move....

          //select 4th jet and match di-jets
          //-----------------------------------------
          float mAve;
          // with RL method
          //----------------
          if(finalIndex == 1){
            //read matrix and assign max R
            float Rmax = 0; 
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
              float thisR = R[ipT][iCSV][ieta];
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
            if(mAve == -1){ ++nCut4b; continue;}; //4th jet CSV<0
            
            h_M_matrix->Fill(mAve);
  	    h_HH_mass_matr->Fill(HH_P.M());
          }

          // with MC truth
          // ----------------------------------       
          if(!isData){
	    jet4index = InTrue;
            nMCTruth++; //1 truth per event..
            mAve = selectBestDiJets(jet4index); //choose best di-jets combination
            if(mAve == -1){ ++nCut4b; continue;}; //4th jet CSV<0
            h_M_true->Fill(mAve);
  	    h_HH_mass_tr->Fill(HH_P.M());
          }

          // with CSV
          // ----------------------------------
          jet4index=3;
          if(3==InTrue)nCSVOk++;
          mAve = selectBestDiJets(jet4index); //choose best di-jets combination
          if(mAve == -1){ ++nCut4b; continue;}; //4th jet CSV<0
          h_M_csv->Fill(mAve);
          h_HH_mass_csv->Fill(HH_P.M());

          // with minMass
          // ----------------------------------
          jet4index=999;
             //to be implemented....

          //now do it with the selected method (repetition..)        
          foundHH = true;
          HHf++;                  
          jet4SelectionMethod(finalIndex); //to initialize jet4index
          selectBestDiJets(jet4index);     //to chose best di-jets
          //if(mAve == -1){ ++nCut4b; continue;}; //4th jet CSV<0 -- not needed now
          //initialize higgs vectors
          for(int i=0;i<4;i++) {
            fJets.push_back(jets_inAcc[fJetsIndex[i]]);
            if(i<3) fJets_CSV.push_back(jets_inAcc[i]);
            else fJets_CSV.push_back(jets_inAcc[fJetsIndex[i]]);
          }
        //debug - make it better
          fJet1_CSV = fJets_CSV.at(0);
          fJet2_CSV = fJets_CSV.at(1);
          fJet3_CSV = fJets_CSV.at(2);
          fJet4_CSV = fJets_CSV.at(3);
          for(int i=0;i<4;i++) fJets_P.push_back(jets_inAcc_P[fJetsIndex[i]]);
          for(std::vector<hh4b_kinSel::Jet>::iterator it = jets_inAcc.begin()+4 ; it != jets_inAcc.end(); ++it){ //additional jets
              aJets.push_back(*it);
          }
          if(jets_inAcc_P[1].M() != fJets_P[1].M()) nMixJets++; //debug - to check
          fJets3avgCSV = (fJets_CSV[0].CSV+fJets_CSV[1].CSV+fJets_CSV[2].CSV)/3;
          fJets3minCSV = std::min({fJets_CSV[0].CSV,fJets_CSV[1].CSV,fJets_CSV[2].CSV});
          fJets4avgCSV = (fJets_CSV[0].CSV+fJets_CSV[1].CSV+fJets_CSV[2].CSV+fJets_CSV[3].CSV)/4;
          Centr = ((fJets_P[0].Pt()/fJets_P[0].E())+(fJets_P[1].Pt()/fJets_P[1].E())+(fJets_P[2].Pt()/fJets_P[2].E())+(fJets_P[3].Pt()/fJets_P[3].E()))/4; //Centrality

          //angles computation:     
          anglesComputation();

          // Fill histos/tree with final variables
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
              fillHistos();
	   
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
              if(withinElipse(H1_P.M(), H2_P.M())){
	        nCut5b++;
                h_H1_H2_massInReg2->Fill(H1_P.M(), H2_P.M());
              }
              else h_H1_H2_massInReg3->Fill(H1_P.M(), H2_P.M());
	      // windows

	      outtree->Fill();	

            }//if FOUND
    } // sample splitting
  } // Event loop

  // Some output (shell & log)  -- debug - put it into a function...
  std::string logfilename;
  if(maxEvents!=0) logfilename = plotsFld+"log_"+sample+"_"+opt+"_"+std::to_string(maxEvents)+".dat";
  else logfilename = plotsFld+"log_"+sample+"_"+opt+".dat";
//  printOutput(logfilename, isData);

// on shell
  //-----------
 // std::cout<<std::setprecision(2);
 // std::cout<<std::fixed;  
  std::cout<< sample << " " << opt << " 4thJet:" << jet4SelectionMethod(finalIndex) << std::endl;
  std::cout<<"0. Input = "<< nCut0 << std::endl;
  std::cout<<"1. Trigger   = "<< nCut1 << " ,  " <<(float)nCut2/nCut0*100<<"%  " <<(float)nCut1/nCut0*100<<"%"<<std::endl;
  //  std::cout<<"2. Number of events after Vtype     = "<< nCut2 << " ,  " <<(float)nCut2/nCut1*100<<"%  " <<(float)nCut2/nCut0*100<<"%"<<std::endl;
  std::cout<<"2. 4 Jets in Acc = "<< nCut3 << " ,  " <<(float)nCut3/nCut2*100<<"%  " <<(float)nCut3/nCut0*100<<"%"<<std::endl;
  std::cout<<"3. 3 Jets CSV = "<< nCut4 << " ,  " <<(float)nCut4/nCut3*100<<"%  " <<(float)nCut4/nCut0*100<<"%"<<std::endl;
  //  std::cout<<"4. 4th CSV>0  = "<< (nCut4-nCut4b) << " ,  " <<(float)(nCut4-nCut4b)/nCut4*100<<"%  " <<(float)(nCut4-nCut4b)/nCut0*100<<"%"<<std::endl;
  std::cout<<"HHfound = "<< HHf << " ,  " <<(float)HHf/nCut4*100<<"%  " <<(float)HHf/nCut0*100<<"%"<<std::endl;
  //  std::cout<<"Number of events with only diJet1 in mass windows  = "<< (nCut5a-nCut5) << " ,  " << (float)(nCut5a-nCut5)/HHf*100<<"%  " << (float)(nCut5a-nCut5)/nCut0*100<<"%"<<std::endl;
  //  std::cout<<"Number of events in mass windows   = "<< nCut5 << " ,  " <<(float)nCut5/(nCut5a-nCut5)*100<<"%  " <<(float)nCut5/nCut0*100<<"%"<<std::endl;
  std::cout<<"In SR (elipse) = "<< nCut5b << " ,  " <<(float)nCut5b/HHf*100<<"%  " <<(float)nCut5b/nCut0*100<<"%"<<std::endl;
  std::cout<<std::endl;

  std::cout<<"# event with 'mixed' jets = " << nMixJets << std::endl;
  std::cout << "# error fJets: " << errfJets << std::endl << std::endl;
  if(!isData){
    std::cout << "'ACCEPTANCE' (Ntrue/Nevents): " << (float)nMCTruth/nEvents*100 << std::endl;
    std::cout << "N true jets " << Ntr[0] << std::endl;
    std::cout << "nMCTruth " << nMCTruth << std::endl;
    std::cout << std::endl;
  } 

/*  std::cout<<"eff B1 for HHfound = "<<effB1<< "  " << (float)effB1/HHf*100 <<"%"<<std::endl;
  std::cout<<"eff B2 for HHfound = "<<effB2<< "  " << (float)effB2/HHf*100 <<"%"<<std::endl;
  std::cout<<"eff Bh1 in windows  = "<<effBH1_wind<< "  " << (float)effBH1_wind/nCut5*100 <<"%"<<std::endl;
  std::cout<<"eff Bh2 in windows  = "<<effBH2_wind<< "  " << (float)effBH2_wind/nCut5*100 <<"%"<<std::endl;*/

  // on file
  //-----------
  std::ofstream ofs (logfilename, std::ofstream::out);
  //ofs<<setprecision(2);
  //ofs<<fixed;  
  ofs<< sample << " " << opt << " 4thJet:" << jet4SelectionMethod(finalIndex) << std::endl;
  ofs<<"Input    "<< nCut0 << std::endl;
  ofs<<"Trigger    "<< nCut1 << "    " <<(float)nCut2/nCut0*100<< "    " <<(float)nCut1/nCut0*100<<std::endl;
  ofs<<"4JetsInAcc    "<< nCut3 << "    " <<(float)nCut3/nCut2*100<< "    " <<(float)nCut3/nCut0*100<<std::endl;
  ofs<<"3JetsCSV    "<< nCut4 << "    " <<(float)nCut4/nCut3*100<< "    " <<(float)nCut4/nCut0*100<<std::endl;
  ofs<<"HH found     "<< HHf   << "    " <<(float)HHf/nCut4*100<< "    " <<(float)HHf/nCut0*100<<std::endl;
  ofs<<"In SR (elipse)    "<< nCut5b << "    " <<(float)nCut5b/HHf*100<< "    " <<(float)nCut5b/nCut0*100<<std::endl;
  ofs<<std::endl;

  if(!isData){
    ofs << "'ACCEPTANCE' (Ntrue/Nevents): " << (float)nMCTruth/nEvents*100 << std::endl;
    ofs << "N true jets " << Ntr[0] <<std::endl;
    ofs << "nMCTruth " << nMCTruth << std::endl;
    ofs << std::endl;
  } 
  ofs.close();
//------------

  //write tree and close
  outtree->Write();
  outfile->Close(); 

  //create Hist file and write histos
  std::string histfilename;
  if(maxEvents!=0) histfilename = plotsFld+"Histograms_"+sample+"_"+opt+"_"+std::to_string(maxEvents)+".root";
  else histfilename = plotsFld+"Histograms_"+sample+"_"+opt+".root";
  writeHistos(histfilename);
  std::cout<<"Wrote output files: "<<std::endl;
  std::cout<<histfilename<<std::endl;
  std::cout<<logfilename<<std::endl;
  std::cout<<outfilename<<std::endl <<std::endl;

}
//-----------

hh4b_kinSel::~hh4b_kinSel(){
}
//---------------

bool hh4b_kinSel::readFiles( TChain * fCh){

  std::string fname, sam_;  
  bool matches = false;
  //inputfilename =  "root://xrootd.ba.infn.it//store/user/arizzi/VHBBHeppyV14/GluGluToBulkGravitonToHHTo4B_M-260_narrow_13TeV-madgraph/VHBB_HEPPY_V14_GluGluToBulkGravitonToHHTo4B_M-260_narrow_13TeV-madgraph__RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/151024_174851/0000/tree_1.root";
  //inputfilename = cms.untracked.vstring("root://xrootd.unl.edu//store/user/arizzi/VHBBHeppyV15/BTagCSV/VHBB_HEPPY_V15_BTagCSV__Run2015D-05Oct2015-v1/151027_145921/0000/tree_1.root");
  if(sample == "BTagCSVRun2015DPromptRecoV4" || sample == "BTagCSVRun2015C-D") {
    sam_ = "BTagCSVRun2015DPromptRecoV4";
    for(int i=1;i<9;i++){ //debug      
      fname = subdataFld+"tree_"+sam_+"_"+std::to_string(i)+".root"; 
      if(i==1) std::cout<<"Adding to Chain files "<< fname << " .." <<std::endl;
      fCh->Add(fname.c_str());
    }  
    matches = true;
  } 
  if(sample == "BTagCSVRun2015DOct05" || sample == "BTagCSVRun2015C-D"){
    sam_ = "BTagCSVRun2015DOct05";
    for(int i=1;i<14;i++){ //debug      
      fname = subdataFld+"tree_"+sam_+"_"+std::to_string(i)+".root"; 
      if(i==1) std::cout<<"Adding to Chain files "<< fname << " .." <<std::endl;
      fCh->Add(fname.c_str());
    }  
    matches = true;
  }
  if(sample == "BTagCSVRun2015C" || sample == "BTagCSVRun2015C-D"){
    sam_ = "BTagCSVRun2015C";
    fname = subdataFld+"tree_"+sam_+".root";       
    fCh->Add(fname.c_str());
    std::cout<<"Adding to Chain file "<< fname << std::endl;
    matches = true;
  }
  if(sample == "M-260"){
    fname = subMCFld+"tree_"+sample+".root"; 
    fCh->Add(fname.c_str());
    std::cout<<"Adding to Chain file "<< fname << std::endl;
    matches = true;
  }
  if(!matches) { std::cout << "ERROR: sample not known" << std::endl; return false; }
  return true;
  //  TFile * f = TFile::Open(inputfilename.c_str());
  //if(!fCh) return; //debug (return..)
  //f->cd();	
}
//---------------

float hh4b_kinSel::selectBestDiJets(int k = 0){ //choose best di-jets combination, compute avg Higgs mass and initialize Higgs vector

  if(jets_inAcc[k].CSV<=0.)  return -1.; //check to avoid underflow in 4th jet.

  fJetsIndex.clear(); //new at every call...
  float M12=0,M13=0,M14=0, M23=0, M24=0, M34=0;
  M12 = (jets_inAcc_P[0] + jets_inAcc_P[1]).M();
  M13 = (jets_inAcc_P[0] + jets_inAcc_P[2]).M();
  M14 = (jets_inAcc_P[0] + jets_inAcc_P[k]).M();
  M23 = (jets_inAcc_P[1] + jets_inAcc_P[2]).M();
  M24 = (jets_inAcc_P[1] + jets_inAcc_P[k]).M();
  M34 = (jets_inAcc_P[2] + jets_inAcc_P[k]).M();
  float DM1, DM2, DM3;
  float Higgs_mAve = 9999;
  DM1= fabs(M12 - M34);
  DM2= fabs(M13 - M24);
  DM3= fabs(M14 - M23);
  if(DM1 < DM2 && DM1<DM3) {
    H1_P = jets_inAcc_P[0]+jets_inAcc_P[1];
    H2_P = jets_inAcc_P[2]+jets_inAcc_P[k];
    fJetsIndex.push_back(0);
    fJetsIndex.push_back(1);
    fJetsIndex.push_back(2);
    fJetsIndex.push_back(k);
    Higgs_mAve=(M12+M34)/2.;
  }
  if(DM2 < DM1 && DM2<DM3) {
    H1_P = jets_inAcc_P[0]+jets_inAcc_P[2];
    H2_P = jets_inAcc_P[1]+jets_inAcc_P[k];
    fJetsIndex.push_back(0);
    fJetsIndex.push_back(2);
    fJetsIndex.push_back(1);
    fJetsIndex.push_back(k);
    Higgs_mAve=(M13+M24)/2.;
  }
  if(DM3 < DM1 && DM3<DM2) {
    H1_P = jets_inAcc_P[0]+jets_inAcc_P[k];
    H2_P = jets_inAcc_P[1]+jets_inAcc_P[2];
    fJetsIndex.push_back(0);
    fJetsIndex.push_back(k);
    fJetsIndex.push_back(1);
    fJetsIndex.push_back(2);
    Higgs_mAve=(M14+M23)/2.;
  }
  //create diJet
  H1 = get_diJet(&H1_P);
  H2 = get_diJet(&H2_P);  
  //create diHiggs
  HH_P = H1_P+H2_P;   
  HH = get_diJet(&HH_P);     

  return Higgs_mAve;
}
//----------------

bool hh4b_kinSel::readMatrix(){
  //read matrix and calculate R - for 4th jet selection
  for(int ipT=0; ipT<binPt; ipT++){
    for(int iCSV=0; iCSV<binCSV; iCSV++){
      for(int ieta=0; ieta<binEta; ieta++){
        nM[ipT][iCSV][ieta] = 0;
        nA[ipT][iCSV][ieta] = 0;
      }
    }
  }
  std::ifstream inmatrix;
  std::string fn = utilsFld+matrixFld+"nMnA_"+MCsample_RL+"_"+opt+".asc";
  inmatrix.open(fn);
  if (inmatrix) {
      int i1, i2, i3;
      float sumnA=0, sumnM=0; 
      for(int ipT=0; ipT<binPt; ipT++){
        for(int iCSV=0; iCSV<binCSV; iCSV++){
          for(int ieta=0; ieta<binEta; ieta++){
            inmatrix >> i1 >> i2 >> i3 >> nM[ipT][iCSV][ieta] >> nA[ipT][iCSV][ieta];
            if(i1!=ipT ||i2!=iCSV ||i3!=ieta) std::cout << "WARNING: matrix mismatch" << std::endl;
            if(nA[ipT][iCSV][ieta]>0) R[ipT][iCSV][ieta] = nM[ipT][iCSV][ieta]/nA[ipT][iCSV][ieta]; 
            else R[ipT][iCSV][ieta] = 1;
            //cout << nM[ipT][iCSV][ieta] << "  " << nA[ipT][iCSV][ieta] <<endl; 
            sumnM+=nM[ipT][iCSV][ieta];
            sumnA+=nA[ipT][iCSV][ieta];
            if(nM[ipT][iCSV][ieta]!=0 || nM[ipT][iCSV][ieta]!=0)  {
              std::cout << ipT << " " << iCSV << " " << ieta << " " << R[ipT][iCSV][ieta] << " +- " << R[ipT][iCSV][ieta]*sqrt(1/nM[ipT][iCSV][ieta]+ 1/nA[ipT][iCSV][ieta]) <<std::endl; 
            }
            else 
              std::cout << "nM or nA null" <<std::endl;
	  }
        }
      }
      Rave=sumnM/sumnA;
      inmatrix.close();
      return true;
  }
  else {
    std::cout << "WARNING: Unable to open matrix file " << fn << std::endl;      
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
//---------------

void hh4b_kinSel::anglesComputation(){
  H1.CosThSt = computeCosThetaStar(H1_P,HH_P);
  H2.CosThSt = computeCosThetaStar(H2_P,HH_P);
  H1.dR = fJets_P[0].DeltaR(fJets_P[1]);
  H2.dR = fJets_P[2].DeltaR(fJets_P[3]);
  H1.dPhi = fJets_P[0].DeltaPhi(fJets_P[1]);
  H2.dPhi = fJets_P[2].DeltaPhi(fJets_P[3]);
  H1.dEta = fJets_P[0].Eta() - fJets_P[1].Eta();
  H2.dEta = fJets_P[2].Eta() - fJets_P[3].Eta();
  H1.dPhi_abs = abs(H1.dPhi);
  H2.dPhi_abs = abs(H2.dPhi);
  H1.dEta_abs = abs(H1.dEta);
  H2.dEta_abs = abs(H2.dEta);
  HH.dR = H1_P.DeltaR(H2_P);  
  HH.dPhi = H1_P.DeltaPhi(H2_P);  
  HH.dEta = H1_P.Eta() - H2_P.Eta();
  HH.dPhi_abs = abs(HH.dPhi);
  HH.dEta_abs = abs(HH.dEta);

}
//---------------

//Retrieve variables
void hh4b_kinSel::setBranches( TChain* tree){
  tree->SetBranchAddress("nprimaryVertices", &(nPV));
  tree->SetBranchAddress("Vtype", &(Vtype_));
  tree->SetBranchAddress("evt",&evt); 
  tree->SetBranchAddress("run",&run);
//  tree->SetBranchAddress("lumi",&lumi);
  tree->SetBranchAddress("nJet",&nJet);              
  tree->SetBranchAddress("Jet_id",&Jet_id);            
  tree->SetBranchAddress("Jet_puId",&Jet_puId);          
  tree->SetBranchAddress("Jet_btagCSV",&Jet_btagCSV);       
  tree->SetBranchAddress("Jet_pt",&Jet_pt);            
  tree->SetBranchAddress("Jet_eta",&Jet_eta);           
  tree->SetBranchAddress("Jet_phi",&Jet_phi);           
  tree->SetBranchAddress("Jet_mass",&Jet_mass);          
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

void hh4b_kinSel::createTree(TTree* outtree){

  outtree->Branch("nprimaryVertices", &(nPV));
  outtree->Branch("Vtype", &(Vtype_));
  outtree->Branch("evt",&evt); 
  outtree->Branch("run",&run);
  outtree->Branch("nJet",&nJet);              
 // outtree->Branch("fJets","std::vector<hh4b_kinSel::Jet>",&fJets);
 // outtree->Branch("aJets","std::vector<hh4b_kinSel::Jet>",&aJets);
//debug
  outtree->Branch("fJet1","hh4b_kinSel::Jet",&fJet1_CSV);
  outtree->Branch("fJet2","hh4b_kinSel::Jet",&fJet2_CSV);
  outtree->Branch("fJet3","hh4b_kinSel::Jet",&fJet3_CSV);
  outtree->Branch("fJet4","hh4b_kinSel::Jet",&fJet4_CSV);
  outtree->Branch("fJets3avgCSV",&fJets3avgCSV);
  outtree->Branch("fJets3minCSV",&fJets3minCSV);
  outtree->Branch("fJets4avgCSV",&fJets4avgCSV);
  outtree->Branch("Centr",&Centr);

  outtree->Branch("H1","diJet",&H1);
  outtree->Branch("H2","diJet",&H2);
  outtree->Branch("HH","diJet",&HH);  //debug
  outtree->Branch("met","TLorentzVector",&met);

  outtree->Branch("fJet3_pT",&fJet3_CSV.pT); //debug
  outtree->Branch("fJet4_CSV",&fJet4_CSV.CSV); //debug
  outtree->Branch("HH_mass",&HH.mass);  //debug
  outtree->Branch("HH_dR",&HH.dR);  //debug
  outtree->Branch("H1_cos",&H1.CosThSt);
  outtree->Branch("H1_pT",&H1.pT);
  outtree->Branch("H1_dR",&H1.dR);
  outtree->Branch("H2_pT",&H2.pT);
  outtree->Branch("H2_dR",&H2.dR);
  outtree->Branch("H1_mass",&H1.mass);
  outtree->Branch("H2_mass",&H2.mass);

  //no Higgs variables (both Higgs together..)
}
//---------------

void hh4b_kinSel::fillHistos(){

   for(std::vector<hh4b_kinSel::Jet>::iterator it = jets_all.begin() ; it != jets_all.end(); ++it){ //additional jets
     h_JetsAll_pT->Fill((*it).pT);
     h_JetsAll_mass->Fill((*it).mass);
     h_JetsAll_eta->Fill((*it).eta);
     h_JetsAll_CSV->Fill((*it).CSV);
   }
   h_nJetsAll->Fill(jets_all.size());

   for(std::vector<hh4b_kinSel::Jet>::iterator it = jets_inAcc.begin() ; it != jets_inAcc.end(); ++it){ //additional jets
     h_JetsAcc_pT->Fill((*it).pT);
     h_JetsAcc_eta->Fill((*it).eta);
     h_JetsAcc_mass->Fill((*it).mass);
     h_JetsAcc_CSV->Fill((*it).CSV);
   }
   h_nJetsAcc->Fill(jets_inAcc.size());

   for(unsigned int i=0;i<fJets.size();i++) {  //4 final jets
     h_fJets_pT->Fill(fJets[i].pT);
     h_fJets_eta->Fill(fJets[i].eta);
     h_fJets_mass->Fill(fJets[i].mass);
     h_fJets_CSV->Fill(fJets[i].CSV);
   }
   h_fJets_Centr->Fill(Centr); //sum over 4 jets
   h_nfJets->Fill(fJets.size());

   for(std::vector<hh4b_kinSel::Jet>::iterator it = aJets.begin() ; it != aJets.end(); ++it){ //additional jets
     h_aJets_pT->Fill((*it).pT);
     h_aJets_eta->Fill((*it).eta);
     h_aJets_mass->Fill((*it).mass);
     h_aJets_CSV->Fill((*it).CSV);
   }
   h_naJets->Fill(aJets.size());
   //final jets sorted in CSV ---
   h_fJet1_mass->Fill(fJets_CSV[0].mass);
    h_fJet2_mass->Fill(fJets_CSV[1].mass);
    h_fJet3_mass->Fill(fJets_CSV[2].mass);
    h_fJet4_mass->Fill(fJets_CSV[3].mass);
    h_fJet1_pT->Fill(fJets_CSV[0].pT);
    h_fJet2_pT->Fill(fJets_CSV[1].pT);
    h_fJet3_pT->Fill(fJets_CSV[2].pT);
    h_fJet4_pT->Fill(fJets_CSV[3].pT);
    h_fJet1_Eta->Fill(fJets_CSV[0].eta);
    h_fJet2_Eta->Fill(fJets_CSV[1].eta);
    h_fJet3_Eta->Fill(fJets_CSV[2].eta);
    h_fJet4_Eta->Fill(fJets_CSV[3].eta);
    h_fJet1_CSV->Fill(fJets_CSV[0].CSV);
    h_fJet2_CSV->Fill(fJets_CSV[1].CSV);
    h_fJet3_CSV->Fill(fJets_CSV[2].CSV);
    h_fJet4_CSV->Fill(fJets_CSV[3].CSV);	
    h_fJets3avg_CSV->Fill(fJets3avgCSV);
    h_fJets3min_CSV->Fill(fJets3minCSV);
    h_fJets4avg_CSV->Fill(fJets4avgCSV);
   h_H1_mass->Fill(H1.mass);
    h_H1_pT->Fill(H1.pT);
    h_H1_Eta->Fill(H1.eta);
    h_H1_Phi->Fill(H1.phi);
    h_H1_CosThSt->Fill(fabs(H1.CosThSt));
    h_H1_deltaR->Fill(H1.dR);
    h_H1_deltaPhi->Fill(H1.dPhi);
    h_H1_deltaEta->Fill(H1.dEta);
    h_H1_deltaPhiVSpT->Fill(H1.dPhi,H1.pT);
   h_H2_mass->Fill(H2.mass);
    h_H2_pT->Fill(H2.pT);
    h_H2_Eta->Fill(H2.eta);
    h_H2_Phi->Fill(H2.phi);
    h_H2_CosThSt->Fill(fabs(H2.CosThSt));
    h_H2_deltaR->Fill(H2.dR);
    h_H2_deltaPhi->Fill(H2.dPhi);
    h_H2_deltaEta->Fill(H2.dEta);
    h_H2_deltaPhiVSpT->Fill(H2.dPhi,H2.pT);
   h_H_mass->Fill(H1.mass);
    h_H_mass->Fill(H2.mass);
    h_H_pT->Fill(H1.pT);
    h_H_pT->Fill(H2.mass);
    h_H_Eta->Fill(H1.eta);
    h_H_Eta->Fill(H2.eta);
    h_H_Phi->Fill(H1.phi);
    h_H_Phi->Fill(H2.phi);
    h_H_CosThSt->Fill(fabs(H1.CosThSt));
    h_H_CosThSt->Fill(fabs(H2.CosThSt));
   h_HH_mass->Fill(HH.mass);
    h_HH_pT->Fill(HH.pT);
    h_HH_Eta->Fill(HH.eta);
    h_HH_Phi->Fill(HH.phi);
    h_H1_H2_mass->Fill(H1.mass, H2.mass);  //leading H first
    h_HH_deltaR->Fill(HH.dR);
    h_HH_deltaPhi->Fill(HH.dPhi);
    h_HH_deltaEta->Fill(HH.dEta);
    //int region=withinRegion(H1_P.mass(), H2_P.mass(), 17.5, 37.5, 125, 125);	
   h_MET->Fill(met.E());     
  
  h_Cuts->Fill(1,nCut0);
  h_Cuts->Fill(3,nCut1);
  //h_Cuts->Fill(5,nCut2);
  h_Cuts->Fill(5,nCut3);
  h_Cuts->Fill(7,nCut4);  
  //h_Cuts->Fill(11,nCut4b);  
  h_Cuts->Fill(9,HHf);
  h_Cuts->Fill(11,nCut5b);
  //  h_Cuts->Fill(17,nCut5);
}
//--------------

void hh4b_kinSel::writeHistos(std::string histfilename){

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
  h_fJets_Centr->Write();
  h_naJets->Write();
  h_aJets_mass->Write();
  h_aJets_pT->Write();
  h_aJets_eta->Write();
  h_aJets_CSV->Write();
  h_fJet1_mass->Write();
  h_fJet2_mass->Write();
  h_fJet3_mass->Write();
  h_fJet4_mass->Write();
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
  h_fJets3avg_CSV->Write();
  h_fJets3min_CSV->Write();
  h_fJets4avg_CSV->Write();
  h_H_mass->Write();
  h_H_pT->Write();
  h_H_Eta->Write();
  h_H_Phi->Write();
  h_H_CosThSt->Write();
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
  h_H2_Eta->Write();
  h_H2_Phi->Write();
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

void hh4b_kinSel::printOutput(std::string logfilename, bool isData){ //improve...

  
}


//get_jetVector --> get TLorentzVector from Jet
TLorentzVector hh4b_kinSel::get_jetVector(Jet* jj){ 
  TLorentzVector jetv;
  jetv.SetPtEtaPhiM(jj->pT,jj->eta,jj->phi,jj->mass);
  return jetv;
}

//cmp_CSV --> Jets sorting by CSV  ---DEBUG!!!!
//bool hh4b_kinSel::cmp_CSV(const Jet& jet1, const Jet& jet2){
  //return jet1.CSV > jet2.CSV;
//}

/*int hh4b_kinSel::dokinSel(std::string sample, bool isData, std::string opt, int finalIndex, int maxEvents , std::string MCsample_RL)
{
   cout << "starting kin selection..." << endl;
   hh4b_kinSel(sample,  isData, opt,  finalIndex,  maxEvents ,  MCsample_RL);
   return 0;
}*/


