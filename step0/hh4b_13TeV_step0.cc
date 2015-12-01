// Code for non resonant HH->bbbb in CMS Run2 
// Step0 - clone of official ntuples -- less branches
//  Author: Martino Dall'Osso
//    Nov 24th, 2015
//     g++  `root-config --libs  --cflags` hh4b_13TeV_step0.cc -o step0
//      ./step0 TTJets 0 0 1
//       doStep0 (std::string sample, bool isData, bool fromLocal, bool isTest) {
//----------------------------------------------------------------

#define hh4b_step0_C

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

// environment
//---------------
//WARNING: CHANGE ME!!
static const std::string dataframVer="V15";
static const std::string MCframVer="V14";

static const std::string dataFld="../data/"+dataframVer+"/";
static const std::string subdataFld=dataFld+"officialNtuples/";
static const std::string MCFld="../data/"+MCframVer+"/";
static const std::string subMCFld=MCFld+"officialNtuples/";

static const std::string step0Fld="../data/Step0/";

//Retrieve variables
void setBranches( TChain* tree, bool isDat, bool isSign){
  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("nprimaryVertices", 1);
  tree->SetBranchStatus("Vtype", 1);
  tree->SetBranchStatus("evt",1);  
  tree->SetBranchStatus("run",1);
  tree->SetBranchStatus("nJet",1);                            
  tree->SetBranchStatus("Jet_id",1);                        
  tree->SetBranchStatus("Jet_puId",1);          
  tree->SetBranchStatus("Jet_btagCSV",1);              
  tree->SetBranchStatus("Jet_pt",1);                        
  tree->SetBranchStatus("Jet_eta",1);                      
  tree->SetBranchStatus("Jet_phi",1);                      
  tree->SetBranchStatus("Jet_mass",1);                    
  tree->SetBranchStatus("Jet_chMult",1);  
  tree->SetBranchStatus("Jet_leadTrackPt",1);  
  tree->SetBranchStatus("met_pt",1);                        
  tree->SetBranchStatus("met_eta",1);                      
  tree->SetBranchStatus("met_phi",1);           
  tree->SetBranchStatus("met_mass",1);              

  if(isDat){
     tree->SetBranchStatus("HLT_BIT_HLT_QuadJet45_TripleBTagCSV0p67_v",1);
     tree->SetBranchStatus("HLT_BIT_HLT_QuadJet45_DoubleBTagCSV0p67_v",1);
     tree->SetBranchStatus("HLT_BIT_HLT_DoubleJet90_Double30_TripleBTagCSV0p67_v",1);
     tree->SetBranchStatus("HLT_BIT_HLT_DoubleJet90_Double30_DoubleBTagCSV0p67_v",1);
     tree->SetBranchStatus("HLT_HH4bAll",1);
  }
  else {
     tree->SetBranchStatus("puWeight", 1);
     tree->SetBranchStatus("HLT_BIT_HLT_QuadJet45_TripleCSV0p5_v",1);
     tree->SetBranchStatus("HLT_BIT_HLT_QuadJet45_DoubleCSV0p5_v",1);
     tree->SetBranchStatus("HLT_BIT_HLT_DoubleJet90_Double30_TripleCSV0p5_v",1);
     tree->SetBranchStatus("HLT_BIT_HLT_DoubleJet90_Double30_DoubleCSV0p5_v",1);
     tree->SetBranchStatus("HLT_HH4bAll",1);
 
     if(isSign){
      tree->SetBranchStatus("nGenHiggsBoson",1);
      tree->SetBranchStatus("GenHiggsBoson_pt",1);
      tree->SetBranchStatus("GenHiggsBoson_mass",1);
      tree->SetBranchStatus("GenHiggsBoson_eta",1);
      tree->SetBranchStatus("GenHiggsBoson_phi",1);
      tree->SetBranchStatus("GenHiggsBoson_status",1);
      tree->SetBranchStatus("nGenBQuarkFromH",1);
      tree->SetBranchStatus("GenBQuarkFromH_pdgId",1);
      tree->SetBranchStatus("GenBQuarkFromH_pt",1);
      tree->SetBranchStatus("GenBQuarkFromH_eta",1);
      tree->SetBranchStatus("GenBQuarkFromH_phi",1);
      tree->SetBranchStatus("GenBQuarkFromH_mass",1);
      tree->SetBranchStatus("GenBQuarkFromH_charge",1);
      tree->SetBranchStatus("GenBQuarkFromH_status",1);
      tree->SetBranchStatus("nGenBQuarkFromHafterISR",1);
      tree->SetBranchStatus("GenBQuarkFromHafterISR_pdgId",1);
      tree->SetBranchStatus("GenBQuarkFromHafterISR_pt",1);
      tree->SetBranchStatus("GenBQuarkFromHafterISR_eta",1);
      tree->SetBranchStatus("GenBQuarkFromHafterISR_phi",1);
      tree->SetBranchStatus("GenBQuarkFromHafterISR_mass",1);
    }
  }
}
//---------------


void getFromT2(std::string , bool fromLoc, bool isDat){
      std::string flistname = subMCFld+"filelist_"+sample;
      std::string in;
      if(!fromLoc){
        if(!isDat) flistname = subMCFld+"filelist_"+sample;
        else flistname = subdataFld+"filelist_"+sample;
      }
      std::ifstream infile;
      std::vector<std::string> rfiles;
      std::string in;
      infile.open(flistname.c_str());
      if(!infile)	{      //check if file exists
	printf( "WARNING: no input file %s \n", flistname.c_str());
        return ;
      }
      while (!infile.eof()) { //debug - check 15
        infile >> in;
        rfiles.push_back(in);
      } 
      for(std::vector<std::string>::iterator it = rfiles.begin() ; it != rfiles.end(); ++it){
        if(isTest && (it - rfiles.begin())==1) break;
        fname = "root://xrootd.ba.infn.it/"+(*it);
        if((it - rfiles.begin()) == 0) std::cout<<"Adding to Chain files "<< fname << " .." <<std::endl;
        fCh->Add(fname.c_str());
      }
      rfiles.clear();
}

//-------------
//   CORE FUNCTION
//-------------
void doStep0 (std::string sample, bool isData, bool fromLocal, bool isTest) {

  bool isSignal= false;
  if(isTest) std::cout << "Processing one file only as Test" << std::endl;

  //read all Files
  // read filenames from file and create chain looking at Local or T2..
  TChain * fCh = new TChain ("tree");

  //data
  //------------
  std::string fname, sam_;  
  bool matches = false;
  if(sample == "BTagCSVRun2015DPromptRecoV4" || sample == "BTagCSVRun2015C-D") {
    sam_ = "BTagCSVRun2015DPromptRecoV4";
    if(fromLocal){
      for(int i=1;i<9;i++){ //debug      
        if(isTest && i>1) break;
        fname = subdataFld+"tree_"+sam_+"_"+std::to_string(i)+".root"; 
        if(i==1) std::cout<<"Adding to Chain files "<< fname << " .." <<std::endl;
        fCh->Add(fname.c_str());
      }  
    }
    else getFromT2(sam_, fromLocal, isData);
    matches = true;
  } 
  if(!isTest){
   if(sample == "BTagCSVRun2015DOct05" || sample == "BTagCSVRun2015C-D"){
    sam_ = "BTagCSVRun2015DOct05";
    if(fromLocal){
      for(int i=1;i<14;i++){ //debug      
      fname = subdataFld+"tree_"+sam_+"_"+std::to_string(i)+".root"; 
      if(i==1) std::cout<<"Adding to Chain files "<< fname << " .." <<std::endl;
      fCh->Add(fname.c_str());
      }  
    }
    else getFromT2(sam_, fromLocal, isData);
    matches = true;
   }
   if(sample == "BTagCSVRun2015C" || sample == "BTagCSVRun2015C-D"){
    sam_ = "BTagCSVRun2015C";
    if(fromLocal){
      fname = subdataFld+"tree_"+sam_+".root";       
      fCh->Add(fname.c_str());
      std::cout<<"Adding to Chain file "<< fname << std::endl;
    }
    else getFromT2(sam_, fromLocal, isData);
    matches = true;
   }  
  }

  //MC
  //------------
  if(sample == "M-260"){
    isSignal = true;
    fname = subMCFld+"tree_"+sample+".root"; 
    fCh->Add(fname.c_str());
    std::cout<<"Adding to Chain file "<< fname << std::endl;
    matches = true;
  }
  else if(sample == "TTJets"){
    if(fromLocal){
      sam_ = "TTJets";
      for(int i=1;i<5;i++){ //debug      
        if(isTest && i>1) break;
        fname = subMCFld+"tree_"+sam_+"_"+std::to_string(i)+".root"; 
        if(i==1) std::cout<<"Adding to Chain files "<< fname << " .." <<std::endl;
        fCh->Add(fname.c_str());
      }
    }
    else getFromT2(sample, fromLocal, isData);
    matches = true;
  }
  else if(strstr(sample.c_str(),"QCD_HT")){
    getFromT2(sample, fromLocal, isData);
    matches = true;
  }
  if(!matches) { std::cout << "ERROR: sample not known" << std::endl; return ; }

  //retrieve variables
  setBranches(fCh, isData, isSignal);

  //output file
  std::string eve = "";
  std::string outfilenameNoSel;
  if(isData){
    if(isTest) outfilenameNoSel=outFld+"tree_Step0_"+dataframVer+"_"+sample+"_Test.root";
    else  outfilenameNoSel=outFld+"tree_Step0_"+dataframVer+"_"+sample+".root";
  }
  else{
    if(isTest) outfilenameNoSel=step0Fld+"tree_Step0_"+MCframVer+"_"+sample+"_Test.root";
    else  outfilenameNoSel=step0Fld+"tree_Step0_"+MCframVer+"_"+sample+".root";
  }

  TFile *outfilenoSel=new TFile(outfilenameNoSel.c_str(), "RECREATE");
  std::cout<<"Start cloning to: "<< outfilenameNoSel << std::endl;
  TTree *outtreeNoSel = fCh->CloneTree();
  outfilenoSel->Close();

  std::cout<<"Wrote output files: "<< outfilenameNoSel<<std::endl <<std::endl;
}
//-----------

int main (int argc, char **argv) {
  
  std::string sample_ ;
  bool isData_, fromLocal_ = true, isTest_ = false; 
  
  for(int i=1; i<argc; i++) {
    if(i==1) sample_ = argv[i];
    if(i==2) isData_ = std::atoi(argv[i]);
    if(i==3) fromLocal_ = std::atoi(argv[i]);
    if(i==4) isTest_ = std::atoi(argv[i]);
  }

  doStep0 (sample_, isData_, fromLocal_, isTest_);
  return  0;
}


