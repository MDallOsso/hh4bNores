
//Struct and Functions for hh4b analysis macros
//---------------------------------------------

float pi=3.14159265;
//for withinRegion
float H_mass = 115.0;
float dijetM_cut_low = 100.;
float dijetM_cut_high = 150.;

//Jet
typedef struct{
  float mass;
  float pT;
  float eta;
  float phi;
  float CosThSt;
  float dR; 
  float dPhi;
  float dEta;
  float dPhi_abs; //modulo
  float dEta_abs; //modulo
} diJet;

//Jet
typedef struct{
  float CSV;
  float mass;
  float pT;
  float eta;
  float phi;
} Jet;

//Jet for plots
typedef struct{
  TH1F* h;
  float norm;
} Jet4Plot;


//get_jetVector --> get TLorentzVector from Jet
diJet get_diJet(TLorentzVector* diJ_P){ 
  diJet diJ;
  diJ.mass = diJ_P->M();
  diJ.pT = diJ_P->Pt();
  diJ.eta = diJ_P->Eta();
  diJ.phi = diJ_P->Phi();
  return diJ;
}

//get_jetVector --> get TLorentzVector from Jet
TLorentzVector get_jetVector(Jet* jj){ 
  TLorentzVector jetv;
  jetv.SetPtEtaPhiM(jj->pT,jj->eta,jj->phi,jj->mass);
  return jetv;
}

//cmp_CSV --> Jets sorting by CSV
bool cmp_CSV(Jet jet1, Jet jet2){
  return jet1.CSV > jet2.CSV;
}

//cmp_CMVA --> Jets sorting by CSV
//bool cmp_CMVA(Jet jet1, Jet jet2){
//  return jet1.CMVA > jet2.CMVA;
//}

//computeCosThetaStar --> star angle computation    
float computeCosThetaStar(TLorentzVector dijet, TLorentzVector dihiggs){     
  dijet.Boost(-dihiggs.BoostVector()); //equal to (-P1.px()/P1.e(),-P1.py()/P1.e(),-P1.pz()/P1.e())
  float costhetast = dijet.CosTheta();
  return costhetast;
}

//withinElipse --> check if two candidates masses are in elipses region
int withinElipse(float mH1, float mH2, float a=30., float b=60., float mH1_c=H_mass, float mH2_c=H_mass)
{
 //45deg antirotation
  float fact = (sqrt(2)/2); //1; //
  float mH1_ = (mH1-mH1_c)*fact - (mH2-mH2_c)*fact;
  float mH2_ = (mH1-mH1_c)*fact + (mH2-mH2_c)*fact;
  float a_ = a*fact;
  float b_ = b*fact ;
  
  float elips = (mH1_/a_)*(mH1_/a_) + (mH2_/b_)*(mH2_/b_);

  int ret=-1;
  if(elips == 1) ret = 1;
  else if(elips > 1) ret = 0;
  return ret;
}

//withinElipse --> check if two candidates masses are in rectangular region
int withinRegion(float mH1, float mH2, float r1=15., float r2=30., float mH1_c=H_mass, float mH2_c=H_mass)
{
  float r=pow(pow(mH1-mH1_c, 2)+pow(mH2-mH2_c, 2), 0.5);
  float angle=atan2(mH2-mH2_c, mH1-mH1_c);
  int ret=-1;
  if (r<r1) ret=0;
  else if (r>r1 && r<r2){
    if (angle>=0 && angle<pi/2.) ret=1;
    else if (angle>=pi/2. && angle<pi) ret=4;
    else if (angle<0 && angle>=-pi/2.) ret=2;
    else if (angle<pi/2.&& angle>=-pi) ret=3;
    else std::cout<<"This is within annulus but not within any CR!"<<std::endl;
  }
  else ret=5;
  return ret;
}

bool checkIfNan(float var){
   if( var != var) return true;
   else return false;
}


