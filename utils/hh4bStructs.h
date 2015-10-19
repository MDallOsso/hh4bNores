
//Struct and Functions for hh4b analysis macros
//---------------------------------------------

double pi=3.14159265;
//for withinRegion
double H_mass = 115.0;
double dijetM_cut_low = 100.;
double dijetM_cut_high = 150.;

//Jet
typedef struct{
  //double CMVA;
  double CSV;
  double mass;
  double pT;
  double eta;
  double phi;
} Jet;

//Jet for plots
typedef struct{
  TH1F* h;
  double norm;
} Jet4Plot;

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
double computeCosThetaStar(TLorentzVector dijet, TLorentzVector dihiggs){     
  dijet.Boost(-dihiggs.BoostVector()); //equal to (-P1.px()/P1.e(),-P1.py()/P1.e(),-P1.pz()/P1.e())
  double costhetast = dijet.CosTheta();
  return costhetast;
}

//withinElipse --> check if two candidates masses are in elipses region
int withinElipse(double mH1, double mH2, double a=30., double b=60., double mH1_c=H_mass, double mH2_c=H_mass)
{
 //45deg antirotation
  double fact = (sqrt(2)/2); //1; //
  double mH1_ = (mH1-mH1_c)*fact - (mH2-mH2_c)*fact;
  double mH2_ = (mH1-mH1_c)*fact + (mH2-mH2_c)*fact;
  double a_ = a*fact;
  double b_ = b*fact ;
  
  double elips = (mH1_/a_)*(mH1_/a_) + (mH2_/b_)*(mH2_/b_);

  int ret=-1;
  if(elips == 1) ret = 1;
  else if(elips > 1) ret = 0;
  return ret;
}

//withinElipse --> check if two candidates masses are in rectangular region
int withinRegion(double mH1, double mH2, double r1=15., double r2=30., double mH1_c=H_mass, double mH2_c=H_mass)
{
  double r=pow(pow(mH1-mH1_c, 2)+pow(mH2-mH2_c, 2), 0.5);
  double angle=atan2(mH2-mH2_c, mH1-mH1_c);
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

