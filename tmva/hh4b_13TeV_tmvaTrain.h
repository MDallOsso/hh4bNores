

class hh4b_tmvaTrain {
  public:
    hh4b_tmvaTrain();
    ~hh4b_tmvaTrain();

    string weightfile = "all";
    string inputSignal = dataFld+"tree_Step1_M-260.root";
    string inputBkg = dataFld+"tree_Step1_BTagCSVRun2015C.root";
    string outTMVAfileName = "TMVA_allTraining.root";
    string appfile;

    void training(string , string , string , string );
    TH1D* BDTapplication(string , string );
    TH1D* Lapplication(string , string );
    TH1F* CutScanL(string );
    TH1F* CutScanBDT(string );
    std::string setTitle(int );

};
