
#define hh4b_tmva_H

class hh4b_tmva {

  std::string datasample, MCsample, opt;
  int appOn, whatToDo;

  private:
    TH1F *bdt_out = new TH1F("bdt_application","bdt_application",100,0,1);
    TH1F *like_out = new TH1F("like_application","like_application",100,0,1);
    TH1F *merit_L = new TH1F("Q_L", "Likelihood Q", 40, 0, 1);
    TH1F *merit_BDT = new TH1F("Q_BDT", "BDT Q", 40, 0, 1);
    std::string outNtuples;

  public:
    hh4b_tmva(std::string , std::string , std::string, int ,  int );
    ~hh4b_tmva();

    void doTmva();
    void training(std::string , std::string , std::string , std::string );
    bool BDTapplication(std::string , std::string );
    bool Lapplication(std::string , std::string );
    bool CutScanL(std::string );
    bool CutScanBDT(std::string );
    std::string setTitle(int );

};
