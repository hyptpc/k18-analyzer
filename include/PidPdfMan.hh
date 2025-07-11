// -*- C++ -*-
#ifndef PIDPDF_HH
#define PIDPDF_HH
#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
#include "THn.h"
#include <TString.h>
#include <TFile.h>
#include <array>
#include <memory>
#include <map>
#include <unordered_map>
#include <vector>
#include "TFitResultPtr.h" 
#include "Kinematics.hh"
#include "DatabasePDG.hh"
#include "TGraphErrors.h"
#include "PidData.hh"
#include "PidCommon.hh"

struct FitResultPars {
  double p[6]{}, perr[6]{};
  double chi2ndf{};
  int    status{};
  bool   valid{false};
};

using FitResultTable =
  std::array<                       // type
  std::array<                       // pid
    std::array<                     // charge
      std::array<                   // beid
	std::array<FitResultPars,pidlikeli::kNmom>,  // momBin
	pidlikeli::kNbe>,	  
      pidlikeli::kNchg>,
    pidlikeli::kNpidAll>,
  pidlikeli::kNtype>;  

class PidPdfMan
{
public:
  static const TString& ClassName();
  static PidPdfMan& GetInstance();
  static bool HasInstance();
  static void Cleanup();  
  ~PidPdfMan();
  
  static constexpr int kNtype = pidlikeli::kNtype;
  static constexpr int kNpid = pidlikeli::kNpid;
  static constexpr int kNpidAll = pidlikeli::kNpidAll;  
  static constexpr int kNchg = pidlikeli::kNchg;
  static constexpr int kNmom = pidlikeli::kNmom;
  static constexpr int kPlus  = pidlikeli::kPlus; 
  static constexpr int kMinus = pidlikeli::kMinus; 
  static constexpr double kDP = pidlikeli::kDP;
  static constexpr double kPmin = pidlikeli::kPmin;
  static constexpr int kNbe = pidlikeli::kNbe;    
  static constexpr double kDBE = pidlikeli::kDBE;
  static constexpr double kBEmin = pidlikeli::kBEmin;
  static constexpr int kPion     = pidlikeli::kPion    ;
  static constexpr int kKaon     = pidlikeli::kKaon    ;
  static constexpr int kProton   = pidlikeli::kProton  ;
  static constexpr int kDeutron  = pidlikeli::kDeutron ;
  static constexpr int kElectron = pidlikeli::kElectron;
  static constexpr int kAllParticles = pidlikeli::kAllParticles;  
  static constexpr int fac_t = pidlikeli::fac_t; // type
  static constexpr int fac_p = pidlikeli::fac_p; // pid
  static constexpr int fac_c = pidlikeli::fac_c; // charge
  static constexpr int fac_b = pidlikeli::fac_b; // be
  static constexpr int fac_m = pidlikeli::fac_m; // mom
  static constexpr double maxpoq = pidlikeli::maxpoq;
  
public:
  // use 1/3,1/3,1/3 prob. when flatPrior=true
  bool Initialize();
  Bool_t Initialize(const TString& file_name);    
  Bool_t IsReady() const { return m_is_ready; }
  TH2D* GetHist(int type,int pid,int charge,int be,int mom) const  ;
  TH2D* GetHistById(long long id) const;
  Bool_t GoodPid(int type, int pid, int chg, int be, int mom) const;
  bool CheckFitFunc() const;
  bool CheckFitFunc(int type, int pid, int chg, int be, int mom) const;  
  void  SetFileName(const TString& file_name) { m_file_name = file_name; }
  // int  MomToBin(double pGeV) const; // convert p [GeV/c] → mom index
  // double BinToMom(int momBin) const;
  // int  ChgToBin(int charge) const;
  // int BEToBin(double beGeV) const;
  // int BEToBinTEMP(double beGeV) const;
  bool GetParam(int pid, int chgmode=-1 /* -1:all, 0:+, 1:- */ );
  bool GetParamKaonUsingOtherParticles();
  bool CreateFitFunction(int it_int, int pid_int, int chg_int, int ibe_int, int m);
  bool CreateFitFunctionFive2DGauss(int chg_int, int ibe_int, int m);  
  bool ExecFitProj(int it_int, int pid_int, int chg_int, int ibe_int, int m);
  bool ExecFit2D(int it_int, int pid_int, int chg_int, int ibe_int, int m);
  bool ExecFitDouble2DGauss(int it_int, int pid_int, int chg_int, int ibe_int, int m);
  bool ExecFitDouble2DGaussAllParticles(int it_int, int pid_int, int chg_int, int ibe_int, int m);
  bool Run();
  bool SetInitFitPars(int pid, int chg, int momid);
  bool SetInitFitParsProj(int type,int pid,int chg,int be,int momid);
  bool SetInitFitPars2D(int type,int pid, int chg,int be,int momid,TH2D*h,TF2*f);  
  bool SetInitFitParsDoubleGauss(int type,int pid,int chg,int be,int momid);
  TFitResultPtr ExecFitFive2DGauss(int chg,int be,int momid);  
  //bool StoreFitResult(int t,int p,int c,int b,int momid,TF2* f);
  bool StoreFitResultProj(int t,int p,int c,int b,int momid);
  bool StoreFitResultKaon(int t,int p,int c,int b,int momid);  
  bool CheckGoodFitM2(int t,int p,int c,int b,int momid);
  bool CheckGoodFitdEdx(int t,int p,int c,int b,int momid);
  bool CheckGoodFitM2dEdx(int t,int p,int c,int b,int momid);
  bool MakeGoodPdf(int pid);
  bool MakeGoodPdfKm();
  void MakeGoodPdfPiKP(int runmode=pidlikeli::kDGeneral);  
  void StoreGoodPdfAll();  
  void FillPriorHist();
  void FillPriorHist(int ip);  
  bool FillPriorHist(int it_int, int ip_int, int ic_int, int ib_int, int im);
  double GetYield(int it, int ip, int ic, int ib, int im);
  double GetTotalYield(int it, int ic, int ib, int im);
  std::vector<double> CalculatePriorProb(const TFitResultPtr& fit_result) const;
  void StorePriorProbs();
  void GetPriorProb();
  double GetPriorProb(int it, int ip, int ic, int ib, int im);
  
  double getFitMuM2(int t,int p,int c,int b,int momid);
  double getFitMudEdx(int t,int p,int c,int b,int momid);
  double getFitSigM2(int t,int p,int c,int b,int momid);
  double getFitSigdEdx(int t,int p,int c,int b,int momid);  
  double getFitAngle(int t,int p,int c,int b,int momid);
  double getFitYield(int t,int p,int c,int b,int momid);
  double getFitErrMuM2(int t,int p,int c,int b,int momid);
  double getFitErrMudEdx(int t,int p,int c,int b,int momid);
  double getFitErrSigM2(int t,int p,int c,int b,int momid);
  double getFitErrSigdEdx(int t,int p,int c,int b,int momid);  
  double getFitErrAngle(int t,int p,int c,int b,int momid);
  double getFitErrYield(int t,int p,int c,int b,int momid);
  double getFitValue(int t,int p,int c,int b,int momid, int paramid);
  double getFitErrValue(int t,int p,int c,int b,int momid, int paramid);
  
  bool SetGEPoint(int t,int p,int c,int b,int momid);
  bool FitGEPoint(int t,int p,int c,int b);    
  std::array<double, pidfunc::kNparamGauss> GetCalcParameters(int t, int pid, int c, int b, int m) const;
  std::array<double, pidfunc::kNparamGauss> GetCalcParametersRef(int pid, int momid) const;
  const FitResultPars& GetStoredFitParamsOfRefPdf(const ParticleFitConfig& config, int momid) const;
  double InterpYield(int t,int p,int c,int b,int momid) const;
   
  Long64_t WriteToRootfile(TFile* fout);
  Bool_t WritePDFToRootfile();
  void DrawResultsToPdf(const std::string& output_filename);
  void DrawResultsToPdf();  
  
  // ---- Getter ----
  TF2*  GetFit2D (int itype, int pid, int chg, int ibe, int mom) const;
  
  // ---- fitting result ----
  // struct FitResultPars{
  //   double p[pidfunc::kNparamGauss]{}, perr[pidfunc::kNparamGauss]{};
  //   double chi2ndf{};
  //   int    status{};
  //   bool   valid{false};
  // };
  bool StoreResult(int type,int pid,int chg,int beid,int momid) const;
  static FitResultPars& at (int type,int pid,int chg,int beid,int momid)
  {
    check(type,pid,chg,beid,momid);
    return tbl_[type][pid][chg][beid][momid];
  }  
  static FitResultTable tbl_;
  
private:  
  PidPdfMan();
  PidPdfMan(const PidPdfMan&);
  PidPdfMan& operator = (const PidPdfMan&);
  void Finalize();  
  static PidPdfMan* s_instance;
  bool m_is_ready = false;
  std::unique_ptr<PidData> m_data;  
  std::array<double, pidfunc::kNparamGauss> GetCalcParametersRef(const ParticleFitConfig& config, int momid) const;

private:
  bool LoadAllHists();
  //  bool FitAllHistsProjection();
  //bool FitAllHists2DRotGauss();
  Bool_t CdMainFile();
  Bool_t ResetParamFile();
  bool OpenParamFile();
  bool CloseCurrentFile();

private:
  // members --------------------------------------------------------------  
  TString    m_file_name;
  std::unique_ptr<TFile> m_file     {};
  std::unique_ptr<TFile> m_mainfile {};
  std::unique_ptr<TFile> m_paramfile{};
  std::unordered_map<long long, std::unique_ptr<TH2D>> m_histMap; // key = full id
  bool m_flat=false;
  
  static void check(int t,int p,int c,int b,int momid)
  {
    if(t<0||t>=kNtype||p<0||p>=kNpidAll||c<0||c>=kNchg||b<0||b>=kNbe||momid<0||momid>=kNmom)
      throw std::out_of_range("PidPdfMan index out of range");
  }  
};

inline const TString&
PidPdfMan::ClassName()
{
  static TString s_name("PidPdfMan");
  return s_name;
}

#endif
