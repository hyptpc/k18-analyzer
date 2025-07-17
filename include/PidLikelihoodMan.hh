// -*- C++ -*-
// // include/PidLikelihoodMan.hh

#ifndef PIDLIKELIHOODMAN_HH
#define PIDLIKELIHOODMAN_HH

#include <vector>
#include <string>
#include <memory>
#include "THn.h"
#include "TF2.h"
#include "PidCommon.hh"
#include "PidData.hh" // for ParticleFitConfig

struct PidResult {
  pidlikeli::Pid pid; // PID
  double postprob; // Posterior
};

struct LHPdfBin {
  //TF2* pdf = nullptr;
  std::unique_ptr<TF2> pdf;
};

// ★★★ 2次元配列 [pid][mom] に変更 ★★★
using LHPdfTable = 
  std::array< // pid
    std::array<LHPdfBin, pidlikeli::kNmom>,
  pidlikeli::kNpid>;


class PidLikelihoodMan
{
public:
  static PidLikelihoodMan& GetInstance();
  ~PidLikelihoodMan();
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  void SetFileName(const TString& file_name) { m_file_name = file_name; }  
  Bool_t IsReady() const { return m_is_ready; }
  static const TString& ClassName();  
  static bool HasInstance();
  std::vector<double> CalculatePosterior(int charge, double momentum, double m2, double dedx) const;
  std::vector<PidResult> GetRankedPid(int charge, double momentum, double m2, double dedx) const;
  void DrawResultsToPdf(const std::string& output_filename);
  void DrawResultsToPdf();    

private:
  PidLikelihoodMan(); // Singleton
  PidLikelihoodMan(const PidLikelihoodMan&) = delete;
  PidLikelihoodMan& operator=(const PidLikelihoodMan&) = delete;
  void Finalize();
  static void Cleanup();
  
  //member
  static PidLikelihoodMan* s_instance;
  bool m_is_ready = false;
  TString  m_file_name;      
  std::unique_ptr<THnD> m_prior_hist;
  //std::array<std::unique_ptr<TF2>, pidlikeli::kNpid> m_ref_pdfs;
  LHPdfTable m_pdf_table;
  
};


inline const TString&
PidLikelihoodMan::ClassName()
{
  static TString s_name("PidLikelihoodMan");
  return s_name;
}
#endif


////old
// #ifndef PIDLIKELIHOOD_HH
// #define PIDLIKELIHOOD_HH
// #include <TF2.h>
// #include <TH1D.h>
// #include <TH2D.h>
// #include <TString.h>
// #include <TFile.h>
// #include <array>
// #include <memory>
// #include <map>
// #include <unordered_map>
// #include "Kinematics.hh"
// #include "DatabasePDG.hh"
// #include "TGraphErrors.h"
// #include "PidData.hh"
// #include "PidCommon.hh"

// // namespace {
// //   static const auto PionMass    = pdg::PionMass();
// //   static const auto KaonMass    = pdg::KaonMass();
// //   static const auto ProtonMass  = pdg::ProtonMass();
// //   static const auto ElectronMass = pdg::ElectronMass();
// //   static const auto DeutronMass = pdg::DeutronMass();
// //   static const auto np = pidlikeli::kNpid;  
// //   static const Double_t mass[np] = {PionMass,KaonMass,ProtonMass,DeutronMass,ElectronMass};      
// // }

// using ArrPiKP = std::array<std::array<double, pidlikeli::kNchg>, pidlikeli::kNpidAll>;

// struct PdfParam {
//   std::unique_ptr<TF2> pdf;
//   double logNorm = -500.;
//   double prior = 1.0 / pidlikeli::kNpid;
//   bool valid = false;
// };

// class PidLikelihoodMan
// {
// public:
//   static const TString& ClassName();
//   static PidLikelihoodMan& GetInstance();
//   static bool HasInstance();
//   static void Cleanup();  
//   ~PidLikelihoodMan();
  
//   static constexpr int kNtype = pidlikeli::kNtype;
//   static constexpr int kNpid = pidlikeli::kNpid;
//   static constexpr int kNpidAll = pidlikeli::kNpidAll;  
//   static constexpr int kNchg = pidlikeli::kNchg;
//   static constexpr int kNmom = pidlikeli::kNmom;
//   static constexpr int kPlus  = pidlikeli::kPlus; 
//   static constexpr int kMinus = pidlikeli::kMinus; 
//   static constexpr double kDP = pidlikeli::kDP;
//   static constexpr double kPmin = pidlikeli::kPmin;
//   static constexpr int kNbe = pidlikeli::kNbe;    
//   static constexpr double kDBE = pidlikeli::kDBE;
//   static constexpr double kBEmin = pidlikeli::kBEmin;
//   static constexpr int kPion     = pidlikeli::kPion    ;
//   static constexpr int kKaon     = pidlikeli::kKaon    ;
//   static constexpr int kProton   = pidlikeli::kProton  ;
//   static constexpr int kDeutron  = pidlikeli::kDeutron ;
//   static constexpr int kElectron = pidlikeli::kElectron;
//   static constexpr int kAllParticles = pidlikeli::kAllParticles;  
//   static constexpr int fac_t = pidlikeli::fac_t;  // type
//   static constexpr int fac_p = pidlikeli::fac_p;  // pid
//   static constexpr int fac_c = pidlikeli::fac_c;  // charge
//   static constexpr int fac_b = pidlikeli::fac_b;  // be
//   static constexpr int fac_m = pidlikeli::fac_m;  // mom
//   static constexpr double maxpoq = pidlikeli::maxpoq;
  
// public:
//   // use 1/3,1/3,1/3 prob. when flatPrior=true
//   bool Initialize();
//   Bool_t Initialize(const TString& file_name);    
//   Bool_t IsReady() const { return m_is_ready; }
//   Bool_t GoodPid(int type, int pid, int chg, int be, int mom) const;
//   void  SetFileName(const TString& file_name) { m_file_name = file_name; }
//   // int  MomToBin(double pGeV) const; // convert p [GeV/c] → mom index
//   // double BinToMom(int momBin) const;
//   // int  ChgToBin(int charge) const;
//   // int BEToBin(double beGeV) const;
//   // int BEToBinTEMP(double beGeV) const;
//   //TF2* GetFunc(int it_int, int pid_int, int chg_int, int ibe_int, int m) const;  
//   TF2* Get2DGauss(int it_int, int pid_int, int chg_int, int ibe_int, int m) const;
//   TF2* Get2DGaussById(long long id) const;  
//   bool GetParam(int pid, int chgmode=-1 /* -1:all, 0:+, 1:- */ );
//   bool CreateFitFunction(int it_int, int pid_int, int chg_int, int ibe_int, int m);
//   bool ExecFitProj(int it_int, int pid_int, int chg_int, int ibe_int, int m);
//   bool GetParamKaonMinus();  
//   bool Run();
//   double GetYield(int it, int ip, int ic, int ib, int im);
//   void GetPriorProb();
//   double GetPriorProb(int it, int ic, int ib, int im);
  
//   double getFitValue(int t,int p,int c,int b,int momid, int paramid);
//   double getFitErrValue(int t,int p,int c,int b,int momid, int paramid);
  
//   std::array<double, 6> GetCalcParameters(int t, int pid, int c, int b, int m) const;
//   double CalcSigM2(int p, double mom) const;
//   double CalcSigdEdx(int p, double mom) const;
//   double InterpYield(int t,int p,int c,int b,int momid) const;
   
//   // log likelihood and posterior
//   bool isGoodBin(int itype,int pid,int chg,double be,double p) const;
//   Bool_t isGoodLogL(double LogL) const;
//   double LogP2D(int itype,int pid,int chg, double be, double p, double invB,double dE) const;
//   double LogNorm(int itype,int pid,int chg, double be, double p, double invB,double dE) const;
//   double LogPrior(int itype,int pid,int chg,double be,double p,double invB,double dE) const;
//   double LogL(int itype,int pid,int chg,double be,double p,double invB,double dE) const;
//   // double LogL(int type,int pid,int charge,double be,double p,
//   //             double invB,double dEdx) const;
//   double Posterior(int type,int pid,int charge,double be,double p,double invB,double dEdx) const;
//   ArrPiKP PriorPiKP(int type,int chg,double be,double p) const;  
//   ArrPiKP PosteriorPiKP(int type,int charge,double be,double p,double invB, double dEdx) const;
//   ArrPiKP LogLPiKP(int type,int chg,double be,double p,double invB,double dEdx) const;
//   Int_t NthBestProbPid(int type,int charge,double be,double p,double invB,double dEdx,int nth) const;
//   Double_t NthBestProb(int type,int charge,double be,double p,double invB,double dEdx,int nth) const;    
//   Long64_t WriteToRootfile(TFile* fout);
//   Bool_t WritePDFToRootfile();
  
//   // likelihood result
//   struct LikeliResult {
//     double logPrior = 0.;    
//     double logNorm = 0.;
//     double logP2D = 0.;    
//     double logL = 0.;
//     double posterior = 0.;
//     bool good_pri = false;
//     bool good_p2d = false;
//     bool good_lnL = false;
//     bool good_norm = false;
//     bool good_pid() const { return good_pri && good_lnL && good_norm; }
//   };  
//   LikeliResult CalcLikeliResult(int type,int pid,int chg,double be,double mom,double invb,double dEdx) const;
  
//   struct LikeliResultPid {
//     std::vector<PidLikelihoodMan::LikeliResult> results;
//     std::vector<int> nthPid;
//     std::vector<double> nthPidProb;
//     std::vector<int> nthPidGood;
//     bool empty() const { return results.empty(); }
//   };
//   LikeliResultPid CalcLikeliResultPid(int type,int chg,double be,double mom,double invb,double dEdx) const;
    
// private:  
//   PidLikelihoodMan();
//   PidLikelihoodMan(const PidLikelihoodMan&);
//   PidLikelihoodMan& operator = (const PidLikelihoodMan&);
//   void Finalize();  
//   static PidLikelihoodMan* s_instance;
//   bool m_is_ready = false;
//   std::unique_ptr<PidData> m_data;

//   using PdfTable = std::array< // DType
//     std::array< // Pid
//       std::array< // Chg
// 	std::array< // BE
// 	  std::array<PdfParam, pidlikeli::kNmom>, 
// 	  pidlikeli::kNbe>,
// 	pidlikeli::kNchg>,
//       pidlikeli::kNpidAll>,
//     pidlikeli::kNtype>;
    
//   PdfTable m_pdf_table;

// private:
//   bool LoadAllPdf();
//   Bool_t CdMainFile();
//   Bool_t ResetParamFile();
//   bool OpenParamFile();
//   bool CloseCurrentFile();

// private:
//   // members --------------------------------------------------------------  
//   TString    m_file_name;
//   std::unique_ptr<TFile> m_file     {};
//   std::unique_ptr<TFile> m_mainfile {};
//   std::unique_ptr<TFile> m_paramfile{};
//   std::unordered_map<long long, std::unique_ptr<TH2D>> m_histMap; // key = full id
//   bool m_flat=false;
  
//   static void check(int t,int p,int c,int b,int momid)
//   {
//     if(t<0||t>=kNtype||p<0||p>=kNpidAll||c<0||c>=kNchg||b<0||b>=kNbe||momid<0||momid>=kNmom)
//       throw std::out_of_range("PidLikelihoodMan index out of range");
//   }  
// };

// inline const TString&
// PidLikelihoodMan::ClassName()
// {
//   static TString s_name("PidLikelihoodMan");
//   return s_name;
// }

// #endif
