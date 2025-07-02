// -*- C++ -*-
#ifndef PIDCOMMON_HH
#define PIDCOMMON_HH
#include <TF2.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TFile.h>
#include <array>
#include <memory>
#include <map>
#include <unordered_map>
#include "Kinematics.hh"
#include "DatabasePDG.hh"
#include "TGraphErrors.h"
//#include "PidData.hh"

namespace pidlikeli {
  enum class DType { General, Lmd, K0, resoKqfK, COUNT }; // { no cut, Lambda->ppi reconstructed, K0->pipi reconstructed, BE region cut for Kstar/quasifreeK-
  enum class Pid { Pi, K, P, D, E, COUNT }; // pi,kaon,proton,deutron,electron
  enum class Chg { Plus, Minus, COUNT }; // charge
  enum class BE { Index0, resoK, qfK, COUNT }; // binding energy region, not prepared yet

  inline constexpr size_t kNtype = static_cast<size_t>(DType::COUNT);
  inline const TString type[kNtype] = {"general","Lmd","K0","rsKqfK"};
  inline constexpr size_t kNpid = static_cast<size_t>(Pid::COUNT);
  inline const TString plist[kNpid] = {"pi","k","p","d","e"};
  inline constexpr size_t kNchg = static_cast<size_t>(Chg::COUNT);
  inline const TString clist[kNchg] = {"+","-"};
  inline constexpr size_t kNbe = static_cast<size_t>(BE::COUNT);
  
  inline constexpr int kPion = 0;
  inline constexpr int kKaon = 1;
  inline constexpr int kProton = 2;
  inline constexpr int kDeutron = 3;
  inline constexpr int kElectron = 4;
  inline constexpr int kAllParticles = 5; 
  inline constexpr int kPlus = 0; 
  inline constexpr int kMinus = 1;
  inline constexpr double kBEmin = -1.0; // GeV
  inline constexpr double kBEmax = 1.0; // GeV
  inline constexpr double kDBE  = 0.2; // GeV, energy step
  inline constexpr double kPmin = 0.0; // GeV/c
  inline constexpr double kPmax = 1.5; // GeV/c
  inline constexpr double kDP  = 0.01; // GeV/c, mom step
  inline constexpr int kNmom = (kPmax-kPmin)/kDP; // should be less than 1000
  // for HistId
  constexpr int fac_t = 10'000'000; // datatype
  constexpr int fac_p = 1'000'000;  // pid
  constexpr int fac_c = 100'000;    // charge
  constexpr int fac_b = 1'000;      // be
  constexpr int fac_m = 1;          // mom
  // fitparam
  inline constexpr int kNfitpar = 6; // for 2D gauss, 1D projected gauss
  // bin
  inline constexpr int nbinpoq = 600;
  inline constexpr double minpoq = -0.9;
  inline constexpr double maxpoq = 0.9; //GeV/c  
  //inline constexpr int nbindedx = 175;
  inline constexpr int nbindedx = 350; 
  inline constexpr double mindedx = 0.;
  inline constexpr double maxdedx = 350.;
  inline constexpr int nbinm2 = 200;
  inline constexpr double minm2 = -1.5;
  inline constexpr double maxm2 = 1.5;
  inline constexpr double binwpoq = (maxpoq-minpoq)/Double_t(nbinpoq);
  inline constexpr double binwdedx = (maxdedx-mindedx)/Double_t(nbindedx);
  inline constexpr double binwm2 = (maxm2-minm2)/Double_t(nbinm2);
  inline constexpr double cutm2min[kNpid] = {-0.1, 0.15, 0.7, 1.0, -0.1};
  inline constexpr double cutm2max[kNpid] = { 0.1 , 0.5,  1.1, 2.0,  0.1};
  inline constexpr double cutdedxmin[kNpid] = {-0.01, 0.15, 0.7, 1.0, -0.01};
  inline constexpr double cutdedxmax[kNpid] = { 0.1 , 0.5,  1.1, 2.0,  0.01};
  inline constexpr double cutm2sigmamax[kNpid] = {0.1, 0.1, 0.4, 0.0,  0.00};
  inline constexpr double cutdedxsigmamax[kNpid] = {20, 0.1, 0.4, 0.0,  0.00};

  static const auto PionMass    = pdg::PionMass();
  static const auto KaonMass    = pdg::KaonMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto ElectronMass = pdg::ElectronMass();
  static const auto DeutronMass = pdg::DeutronMass();
  static const auto np = pidlikeli::kNpid;  
  static const Double_t mass[np] = {PionMass,KaonMass,ProtonMass,DeutronMass,ElectronMass};      
  static constexpr int idPi = kPion;
  static constexpr int idK  = kKaon;
  static constexpr int idP  = kProton;
  static constexpr int idD  = kDeutron;
  static constexpr int idE  = kElectron;    

  // general helper
  template<typename E>
  static constexpr size_t scast(E e) {
    return static_cast<size_t>(e);
  }

  inline Double_t MomToBeta(Double_t mom, Double_t m2){
    double mom2 = mom*mom;
    double ene2 = mom2 + m2;
    double beta = TMath::Sqrt(mom2/ene2);
    return beta;
  }
  inline Double_t MomToBetaPion(Double_t mom){
    double m = PionMass; //GeV/c2
    double m2 = m*m;
    double beta = MomToBeta(mom, m2);
    return beta;
  }
  inline Double_t MomToBetaProton(Double_t mom){
    double m = ProtonMass; // GeV
    double m2 = m*m;
    double beta = MomToBeta(mom, m2);
    return beta;
  }
  inline Double_t MomToBetaKaon(Double_t mom){
    double m = KaonMass; // GeV
    double m2 = m*m;
    double beta = MomToBeta(mom, m2);
    return beta;
  }
  inline Double_t MomToBetaPid(Double_t mom, int pid){
    if(pid==idPi) return MomToBetaPion(mom);
    else if(pid==idK) return MomToBetaKaon(mom);
    else if(pid==idP) return MomToBetaProton(mom);
    else return std::numeric_limits<double>::quiet_NaN();   
  }
  
  inline Double_t dedx_pid(Double_t mom, int pid){
    if(pid==idPi) return Kinematics::HypTPCdEdxPion(mom);
    else if(pid==idK) return Kinematics::HypTPCdEdxKaon(mom);
    else if(pid==idP) return Kinematics::HypTPCdEdxProton(mom);
    else return std::numeric_limits<double>::quiet_NaN();
  }
  inline int MomToBin(double pGeV)
  {
    int idx = int( std::round( (pGeV - kPmin) / kDP ) );
    return (idx < 0 || idx >= kNmom) ? -1 : idx;
  }
  inline double BinToMom(int momBin)
  {
    return kPmin + kDP * momBin;
  }
  inline int ChgToBin(int charge)
  {
    int idx = -2;
    if(charge>0) idx = 0;
    if(charge<0) idx = 1;
    return (idx < 0 || idx >= kNchg) ? -1 : idx;
  }
  inline int BEToBin(double beGeV)
  {
    int idx = int( std::round( (beGeV - kBEmin) / kDBE ) );
    return (idx < 0 || idx >= kNbe) ? -1 : idx;
  }
  inline int BEToBinTEMP(double beGeV)
  {
    int binBE=0;
    //  int idx = int( std::round( (beGeV - kBEmin) / kDBE ) );
    //  return (idx < 0 || idx >= kNbe) ? -1 : idx;
    return binBE;
  }
  
}

namespace pidfunc {
  Double_t RotGauss2D(Double_t* xy, Double_t* par);
  Double_t RotGauss2DFit(Double_t* xy, Double_t* par);  
  Double_t RotGauss2DProjX(Double_t* x, Double_t* par);
  Double_t RotGauss2DProjXFit(Double_t* x, Double_t* par);
  Double_t RotGauss2DProjY(Double_t* x, Double_t* par);
  Double_t RotGauss2DProjYFit(Double_t* x, Double_t* par);
  Double_t RotGauss2DBeta(Double_t* xy, Double_t* par);
  Double_t RotTwoGauss2D(Double_t* xy, Double_t* par);  
  Double_t SigmaDedx(Double_t* x, Double_t* par);  
  Double_t SigmaM2Pid(Double_t* x, int pid);

  Double_t SigmaInvBeta(Double_t* x, Double_t* par);
  Double_t ExpPlusPol3(Double_t *x, Double_t *p);  
  Double_t ExpPlusPol1(Double_t *x, Double_t *p);  
  Double_t Pol1(Double_t *x, Double_t *p);
  Double_t RationalFunc1(Double_t *x, Double_t *p);  
  Double_t LogiFunc1(Double_t *x, Double_t *p);  
}

#endif
