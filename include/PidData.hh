// -*- C++ -*-
#ifndef PIDDATA_HH
#define PIDDATA_HH

#include <array>
#include <memory>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "THn.h"
#include "PidCommon.hh"

struct ParticleFitConfig {
  pidlikeli::Pid pid;
  const char* name;
  pidlikeli::DType source_type;
  pidlikeli::Chg source_chg;
  pidlikeli::BE source_be;
  pidlikeli::Pid sigma_source_pid;
};

struct PdfBin {
  TF2* f2d = nullptr; // 2D gauss
  TF2* f2dc = nullptr;
  TF2* fd2d = nullptr; // double2D gauss for K-pi- (p+pi+ in Kp reaction)
  TF2* f2dkm = nullptr; // K- obtained by fitting with pi- and by selecting BE 
  //  TF2* f2dkmbe = nullptr; // K- obtained
  TF2* ff2d = nullptr; // five 2d gauss
  TF2* pdf = nullptr;
  
  TF1* fx = nullptr;
  TF1* fy = nullptr;

  bool functions_created = false;
  bool five2dgauss_created = false; 
  
  double logNorm = 0.0;
  bool   valid   = false;
  bool   is_goodM2   = false;
  bool   is_gooddEdx   = false;
  bool   is_goodM2dEdx   = false;
};

struct CorrGraph {
  enum class Graph { MeanM2, MeandEdx, SigM2, SigdEdx, RotAngle, Yield,
    MeanM2Mom, MeandEdxMom, SigM2Mom, SigdEdxMom, RotAngleMom, YieldMom, COUNT };
  std::array<TGraphErrors*, static_cast<size_t>(Graph::COUNT)> graphs;
  CorrGraph();
  static const std::array<TString, static_cast<size_t>(Graph::COUNT)> GraphNames;
};

struct CorrFunc {
  enum class Func { MeanM2, MeandEdx, SigM2, SigdEdx, RotAngle, Yield,
    MeanM2Mom, MeandEdxMom, SigM2Mom, SigdEdxMom, RotAngleMom, YieldMom,COUNT };
  std::array<TF1*, static_cast<size_t>(Func::COUNT)> functions;
  CorrFunc();
  static const std::array<TString, static_cast<size_t>(Func::COUNT)> FuncNames;  
};

struct PidData {
  using DTypeArr = std::array<std::array<std::array<std::array<PdfBin, pidlikeli::kNmom>, pidlikeli::kNbe>, pidlikeli::kNchg>, pidlikeli::kNpidAll>;
  std::array<DTypeArr, pidlikeli::kNtype> m_bin;
  
  using CGTypeArr = std::array<std::array<std::array<CorrGraph, pidlikeli::kNbe>, pidlikeli::kNchg>, pidlikeli::kNpidAll>;
  std::array<CGTypeArr, pidlikeli::kNtype> m_ge_init;
  std::array<CGTypeArr, pidlikeli::kNtype> m_ge;
    
  using CFTypeArr = std::array<std::array<std::array<CorrFunc, pidlikeli::kNbe>, pidlikeli::kNchg>, pidlikeli::kNpidAll>;
  std::array<CFTypeArr, pidlikeli::kNtype> m_f_init;
  std::array<CFTypeArr, pidlikeli::kNtype> m_f_good;
  
  std::unique_ptr<THnD> m_h5_priors;
  
  PidData();
  inline PdfBin& at_bin(pidlikeli::DType t, pidlikeli::Pid p, pidlikeli::Chg c, pidlikeli::BE b, int m) {
    return m_bin[pidlikeli::scast(t)][pidlikeli::scast(p)][pidlikeli::scast(c)][pidlikeli::scast(b)][m];
  }
  const inline PdfBin& at_bin(pidlikeli::DType t, pidlikeli::Pid p, pidlikeli::Chg c, pidlikeli::BE b, int m) const {
    return m_bin[pidlikeli::scast(t)][pidlikeli::scast(p)][pidlikeli::scast(c)][pidlikeli::scast(b)][m];
  }
  // inline double& at_prior(pidlikeli::DType t, pidlikeli::Pid p, pidlikeli::Chg c, pidlikeli::BE b, int m) {
  //   return m_prior[pidlikeli::scast(t)][pidlikeli::scast(p)][pidlikeli::scast(c)][pidlikeli::scast(b)][m];
  // }
  
  inline CorrGraph& at_ge(pidlikeli::DType t, pidlikeli::Pid p, pidlikeli::Chg c, pidlikeli::BE b) {
    return m_ge[pidlikeli::scast(t)][pidlikeli::scast(p)][pidlikeli::scast(c)][pidlikeli::scast(b)];
  }
  const inline CorrGraph& at_ge(pidlikeli::DType t, pidlikeli::Pid p, pidlikeli::Chg c, pidlikeli::BE b) const {
    return m_ge[pidlikeli::scast(t)][pidlikeli::scast(p)][pidlikeli::scast(c)][pidlikeli::scast(b)];
  }

  inline CorrGraph& at_ge_init(pidlikeli::DType t, pidlikeli::Pid p, pidlikeli::Chg c, pidlikeli::BE b) {
    return m_ge_init[pidlikeli::scast(t)][pidlikeli::scast(p)][pidlikeli::scast(c)][pidlikeli::scast(b)];
  }

  inline CorrFunc& at_f_good(pidlikeli::DType t, pidlikeli::Pid p, pidlikeli::Chg c, pidlikeli::BE b) {
    return m_f_good[pidlikeli::scast(t)][pidlikeli::scast(p)][pidlikeli::scast(c)][pidlikeli::scast(b)];
  }
};


#endif
