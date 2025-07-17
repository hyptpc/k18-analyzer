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

enum class PdfFuncType {
    F2D,      // original fitting func
    F2DCORR,  // corrected(good) fitting func
    FX_PROJ,  // projectX
    FY_PROJ,  // projectY
    FD2D,     // double gaussian
    FFIVE2D,  // penta gaussian
    F1D_M2,   // 1D Fit function for m2
    F1D_DEDX, // 1D Fit function for dEdx
    PDF,      // final pdf    
    COUNT     // counts of enum elements
};

const char* const PdfFuncName[static_cast<size_t>(PdfFuncType::COUNT)] = {
  "f2d", "f2dc", "fx", "fy", "fd2d", "ff2d", "f1dm2", "f1ddedx", "pdf"
};

struct PdfBin {  
  std::array<TObject*, static_cast<size_t>(PdfFuncType::COUNT)> funcs{};
  template <typename T>
    T* get_func(PdfFuncType type) const {
        return static_cast<T*>(funcs[static_cast<size_t>(type)]);
    }
  template <typename T>
    void set_func(PdfFuncType type, T* func_ptr) {
        funcs[static_cast<size_t>(type)] = func_ptr;
    }
  
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
  std::unique_ptr<THnD> m_h5_good_yields;
  
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
