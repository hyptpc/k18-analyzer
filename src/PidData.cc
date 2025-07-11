// -*- C++ -*-
#include <array>
#include <memory>
#include "TGraphErrors.h"
#include "TF1.h"
#include "TF2.h"
#include "THn.h"
#include "PidCommon.hh"
#include "PidData.hh"

PidData::PidData() {
  const int n_dimensions = pidlikeli::kDimBin;
  int n_bins[n_dimensions] = {pidlikeli::kNtype, pidlikeli::kNpid, pidlikeli::kNchg, pidlikeli::kNbe, pidlikeli::kNmom};
  double x_min[n_dimensions] = {-0.5, -0.5, -0.5, -0.5, -0.5};
  double x_max[n_dimensions] = {pidlikeli::kNtype-0.5, pidlikeli::kNpid-0.5, pidlikeli::kNchg-0.5, pidlikeli::kNbe-0.5, pidlikeli::kNmom-0.5};
  
  m_h5_priors = std::make_unique<THnD>(
				      "h5_priors", "Prior Probabilities",
				      n_dimensions, n_bins, x_min, x_max
				      );
}

const std::array<TString, static_cast<size_t>(CorrGraph::Graph::COUNT)>
CorrGraph::GraphNames = {
  "gMeanM2", "gMeandEdx", "gSigM2", "gSigdEdx", "gRotAngle", "gYield",
  "gMeanM2Mom", "gMeandEdxMom", "gSigM2Mom", "gSigdEdxMom", "gRotAngleMom", "gYieldMom"
};

const std::array<TString, static_cast<size_t>(CorrFunc::Func::COUNT)>
CorrFunc::FuncNames = {
  "fMeanM2", "fMeandEdx", "fSigM2", "fSigdEdx", "fRotAngle", "fYield",
  "fMeanM2Mom", "fMeandEdxMom", "fSigM2Mom", "fSigdEdxMom", "fRotAngleMom", "fYieldMom"
};

// struct CorrFunc
CorrFunc::CorrFunc() {
  using Func = CorrFunc::Func;
  functions[pidlikeli::scast(Func::MeanM2)]
    = new TF1("fMeanM2", pidfunc::Pol1, 0., pidlikeli::maxpoq, 2);
  functions[pidlikeli::scast(Func::MeandEdx)]
    = new TF1("fMeandEdx", pidfunc::Pol1, 0, pidlikeli::maxpoq, 2);
  functions[pidlikeli::scast(Func::SigM2)]
    = new TF1("fSigM2", pidfunc::ExpPlusPol1, 0, pidlikeli::maxpoq, 4);
  functions[pidlikeli::scast(Func::SigdEdx)]
    = new TF1("fSigdEdx", pidfunc::SigmaInvBeta, 0, pidlikeli::maxpoq, 4);
  functions[pidlikeli::scast(Func::Yield)]
    = new TF1("fYield", pidfunc::Pol1, 0, pidlikeli::maxpoq, 2);
}
CorrGraph::CorrGraph() {
  using Graph = CorrGraph::Graph;
  for (size_t i = 0; i < static_cast<size_t>(Graph::COUNT); ++i) {
    graphs[i] = new TGraphErrors();
  }    
}
