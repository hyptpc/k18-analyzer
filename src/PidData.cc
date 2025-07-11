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
