// -*- C++ -*-

#ifndef DC_EXCLUSIVE_PULL_HH
#define DC_EXCLUSIVE_PULL_HH

#include <Rtypes.h>

class DCLocalTrack;

namespace DCExclusivePull
{
  /// Leave-one-out residual [mm]. On failure returns QuietNaN.
  Double_t Residual(const DCLocalTrack& track, Int_t ihit);
  /// Leave-one-out pull = residual / resolution. On failure returns QuietNaN.
  Double_t Pull(const DCLocalTrack& track, Int_t ihit);
}

#endif
