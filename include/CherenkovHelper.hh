// -*- C++ -*-

#ifndef CHERENKOV_HELPER_HH
#define CHERENKOV_HELPER_HH

#include <vector>

#include <Rtypes.h>
#include <TString.h>

class RawData;

namespace CherenkovHelper
{
  /// Offline Npe from raw. BAC: 1 (seg 0–3); KVC: NumOfSegKVC (ch0–3/seg). Other: empty.
  std::vector<Double_t> Compute(const RawData& rawData, const TString& name);
}

#endif
