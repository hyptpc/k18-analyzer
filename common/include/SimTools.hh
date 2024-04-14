// -*- C++ -*-

#ifndef SIM_TOOLS_HH
#define SIM_TOOLS_HH

#include <vector>
#include <string>

#include "ConfMan.hh"
#include "CDCAnalyzer.hh"
#include "DCAnalyzer.hh"
#include "HodoAnalyzer.hh"

#include "RootData.hh"

namespace sim
{
  void Convert(DetectorData*, CDCAnalyzer* cdc, DCAnalyzer* dc, HodoAnalyzer* hodo);
  void Convert(DetectorData*, CDCAnalyzer* cdc);
};

#endif
