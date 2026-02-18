// -*- C++ -*-

#ifndef DC_PARAMETERS_HH
#define DC_PARAMETERS_HH

#include <TString.h>

#include <std_ostream.hh>

//_____________________________________________________________________________
struct DCPairPlaneInfo
{
  Bool_t pair;
  Bool_t honeycomb;
  Bool_t fiber;
  Int_t  id1, id2;
  Double_t CellSize; // sense wire spacing, NOT max drift length
  void Print(const TString& arg="", std::ostream& ost=hddaq::cout) const
  {
    ost << "[DCPairPlaneInfo::Print()] " << arg << std::endl
	<< " pair      : " << pair      << std::endl
	<< " honeycomb : " << honeycomb << std::endl
	<< " fiber     : " << fiber     << std::endl
	<< " layer1    : " << id1       << std::endl
	<< " layer2    : " << id2       << std::endl
	<< " cell size : " << CellSize  << std::endl;
  }
};

extern const DCPairPlaneInfo PPInfoBcIn[], PPInfoBcOut[];
extern const Int_t NPPInfoBcIn, NPPInfoBcOut;

#ifdef DefStatic
const DCPairPlaneInfo PPInfoBcIn[] = {
  // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
  { true, false, false,  0,  1,  8.0 }, { true, false, false,  2,  3,  8.0 },
  { true, false, false,  4,  5,  8.0 }, { true, false, false,  6,  7,  8.0 },
  { true, false, false,  8,  9,  8.0 }, { true, false, false, 10, 11,  8.0 },
  { true, false, false, 12, 13,  8.0 }, { true, false, false, 14, 15,  8.0 }
};

const DCPairPlaneInfo PPInfoBcOut[] = {
  // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
  { true, false, false,  0,  1,  5.0 }, { true, false, false,  2,  3,  5.0 },
  { true, false, false,  4,  5,  5.0 }, { true, false, false,  6,  7,  5.0 },
  { true, false, false,  8,  9,  5.0 }, { true, false, false, 10, 11,  5.0 },
  { true, false, false, 12, 13,  5.0 }, { true, false, false, 14, 15,  5.0 }
};

const Int_t NPPInfoBcIn   = sizeof(PPInfoBcIn)/sizeof(DCPairPlaneInfo);
const Int_t NPPInfoBcOut  = sizeof(PPInfoBcOut)/sizeof(DCPairPlaneInfo);
#endif

// __ Legacy DL Range Arrays (NOT USED IN E72 CODE) _____________________
// Note: These arrays are defined but not currently used in E72 analysis.
// They exist in legacy code (ref/src/DCHit.cc) but are not referenced
// in the current implementation. Kept for reference only.
//
// DL Ranges (BC1&2 for Time range -5 ns <[Time gate]<75 ns)
const Double_t MinDLBc[25] = {
   0.0,
   // BC1
  -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
   // BC2
  -5.0, -5.0, -5.0, -5.0, -5.0, -5.0,
   // BC3
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
   // BC4
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5
};

const Double_t MaxDLBc[25] = {
  0.0,
  // BC1
  75.0, 75.0, 75.0, 75.0, 75.0, 75.0,
  // BC2
  75.0, 75.0, 75.0, 75.0, 75.0, 75.0,
  // BC3
  1.8, 1.8, 1.8, 1.8, 1.8, 1.8,
  // BC4
  1.8, 1.8, 1.8, 1.8, 1.8, 1.8
};
// _______________________________________________________________________


#endif
