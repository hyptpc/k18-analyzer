/*a*
 *  file: DCParameters.hh
 *  date: 2017.04.10
 *
 */

#ifndef DC_PARAMETERS_HH
#define DC_PARAMETERS_HH

#include <std_ostream.hh>

//______________________________________________________________________________
struct DCPairPlaneInfo
{
  bool pair;
  bool honeycomb;
  bool fiber;
  int  id1, id2;
  double CellSize;

  void Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const
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

extern const DCPairPlaneInfo PPInfoBcOut[], PPInfoSdcIn[], PPInfoSdcOut[], PPInfoCFT[], PPInfoCFT16[];
extern const int NPPInfoBcOut, NPPInfoSdcIn, NPPInfoSdcOut, NPPInfoCFT, NPPInfoCFT16;

#ifdef DefStatic
const DCPairPlaneInfo PPInfoBcOut[] = {
  // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
  { true, false, false, 1,  2,  3.0 }, { true, false, false,  3,  4,  3.0 },
  { true, false, false, 5,  6,  3.0 }, { true, false, false,  7,  8,  3.0 },
  { true, false, false, 9, 10,  3.0 }, { true, false, false, 11, 12,  3.0 },
  { false, false, false, 13, 13, 10.0 }
};

const DCPairPlaneInfo PPInfoSdcIn[] = {
  // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
  { true, true, false, 1,  2,  6.0 }, { true, true, false, 3,  4,  6.0 },
  { true, true, false, 5,  6,  6.0 },
  { false, false, true, 7, 7, 3.0 }, { false, false, true, 8, 8, 3.0 }, { false, false, true, 9, 9, 3.0 } // SFT
};
// const DCPairPlaneInfo PPInfoSdcIn[] = {
//   // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
//   { true, true, false, 1,  2,  6.0 }, { true, true, false, 3,  4,  6.0 },
//   { true, true, false, 5,  6,  6.0 }
// };

const DCPairPlaneInfo PPInfoSdcOut[] = {
  // { pair_plane, honeycomb, fiber, id1, id2, CellSize }
  { true, true, false, 1, 2,  9.0 }, { true, true, false, 3, 4,  9.0 },
  { true, true, false, 5, 6, 20.0 }, { true, true, false, 7, 8, 20.0 }, 
  
  { false, false, true, 9, 10, 4.0 }, { false, false, true, 11, 12, 4.0 }, //FBT1
  { false, false, true, 13, 14, 4.0 }, { false, false, true, 15, 16, 4.0 } //FBT2
};

const DCPairPlaneInfo PPInfoCFT[] = {  
  { false, false, true,  0,  0,  3.0 }, { false, false, true,  1,  1,  3.0 },
  { false, false, true,  2,  2,  3.0 }, { false, false, true,  3,  3,  3.0 },
  { false, false, true,  4,  4,  3.0 }, { false, false, true,  5,  5,  3.0 },
  { false, false, true,  6,  6,  3.0 }, { false, false, true,  7,  7,  3.0 },
};
const DCPairPlaneInfo PPInfoCFT16[] = {
  { false, false, true,  0,  0,  3.0 }, { false, false, true,  1,  1,  3.0 },
  { false, false, true,  2,  2,  3.0 }, { false, false, true,  3,  3,  3.0 },
  { false, false, true,  4,  4,  3.0 }, { false, false, true,  5,  5,  3.0 },
  { false, false, true,  6,  6,  3.0 }, { false, false, true,  7,  7,  3.0 },
  { false, false, true,  8,  8,  3.0 }, { false, false, true,  9,  9,  3.0 },
  { false, false, true, 10, 10,  3.0 }, { false, false, true, 11, 11,  3.0 },
  { false, false, true, 12, 12,  3.0 }, { false, false, true, 13, 13,  3.0 },
  { false, false, true, 14, 14,  3.0 }, { false, false, true, 15, 15,  3.0 },
};


const int NPPInfoBcOut  = sizeof(PPInfoBcOut)/sizeof(DCPairPlaneInfo);
const int NPPInfoSdcIn  = sizeof(PPInfoSdcIn)/sizeof(DCPairPlaneInfo);
const int NPPInfoSdcOut = sizeof(PPInfoSdcOut)/sizeof(DCPairPlaneInfo);
const int NPPInfoCFT    = sizeof(PPInfoCFT)   /sizeof(DCPairPlaneInfo);
const int NPPInfoCFT16  = sizeof(PPInfoCFT16)   /sizeof(DCPairPlaneInfo);

#endif

//DL Ranges (BC1&2 for Time range -5 ns <[Time gate]<75 ns)
const double MinDLBc[25] = {
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

const double MaxDLBc[25] = {
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

const double MinDLSdc[] = {
  0.0,
  //SDC1
  -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
  // Dummy Id 07-30
  0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SDC2
  -0.5, -0.5, -0.5, -0.5,
  // SDC3
  -0.5, -0.5, -0.5, -0.5
};

const double MaxDLSdc[] = {
  0.0,
  // SDC1
  10.0, 10.0, 10.0, 10.0, 10.0, 10.0,
  // Dummy Id 07-30
  0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
  // SDC2
  //  12.0, 12.0, 12.0, 12.0,
  4.8, 4.8, 4.8, 4.8,
  // SDC3
  //  25.0, 25.0, 25.0, 25.0
  13.00, 13.00, 13.00, 13.00
};

#endif
