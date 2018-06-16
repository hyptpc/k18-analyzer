/**
 *  file: DetectorID.hh
 *  date: 2017.04.10
 *
 */

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>

#define RunE07 0
#define RunE40 1

// Counters ___________________________________________________________
const int DetIdBH1      =  1;
const int DetIdBH2      =  2;
const int DetIdBAC      =  3;
const int DetIdFBH      =  4;
const int DetIdPVAC     =  5;
const int DetIdFAC      =  6;
const int DetIdSCH      =  7;
const int DetIdTOF      =  8;
const int DetIdSAC      =  9;
const int DetIdLC       = 10;
const int DetIdHtTOF    = 11; // high threshold TOF
const int NumOfSegBH1   = 11;
const int NumOfSegBH2   =  8;
const int NumOfSegBAC   =  2;
const int NumOfSegFBH   = 16;
const int NumOfSegCFBH  = 31;
const int NumOfSegPVAC  =  1;
const int NumOfSegFAC   =  1;
const int NumOfSegSCH   = 64;
const int NumOfSegTOF   = 24;
const int NumOfSegSAC   = 4;
const int NumOfSegLC    = 28;
const int NumOfSegHtTOF = 16;

// Misc _______________________________________________________________
const int DetIdTrig       = 21;
const int DetIdScaler     = 22;
const int DetIdEMC        = 23;
const int DetIdSSDT       = 24;
const int DetIdMsT        = 25;
const int DetIdMtx        = 26;
const int DetIdFpgaBH1    = 27;
const int DetIdFpgaBH2    = 28;
const int DetIdFpgaBH2Mt  = 29;
const int DetIdVmeRm      = 81;
const int DetIdMsTRM      = 82;
const int DetIdHulRM      = 83;
const int NumOfSegTrig    = 32;
const int NumOfSegScaler  = 96;
const int NumOfSegEMC     =  1;
const int NumOfSegSSDT    = 16;
const int SpillEndFlag    = 20; // 0-based
const int NumOfPlaneVmeRm = 2;

#if RunE07
enum eTriggerFlag
  {
    kUB       =  0,
    kUBUB     =  1,
    kKPI      =  2,
    kKK       =  3,
    kUBPS     =  4,
    kUBUBPS   =  5,
    kKPIPS    =  6,
    kKKPS     =  7,
    kMtx1     =  8,
    kMtx2     =  9,
    kMtx3     = 10,
    kMtx4     = 11,
    kKIn      = 12,
    kPiIn     = 13,
    kKOut     = 14,
    kPiOut    = 15,
    kKBeam    = 16,
    kPiBeam   = 17,
    kPBeam    = 18,
    kClock    = 19,
    kSpillEnd = 20
  };
#endif

#if RunE40
enum eTriggerFlag
  {
    kBh21K     =  0,
    kBh22K     =  1,
    kBh23K     =  2,
    kBh24K     =  3,
    kBh25K     =  4,
    kBh26K     =  5,
    kBh27K     =  6,
    kBh28K     =  7,
    kBh2K      =  8,
    kElseOr    =  9,
    kBeam      = 10,
    kBeamTof   = 11,
    kBeamPi    = 12,
    kBeamP     = 13,
    kCoin1     = 14,
    kCoin2     = 15,
    kE03       = 16,
    // kBh2KPs    =
    kBeamPs    = 17,
    kBeamTofPs = 18,
    kBeamPiPs  = 19,
    kBeamPPs   = 10,
    kCoin1Ps   = 21,
    kCoin2Ps   = 22,
    kE03Ps     = 23,
    kClk       = 24,
    kRsv2      = 25,
    kMtx       = 26,
    kSclrTrg   = 27
    // kClkGe   =
  };
#endif

const int DetIdVmeCalib      = 999;
const int NumOfPlaneVmeCalib =   5;
const int NumOfSegVmeCalib   =  32;

// Trackers ___________________________________________________________
const int DetIdBC1  = 101;
const int DetIdBC2  = 102;
const int DetIdBC3  = 103;
const int DetIdBC4  = 104;
const int DetIdSDC1 = 105;
const int DetIdSDC2 = 106;
const int DetIdSDC3 = 107;
const int DetIdBFT  = 110;
const int DetIdSFT  = 111;
const int DetIdSSD1 = 121;
const int DetIdSSD2 = 122;
const int DetIdFBT1 = 131;
const int DetIdFBT2 = 132;

const int PlMinBcIn        =   1;
const int PlMaxBcIn        =  12;
const int PlMinBcOut       =  13;
const int PlMaxBcOut       =  24;
const int PlMinSdcIn       =   1;
// const int PlMaxSdcIn       =  14;
const int PlMaxSdcIn       =   9;
const int PlMinSdcOut      =  31;
const int PlMaxSdcOut      =  38;
const int PlMinSsdIn       =   1;
const int PlMaxSsdIn       =   4;
const int PlMinSsdOut      =   5;
const int PlMaxSsdOut      =   8;
const int PlMinSsdOutSdcIn =   1;
const int PlOffsBc         = 100;
const int PlOffsSdcIn      =   0;
const int PlOffsSft        =   6;
const int PlOffsSdcOut     =  30;
const int PlOffsSsd        =   6;

const int NumOfLayersBc     = 6;
const int NumOfLayersSSD1   = 4;
const int NumOfLayersSSD2   = 4;
const int NumOfLayersSFT    = 3;
const int NumOfLayersSDC1   = 6;
const int NumOfLayersSDC2   = 4;
const int NumOfLayersSDC3   = 4;
const int NumOfLayersBcIn   = PlMaxBcIn   - PlMinBcIn   + 1;
const int NumOfLayersBcOut  = PlMaxBcOut  - PlMinBcOut  + 1;
const int NumOfLayersSsdIn  = PlMaxSsdIn  - PlMinSsdIn  + 1;
const int NumOfLayersSsdOut = PlMaxSsdOut - PlMinSsdOut + 1;
const int NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const int NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;

const int MaxWireBC1      = 256;
const int MaxWireBC2      = 256;
const int MaxWireBC3      =  64;
const int MaxWireBC4      =  64;

const int MaxWireSDC1     =  64;
const int MaxWireSDC2     = 128;
const int MaxWireSDC3X    =  96;
const int MaxWireSDC3Y    =  64;

const int    NumOfSampleSSD      =    8;
const double SamplingIntervalSSD =  25.; //[ns]
const int    NumOfSegSSD0        = 1536;
const int    NumOfSegSSD1        = 1536;
const int    NumOfSegSSD2        = 1536;

const int NumOfPlaneBFT   =   2;
const int NumOfSegBFT     = 160;
// SFT X layer has U D plane.
// SFT UV layers have only U plnane.
// enum SFT_PLANE{ SFT_X1, SFT_X2, SFT_V, SFT_U };
enum SFT_PLANE{ SFT_U, SFT_V, SFT_X1, SFT_X2 };
const int NumOfPlaneSFT   =   4;
const int NumOfSegSFT_X   = 256;
const int NumOfSegSFT_UV  = 320;
const int NumOfSegCSFT    = 48;

// HulRm -----------------------------------------------
const int NumOfHulRm   = 4;

// MsT -------------------------------------------------
enum TypesMst{typeHrTdc, typeLrTdc, typeFlag, NumOfTypesMst};
const int NumOfMstHrTdc = 32;
const int NumOfMstLrTdc = 64;
const int NumOfMstFlag  = 7;
enum dTypesMst
  {
    mstClear,
    mstAccept,
    finalClear,
    cosolationAccept,
    fastClear,
    level2,
    noDecision,
    size_dTypsMsT
  };

// Scaler ----------------------------------------------
const int NumOfScaler  = 2;

#endif
