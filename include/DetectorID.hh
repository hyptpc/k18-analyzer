// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>
#include <vector>

#include <TString.h>

// Counters ___________________________________________________________
const int DetIdBH1      =  1;
const int DetIdBH2      =  2;
const int DetIdBAC      =  3;
const int DetIdSCH      =  6;
const int DetIdTOF      =  7;
const int DetIdHTOF     =  8;
const int DetIdBVH      =  9;
const int DetIdLAC      = 12;
const int DetIdWC       = 13;
const int NumOfSegBH1   = 11;
const int NumOfSegBH2   =  8;
const int NumOfSegBAC   =  2;
const int NumOfSegSCH   = 64;
const int NumOfSegTOF   = 24;
const int NumOfSegHTOF  = 34;
const int NumOfSegBVH   =  4;
const int NumOfSegLAC   = 30;
const int NumOfSegWC    = 20;

// Misc _______________________________________________________________
const int DetIdTrig       = 21;
const int DetIdScaler     = 22;
const int DetIdMsT        = 25;
const int DetIdMtx        = 26;
const int DetIdVmeRm      = 81;
const int DetIdMsTRM      = 82;
const int DetIdHulRM      = 83;
const int NumOfSegScaler  = 96;
const int NumOfPlaneVmeRm = 2;

// Trigger Flag
namespace trigger
{
  enum ETriggerFlag
  {
    kL1SpillOn,
    kL1SpillOff,
    kSpillEnd,
    kSpillOnEnd,
    kTofTiming,
    kMatrix2D1,
    kMatrix2D2,
    kMatrix3D,
    kBeamA,
    kBeamB,
    kBeamC,
    kBeamD,
    kBeamE,
    kBeamF,
    kTrigA,
    kTrigB,
    kTrigC,
    kTrigD,
    kTrigE,
    kTrigF,
    kTrigAPS,
    kTrigBPS,
    kTrigCPS,
    kTrigDPS,
    kTrigEPS,
    kTrigFPS,
    kLevel1A,
    kLevel1B,
    kClockPS,
    kReserve2PS,
    kLevel1OR,
    kEssDischarge,
    NTriggerFlag
  };

  const std::vector<TString> STriggerFlag =
    {
     "L1SpillOn",
     "L1SpillOff",
     "SpillEnd",
     "SpillOnEnd",
     "TofTiming",
     "Matrix2D1",
     "Matrix2D2",
     "Matrix3D",
     "BeamA",
     "BeamB",
     "BeamC",
     "BeamD",
     "BeamE",
     "BeamF",
     "TrigA",
     "TrigB",
     "TrigC",
     "TrigD",
     "TrigE",
     "TrigF",
     "TrigA-PS",
     "TrigB-PS",
     "TrigC-PS",
     "TrigD-PS",
     "TrigE-PS",
     "TrigF-PS",
     "Level1A",
     "Level1B",
     "Clock-PS",
     "Reserve2-PS",
     "Level1OR",
     "EssDischarge",
    };
}
const Int_t NumOfSegTrig = trigger::NTriggerFlag;

const int DetIdVmeCalib      = 999;
const int NumOfPlaneVmeCalib =   5;
const int NumOfSegVmeCalib   =  32;

// Trackers ___________________________________________________________
const int DetIdBC3  = 103;
const int DetIdBC4  = 104;
const int DetIdSDC1 = 105;
const int DetIdSDC2 = 106;
const int DetIdSDC3 = 107;
const int DetIdSDC4 = 108;
const int DetIdBFT  = 110;

const int PlMinBcIn        =   1;
const int PlMaxBcIn        =  12;
const int PlMinBcOut       =  13;
const int PlMaxBcOut       =  24;
const int PlMinSdcIn       =   1;
const int PlMaxSdcIn       =  10;
const int PlMinSdcOut      =  31;
const int PlMaxSdcOut      =  38;
const int PlMinTOF         =  41;
const int PlMaxTOF         =  44;
const int PlOffsBc         = 100;
const int PlOffsSdcIn      =   0;
const int PlOffsSdcOut     =  30;
const int PlOffsTOF        =  40;
const int PlOffsVP         =  20;
const int PlOffsTPCX       = 600;
const int PlOffsTPCY       = 650;

const int NumOfLayersBc     = 6;
const int NumOfLayersSDC1   = 6;
const int NumOfLayersSDC2   = 4;
const int NumOfLayersSDC3   = 4;
const int NumOfLayersSDC4   = 4;
const int NumOfLayersBcIn   = PlMaxBcIn   - PlMinBcIn   + 1;
const int NumOfLayersBcOut  = PlMaxBcOut  - PlMinBcOut  + 1;
const int NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const int NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;
const int NumOfLayersTOF    = PlMaxTOF    - PlMinTOF    + 1;
const int NumOfLayersVP     = 5;
const int NumOfLayersTPC    = 32;
const int NumOfPadTPC       = 5768;
const int NumOfTimeBucket   = 200;

const int MaxWireBC3      =  64;
const int MaxWireBC4      =  64;

const int MaxWireSDC1     =  64;
const int MaxWireSDC2X    =  70;
const int MaxWireSDC2Y    =  40;
const int MaxWireSDC3     = 128;
const int MaxWireSDC4X    =  96;
const int MaxWireSDC4Y    =  64;

const int NumOfPlaneBFT   =   2;
const int NumOfSegBFT     = 160;

// HulRm -----------------------------------------------
const int NumOfHulRm   = 4;

// Matrix ----------------------------------------------
const int NumOfSegSFT_Mtx = 48;

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
