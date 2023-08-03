// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>
#include <map>
#include <vector>
#include <array>
#include <TString.h>

const std::map<TString, std::vector<TString>> DCNameList =
{
  {"BcOut", { "BC3", "BC4" }},
  {"SdcIn", { "SDC1", "SDC2" }},
  {"SdcOut", { "SDC3", "SDC4", "SDC5" }},
};

// Counters ___________________________________________________________
const Int_t NumOfSegBH1   = 11;
const Int_t NumOfSegBH2   =  8;
const Int_t NumOfSegBAC   =  2;
const Int_t NumOfSegSCH   = 64;
const Int_t NumOfSegTOF   = 19;
const Int_t NumOfSegHTOF  = 34;
const Int_t NumOfSegBVH   =  4;
const Int_t NumOfSegAC1   = 30;
const Int_t NumOfSegWC    = 12;
const Int_t NumOfSegSAC3  = 2;
const Int_t NumOfSegSFV   = 6;

// AFT
const Int_t DetIdAFT      = 112;
const Int_t NumOfPlaneAFT = 36;
const Int_t NumOfSegAFTX  = 32;
const Int_t NumOfSegAFTY  = 16;
const Int_t NumOfSegAFT   = 32;
const std::vector<int> NumOfSegAFTarr = { 32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16,
					  32, 32, 16, 16 };
//const Int_t NumOfSegAFT[4]    = {NumOfSegAFTX, NumOfSegAFTX, NumOfSegAFTY, NumOfSegAFTY};

//const Int_t NumOfSegAFT[4]    = {NumOfSegAFTX, NumOfSegAFTX, NumOfSegAFTY, NumOfSegAFTY};

// VMEEASIROC
const Int_t DetIdVMEASIROC = 116;
const Int_t NumOfPlaneVMEEASIROC = 96;
const Int_t NumOfSegVMEEASIROC = 64;

// Misc _______________________________________________________________
const Int_t DetIdTrig       = 21;
const Int_t DetIdScaler     = 22;
const Int_t DetIdMsT        = 25;
const Int_t DetIdMtx        = 26;
const Int_t DetIdVmeRm      = 81;
const Int_t DetIdMsTRM      = 82;
const Int_t DetIdHulRM      = 83;
const Int_t NumOfSegScaler  = 96;
const Int_t NumOfPlaneVmeRm = 2;

// Trigger Flag
namespace trigger
{
  enum ETriggerFlag
  {
    kL1SpillOn,
    kL1SpillOff,
    kSpillOnEnd,
    kSpillOffEnd,
    kCommonStopSdcOut,
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
     "CommonStopSdcOut",
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

const Int_t DetIdVmeCalib      = 999;
const Int_t NumOfPlaneVmeCalib =   5;
const Int_t NumOfSegVmeCalib   =  32;

// Trackers ___________________________________________________________
// const Int_t DetIdBC3  = 103;
// const Int_t DetIdBC4  = 104;
// const Int_t DetIdSDC1 = 105;
// const Int_t DetIdSDC2 = 106;
// const Int_t DetIdSDC3 = 107;
// const Int_t DetIdSDC4 = 108;
// const Int_t DetIdSDC5 = 109;
// const Int_t DetIdBFT  = 110;

const Int_t PlMinBcIn        =   1;
const Int_t PlMaxBcIn        =  12;
const Int_t PlMinBcOut       = 113;
const Int_t PlMaxBcOut       = 124;
const Int_t PlMinSdcIn       =   1;
const Int_t PlMaxSdcIn       =  10;
const Int_t PlMinSdcOut      =  31;
const Int_t PlMaxSdcOut      =  42;
const Int_t PlMinTOF         =  51; // need to change
const Int_t PlMaxTOF         =  54; // need to change
const Int_t PlMinVP          =  21;
const Int_t PlMaxVP          =  30;
const Int_t PlOffsBc         = 100;
const Int_t PlOffsSdcIn      =   0;
const Int_t PlOffsSdcOut     =  30;
const Int_t PlOffsTOF        =  50;
const Int_t PlOffsVP         =  20;
// const Int_t PlOffsTPCX       = 600;
// const Int_t PlOffsTPCY       = 650;

const Int_t NumOfLayersBc     = 6;
const Int_t NumOfLayersSDC1   = 6;
const Int_t NumOfLayersSDC2   = 4;
const Int_t NumOfLayersSDC3   = 4;
const Int_t NumOfLayersSDC4   = 4;
const Int_t NumOfLayersSDC5   = 4;
const Int_t NumOfLayersBcIn   = PlMaxBcIn   - PlMinBcIn   + 1;
const Int_t NumOfLayersBcOut  = PlMaxBcOut  - PlMinBcOut  + 1;
const Int_t NumOfLayersSdcIn  = PlMaxSdcIn  - PlMinSdcIn  + 1;
const Int_t NumOfLayersSdcOut = PlMaxSdcOut - PlMinSdcOut + 1;
const Int_t NumOfLayersTOF    = PlMaxTOF    - PlMinTOF    + 1;
const Int_t NumOfLayersVP     = PlMaxVP     - PlMinVP     + 1;
// const Int_t NumOfLayersTPC    = 32;
// const Int_t NumOfPadTPC       = 5768;
// const Int_t NumOfTimeBucket   = 170;

const Int_t MaxWireBC3      =  64;
const Int_t MaxWireBC4      =  64;

const Int_t MaxWireSDC1     =  64;
const Int_t MaxWireSDC2     =  44;
const Int_t MaxWireSDC3     = 128;
const Int_t MaxWireSDC4     = 128;
const Int_t MaxWireSDC5X    = 128;
const Int_t MaxWireSDC5Y    =  96;

// MaxDriftLength = CellSize/2
const Double_t CellSizeBC3 = 3.0;
const Double_t CellSizeBC4 = 3.0;
const Double_t CellSizeSDC1 =  6.0;
const Double_t CellSizeSDC2 = 10.0;
const Double_t CellSizeSDC3 =  9.0;
const Double_t CellSizeSDC4 =  9.0;
const Double_t CellSizeSDC5 =  9.0;

const Int_t NumOfPlaneBFT   =   2;
const Int_t NumOfSegBFT     = 160;

// HulRm -----------------------------------------------
const Int_t NumOfHulRm   = 4;

// Matrix ----------------------------------------------
const Int_t NumOfSegSFT_Mtx = 48;

// MsT -------------------------------------------------
enum TypesMst{typeHrTdc, typeLrTdc, typeFlag, NumOfTypesMst};
const Int_t NumOfMstHrTdc = 32;
const Int_t NumOfMstLrTdc = 64;
const Int_t NumOfMstFlag  = 7;
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
const Int_t NumOfScaler  = 2;

// Parasite ___________________________________________________________
const Int_t DetIdE72BAC      =  501;
const Int_t DetIdE90SAC      =  502;
const Int_t DetIdE72KVC      =  503;
const Int_t DetIdE42BH2      =  504;
const Int_t DetIdT1          =  505;
const Int_t DetIdT2          =  506;
const Int_t NumOfSegE72BAC   =  1;
const Int_t NumOfSegE90SAC   =  2;
const Int_t NumOfSegE72KVC   =  4;
const Int_t NumOfSegE42BH2   =  8;
const Int_t NumOfSegT1       =  1;
const Int_t NumOfSegT2       =  1;

#endif
