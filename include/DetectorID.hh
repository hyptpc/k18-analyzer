// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <iostream>
#include <map>

#include <Rtypes.h>
#include <TString.h>

static const Char_t* const UorD[2] = {"U", "D"};

const std::map<TString, std::vector<TString>> DCNameList =
{
  {"BcIn", { "BLC1a", "BLC1b" }},
  {"BcOut", { "BLC2a", "BLC2b" }},
};

// Counters ___________________________________________________________
const Int_t DetIdBHT  =  1;
const Int_t DetIdT0   =  2;
const Int_t DetIdBAC  =  3;
const Int_t DetIdKVC1 =  4;
const Int_t DetIdSAC  =  5;
const Int_t DetIdBH2  =  6;
const Int_t DetIdHTOF =  7;
const Int_t DetIdKVC2 =  8;
const Int_t DetIdBVH  = 20;
const Int_t DetIdT1   = 21;
const Int_t DetIdT2   = 22;
const Int_t NumOfSegBHT  = 63;
const Int_t NumOfSegT0   =  5;
const Int_t NumOfSegBAC  =  5;
const Int_t NumOfSegKVC1 =  5;
const Int_t NumOfSegSAC  =  9;
const Int_t NumOfSegBH2  = 16;
const Int_t NumOfSegHTOF = 34;
const Int_t NumOfSegKVC2 =  5;
const Int_t NumOfSegBVH  =  4;
const Int_t NumOfSegT1   =  1;
const Int_t NumOfSegT2   =  1;

const int DetIdVmeRm     =  81;
const int DetIdScaler    =  91;
const int DetIdTrigFlag  =  99;
const Int_t NumOfSegTrigFlag = 32;

#if defined E72
enum EHodoscope {
  kBHT, kT0, kBAC,
  kKVC1, kSAC, kBH2,
  kHTOF,
  kBVH, kT1, kT2,
  // kCVC, kNC,
  kNumHodo
};

const Int_t DetIdHodo[kNumHodo] = {
  DetIdBHT, DetIdT0, DetIdBAC,
  DetIdKVC1, DetIdSAC, DetIdBH2,
  DetIdHTOF,
  DetIdBVH, DetIdT1, DetIdT2,
  // DetIdCVC, DetIdNC,
};

const std::vector<TString> NameHodo = {
  "BHT", "T0", "BAC",
  "KVC1", "SAC", "BH2", "HTOF",
  "BVH", "T1", "T2",
  // "CVC", "NC",
};

const double  NumOfSegHodo[kNumHodo] = {
  NumOfSegBHT, NumOfSegT0, NumOfSegBAC,
  NumOfSegKVC1, NumOfSegSAC, NumOfSegBH2,
  NumOfSegHTOF,
  NumOfSegBVH, NumOfSegT1, NumOfSegT2
  // NumOfSegCVC, NumOfSegNC,
};

#elif defined E73_2024
enum EHodoscope {
  kBHT, kT1, kT0, kDEF,
  // kVeto, kBTC,
  kCVC, kNC,
  // kRC,
  kNumHodo
};

const Int_t DetIdHodo[kNumHodo] = {
  DetIdBHT, DetIdT1, DetIdT0, DetIdDEF,
  // DetIdVeto, DetIdBTC,
  DetIdCVC, DetIdNC,
  // DetIdRC
};

const std::vector<TString> NameHodo = {
  "BHT", "T1", "T0", "DEF",
  // "Veto", "BTC",
  "CVC", "NC",
  // "RC"
};

const double  NumOfSegHodo[kNumHodo] = {
  63, 1, 5, 5,
  // 4, 4,
  35, 6,
};
#endif

// Chambers
const int DetIdCDC    = 100;
const int DetIdBLC1a  = 101;
const int DetIdBLC1b  = 102;
const int DetIdBLC2a  = 103;
const int DetIdBLC2b  = 104;
const int DetIdBPC    = 105;
const int DetIdBPC2   = 105;
const int DetIdBPC1   = 106;
const int DetIdBPCmini= 106;
const int DetIdSDC    = 106;
const int DetIdFDC    = 107;
const int DetIdBLC1   = 111;
const int DetIdBLC2   = 112;
const int DetIdBPC0   = 113;

enum EDC {
  kBLC1a, kBLC1b, kBLC2a, kBLC2b,
  kBPC1, kBPC2,
  // kVFT,
  kBLC1, kBLC2, kBPC0, kNumDC
};

const Int_t DetIdDC[kNumDC] = {
  DetIdBLC1a, DetIdBLC1b, DetIdBLC2a, DetIdBLC2b,
  DetIdBPC1, DetIdBPC2,
  // DetIdVFT,
  DetIdBLC1, DetIdBLC2, DetIdBPC0
};

// Hall etc.
const int DetIdHall       	= 200;
const int DetIdFloor		= 201;
// const int DetIdBeamDump		= 110;
// const int DetIdSideDump		= 111;
// const int DetIdNShield		= 112;
// const int DetIdSideCon		= 113;
// const int DetIdDoorCon		= 114;
// const int DetIdDumpCon		= 115;
const int DetIdDoraemon		= 220;
// const int DetIdUSWK		= 121;
// const int DetIdCDSBobbin	= 125;
// const int DetIdCDSCoil		= 126;
const int DetIdCDCCFRP		= 230;
const int DetIdCDCMylar		= 231;
// const int DetIdCDCEndCap	= 132;
// Target etc.
const int DetIdTarSys		= 240;
const int DetIdRadS		= 241;
const int DetIdTarCFRP		= 242;
const int DetIdTarCap		= 243;
const int DetIdTarRing		= 244;
// const int DetIdTarChm		= 145;
// const int DetIdBeamWindow	= 146;
// const int DetIdScatWindow	= 147;
// const int DetIdMagShield	= 148;

const int DetIdTarCell		= 250;
const int DetIdTarget		= 251;
const int DetIdCellTube		= 252;
const int DetIdCellFlange	= 253;
// const int DetIdBShield		= 154;
// const int DetIdCellWindow	= 154;
// const int DetIdBFrange		= 155;
// const int DetIdCellRing		= 156;
const int DetIdFiducial		= 260;

const int NumOfPlaneVmeRm=3;

namespace beam
{
enum EBeamFlag {
  kAll, kPion, kKaon, /* kProton, */ kUnknown, kBeamFlag
};
const std::vector<TString> BeamFlagList{"", "_Pi", "_K", /* "_P" */};
}

namespace trigger
{
enum ETriggerFlag {
  kStart, kStop, kBeam, kPion, kKaon2, kKaon3, //5
  kKCDH1, kKCDH2, kKCDH3, //8
  kKCDH1Gamma, kKGamma, kPiCDH, //11
  kBeamPbF2, kPiPbF2=kBeamPbF2, //12
  kCDHCosmic, kElectronPbF2=kCDHCosmic, //13
  kPbF2Cosmic, kCosmic=kPbF2Cosmic, //14
  kProtonTrig, kDeuteron=kProtonTrig, //15
  kNumTrig
};

}

// enum eTriggerFlag
//   {
//     kSpillStart =  1,
//     kSpillEnd   =  2,
//     kBeam1      =  3,
//     kBeam2      =  4,
//     kBeam_f     =  5,
//     kPion       =  6,
//     kPion_f     =  7,
//     kKaon1      =  8,
//     kKaon2      =  9,
//     kKaon3      =  10,
//     kKaon1_f       = 11,
//     kKaon2_f       = 12,
//     kKaonStart     = 13,
//     kKaonStartStop = 14,
//     kStartStop     = 15,
//     kMisc          = 16
//   };

// for compatibility
const Int_t LayerMinBcIn        =   1;
const Int_t LayerMaxBcIn        =  16;
const Int_t LayerMinBcOut       = 1001;
const Int_t LayerMaxBcOut       = 1016;
const Int_t LayerMinSdcIn       =   1;
const Int_t LayerMaxSdcIn       =  10;
const Int_t LayerMinSdcOut      =  31;
const Int_t LayerMaxSdcOut      =  42;
const Int_t LayerMinTOF         =  51; // need to change
const Int_t LayerMaxTOF         =  54; // need to change
const Int_t LayerMinVP          =  16;
const Int_t LayerMaxVP          =  26;
const Int_t PlOffsBc         = 100;
const Int_t PlOffsSdcIn      =   0;
const Int_t PlOffsSdcOut     =  30;
const Int_t PlOffsTOF        =  50;
const Int_t PlOffsVP         =  15;

// const Int_t NumOfLayersBc     = 8;
const Int_t NumOfLayersSDC1   = 6;
const Int_t NumOfLayersSDC2   = 4;
const Int_t NumOfLayersSDC3   = 4;
const Int_t NumOfLayersSDC4   = 4;
const Int_t NumOfLayersSDC5   = 4;
const Int_t NumOfLayersBcIn   = LayerMaxBcIn   - LayerMinBcIn   + 1;
const Int_t NumOfLayersBcOut  = LayerMaxBcOut  - LayerMinBcOut  + 1;
const Int_t NumOfLayersSdcIn  = LayerMaxSdcIn  - LayerMinSdcIn  + 1;
const Int_t NumOfLayersSdcOut = LayerMaxSdcOut - LayerMinSdcOut + 1;
const Int_t NumOfLayersTOF    = LayerMaxTOF    - LayerMinTOF    + 1;
const Int_t NumOfLayersVP     = LayerMaxVP     - LayerMinVP     + 1;

#endif
