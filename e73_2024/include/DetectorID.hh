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
const int DetIdCDH    =   0;
const int DetIdBHD    =   1;
const int DetIdBHT    =   1;
const int DetIdT0     =   2;
const int DetIdAC     =   3;
const int DetIdT1     =   4;
const int DetIdE0     =   5;
const int DetIdDEF    =   6;
const int DetIdStart  =  11;
const int DetIdStop   =  12;
const int DetIdCVC    =  14;
const int DetIdSDD    =  20;
const int DetIdSDDVeto=  25;
const int DetIdNC     =  32;
const int DetIdBVC    =  33;
const int DetIdPC     =  35;
const int DetIdIH     =  42;
const int DetIdVFT    =  53;

const Int_t NumOfSegBH2   =  8;
const Int_t NumOfSegTOF   = 19;


const Int_t NumOfSegBHT = 63;
const Int_t NumOfSegAC  = 5;

const int DetIdPbG    =  29;
const int DetIdPbF2   =  30;
const int DetIdVeto   =  31;

const int DetIdBTC    =  33;
const int DetIdCHCbarrel  =  61;
const int DetIdCHCcapF    =  62;
const int DetIdCHCcapB    =  63;
const int DetIdNCbarrel  =  66;
#if E73
const int DetIdVeto1  =  31;
const int DetIdVeto0  =  32;
const int DetIdFinger =  32;
const int DetIdLeak   =  34;
#endif
#if T98
const int DetIdRC     =  98;
#endif
#if E15
const int DetIdNC     =  32;
const int DetIdLB     =  36;
const int DetIdWVC    =  37;
const int DetIdBPD    =  41;
const int DetIdBD     =  90;
#endif
// const int NumOfSegBHD    =  11;
// const int NumOfSegT0     =  5;
const int DetIdVmeRm     =  81;
const int DetIdScaler    =  91;
const int DetIdTrigFlag  =  99;
const Int_t NumOfSegTrigFlag = 32;

enum EHodoscope {
  kBHT, kT1, kT0, kDEF,
  kVeto, kBTC,
#ifdef CDS
  kCDH, kPbG, kPbF2,
#endif
  kCVC, kNC,
  // kRC,
  kNumHodo
};

const Int_t DetIdHodo[kNumHodo] = {
  DetIdBHT, DetIdT1, DetIdT0, DetIdDEF,
  DetIdVeto, DetIdBTC,
#ifdef CDS
  DetIdCDH, DetIdPbG, DetIdPbF2,
#endif
  DetIdCVC, DetIdNC,
  // DetIdRC
};

const std::vector<TString> NameHodo = {
  "BHT", "T1", "T0", "DEF",
  "Veto", "BTC",
#ifdef CDS
  "CDH", "PbG", "PbF2",
#endif
  "CVC", "NC",
  // "RC"
};

const double  NumOfSegHodo[kNumHodo] = {
  63, 1, 5, 5,
  4, 4,
#ifdef CDS
  36, 40, 40,
#endif
  35, 6,
};

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
  kBPC1, kBPC2, kVFT,
#ifdef CDS
  kCDC,
#endif
  kBLC1, kBLC2, kBPC0, kNumDC
};

const Int_t DetIdDC[kNumDC] = {
  DetIdBLC1a, DetIdBLC1b, DetIdBLC2a, DetIdBLC2b,
  DetIdBPC1, DetIdBPC2,
  DetIdVFT,
#ifdef CDS
  DetIdCDC,
#endif
  DetIdBLC1, DetIdBLC2, DetIdBPC0
};

const TString NameDC[kNumDC] = {
  "BLC1a", "BLC1b", "BLC2a", "BLC2b",
  "BPC1", "BPC2",
  "VFT",
#ifdef CDS
  "CDC",
#endif
  "BLC1", "BLC2", "BPC0"
};

const int  NumOfLayerDC[kNumDC] = {
  8, 8, 8, 8,
  8, 8, 14,
#ifdef CDS
  118,
#endif
  16, 16, 16
};

const int  NumOfWireDC[kNumDC] = {
  32, 32, 32, 32,
  15, 32, 64,
#ifdef CDS
  16,
#endif
  32, 32, 32
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
const Int_t PlMinVP          =  16;
const Int_t PlMaxVP          =  26;
const Int_t PlOffsBc         = 100;
const Int_t PlOffsSdcIn      =   0;
const Int_t PlOffsSdcOut     =  30;
const Int_t PlOffsTOF        =  50;
const Int_t PlOffsVP         =  15;

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

#endif
