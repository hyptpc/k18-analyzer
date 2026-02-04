// -*- C++ -*-

#ifndef DETECTOR_ID_HH
#define DETECTOR_ID_HH

#include <initializer_list>
#include <map>
#include <set>

#include <Rtypes.h>
#include <TString.h>

static const Char_t* const UorD[2] = {"U", "D"};

inline const std::map<TString, std::vector<TString>> DCNameList =
{
  {"BcIn", { "BLC1a", "BLC1b" }},
  {"BcOut", { "BLC2a", "BLC2b" }},
};

// __ Counters ___________________________________________________________
const Int_t DetIdBHT  =  1;
const Int_t DetIdT0   =  2;
const Int_t DetIdBH2  =  3;
const Int_t DetIdBAC  =  4;
const Int_t DetIdHTOF =  5;
const Int_t DetIdKVC  =  6;
const Int_t DetIdT1   =  7;
const Int_t DetIdCVC  =  8;
const Int_t DetIdSAC3 =  9;
const Int_t DetIdSFV  = 10;
const Int_t DetIdCOBO = 11;

const Int_t NumOfSegBHT  = 63;
const Int_t NumOfSegT0   =  5;
const Int_t NumOfSegBH2  = 15;
const Int_t NumOfSegBAC  =  5;
const Int_t NumOfSegHTOF = 34;
const Int_t NumOfSegKVC  =  8;
const Int_t NumOfSegT1   =  1;
const Int_t NumOfSegCVC  =  8;
const Int_t NumOfSegSAC3 =  1;
const Int_t NumOfSegSFV  =  1;
const Int_t NumOfSegCOBO =  8;

const Int_t DetIdVmeRm     =  81;
const Int_t DetIdScaler    =  91;
const Int_t DetIdTrigFlag  =  99;
const Int_t NumOfSegTrigFlag = 32;

enum EHodoscope {
  kBHT, kBH2, kBAC,
  kHTOF, kKVC, kT1,
  kCVC, kSAC3, kSFV,
  kCOBO, kNumHodo
};

const Int_t DetIdHodo[kNumHodo] = {
  DetIdBHT, DetIdBH2, DetIdBAC,
  DetIdHTOF, DetIdKVC, DetIdT1,
  DetIdCVC, DetIdSAC3, DetIdSFV,
  DetIdCOBO,
};

inline const std::vector<TString> NameHodo = {
  "BHT", "BH2", "BAC",
  "HTOF", "KVC", "T1",
  "CVC", "SAC3", "SFV",
  "COBO"
};

const Double_t NumOfSegHodo[kNumHodo] = {
  NumOfSegBHT, NumOfSegBH2, NumOfSegBAC,
  NumOfSegHTOF, NumOfSegKVC, NumOfSegT1,
  NumOfSegCVC, NumOfSegSAC3, NumOfSegSFV,
  NumOfSegCOBO
};

enum class HodoGroup : UInt_t { 
  None           = 0,
  NoADC          = 1u << 0,
  NoCluster      = 1u << 1,
  OneSideReadout = 1u << 2,
  Cherenkov      = 1u << 3,
  Ftof           = 1u << 4,
};

constexpr Bool_t HasHodoGroup(UInt_t mask, HodoGroup g)
{
  return (mask & static_cast<UInt_t>(g)) != 0;
}

constexpr UInt_t MakeHodoMask(std::initializer_list<HodoGroup> gArray)
{
  UInt_t m = 0;
  for (auto g : gArray)
    m |= static_cast<UInt_t>(g);
  return m;
}

inline constexpr UInt_t HodoGroupMask[kNumHodo] = {
  0, // kBHT
  0, // kBH2
  MakeHodoMask({HodoGroup::NoCluster, HodoGroup::OneSideReadout, HodoGroup::Cherenkov}), // kBAC
  0, // kHTOF
  MakeHodoMask({HodoGroup::NoCluster, HodoGroup::Cherenkov}), // kKVC
  MakeHodoMask({HodoGroup::NoCluster, HodoGroup::OneSideReadout}), // kT1
  MakeHodoMask({HodoGroup::Ftof}), // kCVC
  MakeHodoMask({HodoGroup::NoCluster, HodoGroup::OneSideReadout, HodoGroup::Cherenkov, HodoGroup::Ftof}), // kSAC3
  MakeHodoMask({HodoGroup::NoADC, HodoGroup::NoCluster, HodoGroup::OneSideReadout, HodoGroup::Ftof}), // kSFV
  MakeHodoMask({HodoGroup::NoADC, HodoGroup::NoCluster, HodoGroup::OneSideReadout}), // kCOBO
};


// __ Chambers ___________________________________________________________
const Int_t DetIdCDC    = 100;
const Int_t DetIdBLC1a  = 101;
const Int_t DetIdBLC1b  = 102;
const Int_t DetIdBLC2a  = 103;
const Int_t DetIdBLC2b  = 104;
const Int_t DetIdBPC    = 105;
const Int_t DetIdBPC2   = 105;
const Int_t DetIdBPC1   = 106;
const Int_t DetIdBPCmini= 106;
const Int_t DetIdSDC    = 106;
const Int_t DetIdFDC    = 107;
const Int_t DetIdBLC1   = 111;
const Int_t DetIdBLC2   = 112;
const Int_t DetIdBPC0   = 113;

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

const Int_t NumOfPlaneVmeRm=3;

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

// for compatibility
const Int_t LayerMinBcIn      =   1;
const Int_t LayerMaxBcIn      =  16;
const Int_t LayerMinBcOut     = 1001;
const Int_t LayerMaxBcOut     = 1016;
const Int_t LayerMinSdcIn     =   1;
const Int_t LayerMaxSdcIn     =  10;
const Int_t LayerMinSdcOut    =  31;
const Int_t LayerMaxSdcOut    =  42;
const Int_t LayerMinTOF       =  51; // need to change
const Int_t LayerMaxTOF       =  54; // need to change
const Int_t LayerMinVP        =  16;
const Int_t LayerMaxVP        =  26;
const Int_t PlOffsBc          = 100;
const Int_t PlOffsBcOut       = 1000;
const Int_t PlOffsSdcIn       =   0;
const Int_t PlOffsSdcOut      =  30;
const Int_t PlOffsTOF         =  50;
const Int_t PlOffsVP          =  15;
const Int_t PlOffsVPHS        = 207; //K1.8 Beam Tracking
const Int_t PlOffsTPCHit      = 700; //K1.8 w/ TPC Tracking

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

// __ TPC ___________________________________________________________
const Int_t NumOfLayersTPC    = 32;
const Int_t NumOfPadTPC       = 5768;
const Int_t NumOfAsadTPC      = 31;
const Int_t NumOfLayersVPTPC  = 5;
const Int_t NumOfLayersVPHS   = 4;

#endif
