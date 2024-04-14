#ifndef SETUP_HH
#define SETUP_HH 1
//#include "TString.h"

namespace e73_2024
{
enum param_hodo{kT1,kBHT,kT0,kT0new,kDEF,kVeto,kBTC,kCDH,kPbG,kPbF2,kCVC,kNC,kRC,kNumHodo};
  const int     hodoid[kNumHodo]={DetIdT1,DetIdBHT,DetIdT0,DetIdDEF,DetIdVeto,DetIdBTC,DetIdCDH,DetIdPbG,DetIdPbF2,DetIdCVC,DetIdNC};
  const TString hodoname[kNumHodo]={"T1","BHT","T0","DEF","Veto","BTC","CDH","PbG","PbF2","CVC","NC"};
  const double  nsegs[kNumHodo]={1,63,5,4,4,2,36,40,36,10,6};

  enum param_chm{kBLC1a,kBLC1b,kBLC2a,kBLC2b,kBPC1,kBPC2,kVFT,kCDC,kBLC1,kBLC2,kBPC0,kNumChm};
  const int     chmid[kNumChm]={DetIdBLC1a,DetIdBLC1b,DetIdBLC2a,DetIdBLC2b,DetIdBPC1,DetIdBPC2,DetIdVFT,DetIdCDC,DetIdBLC1,DetIdBLC2,DetIdBPC0};
  const TString chmname[kNumChm]={"BLC1a", "BLC1b", "BLC2a", "BLC2b","BPC1","BPC2","VFT","CDC","BLC1","BLC2","BPC0"};
  const int  nlayers[kNumChm]={8,8,8,8,8,8,14,118,16,16,16};
  const int  nwires[kNumChm]={32,32,32,32,15,32,64,16,32,32,32};

  enum param_flag{kStart,kStop,kBeam,kPion,kKaon2,kKaon3, //5
		  kKCDH1,kKCDH2,kKCDH3, //8
		  kKCDH1Gamma,kKGamma,kPiCDH, //11
		  kBeamPbF2,kPiPbF2=kBeamPbF2, //12
		  kCDHCosmic,kElectronPbF2=kCDHCosmic, //13
		  kPbF2Cosmic,kCosmic=kPbF2Cosmic, //14
		  kProtonTrig,kDeuteron=kProtonTrig, //15
		  kNumTrig};
};
#endif
