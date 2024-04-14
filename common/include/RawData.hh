#ifndef RAW_DATA_HH
#define RAW_DATA_HH

#include "DetectorID.hh"
#include <vector>

class SDDRawHit;
class HodoRawHit;
class DCRawHit;
class MTDCRawHit;

typedef std::vector<MTDCRawHit*> MTDCRHitContainer;
typedef std::vector<SDDRawHit*>  SDDRHitContainer;
typedef std::vector<HodoRawHit*> HodoRHitContainer;
typedef std::vector<DCRawHit*>   DCRHitContainer;

//______________________________________________________________________________
class RawData
{
public:
  RawData( void );
  virtual ~RawData( void );

// private:
//   RawData( const RawData& );
//   RawData& operator=( const RawData& );

protected:
  bool              m_is_decoded;
  HodoRHitContainer m_BHDRawHC;
  HodoRHitContainer m_T0RawHC;
  HodoRHitContainer m_T1RawHC;
  HodoRHitContainer m_E0RawHC;
  HodoRHitContainer m_DEFRawHC;
  HodoRHitContainer m_ACRawHC;
  SDDRHitContainer  m_SDDRawHC;
  MTDCRHitContainer m_TrigRawHC;

  std::vector<DCRHitContainer> m_BLC1aRawHC;
  std::vector<DCRHitContainer> m_BLC1bRawHC;
  std::vector<DCRHitContainer> m_BLC2aRawHC;
  std::vector<DCRHitContainer> m_BLC2bRawHC;
  std::vector<DCRHitContainer> m_BPC2RawHC;
#ifdef CDS
  HodoRHitContainer m_CDHRawHC;
  std::vector<DCRHitContainer> m_CDCRawHC;
#endif
#ifdef E15
  std::vector<DCRHitContainer> m_FDCRawHC;
  HodoRHitContainer m_BVCRawHC;
  HodoRHitContainer m_CVCRawHC;
  HodoRHitContainer m_NCRawHC;
  HodoRHitContainer m_PCRawHC;
  HodoRHitContainer m_LBRawHC;
  HodoRHitContainer m_WVCRawHC;
  HodoRHitContainer m_BDRawHC;
  HodoRHitContainer m_BPDRawHC;
  HodoRHitContainer m_IHRawHC;
#elif E62
  HodoRHitContainer m_StartRawHC;
  HodoRHitContainer m_StopRawHC;
  std::vector<DCRHitContainer> m_SDCRawHC;
  std::vector<DCRHitContainer> m_FDCRawHC;
#elif E73
  HodoRHitContainer m_PbF2RawHC;
  HodoRHitContainer m_Veto0RawHC;
  HodoRHitContainer m_VetoRawHC;
  HodoRHitContainer m_BTCRawHC;
  HodoRHitContainer m_LeakRawHC;
  //  HodoRHitContainer m_FingerRawHC;
#elif E73_2024
  HodoRHitContainer m_PbGRawHC;
  HodoRHitContainer m_PbF2RawHC;
  HodoRHitContainer m_VetoRawHC;
  HodoRHitContainer m_BTCRawHC;
  HodoRHitContainer m_CVCRawHC;
  HodoRHitContainer m_NCRawHC;
  std::vector<DCRHitContainer> m_BPC1RawHC;
  std::vector<DCRHitContainer> m_VFTRawHC;
#elif T98
  HodoRHitContainer m_PbGRawHC;
  HodoRHitContainer m_PbF2RawHC;
  HodoRHitContainer m_VetoRawHC;
  HodoRHitContainer m_BTCRawHC;
  HodoRHitContainer m_RCRawHC;
  std::vector<DCRHitContainer> m_VFTRawHC;
#endif
  HodoRHitContainer m_VmeCalibRawHC;

public:
  virtual void             ClearAll( void );
  virtual bool             DecodeHits( void );
  virtual bool             DecodeHodoHits( void );
  virtual bool             DecodeDCHits( void );
  virtual bool             DecodeCDCHits( void );
  virtual bool             DecodeMTDCHits( void );
  virtual bool             DecodeCalibHits( void );
  virtual void             PrintHodo( const int &detid );

  SDDRHitContainer&  GetSDDRawHC(  const int &detid );
  HodoRHitContainer& GetHodoRawHC( const int &detid );
  MTDCRHitContainer& GetMTDCRawHC( const int &detid );
  DCRHitContainer&   GetDCRawHC(   const int &detid, int layer );
  HodoRHitContainer& GetVmeCalibRawHC( void );
};

//______________________________________________________________________________
inline SDDRHitContainer&
RawData::GetSDDRawHC( const int &detid )
{
  switch(detid){
  case DetIdSDD:
    return m_SDDRawHC;
  default:
    std::cout<<"E# invalid detector id "<< detid<<std::endl;
  }
  return m_SDDRawHC;
}
//______________________________________________________________________________
inline MTDCRHitContainer&
RawData::GetMTDCRawHC( const int &detid )
{
  switch(detid){
  case DetIdTrigFlag:
    return m_TrigRawHC;
  default:
    std::cout<<"E# invalid detector id "<< detid<<std::endl;
  }
  return m_TrigRawHC;
}
//______________________________________________________________________________
inline HodoRHitContainer&
RawData::GetHodoRawHC( const int &detid )
{
  switch(detid){
  case DetIdBHD:    return m_BHDRawHC;
  case DetIdT0:     return m_T0RawHC;
  case DetIdT1:     return m_T1RawHC;
  case DetIdE0:     return m_E0RawHC;
  case DetIdDEF:    return m_DEFRawHC;
  case DetIdAC:     return m_ACRawHC;
#ifdef E15
  case DetIdBVC:    return m_BVCRawHC;
  case DetIdCVC:    return m_CVCRawHC;
  case DetIdNC:     return m_NCRawHC;
  case DetIdPC:     return m_PCRawHC;
  case DetIdLB:     return m_LBRawHC;
  case DetIdWVC:    return m_WVCRawHC;
  case DetIdBD:     return m_BDRawHC;
  case DetIdBPD:    return m_BPDRawHC;
  case DetIdIH:     return m_IHRawHC;
#endif
#ifdef E62
  case DetIdStart:  return m_StartRawHC;
  case DetIdStop:   return m_StopRawHC;
#endif
#ifdef CDS
  case DetIdCDH:    return m_CDHRawHC;
#endif
#ifdef E73
  case DetIdPbF2:   return m_PbF2RawHC;
  case DetIdVeto:   return m_VetoRawHC;
  case DetIdVeto0:  return m_Veto0RawHC;
  case DetIdBTC:    return m_BTCRawHC;
  case DetIdLeak:   return m_LeakRawHC;
#endif
#ifdef T98
  case DetIdPbG:    return m_PbGRawHC;
  case DetIdPbF2:   return m_PbF2RawHC;
  case DetIdVeto:   return m_VetoRawHC;
  case DetIdBTC:    return m_BTCRawHC;
  case DetIdRC:     return m_RCRawHC;
#endif
#ifdef E73_2024
  case DetIdPbG:    return m_PbGRawHC;
  case DetIdPbF2:   return m_PbF2RawHC;
  case DetIdVeto:   return m_VetoRawHC;
  case DetIdBTC:    return m_BTCRawHC;
  case DetIdCVC:    return m_CVCRawHC;
  case DetIdNC:    return m_NCRawHC;
#endif

  default:
    std::cout<<"E# invalid detector id "<< detid<<std::endl;
  }
  return m_BHDRawHC;
}

//______________________________________________________________________________
inline DCRHitContainer&
RawData::GetDCRawHC( const int &detid, int layer ) 
{
  if( detid==DetIdCDC ){
    if( layer<0 || layer>118 ) layer = 0;
  }else if( detid==DetIdFDC ){
    if( layer<0 || layer>5 ) layer = 0;
  }else if( detid==DetIdVFT ){
    if( layer<0 || layer>14 ) layer = 0;
  }
  else if( layer<0 || layer>7 ) layer = 0;
  
  switch(detid){
#ifdef E15
  case DetIdFDC:    return m_FDCRawHC[layer];
#elif E62
  case DetIdSDC:    return m_SDCRawHC[layer];
  case DetIdFDC:    return m_FDCRawHC[layer];
#elif T98
  case DetIdVFT:    return m_VFTRawHC[layer];
#elif E73_2024
  case DetIdVFT:    return m_VFTRawHC[layer];
  case DetIdBPC1:    return m_BPC1RawHC[layer];
#endif
#ifdef CDS
  case DetIdCDC:    return m_CDCRawHC[layer];
#endif
  case DetIdBLC1a:  return m_BLC1aRawHC[layer];
  case DetIdBLC1b:  return m_BLC1bRawHC[layer];
  case DetIdBLC2a:  return m_BLC2aRawHC[layer];
  case DetIdBLC2b:  return m_BLC2bRawHC[layer];
  case DetIdBPC2:   return m_BPC2RawHC[layer];
  default:
    std::cout<<"E# invalid detector id "<< detid<<std::endl;
  }
  return m_BLC1aRawHC[0];
}

inline HodoRHitContainer&
RawData::GetVmeCalibRawHC( void ) 
{
  return m_VmeCalibRawHC;
}

#endif
