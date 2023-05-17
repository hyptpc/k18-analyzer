// -*- C++ -*-

#ifndef RAW_DATA_HH
#define RAW_DATA_HH

#include <deque>
#include <map>
#include <vector>

#include <TString.h>

#include "DetectorID.hh"

class HodoRawHit;
class DCRawHit;

using HodoRawHitContainer = std::vector<HodoRawHit*>;
using DCRawHitContainer = std::vector<DCRawHit*>;

//_____________________________________________________________________________
class RawData
{
public:
  static TString ClassName();
  RawData();
  ~RawData();

private:
  RawData(const RawData&);
  RawData& operator=(const RawData&);

private:
  template <typename T> using map_t = std::map<TString, T>;

  map_t<Bool_t>              m_is_decoded;
  map_t<HodoRawHitContainer> m_hodo_raw_hit_collection;
  map_t<DCRawHitContainer>   m_dc_raw_hit_collection;

  std::vector<HodoRawHitContainer> m_BFTRawHC;
  std::vector<DCRawHitContainer>   m_BcInRawHC;
  std::vector<DCRawHitContainer>   m_BcOutRawHC;
  std::vector<DCRawHitContainer>   m_SdcInRawHC;
  std::vector<DCRawHitContainer>   m_SdcOutRawHC;
  HodoRawHitContainer              m_ScalerRawHC;

public:
  void                     Clear(const TString& name="");
  void                     ClearAll();
  Bool_t                   DecodeHits(const TString& name="");
  Bool_t                   DecodeCalibHits();
  const map_t<HodoRawHitContainer>& GetHodoRawHitCollection() const
    { return m_hodo_raw_hit_collection; }
  const HodoRawHitContainer& GetHodoRawHitContainer(const TString& name) const;
  const DCRawHitContainer&   GetBcInRawHC(Int_t layer) const;
  const DCRawHitContainer&   GetBcOutRawHC(Int_t layer) const;
  const DCRawHitContainer&   GetSdcInRawHC(Int_t layer) const;
  const DCRawHitContainer&   GetSdcOutRawHC(Int_t layer) const;

private:
  enum EDCDataType { kDcLeading, kDcTrailing, kDcOverflow, kDcNDataType };
  void AddRawHit(const TString& name, Int_t plane, Int_t seg,
                 Int_t ch, Int_t data, Double_t val);
  Bool_t AddHodoRawHit(HodoRawHitContainer& cont,
                       const TString& name, Int_t plane, Int_t seg,
                       Int_t UorD, Int_t data, Double_t val);
  Bool_t AddDCRawHit(DCRawHitContainer& cont,
                     Int_t plane, Int_t wire, Int_t data,
                     Int_t type=kDcLeading);
};

//_____________________________________________________________________________
inline TString
RawData::ClassName()
{
  static TString s_name("RawData");
  return s_name;
}

//_____________________________________________________________________________
inline const DCRawHitContainer&
RawData::GetBcInRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersBcIn) layer = 0;
  return m_BcInRawHC[layer];
}

//_____________________________________________________________________________
inline const DCRawHitContainer&
RawData::GetBcOutRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersBcOut) layer = 0;
  return m_BcOutRawHC[layer];
}

//_____________________________________________________________________________
inline const DCRawHitContainer&
RawData::GetSdcInRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersSdcIn) layer = 0;
  return m_SdcInRawHC[layer];
}

//_____________________________________________________________________________
inline const DCRawHitContainer&
RawData::GetSdcOutRawHC(Int_t layer) const
{
  if(layer<0 || layer>NumOfLayersSdcOut) layer = 0;
  return m_SdcOutRawHC[layer];
}

#endif
