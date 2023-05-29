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

  HodoRawHitContainer              m_ScalerRawHC;

public:
  void                     Clear(const TString& name="");
  void                     ClearAll();
  Bool_t                   DecodeHits(const TString& name="");
  Bool_t                   DecodeCalibHits();
  const auto& GetHodoRawHitCollection() const
    { return m_hodo_raw_hit_collection; }
  const auto& GetDCRawHitCollection() const
    { return m_dc_raw_hit_collection; }
  const HodoRawHitContainer& GetHodoRawHitContainer(const TString& name) const;
  const DCRawHitContainer&   GetDCRawHitContainer(const TString& name) const;
  void        Print(Option_t* arg="") const;

private:
  enum EDCDataType { kDcLeading, kDcTrailing, kDcOverflow, kDcNDataType };
  void AddRawHit(const TString& name, Int_t plane, Int_t seg,
                 Int_t ch, Int_t data, Double_t val);
  Bool_t AddHodoRawHit(HodoRawHitContainer& cont,
                       const TString& name, Int_t plane, Int_t seg,
                       Int_t UorD, Int_t data, Double_t val);
  Bool_t AddDCRawHit(DCRawHitContainer& cont,
                       const TString& name, Int_t plane, Int_t seg,
                       Int_t UorD, Int_t data, Double_t val);
};

//_____________________________________________________________________________
inline TString
RawData::ClassName()
{
  static TString s_name("RawData");
  return s_name;
}

#endif
