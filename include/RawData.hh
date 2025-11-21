// -*- C++ -*-

#ifndef RAW_DATA_HH
#define RAW_DATA_HH

#include <map>
#include <vector>

#include <TString.h>

#include "DetectorID.hh"

class HodoRawHit;
class DCRawHit;
class TPCRawHit;

using HodoRHC = std::vector<HodoRawHit*>;
using DCRHC = std::vector<DCRawHit*>;
using TPCRHC = std::vector<TPCRawHit*>;

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

  map_t<Bool_t>  m_is_decoded;
  map_t<HodoRHC> m_hodo_raw_hit_collection;
  map_t<DCRHC>   m_dc_raw_hit_collection;
  map_t<TPCRHC>  m_tpc_raw_hit_collection;
  TPCRawHit*     m_baseline; //TPC FADC baseline

public:
  void           Clear(const TString& name="");
  Bool_t         DecodeHits(const TString& name="");
  Bool_t         DecodeTPCHits();
  const HodoRHC& GetHodoRawHitContainer(const TString& name) const;
  const HodoRHC& GetHodoRawHitContainer(Int_t det_id) const;
  const DCRHC&   GetDCRawHitContainer(const TString& name) const;
  const DCRHC&   GetDCRawHitContainer(Int_t det_id) const;
  const DCRHC&   GetDCRawHitContainer(Int_t det_id, Int_t plane) const;
  //const TPCRHC&  GetTPCRawHitContainer(Int_t det_id) const;
  const TPCRHC&  GetTPCRawHitContainer(const TString& name) const; //All layers' TPC HC
  const TPCRHC&  GetTPCRawHitContainer(Int_t layer) const; //specific layer's TPC HC
  const TPCRHC&  GetTPCCorHitContainer(Int_t layer) const;
  Bool_t         CorrectBaselineTPC();
  const TPCRawHit* const   GetBaselineTPC() const { return m_baseline; }
  void           Print(Option_t* arg=nullptr) const;

  // aliases
  const HodoRHC& GetHodoRawHC(const TString& name) const
  { return GetHodoRawHitContainer(name); }
  const HodoRHC& GetHodoRawHC(Int_t det_id) const
  { return GetHodoRawHitContainer(det_id); }
  const DCRHC&   GetDCRawHC(const TString& name) const
  { return GetDCRawHitContainer(name); }
  const DCRHC&   GetDCRawHC(Int_t det_id) const
  { return GetDCRawHitContainer(det_id); }
  const DCRHC&   GetDCRawHC(Int_t det_id, Int_t plane) const
  { return GetDCRawHitContainer(det_id, plane); }
  const HodoRHC& GetHodoRawHits(const TString& name) const
  { return GetHodoRawHitContainer(name); }
  const DCRHC&   GetDCRawHits(const TString& name) const
  { return GetDCRawHitContainer(name); }
  const TPCRHC&  GetTPCRawHits(Int_t layer) const
  { return GetTPCRawHitContainer(layer); }
  const TPCRHC&  GetTPCCorHits(Int_t layer) const
  { return GetTPCCorHitContainer(layer); }

  // templates
  template <typename T> Int_t GetEntries(const TString& name) const;
  template <typename T> const T* Get(const TString& name, Int_t i) const;

private:
  Bool_t AddHodoRawHit(const TString& name, Int_t plane, Int_t seg,
                       Int_t UorD, Int_t data, Double_t val);
  Bool_t AddFiberRawHit(const TString& name, Int_t plane, Int_t seg,
                        Int_t UorD, Int_t data, Double_t val);
  Bool_t AddDCRawHit(const TString& name, Int_t plane, Int_t seg,
		     Int_t ch, Int_t data, Double_t val);
		     //Int_t UorD, Int_t data, Double_t val);
  Bool_t AddTPCRawHit(const TString& name, Int_t plane, Int_t seg,
		      Int_t ch, Int_t data, Double_t val, Double_t* par=nullptr, Double_t raw_rms=0); //plane: layer, ch: row

};

//_____________________________________________________________________________
inline TString
RawData::ClassName()
{
  static TString s_name("RawData");
  return s_name;
}

//_____________________________________________________________________________
template <>
inline Int_t
RawData::GetEntries<HodoRawHit>(const TString& name) const
{
  return m_hodo_raw_hit_collection.at(name).size();
}

//_____________________________________________________________________________
template <>
inline Int_t
RawData::GetEntries<DCRawHit>(const TString& name) const
{
  return m_dc_raw_hit_collection.at(name).size();
}

//_____________________________________________________________________________
template <>
inline Int_t
RawData::GetEntries<TPCRawHit>(const TString& name) const
{
  return m_tpc_raw_hit_collection.at(name).size();
}

//_____________________________________________________________________________
template <>
inline const HodoRawHit*
RawData::Get<HodoRawHit>(const TString& name, Int_t i) const
{
  return m_hodo_raw_hit_collection.at(name).at(i);
}

//_____________________________________________________________________________
template <>
inline const DCRawHit*
RawData::Get<DCRawHit>(const TString& name, Int_t i) const
{
  return m_dc_raw_hit_collection.at(name).at(i);
}

//_____________________________________________________________________________
template <>
inline const TPCRawHit*
RawData::Get<TPCRawHit>(const TString& name, Int_t i) const
{
  return m_tpc_raw_hit_collection.at(name).at(i);
}

#endif
