// -*- C++ -*-

#ifndef TPC_RAW_DATA_HH
#define TPC_RAW_DATA_HH

#include <map>
#include <vector>

#include <TString.h>

#include "DetectorID.hh"


class TPCRawHit;

using TPCRHC = std::vector<TPCRawHit*>;

//_____________________________________________________________________________
class TPCRawData
{
public:
  static TString ClassName();
  TPCRawData();
  ~TPCRawData();

private:
  TPCRawData(const TPCRawData&);
  TPCRawData& operator=(const TPCRawData&);

private:
  template <typename T> using map_t = std::map<TString, T>;

  map_t<Bool_t>  m_is_decoded;
  map_t<TPCRHC>  m_tpc_raw_hit_collection;
  TPCRawHit*     m_baseline; //TPC FADC baseline

public:
  void           Clear(const TString& name="");
  Bool_t         DecodeHits(const TString& name="");
  Bool_t         DecodeTPCHits();
  //const TPCRHC&  GetTPCRawHitContainer(Int_t det_id) const;
  const TPCRHC&  GetTPCRawHitContainer(const TString& name) const; //All layers' TPC HC
  const TPCRHC&  GetTPCRawHitContainer(Int_t layer) const; //specific layer's TPC HC
  const TPCRHC&  GetTPCCorHitContainer(Int_t layer) const;
  Bool_t         CorrectBaselineTPC();
  const TPCRawHit* const   GetBaselineTPC() const { return m_baseline; }
  void           Print(Option_t* arg=nullptr) const;

  // aliases
  const TPCRHC&  GetTPCRawHits(Int_t layer) const
  { return GetTPCRawHitContainer(layer); }
  const TPCRHC&  GetTPCCorHits(Int_t layer) const
  { return GetTPCCorHitContainer(layer); }

  // templates
  template <typename T> Int_t GetEntries(const TString& name) const;
  template <typename T> const T* Get(const TString& name, Int_t i) const;

private:
  Bool_t AddTPCRawHit(const TString& name, Int_t plane, Int_t seg,
		      Int_t ch, Int_t data, Double_t val, Double_t* par=nullptr, Double_t raw_rms=0); //plane: layer, ch: row
};

//_____________________________________________________________________________
inline TString
TPCRawData::ClassName()
{
  static TString s_name("TPCRawData");
  return s_name;
}

//_____________________________________________________________________________
template <>
inline Int_t
TPCRawData::GetEntries<TPCRawHit>(const TString& name) const
{
  return m_tpc_raw_hit_collection.at(name).size();
}

//_____________________________________________________________________________
template <>
inline const TPCRawHit*
TPCRawData::Get<TPCRawHit>(const TString& name, Int_t i) const
{
  return m_tpc_raw_hit_collection.at(name).at(i);
}

#endif
