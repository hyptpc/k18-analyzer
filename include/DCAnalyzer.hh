// -*- C++ -*-

#ifndef DC_ANALYZER_HH
#define DC_ANALYZER_HH

#include <map>
#include <vector>

#include <TString.h>

#include "DetectorID.hh"

class DCHit;
class DCLocalTrack;
class RawData;

using DCHC = std::vector<DCHit*>;
using DCLocalTC = std::vector<DCLocalTrack*>;

//_____________________________________________________________________________
class DCAnalyzer
{
public:
  static const TString& ClassName();
  DCAnalyzer();
  DCAnalyzer(const RawData& raw_data);
  ~DCAnalyzer();

private:
  DCAnalyzer(const DCAnalyzer&);
  DCAnalyzer& operator =(const DCAnalyzer&);

private:
  template <typename T> using map_t = std::map<TString, T>;
  const RawData*         m_raw_data;
  map_t<DCHC>            m_dc_hit_collection;
  std::vector<DCHC>      m_TempBcInHC;
  std::vector<DCHC>      m_BcInHC;
  std::vector<DCHC>      m_BcOutHC;
  DCLocalTC              m_BcInTC;
  DCLocalTC              m_BcOutTC;

public:
  // __ Hit Decoding ________________________________________________________
  void   DecodeHits(const TString& name);
  Bool_t DecodeRawHits();
  Bool_t DecodeBcInHits();
  Bool_t DecodeBcOutHits();

  // __ Hit Container Access ________________________________________________
  const DCHC& GetTempBcInHC(Int_t l) const { return m_TempBcInHC.at(l); }
  const DCHC& GetBcInHC(Int_t l)     const { return m_BcInHC.at(l); }
  const DCHC& GetBcOutHC(Int_t l)    const { return m_BcOutHC.at(l); }

  // __ Track Searching _____________________________________________________
  Bool_t TrackSearchBcIn(Int_t T0Seg=-1);
  Bool_t TrackSearchBcOut(Int_t T0Seg=-1);

  // __ Track Container Access ______________________________________________
  const DCLocalTC& GetBcInTrackContainer()  const { return m_BcInTC; }
  const DCLocalTC& GetBcOutTrackContainer() const { return m_BcOutTC; }

  Int_t GetNtracksBcIn() const { return m_BcInTC.size(); }
  Int_t GetNtracksBcOut() const { return m_BcOutTC.size(); }

  const DCLocalTrack* GetTrackBcIn(Int_t l) const { return m_BcInTC.at(l); }
  const DCLocalTrack* GetTrackBcOut(Int_t l) const { return m_BcOutTC.at(l); }

  // __ Cuts and Filters ____________________________________________________
  void ChiSqrCutBcOut(Double_t chisqr);
  
  void TotCutBCOut(Double_t min_tot);
  void TotCut(const TString& name);
  
  void DriftTimeCut(const TString& name);

  // __ Recalculation _______________________________________________________
  Bool_t ReCalcDCHits(std::vector<DCHC>& cont,
                      Bool_t applyRecursively=false);
  Bool_t ReCalcDCHits(Bool_t applyRecursively=false);
  Bool_t ReCalcTrackBcIn(Bool_t applyRecursively=false);
  Bool_t ReCalcTrackBcOut(Bool_t applyRecursively=false);
  Bool_t ReCalcAll();

  // __ Clear and Reset _____________________________________________________
  void ResetTracksBcIn()  { ClearTracksBcIn(); }
  void ResetTracksBcOut() { ClearTracksBcOut();}

private:
  // __ Internal Helpers ____________________________________________________
  void ChiSqrCut(DCLocalTC& cont, Double_t chisqr);
  void EraseEmptyHits(const TString& name);
  void EraseEmptyHits(std::vector<DCHC>& HitCont);
  void TotCut(const TString& name, Double_t min_tot, Bool_t keep_nan=true);
  void DriftTimeCut(const TString& name, Double_t min_dt, Double_t max_dt, Bool_t select_1st);
  Bool_t ReCalcTrack(DCLocalTC& cont, Bool_t applyRecursively=false);

  // __ Clear (Internal) ____________________________________________________
  void ClearDCHits();
  void ClearBcInHits();
  void ClearBcOutHits();

  void ClearTracksBcIn();
  void ClearTracksBcOut();
};



//_____________________________________________________________________________
inline const TString&
DCAnalyzer::ClassName()
{
  static TString s_name("DCAnalyzer");
  return s_name;
}

#endif
