// -*- C++ -*-

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"
#include "BHTHit.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"

using Hodo1HitContainer = std::vector<Hodo1Hit*>;
using Hodo2HitContainer = std::vector<Hodo2Hit*>;
using BHTHitContainer = std::vector<BHTHit*>;

//_____________________________________________________________________________
class HodoAnalyzer
{
public:
  HodoAnalyzer(const RawData& rawData);
  ~HodoAnalyzer();

  static const TString& ClassName();

private:
  HodoAnalyzer(const HodoAnalyzer&);
  HodoAnalyzer& operator =(const HodoAnalyzer&);

private:
  const RawData* m_raw_data;

#if T98
  BHTHitContainer         m_BHDCont;
#elif E73_2024
  BHTHitContainer         m_BHDCont;
#else
  Hodo2HitContainer       m_BHDCont;
#endif
  Hodo2HitContainer       m_T0Cont;
  Hodo2HitContainer       m_T0newCont;
  Hodo2HitContainer       m_E0Cont;
  Hodo2HitContainer       m_DEFCont;
  Hodo2HitContainer       m_StartCont;
  Hodo2HitContainer       m_StopCont;
  Hodo2HitContainer       m_CDHCont;
#if E73
  Hodo2HitContainer       m_PbF2Cont;
  Hodo2HitContainer       m_Veto1Cont;
  Hodo2HitContainer       m_Veto0Cont;
  Hodo2HitContainer       m_BTCCont;
#elif E73_2024
  Hodo2HitContainer       m_PbGCont;
  Hodo2HitContainer       m_PbF2Cont;
  Hodo2HitContainer       m_VetoCont;
  Hodo2HitContainer       m_BTCCont;
#elif T98
  Hodo2HitContainer       m_PbGCont;
  Hodo2HitContainer       m_PbF2Cont;
  Hodo2HitContainer       m_VetoCont;
  Hodo2HitContainer       m_BTCCont;
  Hodo2HitContainer       m_RCCont;
#elif E15
  Hodo2HitContainer       m_BVCCont;
  Hodo2HitContainer       m_CVCCont;
  Hodo2HitContainer       m_NCCont;
  Hodo2HitContainer       m_PCCont;
  Hodo2HitContainer       m_LBCont;
  Hodo2HitContainer       m_WVCCont;
  Hodo2HitContainer       m_BDCont;
  Hodo2HitContainer       m_BPDCont;
  Hodo2HitContainer       m_IHCont;
#endif

public:
  Bool_t DecodeRawHits();
  Bool_t DecodeHodoHits(const Int_t &detid, Hodo2HitContainer &m_Cont);
  Bool_t DecodeBHTHits(const Int_t &detid, BHTHitContainer &m_Cont);
  Bool_t DecodeHodo1Hits(const Int_t &detid, Hodo2HitContainer &m_Cont);

  inline Int_t  GetNHits(Int_t detID) const;
  inline Bool_t AddHit(Int_t detID, Hodo2Hit* hp);
  inline Bool_t AddHit(Int_t detID, BHTHit* hp);

  inline Hodo1Hit* Get1Hit(Int_t detID, Int_t i) const;
  inline Hodo2Hit* GetHit(Int_t detID, Int_t i) const;
  Bool_t ReCalcAll();

};

//_____________________________________________________________________________
inline const TString&
HodoAnalyzer::ClassName()
{
  static const TString s_name("HodoAnalyzer");
  return s_name;
}

//_____________________________________________________________________________
inline Int_t
HodoAnalyzer::GetNHits(Int_t detID) const
{
  switch(detID){
  case DetIdBHD:    return m_BHDCont.size();
  case DetIdT0:     return m_T0Cont.size();
  case DetIdT0new:  return m_T0newCont.size();
  case DetIdE0:     return m_E0Cont.size();
  case DetIdDEF:    return m_DEFCont.size();
  case DetIdStart:  return m_StartCont.size();
  case DetIdStop:   return m_StopCont.size();
  case DetIdCDH:    return m_CDHCont.size();
#if E15
  case DetIdBVC:    return m_BVCCont.size();
  case DetIdCVC:    return m_CVCCont.size();
  case DetIdNC:     return m_NCCont.size();
  case DetIdPC:     return m_PCCont.size();
  case DetIdLB:     return m_LBCont.size();
  case DetIdWVC:    return m_WVCCont.size();
  case DetIdBD:     return m_BDCont.size();
  case DetIdBPD:    return m_BPDCont.size();
  case DetIdIH:     return m_IHCont.size();
#endif
#if E73
  case DetIdPbF2:   return m_PbF2Cont.size();
  case DetIdVeto1:  return m_Veto1Cont.size();
  case DetIdVeto0:  return m_Veto0Cont.size();
  case DetIdBTC:    return m_BTCCont.size();
#endif
#if E73_2024
  case DetIdPbG:    return m_PbGCont.size();
  case DetIdPbF2:   return m_PbF2Cont.size();
  case DetIdVeto:   return m_VetoCont.size();
  case DetIdBTC:    return m_BTCCont.size();
#endif
#if T98
  case DetIdPbG:    return m_PbGCont.size();
  case DetIdPbF2:   return m_PbF2Cont.size();
  case DetIdVeto:   return m_VetoCont.size();
  case DetIdBTC:    return m_BTCCont.size();
  case DetIdRC:     return m_RCCont.size();
#endif
  default:
    return 0;
  }
}

//_____________________________________________________________________________
inline Bool_t
HodoAnalyzer::AddHit(Int_t detID, Hodo2Hit* hp)
{
  switch(detID){
#ifndef T98
#ifndef E73_2024
  case DetIdBHD:    m_BHDCont.push_back(hp);   return true;
#endif
#endif
  case DetIdT0:     m_T0Cont.push_back(hp);    return true;
  case DetIdT0new:  m_T0newCont.push_back(hp); return true;
  case DetIdE0:     m_E0Cont.push_back(hp);    return true;
  case DetIdDEF:    m_DEFCont.push_back(hp);   return true;
  case DetIdStart:  m_StartCont.push_back(hp); return true;
  case DetIdStop:   m_StopCont.push_back(hp);  return true;
  case DetIdCDH:    m_CDHCont.push_back(hp);   return true;
#if E15
  case DetIdBVC:    m_BVCCont.push_back(hp);   return true;
  case DetIdCVC:    m_CVCCont.push_back(hp);   return true;
  case DetIdNC:     m_NCCont.push_back(hp);    return true;
  case DetIdPC:     m_PCCont.push_back(hp);    return true;
  case DetIdLB:     m_LBCont.push_back(hp);    return true;
  case DetIdWVC:    m_WVCCont.push_back(hp);   return true;
  case DetIdBD:     m_BDCont.push_back(hp);    return true;
  case DetIdBPD:    m_BPDCont.push_back(hp);   return true;
  case DetIdIH:     m_IHCont.push_back(hp);    return true;
#endif
#if E73
  case DetIdPbF2:   m_PbF2Cont.push_back(hp);  return true;
  case DetIdVeto1:  m_Veto1Cont.push_back(hp); return true;
  case DetIdVeto0:  m_Veto0Cont.push_back(hp); return true;
  case DetIdBTC:    m_BTCCont.push_back(hp);   return true;
#endif
#if T98
  case DetIdPbG:    m_PbGCont.push_back(hp);   return true;
  case DetIdPbF2:   m_PbF2Cont.push_back(hp);  return true;
  case DetIdVeto:   m_VetoCont.push_back(hp);  return true;
  case DetIdBTC:    m_BTCCont.push_back(hp);   return true;
  case DetIdRC:     m_RCCont.push_back(hp);    return true;
#endif
#if E73_2024
  case DetIdPbG:    m_PbGCont.push_back(hp);   return true;
  case DetIdPbF2:   m_PbF2Cont.push_back(hp);  return true;
  case DetIdVeto:   m_VetoCont.push_back(hp);  return true;
  case DetIdBTC:    m_BTCCont.push_back(hp);   return true;
#endif
  default:
    return false;
  }
}

//_____________________________________________________________________________
inline Bool_t
HodoAnalyzer::AddHit(Int_t detID, BHTHit* hp)
{
  switch(detID){
#if T98
  case DetIdBHD:    m_BHDCont.push_back(hp);   return true;
#elif E73_2024
  case DetIdBHD:    m_BHDCont.push_back(hp);   return true;
#endif
  default:
    return false;
  }
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::Get1Hit(Int_t detID, Int_t i) const
{
  if(GetNHits(detID)<=i) return 0;
  switch(detID){
#if E15
  case DetIdWVC:    return (Hodo1Hit*)m_WVCCont.at(i);
  case DetIdIH:     return (Hodo1Hit*)m_IHCont.at(i);
#endif
#if E73
  case DetIdPbF2:   return (Hodo1Hit*)m_PbF2Cont.at(i);
  case DetIdVeto0:  return (Hodo1Hit*)m_Veto0Cont.at(i);
  case DetIdBTC:    return (Hodo1Hit*)m_BTCCont.at(i);
#endif
#if T98
  case DetIdPbG:    return (Hodo1Hit*)m_PbGCont.at(i);
  case DetIdPbF2:   return (Hodo1Hit*)m_PbF2Cont.at(i);
#endif
#if E73_2024
  case DetIdPbG:    return (Hodo1Hit*)m_PbGCont.at(i);
  case DetIdPbF2:   return (Hodo1Hit*)m_PbF2Cont.at(i);
#endif
  default:
    return 0;
  }
}

//_____________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHit(Int_t detID, Int_t i) const
{
  if(GetNHits(detID)<=i) return 0;
  switch(detID){
  case DetIdBHD:    return m_BHDCont.at(i);
  case DetIdT0:     return m_T0Cont.at(i);
  case DetIdT0new:  return m_T0newCont.at(i);
  case DetIdE0:     return m_E0Cont.at(i);
  case DetIdDEF:    return m_DEFCont.at(i);
  case DetIdStart:  return m_StartCont.at(i);
  case DetIdStop:   return m_StopCont.at(i);
  case DetIdCDH:    return m_CDHCont.at(i);
#if E15
  case DetIdBVC:    return m_BVCCont.at(i);
  case DetIdCVC:    return m_CVCCont.at(i);
  case DetIdNC:     return m_NCCont.at(i);
  case DetIdPC:     return m_PCCont.at(i);
  case DetIdLB:     return m_LBCont.at(i);
  case DetIdBD:     return m_BDCont.at(i);
  case DetIdBPD:    return m_BPDCont.at(i);
#endif
#if E73
  case DetIdVeto1:  return m_Veto1Cont.at(i);
#endif
#if T98
  case DetIdVeto:   return m_VetoCont.at(i);
  case DetIdBTC:    return m_BTCCont.at(i);
  case DetIdRC:     return m_RCCont.at(i);
#endif
#if E73_2024
  case DetIdVeto:   return m_VetoCont.at(i);
  case DetIdBTC:    return m_BTCCont.at(i);
#endif
  default:
    return 0;
  }
}

#endif
