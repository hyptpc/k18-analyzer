// -*- C++ -*-

#ifndef DC_LTRACK_HIT_HH
#define DC_LTRACK_HIT_HH

#include <TString.h>

#include "DCHit.hh"
#include "DCGeomMan.hh"

class DCAnalyzer;

//_____________________________________________________________________________
class DCLTrackHit
{
public:
  static const TString& ClassName();
  DCLTrackHit(DCHit* hit, Double_t pos, Int_t nh);
  DCLTrackHit(const DCLTrackHit& right);
  ~DCLTrackHit();

private:
  DCHit*   m_hit;
  Int_t    m_nth_hit;
  Double_t m_local_hit_pos;
  Double_t m_cal_pos;
  Double_t m_xcal;
  Double_t m_ycal;
  Double_t m_ucal;
  Double_t m_vcal;
  Bool_t   m_honeycomb;

public:
  void     SetLocalHitPos(Double_t xl) { m_local_hit_pos = xl; }
  void     SetCalPosition(Double_t x, Double_t y) { m_xcal = x; m_ycal = y; }
  void     SetCalUV(Double_t u, Double_t v) { m_ucal = u; m_vcal = v; }
  void     SetHoneycomb(Bool_t flag=true) { m_honeycomb = flag; }
  Bool_t   IsHoneycomb() const { return m_honeycomb; }
  Int_t    GetIndex() const { return m_nth_hit; }
  Int_t    GetLayer() const { return m_hit->GetLayer(); }
  Int_t    LayerId() const { return m_hit->LayerId(); }
  Double_t GetWire() const { return m_hit->GetWire(); }
  Int_t    GetTdcVal() const { return m_hit->GetTdcVal(m_nth_hit); }
  Int_t    GetTdcSize() const { return m_hit->GetTdcSize(); }
  Double_t GetDriftTime() const { return m_hit->DT(m_nth_hit); }
  Double_t GetDriftLength() const { return m_hit->DL(m_nth_hit); }
  Double_t GetTrailingTime() const
  { return m_hit->TrailingTime(m_nth_hit); }
  Double_t TimeOverThreshold() const { return m_hit->TOT(m_nth_hit); }
  DCHit*   GetHit() const { return m_hit; }
  Double_t GetTiltAngle() const { return m_hit->GetTiltAngle(); }
  Double_t GetWirePosition() const { return m_hit->GetWirePosition(); }
  Double_t GetLocalHitPos() const { return m_local_hit_pos; }
  Double_t GetLocalCalPos() const;
  Double_t GetXcal() const { return m_xcal; }
  Double_t GetYcal() const { return m_ycal; }
  Double_t GetUcal() const { return m_ucal; }
  Double_t GetVcal() const { return m_vcal; }
  Double_t GetResidual() const;
  Double_t GetResolution() const { return m_hit->GetResolution(); }

  // aliases
  Double_t DriftTime() const { return GetDriftTime(); }
  Double_t DriftLength() const { return GetDriftLength(); }
  Double_t TOT() const { return TimeOverThreshold(); }
  Double_t GetTot() const { return TimeOverThreshold(); }

  ///// for XUV tracking
  void     SetLocalCalPosVXU(Double_t xcl) { m_cal_pos=xcl; }
  Double_t GetLocalCalPosVXU() const { return m_cal_pos; }
  Double_t GetResidualVXU() const { return m_local_hit_pos-m_cal_pos; }

  Double_t GetZ() const { return m_hit->GetZ(); }

  void   JoinTrack() { m_hit->JoinTrack(m_nth_hit); }
  void   QuitTrack() { m_hit->QuitTrack(m_nth_hit); }
  Bool_t BelongToTrack() const { return m_hit->BelongToTrack(m_nth_hit); }

  void Print(const TString& arg="") const;

  Bool_t ReCalc(Bool_t applyRecursively=false);

  friend class DCHit;
};

//_____________________________________________________________________________
inline const TString&
DCLTrackHit::ClassName()
{
  static TString s_name("DCLTrackHit");
  return s_name;
}

#endif
