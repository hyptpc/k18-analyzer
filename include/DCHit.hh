// -*- C++ -*-

#ifndef DC_HIT_HH
#define DC_HIT_HH

#include <cmath>
#include <deque>
#include <numeric>
#include <vector>

#include <TString.h>
#include <TVector3.h>

#include "std_ostream.hh"

class DCRawHit;
class DCLTrackHit;

//_____________________________________________________________________________
class DCHit
{
public:
  static const TString& ClassName();
  DCHit(const DCRawHit* rhit);
  DCHit(Int_t layer);
  DCHit(Int_t layer, Double_t wire);
  ~DCHit();

private:
  DCHit(const DCHit&);
  DCHit& operator =(const DCHit&);

protected:
  const DCRawHit* m_raw_hit;
  Int_t     m_plane; // by DIGIT
  Int_t     m_layer; // by DCGEO
  Double_t  m_wire;  // 0-based (= ch)

  using data_t = std::vector<Double_t>;
  using flag_t = std::deque<Bool_t>;

  data_t m_tdc;
  data_t m_adc;
  data_t m_trailing;

  // DCData normalized
  data_t m_drift_time;
  data_t m_drift_length;
  data_t m_tot;
  flag_t m_belong_to_track;
  flag_t m_is_good;

  Double_t  m_wpos;
  Double_t  m_angle;
  Double_t  m_z;

  std::vector<DCLTrackHit*> m_register_container;

public:
  const DCRawHit* GetRawHit() const { return m_raw_hit; }
  Int_t           PlaneId() const { return m_plane; }
  Int_t           LayerId() const { return m_layer; }
  Int_t           GetLayer() const { return LayerId(); }
  Double_t        WireId() const { return m_wire; }
  Double_t        GetWire() const { return m_wire; }
  Bool_t CalcDCObservables();
  Bool_t Calculate(){ return CalcDCObservables(); }

  void ClearDCData();
  void EraseDCData(Int_t i);

  Int_t GetTdcSize() const { return m_tdc.size(); }
  Int_t GetAdcSize() const { return m_adc.size(); }
  Int_t GetTdcVal(Int_t i=0) const { return m_tdc[i]; }
  Int_t GetAdcVal(Int_t i=0) const { return m_adc[i]; }
  Int_t GetTdcTrailing(Int_t i=0) const { return m_trailing[i]; }
  Int_t GetTdcTrailingSize() const { return m_trailing.size(); }
  Int_t GetTdc1st() const;
  Double_t GetTiltAngle() const { return m_angle; }
  Double_t GetWirePosition() const { return m_wpos; }
  Double_t GetZ() const { return m_z; }
  Double_t GetResolution() const;

  Int_t GetEntries() const { return m_drift_time.size(); }
  Double_t GetDriftTime(Int_t i) const { return m_drift_time.at(i); }
  Double_t GetDriftLength(Int_t i) const { return m_drift_length.at(i); }
  Double_t TimeOverThreshold(Int_t i=0) const { return m_tot.at(i); }
  void JoinTrack(Int_t i=0) { m_belong_to_track.at(i) = true; }
  void QuitTrack(Int_t i=0) { m_belong_to_track.at(i) = false; }
  Bool_t BelongToTrack(Int_t i=0) const { return m_belong_to_track.at(i); }
  Bool_t IsGood(Int_t i=0) const { return m_is_good.at(i); }
  Bool_t IsEmpty() const { return GetEntries() == 0; }

  void SetLayer(Int_t layer){ m_layer = layer; }
  void SetWire(Double_t wire){ m_wire  = wire; }
  void SetTdcVal(Int_t tdc){ m_tdc.push_back(tdc); }
  void SetTdc(Int_t tdc){ m_tdc.push_back(tdc); }
  void SetAdcVal(Int_t adc){ m_adc.push_back(adc); }
  void SetTrailing(Int_t tdc){ m_trailing.push_back(tdc); }
  void SetTdcTrailing(Int_t tdc){ m_trailing.push_back(tdc); }
  void SetDCData(Double_t dt=0., Double_t dl=0.,
                 Double_t tot=TMath::QuietNaN(),
                 Bool_t belong_to_track=false,
                 Bool_t is_good=true);
  void SetDriftLength(Double_t dl){ m_drift_length.push_back(dl); }
  void SetDriftTime(Double_t dt){ m_drift_time.push_back(dt); }
  void SetTiltAngle(Double_t angleDegree){ m_angle = angleDegree; }
  void SetWirePosition(Double_t wpos){ m_wpos = wpos; }

  // aliases
  Int_t GetTdc(Int_t i=0) const { return GetTdcVal(i); }
  Int_t Tdc(Int_t i=0) const { return GetTdcVal(i); }
  Int_t Tdc1st() const { return GetTdc1st(); }
  Int_t GetTrailing(Int_t i=0) const { return GetTdcTrailing(i); }
  Int_t Trailing(Int_t i=0) const { return GetTdcTrailing(i); }
  Double_t DT(Int_t i=0) const { return DriftTime(i); }
  Double_t DL(Int_t i=0) const { return DriftLength(i); }
  Double_t DriftTime(Int_t i=0) const { return GetDriftTime(i); }
  Double_t DriftLength(Int_t i=0) const { return GetDriftLength(i); }
  Double_t TOT(Int_t i=0) const { return TimeOverThreshold(i); }
  Double_t GetTot(Int_t i=0) const { return TimeOverThreshold(i); }
  Double_t TrailingTime(Int_t i=0) const { return DT(i) + TOT(i); }
  Double_t GetDriftTimeSize() const { return GetEntries(); }
  Double_t GetDriftLengthSize() const { return GetEntries(); }

  void DriftTimeCut(Double_t min, Double_t max, Bool_t select_1st);
  void TotCut(Double_t min, Bool_t keep_nan);
  void Print(Option_t* arg="") const;

  void RegisterHits(DCLTrackHit *hit)
    { m_register_container.push_back(hit); }

  Bool_t ReCalcDC(Bool_t applyRecursively=false) { return CalcDCObservables(); }
  static Bool_t Compare(const DCHit* left, const DCHit* right);

protected:
  void ClearRegisteredHits();
};

//_____________________________________________________________________
inline const TString&
DCHit::ClassName()
{
  static TString s_name("DCHit");
  return s_name;
}

//_____________________________________________________________________________
inline Bool_t
DCHit::Compare(const DCHit* left, const DCHit* right)
{
  if(left->LayerId() == right->LayerId())
    return left->WireId() < right->WireId();
  else
    return left->LayerId() < right->LayerId();
}

#endif
