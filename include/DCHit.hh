// -*- C++ -*-

#ifndef DC_HIT_HH
#define DC_HIT_HH

#include <cmath>
#include <vector>
#include <deque>
#include <numeric>

#include <TString.h>
#include <TVector3.h>

#include <std_ostream.hh>

#include "DebugCounter.hh"

using BVec_t = std::deque<Bool_t>;
using IVec_t = std::vector<Int_t>;
using DVec_t = std::vector<Double_t>;

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
  Double_t  m_wire;  // 1-based (ch + 1)
  IVec_t    m_tdc;
  IVec_t    m_adc;
  IVec_t    m_trailing;

  // For DC with HUL MH-TDC
  struct data_t
  {
    Double_t drift_time;
    Double_t drift_length;
    Double_t tot; // [ch]
    Bool_t   belong_to_track;
    Bool_t   dl_is_good;
  };
  std::vector<data_t> m_data;

  Double_t  m_wpos;
  Double_t  m_angle;
  //  BoolVec m_belong_track;
  //  BoolVec m_dl_range;

  ///// for MWPC
  Int_t    m_cluster_size;
  Bool_t   m_mwpc_flag;
  Double_t m_mwpc_wire;
  Double_t m_mwpc_wpos;

  ///// for TOF
  Double_t m_z;

  ///// For E40 Acrylic TOF
  Double_t m_ofs_dt;

  ///// for CFT
  Double_t m_meanseg;
  Double_t m_maxseg;
  Double_t m_adc_low;
  Double_t m_mip_low;
  Double_t m_dE_low;
  Double_t max_adc_low;
  Double_t max_mip_low;
  Double_t max_dE_low;
  Double_t m_r;
  Double_t m_phi;
  // BVec_t   m_belong_track;
  TVector3 m_vtx;
  Double_t m_pos_phi;
  Double_t m_pos_z;
  Double_t m_pos_r;
  Double_t m_time;

  ///// for TPC
  Int_t    m_hitnum;

  std::vector<DCLTrackHit*> m_register_container;

public:
  const DCRawHit* GetRawHit() const { return m_raw_hit; }
  Int_t           PlaneId() const { return m_plane; }
  Int_t           LayerId() const { return m_layer; }
  Int_t           GetLayer() const { return LayerId(); }
  Double_t        WireId() const { return m_wire; }
  Double_t        GetWire() const {
    if(m_mwpc_flag) return m_mwpc_wire;
    else return m_wire;
  }
  Bool_t Calculate();
  Bool_t CalcDCObservables();
  // Bool_t CalcMWPCObservables();
  Bool_t CalcFiberObservables();
  Bool_t CalcCFTObservables();
  //  Bool_t CalcObservablesSimulation(Double_t dlength);

  void SetLayer(Int_t layer){ m_layer = layer; }
  void SetWire(Double_t wire){ m_wire  = wire; }
  void SetTdcVal(Int_t tdc){ m_tdc.push_back(tdc); }
  void SetTdc(Int_t tdc){ m_tdc.push_back(tdc); }
  void SetAdcVal(Int_t adc){ m_adc.push_back(adc); }
  void SetTrailing(Int_t tdc){ m_trailing.push_back(tdc); }
  void SetTdcTrailing(Int_t tdc){ m_trailing.push_back(tdc); }
  void SetDummyData();
  void SetDriftLength(Int_t i, Double_t dl){ m_data.at(i).drift_length = dl; }
  void SetDriftTime(Int_t i, Double_t dt){ m_data.at(i).drift_time = dt; }
  void SetTiltAngle(Double_t angleDegree){ m_angle = angleDegree; }
  void SetClusterSize(Int_t size){ m_cluster_size = size; }
  void SetMWPCFlag(Bool_t flag){ m_mwpc_flag = flag; }
  void SetMeanWire(Double_t mwire){ m_mwpc_wire = mwire; }
  void SetMeanWirePosition(Double_t mwpos){ m_mwpc_wpos = mwpos; }
  void SetWirePosition(Double_t wpos){ m_wpos = wpos; }

  ///// for TOF
  void SetZ(Double_t z) { m_z = z; }

  ///// for TPC
  void SetHitNum(Int_t hitnum) { m_hitnum = hitnum; }

  ///// For E40 Acrylic TOF
  void SetOfsdT(Double_t ofs) { m_ofs_dt = ofs; }

  ///// for CFT
  void SetMeanSeg(Double_t seg) { m_meanseg = seg; }
  void SetMaxSeg(Double_t seg) { m_maxseg = seg; }
  void SetAdcLow(Double_t adc) { m_adc_low = adc; }
  void SetMIPLow(Double_t mip) { m_mip_low = mip; }
  void SetdELow(Double_t dE) { m_dE_low = dE; }
  void SetMaxAdcLow(Double_t adc) { max_adc_low = adc; }
  void SetMaxMIPLow(Double_t mip) { max_mip_low = mip; }
  void SetMaxdELow(Double_t dE) { max_dE_low = dE; }
  void SetPositionR(Double_t r) { m_r = r; }
  void SetPositionPhi(Double_t phi) { m_phi = phi; }
  void SetPosPhi(Double_t phi) { m_pos_phi  = phi; }
  void SetPosZ(Double_t z) { m_pos_z = z; }
  void SetPosR(Double_t r) { m_pos_r = r; }
  void SetTdcCFT(Int_t tdc); // not used
  void SetTime(Double_t time) { m_time = time; }
  void SetVtx(const TVector3& vtx) { m_vtx = vtx; }

  void GateDriftTime(Double_t min, Double_t max, Bool_t select_1st);

  Int_t GetTdcSize() const { return m_tdc.size(); }
  Int_t GetAdcSize() const { return m_adc.size(); }
  Int_t GetTdcVal(Int_t i=0) const { return m_tdc[i]; }
  Int_t GetAdcVal(Int_t i=0) const { return m_adc[i]; }
  Int_t GetTdcTrailing(Int_t i=0) const { return m_trailing[i]; }
  Int_t GetTdcTrailingSize() const { return m_trailing.size(); }
  Int_t GetTdc1st() const;
  Int_t GetEntries() const { return m_data.size(); }

  Double_t GetResolution()       const;

  Double_t GetDriftTime(Int_t i=0) const { return m_data.at(i).drift_time; }
  Double_t GetDriftLength(Int_t i=0) const { return m_data.at(i).drift_length; }
  Double_t TimeOverThreshold(Int_t i=0) const { return m_data.at(i).tot; }

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

  Double_t GetTiltAngle() const { return m_angle; }
  Double_t GetWirePosition() const {
    if(m_mwpc_flag) return m_mwpc_wpos;
    else return m_wpos;
  }

  Int_t GetClusterSize() const { return m_cluster_size; }
  Double_t GetMeamWire() const { return m_mwpc_wire; }
  Double_t GetMeamWirePosition() const { return m_mwpc_wpos; }

  ///// for TOF
  Double_t GetZ() const { return m_z; }

  ///// for TPC
  Double_t GetHitNum() { return m_hitnum; }

  ///// for CFT
  Double_t GetMeanSeg() const { return m_meanseg; }
  Double_t GetMaxSeg() const { return m_maxseg;  }
  Double_t SetMaxSeg() const { return m_maxseg;  }
  Double_t GetAdcLow() const { return m_adc_low; }
  Double_t GetMIPLow() const { return m_mip_low; }
  Double_t GetdELow() const { return m_dE_low ; }
  Double_t GetMaxAdcLow() const { return max_adc_low; }
  Double_t GetMaxMIPLow() const { return max_mip_low; }
  Double_t GetMaxdELow() const { return max_dE_low ; }

  Double_t GetPositionR() const { return m_r;        }
  Double_t GetPositionPhi() const { return m_phi;      }
  Double_t GetPosPhi() const { return m_pos_phi;  }
  Double_t GetPosZ() const { return m_pos_z;    }
  Double_t GetPosR() const { return m_pos_r;    }
  Double_t GetTime() const { return m_time;     }
  TVector3 GetVtx() const { return m_vtx;      }

  void JoinTrack(Int_t i=0) { m_data.at(i).belong_to_track = true; }
  void QuitTrack(Int_t i=0) { m_data.at(i).belong_to_track = false; }
  Bool_t BelongToTrack(Int_t i=0) const { return m_data.at(i).belong_to_track; }
  Bool_t IsWithinRange(Int_t i=0) const { return m_data.at(i).dl_is_good; }

  void RegisterHits(DCLTrackHit *hit)
    { m_register_container.push_back(hit); }

  Bool_t ReCalcDC(Bool_t applyRecursively=false) { return CalcDCObservables(); }
  // Bool_t ReCalcMWPC(Bool_t applyRecursively=false) { return CalcMWPCObservables(); }

  void TotCut(Double_t min_tot, Bool_t adopt_nan);

  void Print(const TString& arg="", std::ostream& ost=hddaq::cout) const;

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

//_____________________________________________________________________
inline std::ostream&
operator <<(std::ostream& ost, const DCHit& hit)
{
  hit.Print("", ost);
  return ost;
}

#endif
