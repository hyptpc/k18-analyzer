// -*- C++ -*-

#ifndef TPC_LTRACK_HIT_HH
#define TPC_LTRACK_HIT_HH

#include <TVector3.h>

#include "DCHit.hh"
#include "MathTools.hh"
#include "ThreeVector.hh"
#include "TPCHit.hh"

class DCAnalyzer;

//_____________________________________________________________________________
class TPCLTrackHit
{
public:
  static const TString& ClassName();
  TPCLTrackHit(TPCHit *hit);
  TPCLTrackHit(const TPCLTrackHit& right);
  ~TPCLTrackHit();

private:
  TPCHit*  m_hit;
  TVector3 m_local_hit_pos;
  TVector3 m_cal_pos;
  Double_t  m_x0;
  Double_t  m_y0;
  Double_t  m_u0;
  Double_t  m_v0;
  Double_t  m_cx;
  Double_t  m_cy;
  Double_t  m_z0;
  Double_t  m_r;
  Double_t  m_dz;

  TVector3  m_res;

public:
  void   SetLocalHitPos(TVector3 xl)        { m_local_hit_pos = xl; }
  void   SetCalPosition(TVector3 cal_pos)        { m_cal_pos = cal_pos; }
  void   SetCalX0Y0(Double_t x0, Double_t y0)          { m_x0 = x0; m_y0 = y0; }
  void   SetCalUV(Double_t u0, Double_t v0)            { m_u0 = u0; m_v0 = v0; }
  void   SetCalHelix(Double_t cx, Double_t cy, Double_t z0, Double_t r, Double_t dz)
  {m_cx = cx; m_cy = cy; m_z0 = z0; m_r = r, m_dz = dz;}
  int    GetLayer()                 const { return m_hit->GetLayer(); }
  void   SetResolution(TVector3 res)        { m_res = res; }

  TPCHit* GetHit() const { return m_hit; }

  const TVector3& GetLocalHitPos()  const { return m_local_hit_pos; }
  TVector3 GetLocalCalPos()  const;
  TVector3 GetLocalCalPosHelix()  const;
  TVector3 GetHelixPosition(Double_t par[5], Double_t t)  const;
  TVector3 GetMomentumHelix()  const;

  Double_t GetTcal(void)  const;

  Double_t GetXcal()     const { return m_cal_pos.x(); }
  Double_t GetYcal()     const { return m_cal_pos.y(); }
  Double_t GetZcal()     const { return m_cal_pos.z(); }

  Double_t GetUcal()     const { return m_u0; }
  Double_t GetVcal()     const { return m_v0; }

  Double_t Getcx()     const { return m_cx; }
  Double_t Getcy()     const { return m_cy; }
  Double_t Getz0()      const { return m_z0; }
  Double_t Getr()       const { return m_r; }
  Double_t Getdz()     const { return m_dz; }

  TVector3 GetResidualVect() const;
  Double_t GetResidual() const;
  Double_t GetResolution() const { return m_hit->GetResolution(); }
  const TVector3& GetResolutionVect() const { return m_res; }
  Bool_t ResidualCut() const;

  void   JoinTrack() { m_hit->JoinTrack(); }
  void   QuitTrack() { m_hit->QuitTrack(); }
  Bool_t BelongToTrack() const { return m_hit->BelongToTrack(); }

  void Print(const std::string& arg="") const;

  friend class TPCHit;
};

//_____________________________________________________________________________
inline const TString&
TPCLTrackHit::ClassName()
{
  static TString s_name("TPCLTrackHit");
  return s_name;
}

#endif
