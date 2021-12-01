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
  TPCHit  *m_hit;
  TVector3  m_local_hit_pos;
  TVector3  m_cal_pos;
  double  m_x0;
  double  m_y0;
  double  m_u0;
  double  m_v0;
  double  m_cx;
  double  m_cy;
  double  m_z0;
  double  m_r;
  double  m_dz;

  TVector3  m_res;

public:
  void   SetLocalHitPos(TVector3 xl)        { m_local_hit_pos = xl; }
  void   SetCalPosition(TVector3 cal_pos)        { m_cal_pos = cal_pos; }
  void   SetCalX0Y0(double x0, double y0)          { m_x0 = x0; m_y0 = y0; }
  void   SetCalUV(double u0, double v0)            { m_u0 = u0; m_v0 = v0; }
  void   SetCalHelix(double cx, double cy, double z0, double r, double dz)
  {m_cx = cx; m_cy = cy; m_z0 = z0; m_r = r, m_dz = dz;}
  int    GetLayer()                 const { return m_hit->GetLayer(); }
  void   SetResolution(TVector3 res)        { m_res = res; }

  TPCHit* GetHit() const { return m_hit; }

  const TVector3& GetLocalHitPos()  const { return m_local_hit_pos; }
  TVector3 GetLocalCalPos()  const;
  TVector3 GetLocalCalPosHelix()  const;
  TVector3 GetHelixPosition(double par[5], double t)  const;
  TVector3 GetMomentumHelix()  const;

  double GetTcal(void)  const;

  double GetXcal()     const { return m_cal_pos.x(); }
  double GetYcal()     const { return m_cal_pos.y(); }
  double GetZcal()     const { return m_cal_pos.z(); }

  double GetUcal()     const { return m_u0; }
  double GetVcal()     const { return m_v0; }

  double Getcx()     const { return m_cx; }
  double Getcy()     const { return m_cy; }
  double Getz0()      const { return m_z0; }
  double Getr()       const { return m_r; }
  double Getdz()     const { return m_dz; }

  TVector3 GetResidualVect() const;
  double GetResidual() const;
  double GetResolution() const { return m_hit->GetResolution(); }
  const TVector3& GetResolutionVect() const { return m_res; }
  bool ResidualCut() const;


  void JoinTrack() { m_hit->JoinTrack(); }
  void QuitTrack() { m_hit->QuitTrack(); }
  bool BelongToTrack() const { return m_hit->BelongToTrack(); }

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
