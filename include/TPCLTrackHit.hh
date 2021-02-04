/**
 *  file: TPCLTrackHit.hh
 *  date: 2020.01.03
 *  based onf DCLTrackHit.hh
 */

#ifndef TPC_LTRACK_HIT_HH
#define TPC_LTRACK_HIT_HH

#include "DCHit.hh"
#include "TPCHit.hh"

#include "MathTools.hh"
#include "ThreeVector.hh"
#include <TVector3.h>

class DCAnalyzer;

//______________________________________________________________________________
class TPCLTrackHit
{
public:
  TPCLTrackHit( TPCHit *hit );
  TPCLTrackHit( const TPCLTrackHit& right );

private:
  ~TPCLTrackHit( void );

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
  void   SetLocalHitPos( TVector3 xl )        { m_local_hit_pos = xl; }
  void   SetCalPosition( TVector3 cal_pos )        { m_cal_pos = cal_pos; }
  void   SetCalX0Y0( double x0, double y0 )          { m_x0 = x0; m_y0 = y0; }
  void   SetCalUV( double u0, double v0 )            { m_u0 = u0; m_v0 = v0; }
  void   SetCalHelix( double cx, double cy, double z0, double r, double dz)
  {m_cx = cx; m_cy = cy; m_z0 = z0; m_r = r, m_dz = dz;}
  int    GetLayer( void )                 const { return m_hit->GetLayer(); }
  void   SetResolution( TVector3 res )        { m_res = res; }

  TPCHit* GetHit( void ) const { return m_hit; }

  TVector3 GetLocalHitPos( void )  const { return m_local_hit_pos; }
  TVector3 GetLocalCalPos( void )  const;
  TVector3 GetLocalCalPos_Helix( void )  const;
  TVector3 GetHelixPosition( double par[5], double t )  const;
  TVector3 GetMomentum_Helix( void )  const;

  double GetXcal( void )     const { return m_cal_pos.x(); }
  double GetYcal( void )     const { return m_cal_pos.y(); }
  double GetZcal( void )     const { return m_cal_pos.z(); }

  double GetUcal( void )     const { return m_u0; }
  double GetVcal( void )     const { return m_v0; }

  double Getcx( void )     const { return m_cx; }
  double Getcy( void )     const { return m_cy; }
  double Getz0( void )      const { return m_z0; }
  double Getr( void )       const { return m_r; }
  double Getdz( void )     const { return m_dz; }

  TVector3 GetResidualVect( void ) const;
  double GetResidual( void ) const;
  double GetResolution( void ) const { return m_hit->GetResolution(); }
  TVector3 GetResolutionVect( void ) const { return m_res; }
  bool ResidualCut( void ) const;


  void JoinTrack( void ) { m_hit->JoinTrack(); }
  void QuitTrack( void ) { m_hit->QuitTrack(); }
  bool BelongToTrack( void ) const { return m_hit->BelongToTrack(); }

  void Print( const std::string& arg="" ) const;


  friend class TPCHit;
};

#endif
