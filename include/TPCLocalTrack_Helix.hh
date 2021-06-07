// -*- C++ -*-

#ifndef TPC_LOCAL_TRACK_HELIX_HH
#define TPC_LOCAL_TRACK_HELIX_HH

#include <vector>
#include <functional>

#include <std_ostream.hh>

//#include "TMinuit.h"
#include "TVector3.h"
#include "ThreeVector.hh"
#include "DetectorID.hh"
#include "TPCHit.hh"
#include "TPCLTrackHit.hh"



class TPCHit;
class DCAnalyzer;

//______________________________________________________________________________
class TPCLocalTrack_Helix
{
public:
  explicit TPCLocalTrack_Helix( void );
  ~TPCLocalTrack_Helix( void );

private:
  TPCLocalTrack_Helix( const TPCLocalTrack_Helix & );
  TPCLocalTrack_Helix & operator =( const TPCLocalTrack_Helix & );
  //  TMinuit *minuit;

private:
  bool   m_is_fitted;     // flag of DoFit()
  bool   m_is_calculated; // flag of Calculate()
  std::vector<TPCLTrackHit*> m_hit_array;

  //equation of Helix
  //x = -X;
  //y = Z - Tgtz;
  //z = Y;
  //x = p[0] + p[3]*cos(t+theta0);
  //y = p[1] + p[3]*sin(t+theta0);
  //z = p[2] + (p[4]*p[3]*t);
  double m_cx;
  double m_cy;
  double m_z0;
  double m_r;
  double m_dz;
  double m_Acx;
  double m_Acy;
  double m_Az0;
  double m_Ar;
  double m_Adz;
  double m_chisqr;
  bool   m_good_for_tracking;
  double m_n_iteration;
  TVector3 m_mom0;
  double m_min_t;
  double m_max_t;
  double m_path;
  int m_isbeam;

public:
  void         AddTPCHit( TPCLTrackHit *hit );
  void         ClearHits( void );
  void         Calculate( void );
  void         DeleteNullHit( void );
  bool         DoHelixFit( int MinHits, int IsBeam );
  //  bool         DoHelixFit( void );
  bool         DoFit( int MinHits, int IsBeam );
  int          GetNDF( void ) const;

  TVector3     GetPosition( double par[5], double t ) const;
  TVector3     CalcHelixMom( double par[5], double y) const;
  TVector3     CalcHelixMom_t( double par[5], double t) const;
  int          GetNHit( void ) const { return m_hit_array.size();  }
  TPCLTrackHit* GetHit( std::size_t nth ) const;
  bool         IsFitted( void ) const { return m_is_fitted; }
  bool         IsCalculated( void ) const { return m_is_calculated; }
  bool         Residual_check( TVector3 pos, TVector3 Res, double resi);
  double       GetTcal( TVector3 pos);
  void         CalcChi2( void);
  double       CalcChi2_circle( double par[3]);
  int     GetHTOFSeg( double min_layer_t, double max_layer_t, double max_layer_y );
  int     GetIsBeam(void) const {return m_isbeam;}
  // void SetMint( double min_t ) { m_min_t = min_t; }
  // void SetMaxt( double max_t ) { m_max_t = max_t; }
  // void SetPath( double path ) { m_path = path; }
  double GetMint(void) const {return m_min_t; }
  double GetMaxt(void) const {return m_max_t; }
  double GetPath(void) const {return m_path; }


  void SetAcx( double Acx ) { m_Acx = Acx; }
  void SetAcy( double Acy ) { m_Acy = Acy; }
  void SetAz0( double Az0 ) { m_Az0 = Az0; }
  void SetAr( double Ar )  { m_Ar = Ar; }
  void SetAdz( double Adz ){  m_Adz = Adz; }
  double Getcx( void ) const { return m_cx; }
  double Getcy( void ) const { return m_cy; }
  double Getz0( void ) const { return m_z0; }
  double Getr( void ) const { return m_r; }
  double Getdz( void ) const { return m_dz; }

  TVector3 GetMom0( void ) const { return m_mom0; }// Momentum at Y = 0


  double GetAcx( void ) const { return m_Acx; }
  double GetAcy( void ) const { return m_Acy; }
  double GetAz0( void ) const { return m_Az0; }
  double GetAr( void ) const { return m_Ar; }
  double GetAdz( void ) const { return m_Adz; }

  double GetChiSquare( void ) const { return m_chisqr; }
  // double GetChiX( void ) const { return m_Chix; }
  // double GetChiY( void ) const { return m_Chiy; }
  // double GetChiU( void ) const { return m_Chiu; }
  // double GetChiV( void ) const { return m_Chiv; }

  // double GetX( double z ) const { return m_x0+m_u0*z; }
  // double GetY( double z ) const { return m_y0+m_v0*z; }
  // double GetS( double z, double tilt ) const { return GetX(z)*std::cos(tilt)+GetY(z)*std::sin(tilt); }
  int    GetNIteration( void ) const { return m_n_iteration; }
  //double GetTheta( void ) const;
  bool   GoodForTracking( void ) const { return m_good_for_tracking; }
  bool   GoodForTracking( bool status )
  { bool ret = m_good_for_tracking; m_good_for_tracking = status; return ret; }
};


//______________________________________________________________________________
// inline
// std::ostream&
// operator <<( std::ostream& ost,
// 	     const TPCLocalTrack_Helix& track )
// {
//   track.Print( "", ost );
//   return ost;
// }


#endif
