// -*- C++ -*-

#ifndef TPC_LOCAL_TRACK_HELIX_HH
#define TPC_LOCAL_TRACK_HELIX_HH

#include <vector>
#include <functional>

#include <std_ostream.hh>

#include "TMinuit.h"
#include "TVector3.h"
#include "ThreeVector.hh"
#include "DetectorID.hh"
#include "TPCHit.hh"
#include "TPCLTrackHit.hh"
#include "TPCCluster.hh"
#include <TH2D.h>

class TPCHit;
class TPCCluster;
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
  TMinuit *minuit;

private:
  bool   m_is_fitted;     // flag of DoFit()
  bool   m_is_calculated; // flag of Calculate()
  std::vector<TPCLTrackHit*> m_hit_array;
  std::vector<TPCCluster*> m_cluster_array;
  // TH2D *hist;

  //equation of Helix
  //x = -X;
  //y = Z - Tgtz;
  //z = Y;
  //x = (1./rho + drho)*cos(phi0) - (1./rho)*cos(phi0 + helixphi)
  //y = (1./rho + drho)*sin(phi0) - (1./rho)*sin(phi0 + helixphi)
  //z = dz - (1./rho)*tanL*helixphi;
  //helixphi = atan2(y-cy, x-cx)

  double m_Adrho;
  double m_Aphi0;
  double m_Arho;
  double m_Adz;
  double m_AtanL;

  double m_drho;
  double m_phi0;
  double m_rho;
  double m_dz;
  double m_tanL;

  double m_chisqr;
  bool   m_good_for_tracking;
  double m_n_iteration;
  // for SSD
  double m_de; // not use now

  TVector3 mom0;

public:
  void         AddTPCHit( TPCLTrackHit *hit );
  void         AddTPCCluster( TPCCluster *cluster );//not supported
  void         ClearHits( void );
  void         Calculate( void );
  void         DeleteNullHit( void );
  bool         DoHelixFit( int MinHits );
  //  bool         DoHelixFit( void );
  bool         DoFit( int MinHits );
  int          GetNDF( void ) const;
  int          GetNHit( void ) const { return m_hit_array.size();  }
  TPCLTrackHit* GetHit( std::size_t nth ) const;
  bool         IsFitted( void ) const { return m_is_fitted; }
  bool         IsCalculated( void ) const { return m_is_calculated; }
  bool         Residual_check( TVector3 pos, TVector3 Res);
  void         CalcChi2( void);
  
  void SetAdrho( double Adrho ) { m_Adrho = Adrho; }
  void SetAphi0( double Aphi0 ) { m_Aphi0 = Aphi0; }
  void SetArho( double Arho ) { m_Arho = Arho; }
  void SetAdz( double Adz )  { m_Adz = Adz; }
  void SetAtanL( double AtanL ){  m_AtanL = AtanL; }
 
  void SetDe( double de ) { m_de = de; }

  double Getdrho( void ) const { return m_drho; }
  double Getphi0( void ) const { return m_phi0; }
  double Getrho( void ) const { return m_rho; }
  double Getdz( void ) const { return m_dz; }
  double GettanL( void ) const { return m_tanL; }
 
  TVector3 GetMom0( void ) const { return mom0; }// Momentum at Y = 0


  double GetAdrho( void ) const { return m_Adrho; }
  double GetAphi0( void ) const { return m_Aphi0; }
  double GetArho( void ) const { return m_Arho; }
  double GetAdz( void ) const { return m_Adz; }
  double GetAtanL( void ) const { return m_AtanL; }

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
  double GetDe( void ) const { return m_de; }
  void   Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const;
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
