// -*- C++ -*-

#include "TPCLocalTrack.hh"

#include <string>
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>

#include <std_ostream.hh>

#include "DCAnalyzer.hh"
#include "DCLTrackHit.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "HodoParamMan.hh"
#include "UserParamMan.hh"

#define HSMagnetON 1

namespace
{
  const std::string& class_name("TPCLocalTrack");
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const double& zK18tgt = gGeom.LocalZ("K18Target");
  const double& zTgt    = gGeom.LocalZ("Target");

  const HodoParamMan& gHodo = HodoParamMan::GetInstance();

  const int ReservedNumOfHits  = 16;
  const int DCLocalMinNHits    =  4;
  const int DCLocalMinNHitsVXU =  2;// for SSD
  const int MaxIteration       = 100;// for Honeycomb
  const double MaxChisqrDiff   = 1.0e-3;
  const int CFTLocalMinNHits   =  3;
}

//______________________________________________________________________________
TPCLocalTrack::TPCLocalTrack( void )
  : m_is_fitted(false),
    m_is_calculated(false),
    m_Ax(0.), m_Ay(0.), m_Au(0.), m_Av(0.),
    m_Chix(0.), m_Chiy(0.), m_Chiu(0.), m_Chiv(0.),
    m_x0(0.), m_y0(0.),
    m_u0(0.), m_v0(0.),
    m_a(0.),  m_b(0.),
    m_chisqr(1.e+10),
    m_good_for_tracking(true),
    m_n_iteration(0),
    m_de(0.)
{
  m_hit_array.reserve( ReservedNumOfHits );
  debug::ObjectCounter::increase(class_name);

}

//______________________________________________________________________________
TPCLocalTrack::~TPCLocalTrack( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
TPCLocalTrack::AddTPCHit( TPCHit *hit )
{
  if( hit )
    m_hit_array.push_back( hit );
}

//______________________________________________________________________________
void
TPCLocalTrack::AddTPCCluster( TPCCluster *cluster )
{
  if( cluster )
    m_cluster_array.push_back( cluster );
}

//______________________________________________________________________________
void
TPCLocalTrack::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsCalculated() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return;
  }

//  const std::size_t n = m_hit_array.size();
//  for( std::size_t i=0; i<n; ++i ){
//    DCLTrackHit *hitp = m_hit_array[i];
//    double z0 = hitp->GetZ();
//    hitp->SetCalPosition( GetX(z0), GetY(z0) );
//    hitp->SetCalUV( m_u0, m_v0 );
//    // for Honeycomb
//    if( !hitp->IsHoneycomb() )
//      continue;
//    double scal = hitp->GetLocalCalPos();
//    double wp   = hitp->GetWirePosition();
//    double dl   = hitp->GetDriftLength();
//    double ss   = scal-wp>0 ? wp+dl : wp-dl;
//    hitp->SetLocalHitPos( ss );
//  }
  m_is_calculated = true;
}

//______________________________________________________________________________
int
TPCLocalTrack::GetNDF( void ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const std::size_t n = m_hit_array.size();
  int ndf = 0;
  for( std::size_t i=0; i<n; ++i ){
    if( m_hit_array[i] ) ++ndf;
  }
  return ndf-4;
}

//______________________________________________________________________________
TPCHit*
TPCLocalTrack::GetHit( std::size_t nth ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( nth<m_hit_array.size() )
    return m_hit_array[nth];
  else
    return 0;
}

//______________________________________________________________________________
void
TPCLocalTrack::DeleteNullHit( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  for( std::size_t i=0; i<m_hit_array.size(); ++i ){
    TPCHit *hit = m_hit_array[i];
    if( !hit ){
      hddaq::cout << func_name << " "
		  << "null hit is deleted" << std::endl;
      m_hit_array.erase( m_hit_array.begin()+i );
      --i;
    }
  }
}
//______________________________________________________________________________
bool
TPCLocalTrack::DoFit( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsFitted() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return false;
  }

#if HSMagnetON
  DoHelixFit();
#else
  DoLinearFit();
#endif
  m_is_fitted = true;
  return true;
}


//______________________________________________________________________________
bool
TPCLocalTrack::DoLinearFit( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  return true;
}

//______________________________________________________________________________
bool
TPCLocalTrack::DoHelixFit( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  return true;
}


//______________________________________________________________________________
double
TPCLocalTrack::GetTheta( void ) const
{
  double cost = 1./std::sqrt(1.+m_u0*m_u0+m_v0*m_v0);
  return std::acos(cost)*math::Rad2Deg();
}

//______________________________________________________________________________
bool
TPCLocalTrack::ReCalc( bool applyRecursively )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

//  std::size_t n = m_hit_array.size();
//  for( std::size_t i=0; i<n; ++i ){
//    TPCHit *hit = m_hit_array[i];
//    if( hit ) hit->ReCalc( applyRecursively );
//  }
//
  bool status = DoFit();
  if( !status ){
    hddaq::cerr << "#W " << func_name << " "
		<< "Recalculation fails" << std::endl;
  }

  return status;
}

//______________________________________________________________________________
void
TPCLocalTrack::Print( const std::string& arg, std::ostream& ost ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  PrintHelper helper( 3, std::ios::fixed, ost );

//  const int w = 8;
//  ost << func_name << " " << arg << std::endl
//      << " X0 : " << std::setw(w) << std::left << m_x0
//      << " Y0 : " << std::setw(w) << std::left << m_y0
//      << " U0 : " << std::setw(w) << std::left << m_u0
//      << " V0 : " << std::setw(w) << std::left << m_v0;
//  helper.setf( std::ios::scientific );
//  ost << " Chisqr : " << std::setw(w) << m_chisqr << std::endl;
//  helper.setf( std::ios::fixed );
//  const std::size_t n = m_hit_array.size();
//  for( std::size_t i=0; i<n; ++i ){
//    DCLTrackHit *hitp = m_hit_array[i];
//    if( !hitp ) continue;
//    int lnum = hitp->GetLayer();
//    double zz = hitp->GetZ();
//    double s  = hitp->GetLocalHitPos();
//    double res = hitp->GetResidual();
//    // double aa = hitp->GetTiltAngle()*math::Deg2Rad();
//    // double scal=GetX(zz)*cos(aa)+GetY(zz)*sin(aa);
//    const std::string& h = hitp->IsHoneycomb() ? "+" : "-";
//    ost << "[" << std::setw(2) << i << "]"
//	<< " #"  << std::setw(2) << lnum << h
//	<< " S " << std::setw(w) << s
//	<< " ( " << std::setw(w) << GetX(zz)
//	<< ", "  << std::setw(w) << GetY(zz)
//	<< ", "  << std::setw(w) << zz
//	<< " )"
//	<< " " << std::setw(w) << s
//	<< " -> " << std::setw(w) << res << std::endl;
//	// << " -> " << std::setw(w) << s-scal << std::endl;
//  }
//  ost << std::endl;
}
