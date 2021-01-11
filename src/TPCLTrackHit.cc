/**
 *  file: TPCLTrackHit.cc
 *  date: 2020.01.03
 *  (ref: DCLTrackHit.cc)
 */

//#include "DCLTrackHit.hh"
#include "TPCLTrackHit.hh"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <std_ostream.hh>

#include "DCAnalyzer.hh"
#include "MathTools.hh"

namespace
{
  const std::string& class_name("TPCLTrackHit");
}

//______________________________________________________________________________
TPCLTrackHit::TPCLTrackHit( TPCHit *hit )
  : m_hit(hit),
    m_x0(-9999.),
    m_y0(-9999.),
    m_u0(-9999.),
    m_v0(-9999.)
{
  m_local_hit_pos = hit->GetPos();
  m_cal_pos = TVector3(0.,0.,0.);
  m_res = TVector3(hit->GetResolutionX(),
		   hit->GetResolutionY(),
		   hit->GetResolutionZ());
  debug::ObjectCounter::increase(class_name);
  m_hit->RegisterHits(this);
}

//______________________________________________________________________________
TPCLTrackHit::TPCLTrackHit( const TPCLTrackHit& right )
  : m_hit(right.m_hit),
    m_x0(right.m_x0),
    m_y0(right.m_y0),
    m_u0(right.m_u0),
    m_v0(right.m_v0)
{
  m_local_hit_pos = right.m_local_hit_pos;
  m_cal_pos = right.m_cal_pos;
  m_res = right.m_res;
  m_hit->RegisterHits(this);
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
TPCLTrackHit::~TPCLTrackHit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetLocalCalPos( void ) const
{
  TVector3 pos = m_local_hit_pos;
  // TVector3 x0(m_x0, m_y0, 0. );
  // TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, 1. );
  //temp
  double zTgtTPC = -143.;
  TVector3 x0(m_x0, m_y0, zTgtTPC );
  TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, zTgtTPC+1. );
  TVector3 u = (x1-x0).Unit();
  TVector3 AP = pos-x0;
  double dist_AX = u.Dot(AP);
  TVector3 AI(x0.x()+(u.x()*dist_AX),
	      x0.y()+(u.y()*dist_AX),
	      x0.z()+(u.z()*dist_AX));
  return AI;
}

//______________________________________________________________________________
TVector3
TPCLTrackHit::GetResidualVect( void ) const
{
  TVector3 Res = m_cal_pos - m_local_hit_pos;
  return Res;
}

//______________________________________________________________________________
double
TPCLTrackHit::GetResidual( void ) const
{
  TVector3 Res = m_cal_pos - m_local_hit_pos;
  return Res.Mag();
}

//______________________________________________________________________________
bool
TPCLTrackHit::ResidualCut( void ) const
{
  bool status = false;
  TVector3 Res = m_cal_pos - m_local_hit_pos;
  double resolution = m_hit->GetResolution();
  if(Res.Mag()<resolution*5.)
    status = true;
  return status;
}

//______________________________________________________________________________
void
TPCLTrackHit::Print( const std::string& arg ) const
{
  m_hit->Print( arg );
  hddaq::cout << "local_hit_pos " << m_local_hit_pos.x() 
	      <<", "<<m_local_hit_pos.y() 
	      <<", "<<m_local_hit_pos.z() <<std::endl
	      << "residual " << GetResidual() << std::endl;
}

