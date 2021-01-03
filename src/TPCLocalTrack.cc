// -*- C++ -*-
//Comment by Ichikawa
//TPCLocalTrack.cc is for lenear fit 
//TPCLocalTrack_Helix.cc will be prepared for Helix tracking 
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
#include "TPCLTrackHit.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "HodoParamMan.hh"
#include "UserParamMan.hh"

#include "TMath.h"
#include "TROOT.h"

//#define HSMagnetON 1 

static int gNumOfHits;
static TVector3 gHitPos[300];
static TVector3 gRes[300];

namespace
{
  const std::string& class_name("TPCLocalTrack");
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const double& zK18tgt = gGeom.LocalZ("K18Target");
  const double& zTgt    = gGeom.LocalZ("Target");


  const int ReservedNumOfHits  = 32;
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  static const double  FitStep[5] = { 1.0e-10, 1.0e-11, 1.0e-10, 1.0e-11};
  static const double  LowLimit[5] = { -300., -100, -300, -100. };
  static const double  UpLimit[5] = {  300., 100, 300, 100. };
  

}

//______________________________________________________________________________
TPCLocalTrack::TPCLocalTrack( void )
  : m_is_fitted(false),
    m_is_calculated(false),
    m_Ax(0.), m_Ay(0.), m_Au(0.), m_Av(0.),
    //m_Chix(0.), m_Chiy(0.), m_Chiu(0.), m_Chiv(0.),
    m_x0(0.), m_y0(0.),
    m_u0(0.), m_v0(0.),
    m_a(0.),  m_b(0.),
    m_chisqr(1.e+10),
    m_good_for_tracking(true),
    m_n_iteration(0),
    m_de(0.)
{
  m_hit_array.reserve( ReservedNumOfHits );
  m_cluster_array.reserve( ReservedNumOfHits );
  debug::ObjectCounter::increase(class_name);
  minuit = new TMinuit(4);
  TROOT minexam("LinearFit","linear fit using TMinuit");
}

//______________________________________________________________________________
TPCLocalTrack::~TPCLocalTrack( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
TPCLocalTrack::AddTPCHit( TPCLTrackHit *hit )
{
  if( hit )
    m_hit_array.push_back( hit );
}

//______________________________________________________________________________
void
TPCLocalTrack::AddTPCCluster( TPCCluster *cluster )
{
  //not supported
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
  
  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    hitp->SetCalX0Y0(m_x0, m_y0);
    hitp->SetCalUV( m_u0, m_v0 );
    hitp->SetCalPosition( hitp->GetLocalCalPos() );
  }
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
static void fcn2(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;



  for( int i=0; i<gNumOfHits; ++i ){
    TVector3 pos = gHitPos[i];
    TVector3 Res = gRes[i];
 
    TVector3 x0(par[0], par[2], 0. );
    TVector3 x1(par[0] + par[1], par[2] + par[3], 1. );
    TVector3 u = (x1-x0).Unit();
    //    TVector3 d = (pos-x0).Cross(u);
    TVector3 AP = pos-x0;
    double dist_AX = u.Dot(AP);
    TVector3 AI(x0.x()+(u.x()*dist_AX),
		x0.y()+(u.y()*dist_AX),
		x0.z()+(u.z()*dist_AX));
    TVector3 d = pos-AI; 

    chisqr += pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    dof++;
  }

  f = chisqr/(double)(dof-4);
}



//______________________________________________________________________________
TPCLTrackHit*
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
    TPCLTrackHit *hit = m_hit_array[i];
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

  DoLinearFit();

  m_is_fitted = true;
  return true;
}





//______________________________________________________________________________
bool
TPCLocalTrack::DoLinearFit( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  
  double par[4]={m_Ax, m_Au, m_Ay, m_Av};
  double err[4]={0};

  const std::size_t n = m_cluster_array.size();
  gNumOfHits = n;
  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    gHitPos[i] =pos;
    gRes[i] = Res;
  }




  minuit = new TMinuit(4);
  TROOT minexam("LinearFit", "Linear fit using TMinuit");

  minuit->SetPrintLevel(-1);
  minuit->SetFCN (fcn2 );

  int ierflg = 0;
  double arglist[10];
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings                                     
                                                                                                     

  TString name[4] = {"x0", "u0", "y0", "v0"};
  for( int i = 0; i<4; i++ )
    {
      minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
    }

  minuit->Command("SET STRategy 0");
  arglist[0] = 1000.;
  arglist[1] = 1.0;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  minuit->mnstat( amin, edm, errdef, nvpar, nparx, icstat);
  //minuit->mnprin(4, amin);
  int Err;
  double bnd1, bnd2;
  for( int i=0; i<4; i++)
    {
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
      //std::cout<<Par[i]<<"  "<<std::endl;                                                          
    }
  
  m_x0=par[0];
  m_u0=par[1];
  m_y0=par[2];
  m_v0=par[3];
  

  //  FitStat = icstat;
  CalcChi2();
  
  return true;
}

//______________________________________________________________________________
void
TPCLocalTrack::CalcChi2( void )
{
  double chisqr=0.0;
  int dof = 0;

  const std::size_t n = m_cluster_array.size();

  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
 
    TVector3 x0(m_x0, m_y0, 0. );
    TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, 1. );
    TVector3 u = (x1-x0).Unit();
    TVector3 AP = pos-x0;
    double dist_AX = u.Dot(AP);
    TVector3 AI(x0.x()+(u.x()*dist_AX),
		x0.y()+(u.y()*dist_AX),
		x0.z()+(u.z()*dist_AX));
    TVector3 d = pos-AI; 
    //    TVector3 d = (pos-x0).Cross(u);
    
    chisqr += pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    dof++;
  }
  m_chisqr = chisqr/(double)(dof-4);
}

//______________________________________________________________________________
// bool
// TPCLocalTrack::DoHelixFit( void )
// {
//   static const std::string func_name("["+class_name+"::"+__func__+"()]");

//   return true;
// }


//______________________________________________________________________________
double
TPCLocalTrack::GetTheta( void ) const
{
  double cost = 1./std::sqrt(1.+m_u0*m_u0+m_v0*m_v0);
  return std::acos(cost)*math::Rad2Deg();
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
