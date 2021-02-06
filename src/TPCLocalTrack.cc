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
#include <TH2D.h>

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
  // Temporary
  const double& zTgtTPC    = -143.;
  //  const int MaxTry    = 10;

  const int ReservedNumOfHits  = 64;
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  static const double  FitStep[4] = { 1.0e-6, 1.0e-10, 1.0e-6, 1.0e-10};
  //static const double  FitStep[4] = { 1.0e-2, 1.0e-6, 1.0e-2, 1.0e-6};
  static const double  LowLimit[4] = { -200., -0.5*acos(-1.), -200, -0.5*acos(-1.) };
  static const double  UpLimit[4] = { 200., 0.5*acos(-1.), 200, 0.5*acos(-1.) };
  static const double  MaxChisqr = 500.;

  const int    theta_ndiv = 100;
  const double theta_min  =   0;
  const double theta_max  = 180;
  const int    r_ndiv =  100;
  const double r_min  = -500;
  const double r_max  =  500;
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
    m_de(0.),
    m_minuit(new TMinuit(4))
{
  m_hit_array.reserve( ReservedNumOfHits );
  m_cluster_array.reserve( ReservedNumOfHits );
  debug::ObjectCounter::increase(class_name);
  TROOT minexam("LinearFit","linear fit using TMinuit");
}

//______________________________________________________________________________
TPCLocalTrack::~TPCLocalTrack( void )
{
  debug::ObjectCounter::decrease(class_name);
  delete m_minuit;
}

//______________________________________________________________________________
void
TPCLocalTrack::ClearHits( void )
{
  m_hit_array.clear();
  m_cluster_array.clear();
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

    double u0 = tan(par[1]);
    double v0 = tan(par[3]);

    TVector3 x0(par[0], par[2], zTgtTPC );
    TVector3 x1(par[0] + u0, par[2] + v0, zTgtTPC+1. );

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
static void fcn2_rt(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;



  for( int i=0; i<gNumOfHits; ++i ){
    TVector3 pos = gHitPos[i];
    TVector3 Res = gRes[i];

    double x_0 = par[0]/sin(par[1]);
    double y_0 = par[2]/sin(par[3]);
    double u0 = -1.*cos(par[1])/sin(par[1]);
    double v0 = -1.*cos(par[3])/sin(par[3]);

    // double m_Au_theta = atan(-1./u0);
    // double m_Av_theta = atan(-1./v0);
    // double m_Ax_r = x_0*sin(m_Au_theta);
    // double m_Ay_r = y_0*sin(m_Av_theta);

    // if(fabs(m_Au_theta-par[1])>0.01){
    //   std::cout<<"u: "<<m_Au_theta<<", "<<par[1]<<std::endl;
    //   getchar();
    // }
    // if(fabs(m_Av_theta-par[3])>0.01){
    //   std::cout<<"v: "<<m_Av_theta<<", "<<par[3]<<std::endl;
    //   getchar();
    // }
    // if(fabs(m_Ax_r-par[0])>0.01){
    //   std::cout<<"x: "<<m_Ax_r<<", "<<par[0]<<std::endl;
    //   getchar();
    // }
    // if(fabs(m_Ay_r-par[2])>0.01){
    //   std::cout<<"y: "<<m_Ay_r<<", "<<par[2]<<std::endl;
    //   getchar();
    // }

    TVector3 x0(x_0, y_0, zTgtTPC );
    TVector3 x1(x_0 + u0, y_0 + v0, zTgtTPC+1. );

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
static void fcn2_rt2(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;

  for( int i=0; i<gNumOfHits; ++i ){
    TVector3 pos = gHitPos[i];
    TVector3 Res = gRes[i];

    double x_0 = par[0]/cos(par[1]);
    double y_0 = par[2]/cos(par[3]);
    double u0 = tan(par[1]);
    double v0 = tan(par[3]);

    TVector3 x0(x_0, y_0, zTgtTPC );
    TVector3 x1(x_0 + u0, y_0 + v0, zTgtTPC+1. );

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
TPCLocalTrack::DoFit( int MinHits)
{

  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsFitted() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return false;
  }

  bool status_dofit =  DoLinearFit(MinHits);
  m_is_fitted = status_dofit;
  //  m_is_fitted = true;

  return status_dofit;
}





//______________________________________________________________________________
bool
TPCLocalTrack::DoLinearFit( int MinHits )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  DeleteNullHit();

  const std::size_t n = m_hit_array.size();

  if(n<MinHits){
    return false;
  }

  if(GetNDF()<1){
    hddaq::cerr << "#W " << func_name << " "
		<< "Min layer should be > NDF" << std::endl;
    return false;
  }


  gNumOfHits = n;
  // r = x * cos(theta) + y * sin(theta)
  static TH2D hist("hist",";theta (deg.); r (mm)",
		   theta_ndiv, theta_min, theta_max, r_ndiv, r_min, r_max);
  hist.Reset();

  //hough translation for ini-param of Ay and Av
  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    gHitPos[i] =pos;
    gRes[i] = Res;
    for( int ti=0; ti<theta_ndiv; ti++ ){
      double theta = theta_min+ti*(theta_max-theta_min)/theta_ndiv;
      hist.Fill(theta, cos(theta*acos(-1)/180.)*pos.Z()
		 +sin(theta*acos(-1)/180.)*pos.Y());
    }
  }
  int maxbin = hist.GetMaximumBin();
  int mx,my,mz;
  hist.GetBinXYZ( maxbin, mx, my, mz );
  double mtheta = hist.GetXaxis()->GetBinCenter(mx)*acos(-1)/180.;
  double mr = hist.GetYaxis()->GetBinCenter(my);
  double p0 = mr/sin(mtheta);
  double p1 = -cos(mtheta)/sin(mtheta);

  double m_Ay = p0+p1*zTgtTPC;
  double m_Av = p1;

  double m_Au_atan = atan(m_Au);
  double m_Av_atan = atan(m_Av);

  double m_Ax_r = m_Ax*cos(m_Au_atan);
  double m_Ay_r = m_Ay*cos(m_Av_atan);

  double m_Au_theta = atan(-1./m_Au);
  double m_Av_theta = atan(-1./m_Av);
  double m_Ax_r2 = m_Ax*sin(m_Au_theta);
  double m_Ay_r2 = m_Ay*sin(m_Av_theta);

  double par[4]={m_Ax, m_Au_atan, m_Ay, m_Av_atan};
  double par2[4]={m_Ax_r, m_Au_atan, m_Ay_r, m_Av_atan};
  double par3[4]={m_Ax_r2, m_Au_theta, m_Ay_r2, m_Av_theta};
  double err[4]={-999.,-999.,-999.,-999.};


  // m_minuit = new TMinuit(4);
  // TROOT minexam("LinearFit", "Linear fit using TMinuit");

  m_minuit->SetPrintLevel(-1);
  m_minuit->SetFCN (fcn2 );

  int ierflg = 0;
  double arglist[10];
  arglist[0] = 1;
  m_minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings
  TString name[4] = {"x0_r", "atan_u0", "y0_r", "atan_v0"};
  for( int i = 0; i<4; i++ )
    {
      m_minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
    }

  m_minuit->Command("SET STRategy 0");
  arglist[0] = 5000.;
  arglist[1] = 0.01;
  m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  //m_minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //m_minuit->mnexcm("SET ERR", arglist, 2, ierflg);

  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  m_minuit->mnstat( amin, edm, errdef, nvpar, nparx, icstat);
  //m_minuit->mnprin(4, amin);
  int Err;
  double bnd1, bnd2;
  for( int i=0; i<4; i++)
    {
      m_minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
      //std::cout<<Par[i]<<"  "<<std::endl;
    }


  double chisqr1 = 100000., chisqr2=100000., chisqr3=100000.;

  m_x0=par[0];
  m_u0=tan(par[1]);
  m_y0=par[2];
  m_v0=tan(par[3]);
  CalcChi2();

  chisqr1 = m_chisqr;

  m_minuit->SetFCN (fcn2_rt2 );
  for( int i = 0; i<4; i++ )
    {
      m_minuit->mnparm(i, name[i], par2[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
    }
  m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  for( int i=0; i<4; i++)
    {
      m_minuit->mnpout(i, name[i], par2[i], err[i], bnd1, bnd2, Err);
    }
  m_x0=par2[0]/cos(par2[1]);
  m_u0=tan(par2[1]);
  m_y0=par2[2]/cos(par2[3]);
  m_v0=tan(par2[3]);
  CalcChi2();
  chisqr2 = m_chisqr;


  m_minuit->SetFCN (fcn2_rt );
  for( int i = 0; i<4; i++ )
    {
      m_minuit->mnparm(i, name[i], par3[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
    }
  m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  for( int i=0; i<4; i++)
    {
      m_minuit->mnpout(i, name[i], par3[i], err[i], bnd1, bnd2, Err);
    }
  m_x0=par3[0]/sin(par3[1]);
  m_u0=-cos(par3[1])/sin(par3[1]);
  m_y0=par3[2]/sin(par3[3]);
  m_v0=-cos(par3[3])/sin(par3[3]);

  CalcChi2();
  chisqr3 = m_chisqr;

  if(fabs(chisqr1)<fabs(chisqr3)&&chisqr1>0.){
    m_x0=par[0];
    m_u0=tan(par[1]);
    m_y0=par[2];
    m_v0=tan(par[3]);
    CalcChi2();
  }
  if(fabs(chisqr2)<fabs(chisqr3)&&fabs(chisqr1)>fabs(chisqr2)&&chisqr2>0.){
    m_x0=par2[0]/cos(par2[1]);
    m_u0=tan(par2[1]);
    m_y0=par2[2]/cos(par2[3]);
    m_v0=tan(par2[3]);
    CalcChi2();
  }


  //  std::cout<<"m_chisqr="<<m_chisqr<<", chisqr1="<<chisqr1<<", chisqr2="<<chisqr2<<std::endl;
  if(m_chisqr > MaxChisqr||std::isnan(m_chisqr))
    return false;

  int false_layer =0;

  for( std::size_t i=0; i<m_hit_array.size(); ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    if(!Residual_check(pos,Res)){
      // std::cout<<"false layer:"<<i<<", n="<<m_hit_array.size()<<std::endl;
      // std::cout<<"m_u0:"<<m_u0<<", m_v0="<<m_v0<<", m_chisqr="<<m_chisqr<<std::endl;
      m_hit_array.erase( m_hit_array.begin()+i );
      ++false_layer;
      --i;
    }
    if(m_hit_array.size()<MinHits)
      return false;
  }


  if(false_layer ==0){
    return true;
  }
  else
    return DoLinearFit(MinHits);
}

//______________________________________________________________________________
bool
TPCLocalTrack::Residual_check(TVector3 pos, TVector3  Res)
{
  bool status_rescheck=false;
  // TVector3 x0(m_x0, m_y0, 0. );
  // TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, 1. );
  TVector3 x0(m_x0, m_y0, zTgtTPC );
  TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, zTgtTPC+1. );
  TVector3 u = (x1-x0).Unit();
  TVector3 AP = pos-x0;
  double dist_AX = u.Dot(AP);
  TVector3 AI(x0.x()+(u.x()*dist_AX),
	      x0.y()+(u.y()*dist_AX),
	      x0.z()+(u.z()*dist_AX));
  TVector3 d = pos-AI;

  if(d.Mag()<Res.Mag()*5.)
    status_rescheck = true;

  return status_rescheck;
}


//______________________________________________________________________________
void
TPCLocalTrack::CalcChi2( void )
{
  double chisqr=0.0;
  int dof = 0;

  const std::size_t n = m_hit_array.size();

  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();

    // TVector3 x0(m_x0, m_y0, 0. );
    // TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, 1. );
    TVector3 x0(m_x0, m_y0, zTgtTPC );
    TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, zTgtTPC+1. );
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
