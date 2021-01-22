// -*- C++ -*-
//Comment by Ichikawa
//TPCLocalTrack_Helix.cc is for Helix fit
//TPCLocalTrack.cc is for lenear fit 
#include "TPCLocalTrack_Helix.hh"

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
#include "Utility_Helix.hh"

#include "TMath.h"
#include "TROOT.h"

//#define HSMagnetON 1 

static int gNumOfHits;
static TVector3 gHitPos[300];
static TVector3 gRes[300];

namespace
{
  const std::string& class_name("TPCLocalTrack_Helix");
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const double& zK18tgt = gGeom.LocalZ("K18Target");
  const double& zTgt    = gGeom.LocalZ("Target");
  // Temporary
  const double& zTgtTPC    = -143.;  
  //  const int MaxTry    = 10;  

  const int ReservedNumOfHits  = 64;
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  static const double  FitStep[5] = { 1.0e-20, 1.0e-21, 1.0e-22, 1.0e-20, 1.0e-20 };
  static const double  LowLimit[5] = { -500., -2*TMath::Pi(), -0.5, -300., -10. };
  static const double  UpLimit[5] = { 500., 2*TMath::Pi(), 0.5, 300., 10. };
  static const double  MaxChisqr = 500.;

  const int    theta_ndiv = 100;
  const double theta_min  =   0;
  const double theta_max  = 180;
  const int    r_ndiv =  100;
  const double r_min  = -500;
  const double r_max  =  500;
}

//______________________________________________________________________________
TPCLocalTrack_Helix::TPCLocalTrack_Helix( void )
  : m_is_fitted(false),
    m_is_calculated(false),
    m_Adrho(0.), m_Aphi0(0.), m_Arho(0.), m_Adz(0.), m_AtanL(0.), 
    m_drho(0.), m_phi0(0.), m_rho(0.), m_dz(0.), m_tanL(0.), 
    m_chisqr(1.e+10),
    m_good_for_tracking(true),
    m_n_iteration(0),
    minuit(new TMinuit(5)),
    m_de(0.)
{
  m_hit_array.reserve( ReservedNumOfHits );
  m_cluster_array.reserve( ReservedNumOfHits );
  debug::ObjectCounter::increase(class_name);
  TROOT minexam("LinearFit","linear fit using TMinuit");
}

//______________________________________________________________________________
TPCLocalTrack_Helix::~TPCLocalTrack_Helix( void )
{
  debug::ObjectCounter::decrease(class_name);
  delete minuit;
}

//______________________________________________________________________________
void
TPCLocalTrack_Helix::ClearHits( void )
{
  m_hit_array.clear();
  m_cluster_array.clear();
}


//______________________________________________________________________________
void
TPCLocalTrack_Helix::AddTPCHit( TPCLTrackHit *hit )
{
  if( hit )
    m_hit_array.push_back( hit ); 
}

//______________________________________________________________________________
void
TPCLocalTrack_Helix::AddTPCCluster( TPCCluster *cluster )
{
  //not supported
  if( cluster )
    m_cluster_array.push_back( cluster );
}

//______________________________________________________________________________
void
TPCLocalTrack_Helix::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsCalculated() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return;
  }
  ////////// Need to change!!!!//////
  const std::size_t n = m_hit_array.size();
  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    hitp->SetCalHelix(m_drho, m_phi0, m_rho, m_dz, m_tanL);
    hitp->SetCalPosition( hitp->GetLocalCalPos_Helix() );
  }
  m_is_calculated = true;
}


//______________________________________________________________________________
int
TPCLocalTrack_Helix::GetNDF( void ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const std::size_t n = m_hit_array.size();
  int ndf = 0;
  for( std::size_t i=0; i<n; ++i ){
    if( m_hit_array[i] ) ++ndf;
  }
  return ndf-5;
}


//______________________________________________________________________________
static void fcn2(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;

  for( int i=0; i<gNumOfHits; ++i ){
    TVector3 pos(-gHitPos[i].X(), 
		 gHitPos[i].Z()-zTgtTPC, 
		 gHitPos[i].Y());
    
    TVector3 fittmp;
    if( !Utility_Helix::PointToHelix(pos, fittmp, par) )
      {
	chisqr+=9999;
	continue;
      }
    TVector3 fittmp_(-fittmp.X(),
		     fittmp.Z(),
		     fittmp.Y()+zTgtTPC);
    
    //    TVector3 d = pos - fittmp; 
    TVector3 d = gHitPos[i] - fittmp_; 
    TVector3 Res = gRes[i];
    
    chisqr += pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    dof++;
  }

  f = chisqr/(double)(dof-5);
}

//______________________________________________________________________________
TPCLTrackHit*
TPCLocalTrack_Helix::GetHit( std::size_t nth ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( nth<m_hit_array.size() )
    return m_hit_array[nth];
  else
    return 0;
}

//______________________________________________________________________________
void
TPCLocalTrack_Helix::DeleteNullHit( void )
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
TPCLocalTrack_Helix::DoFit( int MinHits)
{

  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( IsFitted() ){
    hddaq::cerr << "#W " << func_name << " "
		<< "already called" << std::endl;
    return false;
  }
  bool status_dofit =  DoHelixFit(MinHits);
  m_is_fitted = status_dofit;
  //  m_is_fitted = true;

  return status_dofit;
}





//______________________________________________________________________________
bool
TPCLocalTrack_Helix::DoHelixFit( int MinHits )
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
  double tmp_par[5]={0};
  tmp_par[0]=m_Adrho;
  tmp_par[1]=m_Aphi0;
  tmp_par[2]=m_Arho;
  
  //hough translation for ini-param of Adz and AtanL
  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    gHitPos[i] =pos;
    gRes[i] = Res;
    for( int ti=0; ti<theta_ndiv; ti++ ){
      double theta = theta_min+ti*(theta_max-theta_min)/theta_ndiv;
      // hist.Fill(theta, cos(theta*acos(-1)/180.)*pos.Z()
      // 		 +sin(theta*acos(-1)/180.)*pos.Y());
      double tmpx = -pos.X();
      double tmpy = pos.Z()-zTgtTPC;
      double tmpz = pos.Y();
      double tmp_helixphi = Utility_Helix::CalcHelixPhi(tmpx, tmpy, tmp_par);
      hist.Fill(theta, cos(theta*acos(-1)/180.)*tmp_helixphi
		+sin(theta*acos(-1)/180.)*tmpz);
   }
  }
  int maxbin = hist.GetMaximumBin();
  int mx,my,mz;
  hist.GetBinXYZ( maxbin, mx, my, mz );
  double mtheta = hist.GetXaxis()->GetBinCenter(mx)*acos(-1)/180.;
  double mr = hist.GetYaxis()->GetBinCenter(my);
  double p0 = mr/sin(mtheta);
  double p1 = -cos(mtheta)/sin(mtheta);  
  
  double m_Adz = p0;
  double m_AtanL = -1.*p1*m_Arho;

  double par[5]={m_Adrho, m_Aphi0, m_Arho, m_Adz, m_AtanL};
  double err[5]={-999.,-999.,-999.,-999.,-999.};


  // minuit = new TMinuit(4);
  // TROOT minexam("LinearFit", "Linear fit using TMinuit");

  minuit->SetPrintLevel(-1);
  minuit->SetFCN (fcn2 );

  int ierflg = 0;
  double arglist[10];
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings                                     
  TString name[5] = {"drho", "phi0", "rho", "dz", "tanL"};
  for( int i = 0; i<5; i++ )
    {
      minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
    }

  minuit->Command("SET STRategy 0");
  arglist[0] = 1000.;
  arglist[1] = 1.0;
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  //minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //minuit->mnexcm("SET ERR", arglist, 2, ierflg);
 
  double amin, edm, errdef;
  int nvpar, nparx, icstat;
  minuit->mnstat( amin, edm, errdef, nvpar, nparx, icstat);
  //minuit->mnprin(4, amin);
  int Err;
  double bnd1, bnd2;
  for( int i=0; i<5; i++)
    {
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
      //std::cout<<Par[i]<<"  "<<std::endl;
    }


  double chisqr1 = 100000.;
 
  m_drho = par[0];
  m_phi0 = par[1];
  m_rho = par[2];
  m_dz = par[3];
  m_tanL = par[4];
  CalcChi2(); 
  
  mom0 = Utility_Helix::CalcHelixMom(par, 0.);

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
    return DoHelixFit(MinHits);
}

//______________________________________________________________________________
bool
TPCLocalTrack_Helix::Residual_check(TVector3 pos, TVector3  Res)
{
  bool status_rescheck=false;
  // TVector3 x0(m_x0, m_y0, 0. );
  // TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, 1. );

  double par[5] = {m_drho, m_phi0, m_rho, m_dz, m_tanL};
  TVector3 pos_(-pos.X(), 
		pos.Z()-zTgtTPC, 
		pos.Y());
  
  TVector3 fittmp;
  if( !Utility_Helix::PointToHelix(pos_, fittmp, par) )
    {
      std::cout<<"Residual_Check: No points"<<std::endl;
      return false;
    }
  TVector3 fittmp_(-fittmp.X(),
		   fittmp.Z(),
		   fittmp.Y()+zTgtTPC);
    
  TVector3 d = pos - fittmp_; 
  
  if(d.Mag()<Res.Mag()*5.)
    status_rescheck = true;
  
  return status_rescheck;
}


//______________________________________________________________________________
void
TPCLocalTrack_Helix::CalcChi2( void )
{
  double chisqr=0.0;
  int dof = 0;
  double par[5] = {m_drho, m_phi0, m_rho, m_dz, m_tanL};

  const std::size_t n = m_hit_array.size();

  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();

    TVector3 pos_(-pos.X(), 
		  pos.Z()-zTgtTPC, 
		  pos.Y());
  
    TVector3 fittmp;
    if( !Utility_Helix::PointToHelix(pos_, fittmp, par) )
      {
	std::cout<<"Residual_Check: No points"<<std::endl;
	chisqr+=9999;
	continue;
      }
    TVector3 fittmp_(-fittmp.X(),
		     fittmp.Z(),
		     fittmp.Y()+zTgtTPC);
    
    TVector3 d = pos - fittmp_; 
    
    
    chisqr += pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    dof++;
  }
  m_chisqr = chisqr/(double)(dof-5);
}

//______________________________________________________________________________
// bool
// TPCLocalTrack_Helix::DoHelixFit( void )
// {
//   static const std::string func_name("["+class_name+"::"+__func__+"()]");

//   return true;
// }



//______________________________________________________________________________
void
TPCLocalTrack_Helix::Print( const std::string& arg, std::ostream& ost ) const
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
