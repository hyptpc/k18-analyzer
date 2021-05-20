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
#include <TF1.h>
#include <std_ostream.hh>

//#include "DCAnalyzer.hh"
#include "DCLTrackHit.hh"
//#include "TPCLTrackHit.hh"
#include "DCGeomMan.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "HodoParamMan.hh"
#include "UserParamMan.hh"
#include "TMinuit.h"

#include "TMath.h"
#include "TROOT.h"

//#define HSMagnetON 1


namespace
{
  static int gNumOfHits;
  static TVector3 gHitPos[1000];
  static TVector3 gRes[1000];
  const std::string& class_name("TPCLocalTrack_Helix");
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const double& zK18tgt = gGeom.LocalZ("K18Target");
  const double& zTgt    = gGeom.LocalZ("Target");
  // Temporary
  const double& zTgtTPC    = -143.;
  const double& HS_field_0 = 0.9860;
  const double& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const double& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");


  const int ReservedNumOfHits  = 64;
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();

  //cx, cy, z0, r, dz
  static const double  FitStep[5] = { 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5 };
  static const double  LowLimit[5] = { -7000., -7000., -7000., 0., -10. };
  static const double  UpLimit[5] = { 7000., 7000., 7000., 7000., 10. };
  //rdiff, theta, z0, r, dz
  static const double  FitStep2[5] = { 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-4, 1.0e-5 };
  static const double  LowLimit2[5] = { -200., -acos(-1), -7000., 0., -10. };
  static const double  UpLimit2[5] = { 200., acos(-1), 7000., 7000., 10. };

  static const double  MaxChisqr = 500.;
  static const int  MaxTryMinuit = 2;

  const int    theta_ndiv = 360;
  const double theta_min  =   0;
  const double theta_max  = 180;
  const int    r_ndiv =  2000;
  const double r_min  = -5000;
  const double r_max  =  5000;

  //for Helix tracking
  //[0]~[4] are the Helix parameters,
  //([5],[6],[7]) = (x, y, z)
  static std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
  //static TF1 fint("fint",s_tmp.c_str(),-10.,10.);
  static TF1 fint("fint",s_tmp.c_str(),-4.,4.);
  int ii =0;
}

//______________________________________________________________________________
TPCLocalTrack_Helix::TPCLocalTrack_Helix( void )
  : m_is_fitted(false),
    m_is_calculated(false),
    m_cx(0.), m_cy(0.), m_z0(0.), m_r(0.), m_dz(0.),
    m_Acx(0.), m_Acy(0.), m_Az0(0.), m_Ar(0.), m_Adz(0.),
    m_chisqr(1.e+10),
    m_good_for_tracking(true),
    m_n_iteration(0),
    //minuit(new TMinuit(5)),
    m_mom0(0.,0.,0.)
{
  m_hit_array.reserve( ReservedNumOfHits );

  debug::ObjectCounter::increase(class_name);
  //TROOT minexam("HelixFit","Helix fit using TMinuit");
}

//______________________________________________________________________________
TPCLocalTrack_Helix::~TPCLocalTrack_Helix( void )
{
  debug::ObjectCounter::decrease(class_name);
  // delete minuit;
}

//______________________________________________________________________________
void
TPCLocalTrack_Helix::ClearHits( void )
{
  m_hit_array.clear();

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
    hitp->SetCalHelix(m_cx, m_cy, m_z0, m_r, m_dz);
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
TVector3
TPCLocalTrack_Helix::GetPosition(double par[5], double t) const
{

  // double  x = p[0] + p[3]*cos(t+theta0);
  // double  y = p[1] + p[3]*sin(t+theta0);
  // double  z = p[2] + (p[4]*p[3]*t);
  //This is the eqation of Helix
  double  x = par[0] + par[3]*cos(t);
  double  y = par[1] + par[3]*sin(t);
  double  z = par[2] + (par[4]*par[3]*t);

  return TVector3(x, y, z);
}


//______________________________________________________________________________
TVector3
TPCLocalTrack_Helix::CalcHelixMom(double par[5], double y) const
{

  const double Const = 0.299792458; // =c/10^9
  //const double dMagneticField = 1.; //T, "-1" is needed. // Should be given by field param
  const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  double t = (y-par[2])/(par[3]*par[4]);
  double pt = fabs(par[3])*(Const*dMagneticField); // MeV/c
  //From here!!!!
  double tmp_px = pt*(-1.*sin(t));
  double tmp_py = pt*(cos(t));
  double tmp_pz = pt*(par[4]);
  double px = -tmp_px;
  double py = tmp_pz;
  double pz = tmp_py;


  return TVector3(px,py,pz);

}


//______________________________________________________________________________
//static void fcn2(int &npar, double *gin, double &f, double *par, int iflag)
static inline void fcn2(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;

  double fpar[8];
  //std::cout<<"paramter in fcn"<<std::endl;
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
    //std::cout<<"par ["<<ip<<"]: "<<par[ip]<<std::endl;
  }

  for( int i=0; i<gNumOfHits; ++i ){
    TVector3 pos(-gHitPos[i].X(),
		 gHitPos[i].Z()-zTgtTPC,
		 gHitPos[i].Y());
    fpar[5] = pos.X();
    fpar[6] = pos.Y();
    fpar[7] = pos.Z();

    fint.SetParameters(fpar);
    double min_t = fint.GetMinimumX();
    double  x = par[0] + par[3]*cos(min_t);
    double  y = par[1] + par[3]*sin(min_t);
    double  z = par[2] + (par[4]*par[3]*min_t);

    TVector3 fittmp(x, y, z);
    TVector3 fittmp_(-1.*fittmp.X(),
		     fittmp.Z(),
		     fittmp.Y()+zTgtTPC);


    double tmp_t = atan2(pos.Y()-par[1], pos.X()-par[0]);
    double  tmpx = par[0] + par[3]*cos(tmp_t);
    double  tmpy = par[1] + par[3]*sin(tmp_t);
    double  tmpz = par[2] + (par[4]*par[3]*tmp_t);


    TVector3 fittmp2_(-tmpx,
    		      tmpz,
    		      tmpy+zTgtTPC);

    TVector3 d = gHitPos[i] - fittmp_;
    TVector3 Res = gRes[i];



    chisqr += pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    dof++;

    //double tmp_chisqr = pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    // if(ii==0&&tmp_chisqr>500.){
    //   std::cout<<"tmp_chisqr: "<<tmp_chisqr<<std::endl;
    //   std::cout<<"gHitPos: "<<gHitPos[i]<<", fittmp_: "<<fittmp_
    // 	       <<", fittmp2_"<<fittmp2_
    // 	       <<", min_t="<<min_t<<", tmp_t="<<tmp_t<<std::endl;
    // }
  }
  f = chisqr/(double)(dof-5);
  //std::cout<<"chisqr in fcn: "<<f<<std::endl;
  // if(ii==0)
  //   std::cout<<"chisqr in fcn: "<<f<<std::endl;
  // //   getchar();
  // ++ii;
}

/*
static inline void fcn2_rd(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;

  double fpar[8];
  //std::cout<<"paramter in fcn"<<std::endl;
  double cx = (par[3]+par[0])*cos(par[1]);
  double cy = (par[3]+par[0])*sin(par[1]);
  par[0] = cx;
  par[1] = cy;
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
    //std::cout<<"par ["<<ip<<"]: "<<par[ip]<<std::endl;
  }

  for( int i=0; i<gNumOfHits; ++i ){
    TVector3 pos(-gHitPos[i].X(),
		 gHitPos[i].Z()-zTgtTPC,
		 gHitPos[i].Y());
    fpar[5] = pos.X();
    fpar[6] = pos.Y();
    fpar[7] = pos.Z();

    fint.SetParameters(fpar);
    double min_t = fint.GetMinimumX();
    double  x = par[0] + par[3]*cos(min_t);
    double  y = par[1] + par[3]*sin(min_t);
    double  z = par[2] + (par[4]*par[3]*min_t);

    TVector3 fittmp(x, y, z);
    TVector3 fittmp_(-1.*fittmp.X(),
		     fittmp.Z(),
		     fittmp.Y()+zTgtTPC);


    double tmp_t = atan2(pos.Y()-par[1], pos.X()-par[0]);
    double  tmpx = par[0] + par[3]*cos(tmp_t);
    double  tmpy = par[1] + par[3]*sin(tmp_t);
    double  tmpz = par[2] + (par[4]*par[3]*tmp_t);


    TVector3 fittmp2_(-tmpx,
    		      tmpz,
    		      tmpy+zTgtTPC);

    TVector3 d = gHitPos[i] - fittmp_;
    TVector3 Res = gRes[i];



    chisqr += pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    dof++;

    //double tmp_chisqr = pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    // if(ii==0&&tmp_chisqr>500.){
    //   std::cout<<"tmp_chisqr: "<<tmp_chisqr<<std::endl;
    //   std::cout<<"gHitPos: "<<gHitPos[i]<<", fittmp_: "<<fittmp_
    // 	       <<", fittmp2_"<<fittmp2_
    // 	       <<", min_t="<<min_t<<", tmp_t="<<tmp_t<<std::endl;
    // }
  }
  f = chisqr/(double)(dof-5);
  //std::cout<<"chisqr in fcn: "<<f<<std::endl;
  // if(ii==0)
  //   std::cout<<"chisqr in fcn: "<<f<<std::endl;
  // //   getchar();
  // ++ii;
}
*/

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
  // static TH2D hist("hist",";theta (deg.); r (mm)",
  //  		   theta_ndiv, theta_min, theta_max, r_ndiv, r_min, r_max);
  // static TH2D hist("hist",";theta (deg.); r (mm)",
  //   		   theta_ndiv, theta_min, theta_max, r_ndiv, r_min, r_max);
  TH2D *hist = new TH2D("hist",";theta (deg.); r (mm)",
			theta_ndiv, theta_min, theta_max, r_ndiv, r_min, r_max);
  //hist.Reset();
  //hough translation for ini-param of Adz and AtanL


  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    gHitPos[i] =pos;
    gRes[i] = Res;

    for( int ti=0; ti<theta_ndiv; ti++ ){
      double theta = theta_min+ti*(theta_max-theta_min)/theta_ndiv;
      double tmpx = -pos.X();
      double tmpy = pos.Z()-zTgtTPC;
      double tmpz = pos.Y();
      double tmp_t = atan2(tmpy-m_Acy, tmpx-m_Acx);
      double tmp_xval = m_Ar * tmp_t;
      hist->Fill(theta, cos(theta*acos(-1.)/180.)*tmp_xval
		 +sin(theta*acos(-1.)/180.)*tmpz);
    }
  }
  int maxbin = hist->GetMaximumBin();
  int mx,my,mz;
  hist->GetBinXYZ( maxbin, mx, my, mz );

  double mtheta = hist->GetXaxis()->GetBinCenter(mx)*acos(-1.)/180.;
  double mr = hist->GetYaxis()->GetBinCenter(my);

  double p0 = mr/sin(mtheta);
  double p1 = -cos(mtheta)/sin(mtheta);
  //  std::cout<<"MaxBin(y,theta): "<<hist.GetMaximum()<<std::endl;
  // std::cout<<"mr: "<<mr<<", mtheta: "<<mtheta<<std::endl;
  // std::cout<<"p0: "<<p0<<", p1: "<<p1<<", m_Ar:"<<m_Ar<<std::endl;
  double m_Az0 = p0;
  double m_Adz = p1;

  delete hist;
  //  std::cout<<"Hough check"<<std::endl;
  // for( std::size_t i=0; i<n; ++i ){
  //   TPCLTrackHit *hitp = m_hit_array[i];
  //   TVector3 pos = hitp->GetLocalHitPos();
  //   TVector3 Res = hitp->GetResolutionVect();
  //   double tmpx = -pos.X();
  //   double tmpy = pos.Z()-zTgtTPC;
  //   double tmpz = pos.Y();
  //   double tmp_t = atan2(tmpy-m_Acy, tmpx-m_Acx);
  //   double calz = p0 + p1*m_Ar*tmp_t;
  //   std::cout<<"posz:"<<tmpz<<", calz:"<<calz<<", tmp_t:"<<tmp_t<<std::endl;
  // }



  double par[5]={m_Acx, m_Acy, m_Az0, m_Ar, m_Adz};
  // double par_rd = sqrt(m_Acx*m_Acx + m_Acy*m_Acy) - m_Ar;
  // double par_theta = atan2(m_Acy, m_Acx);
  // double par2[5]={par_rd, par_theta, m_Az0, m_Ar, m_Adz};
  double err[5]={-999.,-999.,-999.,-999.,-999.};

  TMinuit *minuit = new TMinuit(5);
  //static TMinuit minuit(5);
  //TROOT minexam("HelixFit", "Helix fit using TMinuit");

  ii = 0;

  minuit->SetPrintLevel(-1);
  minuit->SetFCN( fcn2 );

  int ierflg = 0;
  double arglist[10];
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings
  double min_chisqr = 1.e+10;
  double min_par[5]={0};

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};

  for( int i = 0; i<5; i++ )
    {
      minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
    }

  minuit->Command("SET STRategy 0");
  // arglist[0] = 5000.;
  // arglist[1] = 0.01;


  arglist[0] = 1000;
  arglist[1] = 0.1;
  int itry=0;
  while(m_chisqr>1.5){

    if(itry>MaxTryMinuit)
      break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    // double amin, edm, errdef;
    // int nvpar, nparx, icstat;
    // minuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat);
    //minuit->mnprin(4, amin);
    int Err;
    double bnd1, bnd2;
    for( int i=0; i<5; i++)
      {
	minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
	//std::cout<<Par[i]<<"  "<<std::endl;
      }


    m_cx = par[0];
    m_cy = par[1];
    m_z0 = par[2];
    m_r  = par[3];
    m_dz = par[4];
    CalcChi2();
    if(min_chisqr>m_chisqr){
      min_chisqr = m_chisqr;
      min_par[0] = m_cx;
      min_par[1] = m_cy;
      min_par[2] = m_z0;
      min_par[3] = m_r;
      min_par[4] = m_dz;
    }


    // std::cout<<"while, itry: "<<itry
    // 	     <<", arglist0: "<<arglist[0]
    // 	     <<", 1: "<<arglist[1]
    // 	     <<", chisqr: "<<m_chisqr<<std::endl;
    ++itry;
    arglist[0] = arglist[0]*5;
    arglist[1] = arglist[1]/10.;
  }



  /*
  minuit->SetFCN( fcn2_rd );
  TString name2[5] = {"rd", "theta", "z0", "r", "dz"};

  for( int i = 0; i<5; i++ )
    {
      minuit->mnparm(i, name2[i], par2[i], FitStep2[i], LowLimit2[i], UpLimit2[i], ierflg);
    }
  arglist[0] = 1000;
  arglist[1] = 0.1;

  int itry2=0;

  while(m_chisqr>1.5){

    if(itry2>MaxTryMinuit)
      break;

    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    // double amin, edm, errdef;
    // int nvpar, nparx, icstat;
    // minuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat);
    //minuit->mnprin(4, amin);
    int Err;
    double bnd1, bnd2;
    for( int i=0; i<5; i++)
      {
	minuit->mnpout(i, name2[i], par2[i], err[i], bnd1, bnd2, Err);
	//std::cout<<Par[i]<<"  "<<std::endl;
      }


    m_cx = (par2[3]+par2[0])*cos(par2[1]);
    m_cy = (par2[3]+par2[0])*sin(par2[1]);
    m_z0 = par2[2];
    m_r  = par2[3];
    m_dz = par2[4];

    CalcChi2();

    if(min_chisqr>m_chisqr){
      min_chisqr = m_chisqr;
      min_par[0] = m_cx;
      min_par[1] = m_cy;
      min_par[2] = m_z0;
      min_par[3] = m_r;
      min_par[4] = m_dz;
    }
    std::cout<<"while2, itry: "<<itry
	     <<", arglist0: "<<arglist[0]
	     <<", 1: "<<arglist[1]
	     <<", chisqr: "<<m_chisqr<<std::endl;
    ++itry2;
    arglist[0] = arglist[0]*5;
    arglist[1] = arglist[1]/10.;
  }
  */

  m_chisqr = min_chisqr;
  m_cx =   min_par[0];
  m_cy =   min_par[1];
  m_z0 =   min_par[2];
  m_r  =   min_par[3];
  m_dz =   min_par[4];


  m_mom0 = CalcHelixMom(par, 0.);
  //std::cout<<"m_chisqr: "<<m_chisqr<<std::endl;

  delete  minuit;

  if(m_chisqr > MaxChisqr||std::isnan(m_chisqr))
    return false;

  int false_layer =0;

  for( std::size_t i=0; i<m_hit_array.size(); ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    if(!Residual_check(pos,Res)){
      //std::cout<<"false layer:"<<i<<", n="<<m_hit_array.size()<<std::endl;
      m_hit_array.erase( m_hit_array.begin()+i );
      ++false_layer;
      --i;
    }
    if(m_hit_array.size()<MinHits)
      return false;
  }
  //getchar();

  if(false_layer ==0){
    return true;
  }
  else
    return DoHelixFit(MinHits);


  return true;
}

//______________________________________________________________________________
bool
TPCLocalTrack_Helix::Residual_check(TVector3 pos, TVector3  Res)
{

  bool status_rescheck=false;

  double par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 pos_(-pos.X(),
		pos.Z()-zTgtTPC,
		pos.Y());
  double fpar[8];
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
  }
  fpar[5] = pos_.X();
  fpar[6] = pos_.Y();
  fpar[7] = pos_.Z();

  fint.SetParameters(fpar);
  double min_t = fint.GetMinimumX();
  TVector3 fittmp = GetPosition(par, min_t);
  TVector3 fittmp_(-fittmp.X(),
		   fittmp.Z(),
		   fittmp.Y()+zTgtTPC);


  TVector3 d = pos - fittmp_;

  if(d.Mag()<Res.Mag()*5.)
    status_rescheck = true;
  // else
  //   std::cout<<"Residual: "<<d.Mag()<<std::endl;

  return status_rescheck;
}


//______________________________________________________________________________
void
TPCLocalTrack_Helix::CalcChi2( void )
{

  double chisqr=0.0;
  int dof = 0;
  double par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  double fpar[8];
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
  }

  const std::size_t n = m_hit_array.size();

  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    TVector3 pos_(-pos.X(),
		  pos.Z()-zTgtTPC,
		  pos.Y());

    fpar[5] = pos_.X();
    fpar[6] = pos_.Y();
    fpar[7] = pos_.Z();
    fint.SetParameters(fpar);
    double min_t = fint.GetMinimumX();
    TVector3 fittmp = GetPosition(par, min_t);
    TVector3 fittmp_(-fittmp.X(),
		     fittmp.Z(),
		     fittmp.Y()+zTgtTPC);
    TVector3 d = pos - fittmp_;

    chisqr += pow( d.x()/Res.x(), 2) + pow( d.y()/Res.y(), 2) + pow( d.z()/Res.z(), 2);
    dof++;
  }
  m_chisqr = chisqr/(double)(dof-5);

}
