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
#include "TF2.h"

#include "TMath.h"
#include "TROOT.h"

//#define HSMagnetON 1


namespace
{
  static int gNumOfHits;
  // static TVector3 gHitPos[1000];
  // static TVector3 gRes[1000];
  static std::vector<TVector3> gHitPos;
  static std::vector<TVector3> gRes;
  const std::string& class_name("TPCLocalTrack_Helix");
  const UserParamMan& gUser = UserParamMan::GetInstance();
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
  const double& zK18tgt = gGeom.LocalZ("K18Target");
  const double& zTgt    = gGeom.LocalZ("Target");
  // Temporary
  const double& zTgtTPC    = -143.;
  const double& r_HTOF    = 337.;
  const double& HS_field_0 = 0.9860;
  const double& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const double& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

  const int ReservedNumOfHits  = 32*4;
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();

  //cx, cy, z0, r, dz
  static const double  FitStep[5] = { 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5 };
  // static const double  LowLimit[5] = { -7000., -7000., -7000., 0., -10. };
  // static const double  UpLimit[5] = { 7000., 7000., 7000., 7000., 10. };
  static const double  LowLimit[5] = { -20000., -20000., -700., 0., -10. };
  static const double  UpLimit[5] = { 20000., 20000., 700., 20000., 10. };
  //rdiff, theta, z0, r, dz
  static const double  FitStep2[5] = { 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-4, 1.0e-5 };
  static const double  LowLimit2[5] = { -200., -acos(-1), -7000., 0., -10. };
  static const double  UpLimit2[5] = { 200., acos(-1), 7000., 20000., 10. };

  static const double  MaxChisqr = 500.;
  //  static const int  MaxTryMinuit = 3;
  static const int  MaxTryMinuit = 0;

  const int    theta_ndiv = 360;
  const double theta_min  =   0;
  const double theta_max  = 180;
  const int    r_ndiv =  2000;
  const double r_min  = -5000;
  const double r_max  =  5000;

  double min_par[5]={0};
  //for Helix tracking
  //[0]~[4] are the Helix parameters,
  //([5],[6],[7]) = (x, y, z)
  static std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
  //static TF1 fint("fint",s_tmp.c_str(),-10.,10.);
  static TF1 fint("fint",s_tmp.c_str(),-4.,4.);

  //for circle tracking
  //[0]~[2] are the circle parameters,
  //([3],[4]) = (x, y)
  static std::string s_tmp_circ="pow([3]-([0]+([2]*cos(x))),2)+pow([4]-([1]+([2]*sin(x))),2)";
  static TF1 fint_circ("fint_circ",s_tmp_circ.c_str(),-acos(-1.),acos(-1.));

  //for linear
  //[0],[1],[2] are the Linear parameters,
  //([3],[4]) = (t,z)
  static std::string s_tmp_li="pow([3]-x,2)+pow([4]-([0]+([1]*[2]*x)),2)";
  static TF1 fint_li("fint_li",s_tmp_li.c_str(),-4.,4.);

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
  double px = -tmp_px*0.001;
  double py = tmp_pz*0.001;
  double pz = tmp_py*0.001;


  return TVector3(px,py,pz);

}


//______________________________________________________________________________
int
TPCLocalTrack_Helix::GetHTOFSeg(double min_layer_t, double max_layer_t, double max_layer_y){

  //circle function 1 (tracking)
  //x = [0] + [2]*cos(t);
  //y = [1] + [2]*sin(t);

  //circle function 2 (HTOF position)
  //x = [3] + [5]*cos(t);
  //y = [4] + [5]*sin(t);
  
  double min_t=0., max_t=0.;
  if(min_layer_t<max_layer_t){
    min_t = min_layer_t;
    max_t = max_layer_t + 0.5;
  }

  else{
    min_t = max_layer_t-0.5;
    max_t = min_layer_t;
  }
  
  static TF2 fhtof("fhtof",
		   "pow(([0]+[2]*cos(x))-([3]+[5]*cos(y)),2)+pow(([1]+[2]*sin(x))-([4]+[5]*sin(y)),2)",
		   min_t, max_t, 0., 2.*acos(-1.));
  fhtof.SetParameter(0, m_cx);
  fhtof.SetParameter(1, m_cy);
  fhtof.SetParameter(2, m_r);
  fhtof.SetParameter(3, 0.);
  fhtof.SetParameter(4, -1.*zTgtTPC);
  fhtof.SetParameter(5, r_HTOF);
  double close_t1, close_t2;
  fhtof.GetMinimumXY(close_t1, close_t2);
  double seg_deg = 2.*acos(-1.)/32.;
  int htofseg =0;
  if(close_t2<22.*seg_deg){
    for(int ihtof =0; ihtof<22; ++ihtof){
      if((double)ihtof*seg_deg<close_t2&&
	 (double)(ihtof+1)*seg_deg>close_t2)
	htofseg = ihtof+12;
    }
  }
  else if(22.*seg_deg<close_t2&&
	  23.*seg_deg>close_t2)
    htofseg = 0;
  else if(23.*seg_deg<close_t2&&
	  24.*seg_deg>close_t2){
    if(max_layer_y >0.)
      htofseg = 1;
    else
      htofseg = 2;
  }
  else if(24.*seg_deg<close_t2&&
	  25.*seg_deg>close_t2){
    if(max_layer_y >0.)
      htofseg = 3;
    else
      htofseg = 4;
  }
  else{
    for(int ihtof =25; ihtof<32; ++ihtof){
      if((double)ihtof*seg_deg<close_t2&&
	 (double)(ihtof+1)*seg_deg>close_t2)
	htofseg = ihtof-20;
    }
  }
  return htofseg;
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

  }
  f = chisqr/(double)(dof-5);
}

//______________________________________________________________________________
static inline void fcn2_circ(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;

  double fpar[5];
  //std::cout<<"paramter in fcn"<<std::endl;
  for(int ip=0; ip<3; ++ip){
    fpar[ip] = par[ip];
    //std::cout<<"par ["<<ip<<"]: "<<par[ip]<<std::endl;
  }

  for( int i=0; i<gNumOfHits; ++i ){
    TVector2 pos_(gHitPos[i].X(),
		  gHitPos[i].Z());

    TVector2 pos(-gHitPos[i].X(),
		 gHitPos[i].Z()-zTgtTPC);
    fpar[3] = pos.X();
    fpar[4] = pos.Y();

    fint_circ.SetParameters(fpar);
    double min_t = fint_circ.GetMinimumX();
    double  x = par[0] + par[2]*cos(min_t);
    double  y = par[1] + par[2]*sin(min_t);
 
    TVector2 fittmp(x, y);
    TVector2 fittmp_(-1.*fittmp.X(),
		     fittmp.Y()+zTgtTPC);
    // double tmp_t = atan2(pos.Y()-par[1], pos.X()-par[0]);
    // double  tmpx = par[0] + par[3]*cos(tmp_t);
    // double  tmpy = par[1] + par[3]*sin(tmp_t);
    // TVector2 fittmp2_(-tmpx,
    // 		      tmpy+zTgtTPC);
    TVector2 d = pos_ - fittmp_;

    chisqr += pow( d.X()/gRes[i].x(), 2) + pow( d.Y()/gRes[i].z(), 2);
    dof++;
  }
  f = chisqr/(double)(dof-3);
  //  std::cout<<"f_circ:"<<f<<std::endl;
}


//______________________________________________________________________________
static inline void fcn2_li(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;
  
  double fpar[5]={0};
  //std::cout<<"paramter in fcn"<<std::endl;
  for(int ip=0; ip<2; ++ip){
    fpar[ip] = par[ip];
    //std::cout<<"par ["<<ip<<"]: "<<par[ip]<<std::endl;
  }

  for( int i=0; i<gNumOfHits; ++i ){
    TVector3 pos(-gHitPos[i].X(),
		 gHitPos[i].Z()-zTgtTPC,
		 gHitPos[i].Y());
    double tmpx = pos.X();
    double tmpy = pos.Y();
    //      double tmp_t = atan2(tmpy-m_Acy, tmpx-m_Acx);
    double tmp_t = atan2(tmpy - min_par[1], 
			 tmpx - min_par[0]);
    fpar[2] = min_par[3];
    fpar[3] = tmp_t;
    fpar[4] = pos.Z();
    fint_li.SetParameters(fpar);
    
    //double diff_z = fint_li.GetMinimum();
    double diff_z = pow(fpar[4]-(fpar[0]+(fpar[1]*fpar[2]*fpar[3])),2);
    

    chisqr +=  diff_z/pow(gRes[i].y(),2);
    dof++;
  }
  f = chisqr/(double)(dof-2);
  //std::cout<<"f_li:"<<f<<std::endl;
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

  if(n>ReservedNumOfHits)
    return false;
  
  gNumOfHits = n;
  gHitPos.clear();
  gRes.clear();
  
  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    gHitPos.push_back(pos);
    gRes.push_back(Res);
  }

  //pre circle fit
  double par_circ[3]={m_Acx, m_Acy, m_Ar};
  double err_circ[3]={-999.,-999.,-999.};
  TMinuit *minuit_circ = new TMinuit(3);
  minuit_circ->SetPrintLevel(-1);
  minuit_circ->SetFCN( fcn2_circ );

  int ierflg_circ = 0;
  double arglist_circ[10];
  arglist_circ[0] = 3.52;//for 3 parameters
  minuit_circ->mnexcm("SET ERR", arglist_circ,1,ierflg_circ);
  arglist_circ[0] = 1;
  minuit_circ->mnexcm("SET NOW", arglist_circ,1,ierflg_circ); // No warnings
  double min_chisqr_circ = 1.e+10;
  double min_par_circ[3]={0};

  TString name_circ[3] = {"cx", "cy","r"};

  for( int i = 0; i<3; i++ ){
    if(i<2)
      minuit_circ->mnparm(i, name_circ[i], par_circ[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg_circ);
    else
      minuit_circ->mnparm(i, name_circ[i], par_circ[i], FitStep[i+1], LowLimit[i+1], UpLimit[i+1], ierflg_circ);
  }
  

  minuit_circ->Command("SET STRategy 0");
 
  arglist_circ[0] = 1000*5*5*5;
  arglist_circ[1] = 0.1/(10.*10.*10.);
  int Err_circ;
  double bnd1_circ, bnd2_circ;

  minuit_circ->mnexcm("MIGRAD", arglist_circ, 2, ierflg_circ);

  for( int i=0; i<3; i++){
      minuit_circ->mnpout(i, name_circ[i], par_circ[i], err_circ[i], bnd1_circ, bnd2_circ, Err_circ);
      //std::cout<<Par[i]<<"  "<<std::endl;
  }
  min_par_circ[0] = par_circ[0];
  min_par_circ[1] = par_circ[1];
  min_par_circ[2] = par_circ[2];
  min_chisqr_circ = CalcChi2_circle(par_circ);
  
  // for invert charge fit
  double tmpx_min = -gHitPos[0].X();
  double tmpy_min = gHitPos[0].Z()-zTgtTPC;
  double tmpt_min = atan2(tmpy_min - par_circ[1], tmpx_min - par_circ[0]);
  double tmpx_max = -gHitPos[n-1].X();
  double tmpy_max = gHitPos[n-1].Z()-zTgtTPC;
  double tmpt_max = atan2(tmpy_max - par_circ[1], tmpx_max - par_circ[0]);
  
  double mid_t_circ = (tmpt_min + tmpt_max)/2.;
  double mid_x_circ = par_circ[0] + par_circ[2]*cos(mid_t_circ);
  double mid_y_circ = par_circ[1] + par_circ[2]*cos(mid_t_circ);
  
  par_circ[0] = min_par_circ[0] - 2.*(min_par_circ[0] - mid_x_circ);
  par_circ[1] = min_par_circ[1] - 2.*(min_par_circ[1] - mid_y_circ);

  for( int i = 0; i<3; i++ ){
    if(i<2)
      minuit_circ->mnparm(i, name_circ[i], par_circ[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg_circ);
    else
      minuit_circ->mnparm(i, name_circ[i], par_circ[i], FitStep[i+1], LowLimit[i+1], UpLimit[i+1], ierflg_circ);
  }
  minuit_circ->mnexcm("MIGRAD", arglist_circ, 2, ierflg_circ);
  
  for( int i=0; i<3; i++){
    minuit_circ->mnpout(i, name_circ[i], par_circ[i], err_circ[i], bnd1_circ, bnd2_circ, Err_circ);
  }

  if(min_chisqr_circ>CalcChi2_circle(par_circ)){
    min_par_circ[0] = par_circ[0];
    min_par_circ[1] = par_circ[1];
    min_par_circ[2] = par_circ[2];
    min_chisqr_circ = CalcChi2_circle(par_circ);
  }
  delete  minuit_circ;
  m_cx = min_par_circ[0];
  m_cy = min_par_circ[1];
  m_r = min_par_circ[2];
  min_par[0] = m_cx;
  min_par[1] = m_cy;
  min_par[3] = m_r;

  // std::cout<<"comp circle par:"<<"{"<<m_Acx<<", "<<m_Acy<<", "<<m_Ar<<"}, "
  //  	   <<"{"<<m_cx<<", "<<m_cy<<", "<<m_r<<"}"<<std::endl;
  // std::cout<<"chisqr circle: "<<min_chisqr_circ<<std::endl;
  
  TH2D *hist = new TH2D("hist",";theta (deg.); r (mm)",
			theta_ndiv, theta_min, theta_max, r_ndiv, r_min, r_max);

  // r = x * cos(theta) + y * sin(theta)
  // static TH2D hist("hist",";theta (deg.); r (mm)",
  //  		   theta_ndiv, theta_min, theta_max, r_ndiv, r_min, r_max);
  // static TH2D hist("hist",";theta (deg.); r (mm)",
  //   		   theta_ndiv, theta_min, theta_max, r_ndiv, r_min, r_max);
  //hist.Reset();
  //hough translation for ini-param of Adz and AtanL

  for( std::size_t i=0; i<n; ++i ){
    for( int ti=0; ti<theta_ndiv; ti++ ){
      double theta = theta_min+ti*(theta_max-theta_min)/theta_ndiv;
      double tmpx = -gHitPos[i].X();
      double tmpy = gHitPos[i].Z()-zTgtTPC;
      double tmpz = gHitPos[i].Y();
      //      double tmp_t = atan2(tmpy-m_Acy, tmpx-m_Acx);
      double tmp_t = atan2(tmpy - min_par_circ[1], 
			   tmpx - min_par_circ[0]);
      //double tmp_xval = m_Ar * tmp_t;
      double tmp_xval = min_par_circ[2] * tmp_t;
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

  //pre t-y fit
  double par_li[2]={m_Az0, m_Adz};
  double err_li[2]={-999.,-999.};
  TMinuit *minuit_li = new TMinuit(2);
  minuit_li->SetPrintLevel(-1);
  minuit_li->SetFCN( fcn2_li );

  int ierflg_li = 0;
  double arglist_li[10];
  arglist_li[0] = 2.3;
  minuit_li->mnexcm("SET ERR", arglist_li,1,ierflg_li); // No warnings
  arglist_li[0] = 1;
  minuit_li->mnexcm("SET NOW", arglist_li,1,ierflg_li); // No warnings

  TString name_li[2] = {"z0", "dz"};
  minuit_li->mnparm(0, name_li[0], par_li[0], FitStep[2], LowLimit[2], UpLimit[2], ierflg_li);
  minuit_li->mnparm(1, name_li[1], par_li[1], FitStep[4], LowLimit[4], UpLimit[4], ierflg_li);
  
  minuit_li->Command("SET STRategy 0");
 
  arglist_li[0] = 1000*5*5*5;
  arglist_li[1] = 0.1/(10.*10.*10.);
  int Err_li;
  double bnd1_li, bnd2_li;

  minuit_li->mnexcm("MIGRAD", arglist_li, 2, ierflg_li);

  for( int i=0; i<2; i++){
    minuit_li->mnpout(i, name_li[i], par_li[i], err_li[i], bnd1_li, bnd2_li, Err_li);
  }
  delete minuit_li;
  
  double min_chisqr = 1.e+10;
  
  m_z0 = par_li[0];
  m_dz = par_li[1];
  min_par[2] = m_z0;
  min_par[4] = m_dz;  
  CalcChi2();
  min_chisqr = m_chisqr;
  
  // std::cout<<"comp li par:"<<"{"<<m_Az0<<", "<<m_Adz<<"}, "
  //    	   <<"{"<<m_z0<<", "<<m_dz<<"}"<<std::endl;

  //  std::cout<<"chisqr fist: "<<min_chisqr<<std::endl;


  //double par[5]={m_Acx, m_Acy, m_Az0, m_Ar, m_Adz};
  double par[5]={m_cx, m_cy, m_z0, m_r, m_dz};

  double err[5]={-999.,-999.,-999.,-999.,-999.};

  TMinuit *minuit = new TMinuit(5);


  ii = 0;

  minuit->SetPrintLevel(-1);
  minuit->SetFCN( fcn2 );

  int ierflg = 0;
  double arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist,1,ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};

  for( int i = 0; i<5; i++ ){
    minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }

  minuit->Command("SET STRategy 0");
  // arglist[0] = 5000.;
  // arglist[1] = 0.01;


  arglist[0] = 1000;
  arglist[1] = 0.1;

  arglist[0] = arglist[0]*5*5*5;
  arglist[1] = arglist[1]/(10.*10.*10.);
  int itry=0;
  int Err;
  double bnd1, bnd2;


  while(m_chisqr>1.5){

    if(itry>MaxTryMinuit)
      break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    minuit->mnimpr();
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    // double amin, edm, errdef;
    // int nvpar, nparx, icstat;
    // minuit.mnstat( amin, edm, errdef, nvpar, nparx, icstat);
    //minuit->mnprin(4, amin);

    for( int i=0; i<5; i++){
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

  m_chisqr = min_chisqr;
  m_cx =   min_par[0];
  m_cy =   min_par[1];
  m_z0 =   min_par[2];
  m_r  =   min_par[3];
  m_dz =   min_par[4];


  // for invert charge fit
  double min_layer_t=0., max_layer_t=0.;
  int min_layer =100, max_layer = -1;
  for( std::size_t i=0; i<m_hit_array.size(); ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    double t = GetTcal(pos);
    if(min_layer > hitp->GetLayer()){
      min_layer =  hitp->GetLayer();
      min_layer_t = t;
    }
    if(max_layer < hitp->GetLayer()){
      max_layer =  hitp->GetLayer();
      max_layer_t = t;
    }
  }
  double mid_t = (min_layer_t + max_layer_t)/2.;
  TVector3 mid_pos = GetPosition(min_par, mid_t);
  double mid_x = mid_pos.x();
  double mid_y = mid_pos.y();
  

  par[0] = m_cx - 2.*(m_cx - mid_x);
  par[1] = m_cy - 2.*(m_cy - mid_y);
  par[2] = m_z0;
  par[3] = m_r;
  par[4] = m_dz;  

  for( int i = 0; i<5; i++ ){
    minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  
  for( int i=0; i<5; i++){
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
  m_chisqr = min_chisqr;
  m_cx =   min_par[0];
  m_cy =   min_par[1];
  m_z0 =   min_par[2];
  m_r  =   min_par[3];
  m_dz =   min_par[4];
  
  //  std::cout<<"final_chisqr: "<<min_chisqr<<std::endl;
  m_mom0 = CalcHelixMom(min_par, 0.);
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
double
TPCLocalTrack_Helix::GetTcal(TVector3 pos)
{
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

  return min_t;
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
    //    chisqr += pow( d.x()/Res.x(), 2) + pow( d.z()/Res.z(), 2);
    dof++;
  }
  m_chisqr = chisqr/(double)(dof-5);

}



//______________________________________________________________________________
double
TPCLocalTrack_Helix::CalcChi2_circle( double par[3] )  
{
  double chisqr=0.0;
  int dof = 0;
  double fpar[5];
  for(int ip=0; ip<3; ++ip){
    fpar[ip] = par[ip];
  }

  const std::size_t n = m_hit_array.size();

  for( std::size_t i=0; i<n; ++i ){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector2 pos2(pos.X(), pos.Z());
    TVector3 Res = hitp->GetResolutionVect();
    TVector2 pos_(-pos.X(),
		  pos.Z()-zTgtTPC);

    fpar[3] = pos_.X();
    fpar[4] = pos_.Y();
    fint_circ.SetParameters(fpar);
    double min_t = fint_circ.GetMinimumX();
    double  x = par[0] + par[2]*cos(min_t);
    double  y = par[1] + par[2]*sin(min_t);
 
    TVector2 fittmp(x, y);
    TVector2 fittmp_(-1.*fittmp.X(),
		     fittmp.Y()+zTgtTPC);

    TVector2 d = pos2 - fittmp_;

    chisqr += pow( d.X()/Res.x(), 2) + pow( d.Y()/Res.z(), 2);
    dof++;
  }
  chisqr = chisqr/(double)(dof-3);
  return chisqr;
}
