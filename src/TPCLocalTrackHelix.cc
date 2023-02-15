// -*- C++ -*-

//Comment by Ichikawa
//TPCLocalTrackHelix.cc is for Helix fit
//Pre circle fit (TMinuit -> reduced chi2 method)
#include <chrono>
#include "TPCLocalTrackHelix.hh"

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

#include "TPCPadHelper.hh"
#include "HoughTransform.hh"
#include "ConfMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "TMinuit.h"
#include "TF2.h"

#include "TMath.h"
#include "TROOT.h"

#define DebugDisp        0
#define InvertChargeTest 0

namespace
{
  //for Minimization
  static int gNumOfHits;
  static std::vector<TVector3> gHitPos;
  static std::vector<TVector3> gRes;
  static std::vector<double> gTcal;
  static double gPar[5] = {0};
  static double gChisqr = 1.e+10;
  static const double  MaxChisqr = 500.;
  //  static const double  MaxChisqr = 10000.;
  //static const double  MaxChisqr = 100.;
  //static const int  MaxTryMinuit = 3;
  static const int  MaxTryMinuit = 0;
  static const Double_t MaxResidual = 5.;

  // B-field
  const double& HS_field_0 = 0.9860;
  const double& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
  const double& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

  const int ReservedNumOfHits  = 32*10;

  //cx, cy, z0, r, dz
  static const double  FitStep[5] = { 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-4, 1.0e-5 };
  // static const double  LowLimit[5] = { -7000., -7000., -7000., 0., -10. };
  // static const double  UpLimit[5] = { 7000., 7000., 7000., 7000., 10. };
  // static const double  LowLimit[5] = { -20000., -20000., -7000., 0., -10. };
  // static const double  UpLimit[5] = { 20000., 20000., 7000., 20000., 10. };
  static const double  LowLimitBeam[5] = { -40000., -40000., -7000., 5666., -10. };// about 1.7 GeV/c
  static const double  LowLimit[5] = { -40000., -40000., -7000., 0., -10. };
  static const double  UpLimit[5] = { 40000., 40000., 7000., 40000., 10. };
  // static const double  LowLimit[5] = { -4000., -4000., -700., 0., -10. };
  // static const double  UpLimit[5] = { 4000., 4000., 700., 10000., 10. };
  //static const double  LowLimit[5] = { -100000., -100000., -30000., 0., -10. };
  //static const double  UpLimit[5] = { 100000., 100000., 30000., 100000., 10. };
  //rdiff, theta, z0, r, dz
  static const double  FitStep2[5] = { 1.0e-4, 1.0e-5, 1.0e-4, 1.0e-4, 1.0e-5 };
  //  static const double  FitStep2[5] = { 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-5, 1.0e-6 };
  static const double  LowLimit2[5] = { -200., -acos(-1), -7000., 0., -10. };
  static const double  UpLimit2[5] = { 200., acos(-1), 7000., 20000., 10. };
  // static const double  LowLimit2[5] = { -200., -acos(-1), -700., 0., -10. };
  // static const double  UpLimit2[5] = { 200., acos(-1), 700., 10000., 10. };
  //static const double  LowLimit2[5] = { -1000., -acos(-1), -7000., 0., -10. };
  //static const double  UpLimit2[5] = { 1000., acos(-1), 7000., 100000., 10. };

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

}

//______________________________________________________________________________
TPCLocalTrackHelix::TPCLocalTrackHelix()
  : m_is_fitted(false),
    m_is_calculated(false),
    m_flag(0),
    m_hit_order(), m_hit_t(),
    m_Acx(0.), m_Acy(0.), m_Az0(0.), m_Ar(0.), m_Adz(0.),
    m_cx(0.), m_cy(0.), m_z0(0.), m_r(0.), m_dz(0.),
    m_chisqr(1.e+10),
    m_n_iteration(0),
    m_mom0(0.,0.,0.),
    m_mom0_corP(0.,0.,0.),
    m_mom0_corN(0.,0.,0.),
    m_median_t(0.), m_min_t(0.), m_max_t(0.),
    m_path(0.),
    m_transverse_path(0.),
    m_charge(0), m_fitflag(0),
    m_isBeam(0), m_isK18(0), m_isKurama(0), m_isAccidental(0),
    m_searchtime(0), m_fittime(0),
    m_cx_exclusive(), m_cy_exclusive(), m_z0_exclusive(),
    m_r_exclusive(), m_dz_exclusive(), m_chisqr_exclusive(),
    m_vp()
{
  m_hit_array.reserve(ReservedNumOfHits);
  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
TPCLocalTrackHelix::~TPCLocalTrackHelix()
{
  debug::ObjectCounter::decrease(ClassName());
}

//______________________________________________________________________________
TPCLocalTrackHelix::TPCLocalTrackHelix(TPCLocalTrackHelix *init){

  this -> m_is_fitted = false;
  this -> m_is_calculated = false;
  this -> m_flag = 0;
  for(int i=0;i<init -> m_hit_array.size();i++){
    this -> m_hit_array.push_back(new TPCLTrackHit(init -> m_hit_array[i] -> GetHit()));
  }

  this -> m_Acx = init -> m_Acx ;
  this -> m_Acy = init -> m_Acy ;
  this -> m_Az0 = init -> m_Az0 ;
  this -> m_Ar = init -> m_Ar ;
  this -> m_Adz = init -> m_Adz ;
  this -> m_cx = init -> m_cx ;
  this -> m_cy = init -> m_cy ;
  this -> m_z0 = init -> m_z0 ;
  this -> m_r = init -> m_r ;
  this -> m_dz = init -> m_dz ;
  this -> m_chisqr = init -> m_chisqr ;
  this -> m_n_iteration = init -> m_n_iteration ;
  this -> m_mom0 = init -> m_mom0 ;
  this -> m_mom0_corP = init -> m_mom0_corP ;
  this -> m_mom0_corN = init -> m_mom0_corN ;
  this -> m_median_t = init -> m_median_t ;
  this -> m_min_t = init -> m_min_t ;
  this -> m_max_t = init -> m_max_t ;
  this -> m_path = init -> m_path ;
  this -> m_transverse_path = init -> m_transverse_path ;
  this -> m_charge = init -> m_charge ;
  this -> m_fitflag = init -> m_fitflag ;
  this -> m_isBeam = init -> m_isBeam ;
  this -> m_isK18 = init -> m_isK18 ;
  this -> m_isKurama = init -> m_isKurama ;
  this -> m_isAccidental = init -> m_isAccidental ;
  this -> m_searchtime = init -> m_searchtime ; //millisec
  this -> m_fittime = init -> m_fittime ; //millisec

  for(int i=0;i<init -> m_hit_order.size();i++){
    this -> m_hit_order.push_back(init -> m_hit_order[i]);
  }
  for(int i=0;i<init -> m_hit_t.size();i++){
    this -> m_hit_t.push_back(init -> m_hit_t[i]);
  }

  for(int i=0;i<init -> m_cx_exclusive.size();i++){
    this -> m_cx_exclusive.push_back(init -> m_cx_exclusive[i]);
    this -> m_cy_exclusive.push_back(init -> m_cy_exclusive[i]);
    this -> m_z0_exclusive.push_back(init -> m_z0_exclusive[i]);
    this -> m_r_exclusive.push_back(init -> m_r_exclusive[i]);
    this -> m_dz_exclusive.push_back(init -> m_dz_exclusive[i]);
    this -> m_chisqr_exclusive.push_back(init -> m_chisqr_exclusive[i]);
  }

  for(int i=0;i<init -> m_vp.size();i++){
    this -> m_vp.push_back(init -> m_vp[i]);
  }
  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
static inline bool CompareTheta(int a, int b){
  return gTcal[a] < gTcal[b];
}

//______________________________________________________________________________
static inline void fcn2(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;

  double fpar[8];
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
  }

  for(int i=0; i<gNumOfHits; ++i){
    TVector3 pos(-gHitPos[i].X(),
		 gHitPos[i].Z()-tpc::ZTarget,
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
		     fittmp.Y()+tpc::ZTarget);

    double tmp_t = TMath::ATan2(pos.Y()-par[1], pos.X()-par[0]);
    double tmpx = par[0] + par[3]*cos(tmp_t);
    double tmpy = par[1] + par[3]*sin(tmp_t);
    double tmpz = par[2] + (par[4]*par[3]*tmp_t);

    TVector3 fittmp2_(-tmpx,
    		      tmpz,
    		      tmpy+tpc::ZTarget);

    TVector3 d = gHitPos[i] - fittmp_;
    TVector3 Res = gRes[i];
    chisqr += pow(d.x()/Res.x(), 2) + pow(d.y()/Res.y(), 2) + pow(d.z()/Res.z(), 2);

    dof++;
    dof++;
  }
  f = chisqr/(double)(dof-5);
}

//______________________________________________________________________________
[[maybe_unused]] static inline void fcn2_circ(int &npar, double *gin, double &f, double *par, int iflag)
{
  double chisqr=0.0;
  int dof = 0;

  double fpar[5];
  //std::cout<<"paramter in fcn"<<std::endl;
  for(int ip=0; ip<3; ++ip){
    fpar[ip] = par[ip];
    //std::cout<<"par ["<<ip<<"]: "<<par[ip]<<std::endl;
  }

  for(int i=0; i<gNumOfHits; ++i){
    TVector2 pos_(gHitPos[i].X(),
		  gHitPos[i].Z());

    TVector2 pos(-gHitPos[i].X(),
		 gHitPos[i].Z()-tpc::ZTarget);
    fpar[3] = pos.X();
    fpar[4] = pos.Y();

    fint_circ.SetParameters(fpar);
    double min_t = fint_circ.GetMinimumX();
    double  x = par[0] + par[2]*cos(min_t);
    double  y = par[1] + par[2]*sin(min_t);

    TVector2 fittmp(x, y);
    TVector2 fittmp_(-1.*fittmp.X(),
		     fittmp.Y()+tpc::ZTarget);
    // double tmp_t = atan2(pos.Y()-par[1], pos.X()-par[0]);
    // double  tmpx = par[0] + par[3]*cos(tmp_t);
    // double  tmpy = par[1] + par[3]*sin(tmp_t);
    // TVector2 fittmp2_(-tmpx,
    // 		      tmpy+tpc::ZTarget);
    TVector2 d = pos_ - fittmp_;
    chisqr += pow(d.X()/gRes[i].x(), 2) + pow(d.Y()/gRes[i].z(), 2);

    //double dxy = d.Mod();
    //double resxy = sqrt(pow(gRes[i].x(), 2) + pow(gRes[i].z(), 2));
    //chisqr += pow(dxy/resxy,2);

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
  for(int ip=0; ip<2; ++ip){
    fpar[ip] = par[ip];
  }

  for(int i=0; i<gNumOfHits; ++i){
    TVector3 pos(-gHitPos[i].X(),
		 gHitPos[i].Z()-tpc::ZTarget,
		 gHitPos[i].Y());
    double tmpx = pos.X();
    double tmpy = pos.Y();
    double tmp_t = TMath::ATan2(tmpy - gPar[1],
			 tmpx - gPar[0]);
    fpar[2] = gPar[3];
    fpar[3] = tmp_t;
    fpar[4] = pos.Z();

    double diff_z = TMath::Power(fpar[4]-(fpar[0]+(fpar[1]*fpar[2]*fpar[3])),2);

    chisqr += diff_z/pow(gRes[i].y(),2);
    dof++;
  }
  f = chisqr/(double)(dof-2);

}

//______________________________________________________________________________
static inline double CalcChi2(double *HelixPar)
{

  double fpar[8];
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = HelixPar[ip];
  }

  int dof = 0; double chisqr = 0.;
  for(int i=0; i<gNumOfHits; ++i){
    TVector3 pos_(-gHitPos[i].X(),
		  gHitPos[i].Z() - tpc::ZTarget,
		  gHitPos[i].Y());
    fpar[5] = pos_.X();
    fpar[6] = pos_.Y();
    fpar[7] = pos_.Z();

    fint.SetParameters(fpar);
    double min_t = fint.GetMinimumX();
    double calcx = fpar[0] + fpar[3]*TMath::Cos(min_t);
    double calcy = fpar[1] + fpar[3]*TMath::Sin(min_t);
    double calcz = fpar[2] + fpar[4]*fpar[3]*min_t;
    TVector3 fittmp(calcx, calcy, calcz);
    TVector3 fittmp_(-fittmp.X(),
		     fittmp.Z(),
		     fittmp.Y() + tpc::ZTarget);
    TVector3 d = gHitPos[i] - fittmp_;
    chisqr += TMath::Power(d.x()/gRes[i].X(), 2) + TMath::Power(d.y()/gRes[i].Y(), 2) + TMath::Power(d.z()/gRes[i].Z(), 2);

    dof++;
    dof++;
  }

  return chisqr/(double)(dof-5);
}

//______________________________________________________________________________
static inline double CircleFit(const double *mX,const double *mY, const int npoints, double* mXCenter, double* mYCenter, double* mRadius)
{
  double xx, yy, xx2, yy2;
  double f, g, h, p, q, t, g0, g02, a, b, c, d;
  double xroot, ff, fp, xd, yd, g1;
  double dx, dy, dradius2, xnom;

  double xgravity = 0.0;
  double ygravity = 0.0;
  double x2 = 0.0;
  double y2 = 0.0;
  double xy = 0.0;
  double xx2y2 = 0.0;
  double yx2y2 = 0.0;
  double x2y22 = 0.0;
  double radius2 = 0.0;
  double mVariance = 0.0;

  if (npoints <= 3){
    fprintf(stderr,"CircleFit: npoints %d <= 3\n",npoints);
    return -1;
  }else  if (npoints > 499){
    fprintf(stderr,"CircleFit: npoints %d > 499\n",npoints);
    return -1;
  }

  for (int i=0; i<npoints; i++) {
    xgravity += mX[i];
    ygravity += mY[i];
  }
  xgravity /= npoints;
  ygravity /= npoints;

  for (int i=0; i<npoints; i++) {
    xx  = mX[i]-xgravity;
    yy  = mY[i]-ygravity;
    xx2 = xx*xx;
    yy2 = yy*yy;
    x2  += xx2;
    y2  += yy2;
    xy  += xx*yy;
    xx2y2 += xx*(xx2+yy2);
    yx2y2 += yy*(xx2+yy2);
    x2y22 += (xx2+yy2)*(xx2+yy2);
  }
  if (xy == 0.){
    fprintf(stderr,"CircleFit: xy = %f,    grav=%f, %f\n",xy,xgravity,ygravity);
    return -1;
  }

  f = (3.*x2+y2)/npoints;
  g = (x2+3.*y2)/npoints;
  h = 2*xy/npoints;
  p = xx2y2/npoints;
  q = yx2y2/npoints;
  t = x2y22/npoints;
  g0 = (x2+y2)/npoints;
  g02 = g0*g0;
  a = -4.0;
  b = (f*g-t-h*h)/g02;
  c = (t*(f+g)-2.*(p*p+q*q))/(g02*g0);
  d = (t*(h*h-f*g)+2.*(p*p*g+q*q*f)-4.*p*q*h)/(g02*g02);
  xroot = 1.0;
  for (int i=0; i<5; i++) {
    ff = (((xroot+a)*xroot+b)*xroot+c)*xroot+d;
    fp = ((4.*xroot+3.*a)*xroot+2.*b)*xroot+c;
    xroot -= ff/fp;
  }
  g1 = xroot*g0;
  xnom = (g-g1)*(f-g1)-h*h;
  if (xnom == 0.){
    fprintf(stderr,"CircleFit: xnom1 = %f\n",xnom);
    return -1;
  }

  yd = (q*(f-g1)-h*p)/xnom;
  xnom = f-g1;
  if (xnom == 0.){
    fprintf(stderr,"CircleFit: xnom2 = %f\n",xnom);
    return -1;
  }

  xd = (p-h*yd )/xnom;

  radius2 = xd*xd+yd*yd+g1;
  *mXCenter = xd+xgravity;
  *mYCenter = yd+ygravity;
  for (int i=0; i<npoints; i++) {
    dx = mX[i]-(*mXCenter);
    dy = mY[i]-(*mYCenter);
    dradius2 = dx*dx+dy*dy;
    mVariance += dradius2+radius2-2.*sqrt(dradius2*radius2);
  }

  *mRadius  = (double) sqrt(radius2);

  return  mVariance;
}

//______________________________________________________________________________
static inline void LineFit()
{
  //pre t-y fit
  double par_li[2] = {gPar[2], gPar[4]};
  double err_li[2] = {-999., -999.};
  TMinuit *minuit_li = new TMinuit(2);
  minuit_li->SetPrintLevel(-1);
  minuit_li->SetFCN(fcn2_li);

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

  for(int i=0; i<2; i++){
    minuit_li->mnpout(i, name_li[i], par_li[i], err_li[i], bnd1_li, bnd2_li, Err_li);
  }
  delete minuit_li;

  gPar[2] = par_li[0];
  gPar[4] = par_li[1];
}

//______________________________________________________________________________
static inline void HelixFit(int IsBeam)
{

  double par[5] = {gPar[0], gPar[1], gPar[2], gPar[3], gPar[4]};
  double err[5] = {-999., -999., -999., -999., -999.};

  TMinuit *minuit = new TMinuit(5);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn2);

  int ierflg = 0;
  double arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist,1,ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};
  for(int i=0; i<5; i++){
    if(IsBeam==1) minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimitBeam[i], UpLimit[i], ierflg);
    else minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }

  minuit->Command("SET STRategy 0");
  // arglist[0] = 5000.;
  // arglist[1] = 0.01;
  arglist[0] = 1000;
  arglist[1] = 0.1;

  arglist[0] = arglist[0]*5*5*5;
  arglist[1] = arglist[1]/(10.*10.*10.);

  int Err;
  double bnd1, bnd2;

  int itry=0;
  while(gChisqr>1.5){
    if(itry>MaxTryMinuit) break;
    minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    minuit->mnimpr();
    //minuit->mnexcm("MINOS", arglist, 0, ierflg);
    //minuit->mnexcm("SET ERR", arglist, 2, ierflg);

    // double amin, edm, errdef;
    // int nvpar, nparx, icstat;
    // minuit.mnstat(amin, edm, errdef, nvpar, nparx, icstat);
    //minuit->mnprin(4, amin);
    for(int i=0; i<5; i++){
      minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
    }

    double Chisqr = CalcChi2(par);
    if(gChisqr>Chisqr){
      gChisqr = Chisqr;
      gPar[0] = par[0];
      gPar[1] = par[1];
      gPar[2] = par[2];
      gPar[3] = par[3];
      gPar[4] = par[4];
    }

    arglist[0] = arglist[0]*5;
    arglist[1] = arglist[1]/10.;
    ++itry;
  }

  delete minuit;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::ClearHits()
{
  m_hit_array.clear();
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::AddTPCHit(TPCLTrackHit *hit)
{
  if(hit){
    m_hit_order.push_back(m_hit_array.size());
    m_hit_array.push_back(hit);
  }
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::Calculate()
{
  if(IsCalculated()){
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "already called" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    hitp->SetCalHelix(m_cx, m_cy, m_z0, m_r, m_dz);
    hitp->SetTheta(m_hit_t[i]);
    hitp->SetCalPosition(hitp->GetLocalCalPosHelix());
  }
  m_is_calculated = true;

}

//______________________________________________________________________________
void
TPCLocalTrackHelix::CalculateExclusive()
{

  if(!IsCalculated()){
    hddaq::cerr << "#W " << FUNC_NAME << "No inclusive calculation" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    hitp->SetCalHelixExclusive(m_cx_exclusive[i], m_cy_exclusive[i], m_z0_exclusive[i], m_r_exclusive[i], m_dz_exclusive[i]);
    hitp->SetCalPositionExclusive(hitp->GetLocalCalPosHelixExclusive());
  }
}

//______________________________________________________________________________
int
TPCLocalTrackHelix::GetNDF() const
{

  const std::size_t n = m_hit_array.size();
  int ndf = 0;
  for(std::size_t i=0; i<n; ++i){
    if(m_hit_array[i]){
      ++ndf;
      ++ndf;
    }
  }
  return ndf-5;
}

//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::GetPosition(double par[5], double t) const
{

  // double  x = p[0] + p[3]*cos(t+theta0);
  // double  y = p[1] + p[3]*sin(t+theta0);
  // double  z = p[2] + p[3]*p[4]*(t+theta0);
  //This is the eqation of Helix
  double  x = par[0] + par[3]*cos(t);
  double  y = par[1] + par[3]*sin(t);
  double  z = par[2] + (par[4]*par[3]*t);

  return TVector3(x, y, z);
}


//______________________________________________________________________________
TVector3
TPCLocalTrackHelix::CalcHelixMom(double par[5], double y) const
{

  const double Const = 0.299792458; // =c/10^9
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
TVector3
TPCLocalTrackHelix::CalcHelixMom_t(double par[5], double t) const
{

  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  //  double t = (y-par[2])/(par[3]*par[4]);
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

// momentum correction for positive particle
//______________________________________________________________________________
[[maybe_unused]] TVector3
TPCLocalTrackHelix::CalcHelixMom_corP(double par[5], double y) const
{

  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  // obtained by Beam through analysis
  double cor_p1 = 0.6222;
  double cor_p0 = 0.1198*1000.;

  double t = (y-par[2])/(par[3]*par[4]);
  double pt = fabs(par[3])*(Const*dMagneticField); // MeV/c

  pt = pt*cor_p1 + cor_p0;
  if(pt<10.)
    pt =10.;

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
[[maybe_unused]] TVector3
TPCLocalTrackHelix::CalcHelixMom_t_corP(double par[5], double t) const
{

  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  // obtained by Beam through analysis
  double cor_p1 = 0.6222;
  double cor_p0 = 0.1198*1000.;

  //  double t = (y-par[2])/(par[3]*par[4]);
  double pt = fabs(par[3])*(Const*dMagneticField); // MeV/c
  pt = pt*cor_p1 + cor_p0;
  if(pt<10.)
    pt =10.;

  //From here!!!!
  double tmp_px = pt*(-1.*sin(t));
  double tmp_py = pt*(cos(t));
  double tmp_pz = pt*(par[4]);
  double px = -tmp_px*0.001;
  double py = tmp_pz*0.001;
  double pz = tmp_py*0.001;

  return TVector3(px,py,pz);
}

// momentum correction for negative particle
//______________________________________________________________________________
[[maybe_unused]] TVector3
TPCLocalTrackHelix::CalcHelixMom_corN(double par[5], double y) const
{

  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  // obtained by Beam through analysis
  double cor_p1 = 1.07;
  double cor_p0 = -0.0353*1000.;

  double t = (y-par[2])/(par[3]*par[4]);
  double pt = fabs(par[3])*(Const*dMagneticField); // MeV/c
  pt = pt*cor_p1 + cor_p0;
  if(pt<10.)
    pt =10.;

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
[[maybe_unused]] TVector3
TPCLocalTrackHelix::CalcHelixMom_t_corN(double par[5], double t) const
{

  const double Const = 0.299792458; // =c/10^9
  const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  // obtained by Beam through analysis
  double cor_p1 = 1.07;
  double cor_p0 = -0.0353*1000.;

  //  double t = (y-par[2])/(par[3]*par[4]);
  double pt = fabs(par[3])*(Const*dMagneticField); // MeV/c
  pt = pt*cor_p1 + cor_p0;
  if(pt<10.)
    pt =10.;

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
TPCLTrackHit*
TPCLocalTrackHelix::GetHit(std::size_t nth) const
{
  if(nth<m_hit_array.size())
    return m_hit_array[nth];
  else
    return 0;
}

//______________________________________________________________________________
TPCLTrackHit*
TPCLocalTrackHelix::GetHitInOrder(std::size_t nth) const
{
  int order = m_hit_order[nth];
  if(nth<m_hit_array.size())
    return m_hit_array[order];
  else
    return 0;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::DeleteNullHit()
{

  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    if(!hit){
      hddaq::cout << FUNC_NAME << " "
		  << "null hit is deleted" << std::endl;
      m_hit_array.erase(m_hit_array.begin()+i);
      --i;
    }
  }
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SetHoughFlag(int hough_flag)
{
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TPCHit *hit = hitp->GetHit();
    if( !hit ) continue;
    hit->SetHoughFlag(hough_flag);
  }
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SetFlag(int flag)
{
  if((flag&1)==1) m_isBeam=1;
  if((flag&2)==2) m_isK18=1;
  if((flag&4)==4) m_isKurama=1;
  if((flag&8)==8) m_isAccidental=1;
}

/*
//______________________________________________________________________________
bool
TPCLocalTrackHelix::DoFit(int MinHits)
{

  bool status = DoHelixFit(MinHits);
  m_is_fitted = status;
  if(m_chisqr<MaxChisqr) return status;
  else return false;
}

//______________________________________________________________________________
bool
TPCLocalTrackHelix::DoHelixFit(int MinHits)
{

  DeleteNullHit();

  const std::size_t n = m_hit_array.size();
  std::vector<Int_t> pads;
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    pads.push_back(hit -> GetHit() -> GetPad());
  }
  std::sort(pads.begin(),pads.end());
  pads.erase(std::unique(pads.begin(),pads.end()),pads.end());
  if(pads.size()<MinHits) return false;

  if(GetNDF()<1){
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "Min layer should be > NDF" << std::endl;
    return false;
  }

  if(n>ReservedNumOfHits){
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "n > ReservedNumOfHits" << std::endl;
    return false;
  }
  gNumOfHits = n;
  gHitPos.clear();
  gRes.clear();

  //for pre circle fit
  double xp[n];
  double yp[n];
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    gHitPos.push_back(pos);
    gRes.push_back(Res);
    xp[i] = -pos.X();
    yp[i] = pos.Z()-tpc::ZTarget;
  }

  double par_circ[3]={0};
  CircleFit(xp, yp, n, &par_circ[0], &par_circ[1], &par_circ[2]);
  gPar[0] = par_circ[0];
  gPar[1] = par_circ[1];
  gPar[2] = 0.;
  gPar[3] = par_circ[2];
  gPar[4] = 0.;

  //Pre linear fitting
  tpc::HoughTransformLineYPhi(gHitPos, gPar, 1000.);
  LineFit();
  gChisqr = CalcChi2(gPar);

#if DebugDisp
  std::cout<<"Before helix fitting"
	   <<" m_chisqr: "<<gChisqr
	   <<", m_cx: "<<m_Acx
	   <<", m_cy: "<<m_Acy
	   <<", m_z0: "<<m_Az0
	   <<", m_r: "<<m_Ar
	   <<", m_dz: "<<m_Adz
	   <<std::endl;
#endif

  //Helix fitting
  HelixFit(m_isBeam);
  m_chisqr = gChisqr;
  m_cx = gPar[0];
  m_cy = gPar[1];
  m_z0 = gPar[2];
  m_r  = gPar[3];
  m_dz = gPar[4];

#if InvertChargeTest
  InvertChargeCheck();
#endif

  int delete_hit = -1;
  int false_layer = FinalizeTrack(delete_hit);
  if(m_path>500.||fabs(m_min_t-m_max_t)>acos(-1)){
    if(fabs(m_min_t - m_median_t)>fabs(m_max_t - m_median_t)){
#if DebugDisp
      std::cout<<"delete mint, m_median_t="
	       <<m_median_t<<", mint="<<m_min_t<<std::endl;
#endif
      int mint_hit = m_hit_order[0];
      TPCLTrackHit *hitp = m_hit_array[mint_hit];
      TPCHit *hit = hitp->GetHit();
      hit->SetHoughFlag(0);
      m_hit_array.erase(m_hit_array.begin()+mint_hit);
      m_hit_t.erase(m_hit_t.begin()+mint_hit);
      gTcal.erase(gTcal.begin()+mint_hit);
      Sort();
    }
    else{
#if DebugDisp
      std::cout<<"delete maxt, m_median_t="
	       <<m_median_t<<", maxt="<<m_max_t<<std::endl;
#endif
      int maxt_hit = m_hit_order[m_hit_order.size()-1];
      TPCLTrackHit *hitp = m_hit_array[maxt_hit];
      TPCHit *hit = hitp->GetHit();
      hit->SetHoughFlag(0);
      m_hit_array.erase(m_hit_array.begin()+maxt_hit);
      m_hit_t.erase(m_hit_t.begin()+maxt_hit);
      gTcal.erase(gTcal.begin()+maxt_hit);
      Sort();
    }
  }
  else if(false_layer>0){
    TPCLTrackHit *hitp = m_hit_array[delete_hit];
    TPCHit *hit = hitp->GetHit();
    hit->SetHoughFlag(0);
    m_hit_array.erase(m_hit_array.begin()+delete_hit);
    m_hit_t.erase(m_hit_t.begin()+delete_hit);
    gTcal.erase(gTcal.begin()+delete_hit);
    Sort();
#if DebugDisp
    TVector3 pos = hitp->GetLocalHitPos();
    std::cout<<"delete hits ["<<delete_hit<<"]=("
    	     <<pos.x()<<", "
      	     <<pos.y()<<", "
      	     <<pos.z()<<std::endl;
#endif
  }

  m_n_iteration++;
  if(m_hit_array.size()<MinHits) return false;
  if(false_layer==0) return true;
  else return DoHelixFit(MinHits);
}
*/
//______________________________________________________________________________
bool
TPCLocalTrackHelix::DoFit(int MinHits)
{

  bool status = DoHelixFit(); //track chisqr minimization
  m_is_fitted = status;

  // #hits < MinHits
  const std::size_t n = m_hit_array.size();
  std::vector<Int_t> pads;
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    pads.push_back(hit -> GetHit() -> GetPad());
  }
  std::sort(pads.begin(),pads.end());
  pads.erase(std::unique(pads.begin(),pads.end()),pads.end());
  if(pads.size()<MinHits) return false;

  // #hits >= MinHits,
  if(m_chisqr<MaxChisqr) return status;
  else return false;
}

//______________________________________________________________________________
bool
TPCLocalTrackHelix::DoHelixFit()
{

  DeleteNullHit();
  const std::size_t n = m_hit_array.size();
  if(GetNDF()<1){
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "Min layer should be > NDF" << std::endl;
    return false;
  }

  if(n>ReservedNumOfHits){
    hddaq::cerr << "#W " << FUNC_NAME << " "
		<< "n > ReservedNumOfHits" << std::endl;
    return false;
  }
  gNumOfHits = n;
  gHitPos.clear();
  gRes.clear();

  //for pre circle fit
  double xp[n];
  double yp[n];
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    gHitPos.push_back(pos);
    gRes.push_back(Res);
    xp[i] = -pos.X();
    yp[i] = pos.Z()-tpc::ZTarget;
  }

  double par_circ[3]={0};
  CircleFit(xp, yp, n, &par_circ[0], &par_circ[1], &par_circ[2]);
  gPar[0] = par_circ[0];
  gPar[1] = par_circ[1];
  gPar[2] = 0.;
  gPar[3] = par_circ[2];
  gPar[4] = 0.;

  //Pre linear fitting
  tpc::HoughTransformLineYPhi(gHitPos, gPar, 1000.);
  LineFit();
  gChisqr = CalcChi2(gPar);

#if DebugDisp
  std::cout<<"Before helix fitting"
	   <<" m_chisqr: "<<gChisqr
	   <<", m_cx: "<<m_Acx
	   <<", m_cy: "<<m_Acy
	   <<", m_z0: "<<m_Az0
	   <<", m_r: "<<m_Ar
	   <<", m_dz: "<<m_Adz
	   <<std::endl;
#endif

  //Helix fitting
  HelixFit(m_isBeam);
  m_chisqr = gChisqr;
  m_cx = gPar[0];
  m_cy = gPar[1];
  m_z0 = gPar[2];
  m_r  = gPar[3];
  m_dz = gPar[4];

#if InvertChargeTest
  InvertChargeCheck();
#endif

  int delete_hit = -1;
  int false_layer = FinalizeTrack(delete_hit);
  if(m_path>500.||fabs(m_min_t-m_max_t)>acos(-1)){
    if(fabs(m_min_t - m_median_t)>fabs(m_max_t - m_median_t)){
#if DebugDisp
      std::cout<<"delete mint, m_median_t="
	       <<m_median_t<<", mint="<<m_min_t<<std::endl;
#endif
      int mint_hit = m_hit_order[0];
      TPCLTrackHit *hitp = m_hit_array[mint_hit];
      TPCHit *hit = hitp->GetHit();
      hit->SetHoughFlag(0);
      m_hit_array.erase(m_hit_array.begin()+mint_hit);
      m_hit_t.erase(m_hit_t.begin()+mint_hit);
      gTcal.erase(gTcal.begin()+mint_hit);
      Sort();
    }
    else{
#if DebugDisp
      std::cout<<"delete maxt, m_median_t="
	       <<m_median_t<<", maxt="<<m_max_t<<std::endl;
#endif
      int maxt_hit = m_hit_order[m_hit_order.size()-1];
      TPCLTrackHit *hitp = m_hit_array[maxt_hit];
      TPCHit *hit = hitp->GetHit();
      hit->SetHoughFlag(0);
      m_hit_array.erase(m_hit_array.begin()+maxt_hit);
      m_hit_t.erase(m_hit_t.begin()+maxt_hit);
      gTcal.erase(gTcal.begin()+maxt_hit);
      Sort();
    }
  }
  else if(false_layer>0){
    TPCLTrackHit *hitp = m_hit_array[delete_hit];
    TPCHit *hit = hitp->GetHit();
    hit->SetHoughFlag(0);
    m_hit_array.erase(m_hit_array.begin()+delete_hit);
    m_hit_t.erase(m_hit_t.begin()+delete_hit);
    gTcal.erase(gTcal.begin()+delete_hit);
    Sort();
#if DebugDisp
    TVector3 pos = hitp->GetLocalHitPos();
    std::cout<<"delete hits ["<<delete_hit<<"]=("
    	     <<pos.x()<<", "
      	     <<pos.y()<<", "
      	     <<pos.z()<<std::endl;
#endif
  }

  m_n_iteration++;
  if(m_hit_array.size()<3) return false;
  if(false_layer==0) return true;
  else return DoHelixFit();
}

//______________________________________________________________________________
bool
TPCLocalTrackHelix::ResidualCheck(TVector3 pos, TVector3 Res, double &resi)
{

  bool status = false;

  double par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 pos_(-pos.X(),
		pos.Z()-tpc::ZTarget,
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
		   fittmp.Y()+tpc::ZTarget);
  TVector3 d = pos - fittmp_;

  resi = d.Mag();
  if(d.Mag()<Res.Mag()*MaxResidual)
    status = true;

  return status;
}

//______________________________________________________________________________
double
TPCLocalTrackHelix::GetTcal(TVector3 pos)
{
  double par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 pos_(-pos.X(),
		pos.Z()-tpc::ZTarget,
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
TPCLocalTrackHelix::DoHelixFitExclusive()
{
  const std::size_t n = m_hit_array.size();
  m_cx_exclusive.resize(n);
  m_cy_exclusive.resize(n);
  m_z0_exclusive.resize(n);
  m_r_exclusive.resize(n);
  m_dz_exclusive.resize(n);
  m_chisqr_exclusive.resize(n);

  gNumOfHits = n-1;
  for(Int_t ihit=0; ihit<n; ++ihit){
    gHitPos.clear();
    gRes.clear();
    gHitPos.resize(n-1);
    gRes.resize(n-1);

    gChisqr = 1.e+10;
    gPar[0] = m_cx;
    gPar[1] = m_cy;
    gPar[2] = m_z0;
    gPar[3] = m_r;
    gPar[4] = m_dz;

    int flag=0;
    for(Int_t j=0; j<n; ++j){
      if(j==ihit) continue; //exclude ith hit
      TPCLTrackHit *hitp = m_hit_array[j];
      TVector3 pos = hitp->GetLocalHitPos();
      TVector3 res = hitp->GetResolutionVect();
      gHitPos[flag] = pos;
      gRes[flag] = res;
      flag++;
    } //j

    //Helix fitting
    HelixFit(m_isBeam);

    m_chisqr_exclusive[ihit] = gChisqr;
    m_cx_exclusive[ihit] = gPar[0];
    m_cy_exclusive[ihit] = gPar[1];
    m_z0_exclusive[ihit] = gPar[2];
    m_r_exclusive[ihit]  = gPar[3];
    m_dz_exclusive[ihit] = gPar[4];

#if DebugDisp
    std::cout<<"exclusive "<<ihit<<" th chisqr : "<<m_chisqr_exclusive[ihit]<<
      " par : "<<m_cx_exclusive[ihit]<<
      " "<<m_cy_exclusive[ihit]<<
      " "<<m_z0_exclusive[ihit]<<
      " "<<m_r_exclusive[ihit]<<
      " "<<m_dz_exclusive[ihit]<<std::endl;
#endif
  } //ihit
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::Print(const TString& arg, bool print_allhits) const
{

  TString tracksize = Form(" #clusters = %d", (int)m_hit_array.size());
  std::cout<<arg.Data()<<std::endl;
  std::cout<<"Track info : "<<tracksize.Data()<<std::endl;
  if(print_allhits){
    for(std::size_t i=0; i<m_hit_array.size(); ++i){
      TPCLTrackHit *hitp = m_hit_array[i];
      if( !hitp ) continue;
      hitp->Print();
    }
  }
}

//______________________________________________________________________________
double
TPCLocalTrackHelix::GetdEdx(double truncatedMean)
{

  DeleteNullHit();
  std::vector<Double_t> dEdx_vect;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    Double_t clde = hitp->GetDe();
    Double_t pathHit = hitp->GetPathHelix();
    Double_t dEdx_cor = clde/pathHit;
    dEdx_vect.push_back(dEdx_cor);
  }

  double dEdx = 0.;
  std::sort(dEdx_vect.begin(), dEdx_vect.end());
  int n_truncated = (int)(dEdx_vect.size()*truncatedMean);
  for( int ih=0; ih<dEdx_vect.size(); ++ih ){
    if(ih<n_truncated) dEdx += dEdx_vect[ih];
  }
  dEdx /= (double)n_truncated;

  return dEdx;
}

//______________________________________________________________________________
double
TPCLocalTrackHelix::GetTrackdE()
{

  DeleteNullHit();
  double dE = 0;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    dE += hitp->GetDe();
  }
  return dE;
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::SetParamUsingHoughParam()
{

  DeleteNullHit();

  gPar[0] = m_Acx;
  gPar[1] = m_Acy;
  gPar[2] = m_Az0;
  gPar[3] = m_Ar;
  gPar[4] = m_Adz;
  gChisqr = CalcChi2(gPar);

  m_cx = gPar[0];
  m_cy = gPar[1];
  m_z0 = gPar[2];
  m_r  = gPar[3];
  m_dz = gPar[4];
  m_chisqr = gChisqr;

  int dummy;
  FinalizeTrack(dummy);

}

//______________________________________________________________________________
void
TPCLocalTrackHelix::AddVPHit(TVector3 vp)
{
  m_vp.push_back(vp);
}

//______________________________________________________________________________
bool
TPCLocalTrackHelix::DoVPFit()
{

  const std::size_t n = m_vp.size();
  gNumOfHits = n;
  gHitPos.clear();
  gRes.clear();

  //Pre circle fit
  TVector3 res_vp(0.3, 0.55, 0.3); //dummy value
  double xp[n]; double yp[n];
  for(std::size_t i=0; i<n; ++i){
    gHitPos.push_back(m_vp[i]);
    gRes.push_back(res_vp);
    xp[i] = -m_vp[i].X();
    yp[i] = m_vp[i].Z() - tpc::ZTarget;
  }

  double par_circ[3]={0};
  CircleFit(xp, yp, n, &par_circ[0], &par_circ[1], &par_circ[2]);
  m_cx = par_circ[0];
  m_cy = par_circ[1];
  m_r = par_circ[2];

  double par[5] = {m_cx, m_cy, 0, m_r, 0};
  tpc::HoughTransformLineYPhi(gHitPos, par, 100.);

  m_z0 = par[2];
  m_dz = par[4];

  return true;
}

//______________________________________________________________________________
bool
TPCLocalTrackHelix::ResidualCheck(TVector3 pos, double xzwindow, double ywindow, double &resi)
{

  bool status = false;
  double par[5] = {m_cx, m_cy, m_z0, m_r, m_dz};
  TVector3 pos_(-pos.X(),
		pos.Z() - tpc::ZTarget,
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
		   fittmp.Y() + tpc::ZTarget);
  TVector3 d = pos - fittmp_;
  resi = d.Mag();
  double xz_resi = TMath::Sqrt(d.x()*d.x()+d.z()*d.z());
  double y_resi = TMath::Sqrt(d.y()*d.y());
  if(xz_resi<xzwindow && y_resi<ywindow) status = true;
  return status;
}

//______________________________________________________________________________
bool
TPCLocalTrackHelix::ResidualCheck(TVector3 pos, double xzwindow, double ywindow)
{
  double resi;
  return ResidualCheck(pos, xzwindow, ywindow, resi);
}

//______________________________________________________________________________
void
TPCLocalTrackHelix::Sort()
{

  m_hit_order.clear();
  for(std::size_t i=0; i<m_hit_t.size(); ++i){
    m_hit_order.push_back(i);
  }
  for(std::size_t i=0; i<m_hit_order.size(); ++i){
    std::sort(m_hit_order.begin(), m_hit_order.end(), CompareTheta);
  }
}

//______________________________________________________________________________
int
TPCLocalTrackHelix::FinalizeTrack(int &delete_hit)
{

  int false_layer = 0;
  double Max_residual = -100.;
  double minlayer_t = 0., maxlayer_t = 0.;
  int minlayer = 33, maxlayer = -1;

  gTcal.clear();
  m_hit_t.clear(); m_hit_order.clear();
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    m_hit_order.push_back(i);
    TPCLTrackHit *hitp = m_hit_array[i];
    double resi = 0.;
    TVector3 pos = hitp->GetLocalHitPos();
    double t = GetTcal(pos);
    int layer = hitp->GetLayer();
    m_hit_t.push_back(t);
    gTcal.push_back(t);

    TVector3 Res = hitp->GetResolutionVect();
    if(!ResidualCheck(pos,Res,resi)){
      if(Max_residual<resi){
	Max_residual = resi;
	delete_hit = i;
      }
      ++false_layer;
    }

    if(layer<minlayer){
      minlayer = layer;
      minlayer_t = t;
    }
    if(layer>maxlayer){
      maxlayer = layer;
      maxlayer_t = t;
    }
  }
  if(minlayer_t<maxlayer_t) m_charge = 1;
  else m_charge = -1;

for(std::size_t i=0; i<m_hit_order.size(); ++i){
  std::sort(m_hit_order.begin(), m_hit_order.end(), CompareTheta);
 }

 int mint = m_hit_order[0];
 int maxt = m_hit_order[m_hit_order.size()-1];
 m_min_t = m_hit_t[mint];
 m_max_t = m_hit_t[maxt];
 m_median_t = TMath::Median(m_hit_t.size(), m_hit_t.data());
 m_path = (m_max_t-m_min_t)*sqrt(m_r*m_r*(1.+m_dz*m_dz));
 m_transverse_path = (m_max_t-m_min_t)*m_r;
 m_mom0 = CalcHelixMom(gPar, 0.);
 m_mom0_corP = CalcHelixMom_corP(gPar, 0.);
 m_mom0_corN = CalcHelixMom_corN(gPar, 0.);

#if DebugDisp
  std::cout<<"m_chisqr: "<<m_chisqr
	   <<", m_cx: "<<m_cx
	   <<", m_cy: "<<m_cy
	   <<", m_z0: "<<m_z0
	   <<", m_r: "<<m_r
	   <<", m_dz: "<<m_dz
	   <<", m_charge: "<<m_charge
	   <<", m_path: "<<m_path
	   <<", fabs(m_min_t-m_max_t): "<<fabs(m_min_t-m_max_t)<<std::endl;
#endif

  return false_layer;
}


//______________________________________________________________________________
void
TPCLocalTrackHelix::InvertChargeCheck()
{

  // for invert charge fit
  double test_minlayer_t = 0., test_maxlayer_t = 0.;
  int test_minlayer = 33, test_maxlayer = -1;
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    double t = GetTcal(pos);
    int layer = hitp->GetLayer();
    if(test_minlayer > layer){
      test_minlayer = layer;
      test_minlayer_t = t;
    }
    if(test_maxlayer < layer){
      test_maxlayer = layer;
      test_maxlayer_t = t;
    }
  }
  double mid_t = (test_minlayer_t + test_maxlayer_t)/2.;
  TVector3 mid_pos = GetPosition(gPar, mid_t);
  double mid_x = mid_pos.x();
  double mid_y = mid_pos.y();

  double par[5];
  par[0] = m_cx - 2.*(m_cx - mid_x);
  par[1] = m_cy - 2.*(m_cy - mid_y);
  par[2] = m_z0;
  par[3] = m_r;
  par[4] = m_dz;
  double err[5]={-999., -999., -999., -999., -999.};

  TMinuit *minuit = new TMinuit(5);
  minuit->SetPrintLevel(-1);
  minuit->SetFCN(fcn2);

  int ierflg = 0;
  double arglist[10];
  arglist[0] = 5.89;
  minuit->mnexcm("SET ERR", arglist,1,ierflg); //Num of parameter
  arglist[0] = 1;
  minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings

  TString name[5] = {"cx", "cy", "z0", "r", "dz"};
  for(int i = 0; i<5; i++){
    if(GetIsBeam()==1)
      minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimitBeam[i], UpLimit[i], ierflg);
    else
      minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  minuit->mnexcm("MIGRAD", arglist, 2, ierflg);

  int Err;
  double bnd1, bnd2;
  for(int i=0; i<5; i++){
    minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
  }

  double chisqr = CalcChi2(par);
  if(gChisqr > chisqr){
    gChisqr = chisqr;
    gPar[0] = par[0];
    gPar[1] = par[1];
    gPar[2] = par[2];
    gPar[3] = par[3];
    gPar[4] = par[4];
  }
  m_chisqr = gChisqr;
  m_cx = gPar[0];
  m_cy = gPar[1];
  m_z0 = gPar[2];
  m_r  = gPar[3];
  m_dz = gPar[4];

}
