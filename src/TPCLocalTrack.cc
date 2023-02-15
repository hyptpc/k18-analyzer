// -*- C++ -*-

//Comment by Ichikawa
//TPCLocalTrack.cc is for lenear fit

#include "TPCLocalTrack.hh"

#include <string>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <stdexcept>
#include <sstream>

#include <TF2.h>
#include <TGraph2D.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TROOT.h>
#include <Math/Functor.h>
#include <Math/Vector3D.h>
#include <TPolyLine3D.h>
#include <Fit/Fitter.h>
#include <std_ostream.hh>

#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "HoughTransform.hh"
#include "TPCLTrackHit.hh"
#include "TPCPadHelper.hh"
#include "UserParamMan.hh"

#define DebugEvDisp 0

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
#endif

namespace
{
  //for Minimization
  static int gNumOfHits;
  static std::vector<TVector3> gHitPos;
  static std::vector<TVector3> gRes;
  static std::vector<double> gTcal;
  static double gPar[4] = {0};
  static double gChisqr = 1.e+10;

  static const Double_t MaxResidual = 10.;// factor to resolution

  const Int_t ReservedNumOfHits  = 64;
  const Double_t  FitStep[4] = { 1.0e-6, 1.0e-10, 1.0e-6, 1.0e-10};
  //const Double_t  FitStep[4] = { 1.0e-2, 1.0e-6, 1.0e-2, 1.0e-6};
  const Double_t  LowLimit[4] = { -400., -1.0*TMath::ACos(-1.), -400, -1.0*TMath::ACos(-1.) };
  const Double_t  UpLimit[4] = { 400., 1.0*TMath::ACos(-1.), 400, 1.0*TMath::ACos(-1.) };
  static const Double_t  MaxChisqr = 10000.;

  const Int_t    theta_ndiv = 100;
  const Double_t theta_min  =   0;
  const Double_t theta_max  = 180;
  const Int_t    r_ndiv =  100;
  const Double_t r_min  = -600;
  const Double_t r_max  =  600;

}

//______________________________________________________________________________
static inline double CalcChi2(double *par)
{

  double chisqr = 0.;
  for(int i=0; i<gNumOfHits; ++i){
    TVector3 x0(par[0], par[1], tpc::ZTarget);
    TVector3 x1(par[0] + par[2], par[1] + par[3], tpc::ZTarget+1.);
    TVector3 u = (x1-x0).Unit();
    TVector3 AP = gHitPos[i] - x0;
    Double_t dist_AX = u.Dot(AP);
    TVector3 AI(x0.x()+(u.x()*dist_AX),
		x0.y()+(u.y()*dist_AX),
		x0.z()+(u.z()*dist_AX));
    TVector3 d = gHitPos[i]-AI;
    chisqr += TVector3(d.x()/gRes[i].x(), d.y()/gRes[i].y(), d.z()/gRes[i].z()).Mag2();
  }

  return chisqr/(double)(gNumOfHits-4);
}
//_____________________________________________________________________________
static void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t chisqr=0.0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 pos = gHitPos[i];
    TVector3 Res = gRes[i];

    Double_t u0 = tan(par[1]);
    Double_t v0 = tan(par[3]);

    TVector3 x0(par[0], par[2], tpc::ZTarget);
    TVector3 x1(par[0] + u0, par[2] + v0, tpc::ZTarget+1.);

    TVector3 u = (x1-x0).Unit();
    //    TVector3 d = (pos-x0).Cross(u);
    TVector3 AP = pos-x0;
    Double_t dist_AX = u.Dot(AP);
    TVector3 AI(x0.x()+(u.x()*dist_AX),
		x0.y()+(u.y()*dist_AX),
		x0.z()+(u.z()*dist_AX));
    TVector3 d = pos-AI;
    chisqr += TVector3(d.x()/Res.x(), d.y()/Res.y(), d.z()/Res.z()).Mag2();
  }
  f = chisqr/(gNumOfHits - 4);
}

//_____________________________________________________________________________
[[maybe_unused]]
static void fcn2_rt(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t chisqr=0.0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 pos = gHitPos[i];
    TVector3 Res = gRes[i];
    Double_t x_0 = par[0]/TMath::Sin(par[1]);
    Double_t y_0 = par[2]/TMath::Sin(par[3]);
    Double_t u0 = -1.*TMath::Cos(par[1])/TMath::Sin(par[1]);
    Double_t v0 = -1.*TMath::Cos(par[3])/TMath::Sin(par[3]);

    // Double_t m_Au_theta = TMath::ATan(-1./u0);
    // Double_t m_Av_theta = TMath::ATan(-1./v0);
    // Double_t m_Ax_r = x_0*TMath::Sin(m_Au_theta);
    // Double_t m_Ay_r = y_0*TMath::Sin(m_Av_theta);

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

    TVector3 x0(x_0, y_0, tpc::ZTarget);
    TVector3 x1(x_0 + u0, y_0 + v0, tpc::ZTarget+1.);

    TVector3 u = (x1-x0).Unit();
    //    TVector3 d = (pos-x0).Cross(u);
    TVector3 AP = pos-x0;
    Double_t dist_AX = u.Dot(AP);
    TVector3 AI(x0.x()+(u.x()*dist_AX),
		x0.y()+(u.y()*dist_AX),
		x0.z()+(u.z()*dist_AX));
    TVector3 d = pos-AI;
    chisqr += TVector3(d.x()/Res.x(), d.y()/Res.y(), d.z()/Res.z()).Mag2();
  }
  f = chisqr/(gNumOfHits - 4);
}

//_____________________________________________________________________________
[[maybe_unused]]
static void fcn2_rt2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t chisqr=0.0;
  for(Int_t i=0; i<gNumOfHits; ++i){
    TVector3 pos = gHitPos[i];
    TVector3 Res = gRes[i];

    Double_t x_0 = par[0]/TMath::Cos(par[1]);
    Double_t y_0 = par[2]/TMath::Cos(par[3]);
    Double_t u0 = tan(par[1]);
    Double_t v0 = tan(par[3]);

    TVector3 x0(x_0, y_0, tpc::ZTarget);
    TVector3 x1(x_0 + u0, y_0 + v0, tpc::ZTarget+1.);

    TVector3 u = (x1-x0).Unit();
    //    TVector3 d = (pos-x0).Cross(u);
    TVector3 AP = pos-x0;
    Double_t dist_AX = u.Dot(AP);
    TVector3 AI(x0.x()+(u.x()*dist_AX),
		x0.y()+(u.y()*dist_AX),
		x0.z()+(u.z()*dist_AX));
    TVector3 d = pos-AI;
    chisqr += TVector3(d.x()/Res.x(), d.y()/Res.y(), d.z()/Res.z()).Mag2();
  }
  f = chisqr/(gNumOfHits - 4);
}

//______________________________________________________________________________
static inline void LinearFit()
{

  double par[4] = {gPar[0], TMath::ATan(gPar[2]), gPar[1], TMath::ATan(gPar[3])};
  double err[4] = {-999., -999., -999., -999.};

  TMinuit *m_minuit = new TMinuit(4);
  m_minuit->SetPrintLevel(-1);
  m_minuit->SetFCN(fcn2);

  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 1;
  m_minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings
  TString name[4] = {"x0_r", "atan_u0", "y0_r", "atan_v0"};
  for(Int_t i=0; i<4; i++){
    m_minuit->mnparm(i, name[i], par[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }

  m_minuit->Command("SET STRategy 0");
  arglist[0] = 5000.;
  arglist[1] = 0.01;
  // arglist[0] = 10000.;
  // arglist[1] = 0.001;
  m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  //m_minuit->mnexcm("MINOS", arglist, 0, ierflg);
  //m_minuit->mnexcm("SET ERR", arglist, 2, ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  m_minuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);
  Int_t Err;
  Double_t bnd1, bnd2;
  for(Int_t i=0; i<4; i++){
    m_minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
  }
  delete m_minuit;
  double par_linear[4] = {par[0], par[2], TMath::Tan(par[1]), TMath::Tan(par[3])};
  double Chisqr = CalcChi2(par_linear);
  if(gChisqr>Chisqr){
    gChisqr = Chisqr;
    gPar[0] = par_linear[0];
    gPar[1] = par_linear[1];
    gPar[2] = par_linear[2];
    gPar[3] = par_linear[3];
  }
}

//_____________________________________________________________________________
TPCLocalTrack::TPCLocalTrack()
  : m_is_fitted(false),
    m_is_calculated(false),
    m_fitflag(0),
    m_hit_array(),
    m_Ax(0.), m_Ay(0.), m_Au(0.), m_Av(0.),
    m_x0(0.), m_y0(0.), m_u0(0.), m_v0(0.),
    m_chisqr(1.e+10),
    m_n_iteration(0),
    m_searchtime(0), m_fittime(0),
    m_x0_exclusive(), m_y0_exclusive(),
    m_u0_exclusive(), m_v0_exclusive()
{
  m_hit_array.reserve(ReservedNumOfHits);
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCLocalTrack::~TPCLocalTrack()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCLocalTrack::AddTPCHit(TPCLTrackHit *hit)
{
  if(hit) m_hit_array.push_back(hit);
}

//_____________________________________________________________________________
void
TPCLocalTrack::Calculate()
{
  if(IsCalculated()){
    hddaq::cerr << FUNC_NAME << " already called" << std::endl;
    return;
  }

  for(const auto& hit: m_hit_array){
    hit->SetCalX0Y0(m_x0, m_y0);
    hit->SetCalUV(m_u0, m_v0);
    hit->SetCalPosition(hit->GetLocalCalPos());
  }
  m_is_calculated = true;
}

//_____________________________________________________________________________
void
TPCLocalTrack::CalculateExclusive()
{

  if(!IsCalculated()){
    hddaq::cerr << "#W " << FUNC_NAME << " No inclusive calculation" << std::endl;
    return;
  }

  const std::size_t n = m_hit_array.size();
  for(std::size_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    hit->SetCalX0Y0Exclusive(m_x0_exclusive[i], m_y0_exclusive[i]);
    hit->SetCalUVExclusive(m_u0_exclusive[i], m_v0_exclusive[i]);
    hit->SetCalPositionExclusive(hit->GetLocalCalPosExclusive());
  }
}

//_____________________________________________________________________________
void
TPCLocalTrack::ClearHits()
{
  m_hit_array.clear();
}

//_____________________________________________________________________________
void
TPCLocalTrack::DeleteNullHit()
{
  auto itr = m_hit_array.begin();
  while(itr != m_hit_array.end()){
    if(!*itr || !(*itr)->IsGood()){
      itr = m_hit_array.erase(itr);
    }else{
      itr++;
    }
  }
}

//______________________________________________________________________________
void
TPCLocalTrack::SetHoughFlag(int hough_flag)
{
  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TPCHit *hit = hitp->GetHit();
    if( !hit ) continue;
    hit->SetHoughFlag(hough_flag);
  }
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::DoFit(Int_t MinNumOfHits)
{

  Bool_t status = DoLinearFit(MinNumOfHits);
  m_is_fitted = status;
  if(m_chisqr<MaxChisqr) return status;
  else return false;
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::DoLinearFit(Int_t MinNumOfHits)
{

  const Int_t n = m_hit_array.size();
  std::vector<Int_t> pads;
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hit = m_hit_array[i];
    pads.push_back(hit -> GetHit() -> GetLayer());
  }
  std::sort(pads.begin(),pads.end());
  pads.erase(std::unique(pads.begin(),pads.end()),pads.end());
  if(pads.size()<MinNumOfHits || GetNDF()<1) return false;

  gNumOfHits = n;
  gHitPos.clear();
  gRes.clear();

  //hough translation for ini-param of Ay and Av
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 res = hitp->GetResolutionVect();
    gHitPos.push_back(pos);
    gRes.push_back(res);
  }

  //YZ Line Hough-transform
  Double_t LinearPar[2]; Int_t MaxBin[3];
  tpc::HoughTransformLineYZ(gHitPos, MaxBin, LinearPar, MinNumOfHits);
  m_Ay = LinearPar[0] + LinearPar[1]*tpc::ZTarget;
  m_Av = LinearPar[1];

  gPar[0] = m_Ax;
  gPar[1] = m_Ay;
  gPar[2] = m_Au;
  gPar[3] = m_Av;
  gChisqr = CalcChi2(gPar);

#if DebugDisp
  std::cout<<"Before linear fitting"
	   <<" m_chisqr: "<<gChisqr
	   <<", m_x0: "<<m_Ax
	   <<", m_y0: "<<m_Ay
	   <<", m_u0: "<<m_Au
	   <<", m_v0: "<<m_Av
	   <<std::endl;
#endif

  //Track fitting
  LinearFit();

  m_chisqr = gChisqr;
  m_x0 = gPar[0];
  m_y0 = gPar[1];
  m_u0 = gPar[2];
  m_v0 = gPar[3];

#if 0
  Double_t chisqr1 = 100000., chisqr2=100000., chisqr3=100000.;
  chisqr1 = m_chisqr;
  Double_t m_Ax_r = m_Ax*TMath::Cos(m_Au_atan);
  Double_t m_Ay_r = m_Ay*TMath::Cos(m_Av_atan);
  Double_t m_Au_theta = TMath::ATan(-1./m_Au);
  Double_t m_Av_theta = TMath::ATan(-1./m_Av);
  Double_t m_Ax_r2 = m_Ax*TMath::Sin(m_Au_theta);
  Double_t m_Ay_r2 = m_Ay*TMath::Sin(m_Av_theta);
  Double_t par2[4]={m_Ax_r, m_Au_atan, m_Ay_r, m_Av_atan};
  Double_t par3[4]={m_Ax_r2, m_Au_theta, m_Ay_r2, m_Av_theta};

  TMinuit *m_minuit = new TMinuit(4);
  m_minuit->SetFCN(fcn2_rt2);
  for(Int_t i=0; i<4; i++){
    m_minuit->mnparm(i, name[i], par2[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  for(Int_t i=0; i<4; i++){
    m_minuit->mnpout(i, name[i], par2[i], err[i], bnd1, bnd2, Err);
  }
  m_x0=par2[0]/TMath::Cos(par2[1]);
  m_u0=tan(par2[1]);
  m_y0=par2[2]/TMath::Cos(par2[3]);
  m_v0=tan(par2[3]);
  CalcChisquare();
  chisqr2 = m_chisqr;

  m_minuit->SetFCN(fcn2_rt);
  for(Int_t i = 0; i<4; i++){
    m_minuit->mnparm(i, name[i], par3[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  for(Int_t i=0; i<4; i++){
    m_minuit->mnpout(i, name[i], par3[i], err[i], bnd1, bnd2, Err);
  }
  delete m_minuit;

  m_x0=par3[0]/TMath::Sin(par3[1]);
  m_u0=-TMath::Cos(par3[1])/TMath::Sin(par3[1]);
  m_y0=par3[2]/TMath::Sin(par3[3]);
  m_v0=-TMath::Cos(par3[3])/TMath::Sin(par3[3]);

  CalcChisquare();
  chisqr3 = m_chisqr;

  if(fabs(chisqr1)<fabs(chisqr3)&&chisqr1>0.){
    m_x0=par[0];
    m_u0=tan(par[1]);
    m_y0=par[2];
    m_v0=tan(par[3]);
    CalcChisquare();
  }
  if(fabs(chisqr2)<fabs(chisqr3)&&fabs(chisqr1)>fabs(chisqr2)&&chisqr2>0.){
    m_x0=par2[0]/TMath::Cos(par2[1]);
    m_u0=tan(par2[1]);
    m_y0=par2[2]/TMath::Cos(par2[3]);
    m_v0=tan(par2[3]);
    CalcChisquare();
  }
  std::cout<<"chisqr1="<<chisqr1<<", chisqr2="<<chisqr2<<", chisqr3="<<chisqr3<<std::endl;
  std::cout<<"m_chisqr="<<m_chisqr<<std::endl;
#endif

  if(m_chisqr > MaxChisqr || TMath::IsNaN(m_chisqr)) return false;

  Int_t false_layer =0;
  for(Int_t i=0; i<m_hit_array.size(); ++i){
    auto hit = m_hit_array[i];
    TVector3 pos = hit->GetLocalHitPos();
    TVector3 Res = hit->GetResolutionVect();
    if(!ResidualIsWithinResolution(pos,Res)){
      m_hit_array.erase(m_hit_array.begin()+i);
      ++false_layer;
      --i;
    }
    if(m_hit_array.size()<MinNumOfHits) return false;
  }

  m_n_iteration++;
  if(false_layer == 0){
    return true;
  }
  else
    return DoLinearFit(MinNumOfHits);
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::ResidualIsWithinResolution(const TVector3& position,
                                          const TVector3& resolution)
{
  TVector3 x0(m_x0, m_y0, tpc::ZTarget);
  TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, tpc::ZTarget+1.);
  TVector3 u = (x1-x0).Unit();
  TVector3 AP = position-x0;
  Double_t dist_AX = u.Dot(AP);
  TVector3 AI(x0.x()+(u.x()*dist_AX),
	      x0.y()+(u.y()*dist_AX),
	      x0.z()+(u.z()*dist_AX));
  TVector3 residual = position-AI;
  if(residual.Mag() < resolution.Mag()*MaxResidual) return true;
  else return false;
}

//_____________________________________________________________________________
Double_t
TPCLocalTrack::GetTheta() const
{
  Double_t cost = 1./TMath::Sqrt(1.+m_u0*m_u0+m_v0*m_v0);
  return TMath::ACos(cost)*TMath::RadToDeg();
}

//_____________________________________________________________________________
void
TPCLocalTrack::Print(const TString& arg, bool print_allhits) const
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
Double_t
TPCLocalTrack::GetTrackdE()
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
TPCLocalTrack::SetParamUsingHoughParam()
{

  DeleteNullHit();

  gPar[0] = m_Ax;
  gPar[1] = m_Ay;
  gPar[2] = m_Au;
  gPar[3] = m_Av;
  gChisqr = CalcChi2(gPar);

  m_x0 = gPar[0];
  m_y0 = gPar[1];
  m_u0 = gPar[2];
  m_v0 = gPar[3];
  m_chisqr = gChisqr;

#if DebugDisp
  std::cout<<" m_chisqr: "<<m_Chisqr
	   <<", m_x0: "<<m_Ax
	   <<", m_y0: "<<m_Ay
	   <<", m_u0: "<<m_Au
	   <<", m_v0: "<<m_Av
	   <<std::endl;
#endif

}

//_____________________________________________________________________________
void
TPCLocalTrack::DoLinearFitExclusive()
{
  const Int_t n = m_hit_array.size();
  m_x0_exclusive.resize(n);
  m_y0_exclusive.resize(n);
  m_u0_exclusive.resize(n);
  m_v0_exclusive.resize(n);

  gNumOfHits = n-1;
  for(Int_t i=0; i<n; ++i){
    gHitPos.clear();
    gRes.clear();
    gHitPos.resize(n-1);
    gRes.resize(n-1);

    gChisqr = 1.e+10;
    gPar[0] = m_x0;
    gPar[1] = m_y0;
    gPar[2] = m_u0;
    gPar[3] = m_v0;

    int flag=0;
    for(Int_t j=0; j<n; ++j){
      if(j==i) continue; //exclude ith hit

      TPCLTrackHit *hitp = m_hit_array[j];
      TVector3 pos = hitp->GetLocalHitPos();
      TVector3 res = hitp->GetResolutionVect();
      gHitPos[flag] = pos;
      gRes[flag] = res;
      flag++;
    } //j

    //Track fitting
    LinearFit();

    m_chisqr = gChisqr;
    m_x0_exclusive[i] = gPar[0];
    m_y0_exclusive[i] = gPar[1];
    m_u0_exclusive[i] = gPar[2];
    m_v0_exclusive[i] = gPar[3];
  } //i
}
