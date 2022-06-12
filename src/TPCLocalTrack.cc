// -*- C++ -*-

//Comment by Ichikawa
//TPCLocalTrack.cc is for lenear fit
//TPCLocalTrack_Helix.cc will be prepared for Helix tracking

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

#include "DCAnalyzer.hh"
#include "DCLTrackHit.hh"
#include "DCGeomMan.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "HodoParamMan.hh"
#include "TPCLTrackHit.hh"
#include "TPCPadHelper.hh"
#include "UserParamMan.hh"

#define DebugEvDisp 1

#if DebugEvDisp
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TSystem.h>
namespace { TApplication app("DebugApp", nullptr, nullptr); }
#endif


namespace
{
Int_t    gNumOfHits;
std::vector<TVector3> gHitPos;
std::vector<TVector3> gRes;
}

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const auto& gGeom = DCGeomMan::GetInstance();
const auto& gHodo = HodoParamMan::GetInstance();
const Double_t& zK18tgt = gGeom.LocalZ("K18Target");
const Double_t& zTgt    = gGeom.LocalZ("Target");
// Temporary
//  const Int_t MaxTry    = 10;

const Int_t ReservedNumOfHits  = 64;
static const Double_t  FitStep[TPCLocalTrack::NumOfParam] = { 1.0e-6, 1.0e-10, 1.0e-6, 1.0e-10};
//static const Double_t  FitStep[TPCLocalTrack::NumOfParam] = { 1.0e-2, 1.0e-6, 1.0e-2, 1.0e-6};
static const Double_t  LowLimit[TPCLocalTrack::NumOfParam] = { -400., -1.0*TMath::ACos(-1.), -400, -1.0*TMath::ACos(-1.) };
static const Double_t  UpLimit[TPCLocalTrack::NumOfParam] = { 400., 1.0*TMath::ACos(-1.), 400, 1.0*TMath::ACos(-1.) };
//  static const Double_t  MaxChisqr = 500.;
static const Double_t  MaxChisqr = 10000.;

const Int_t    theta_ndiv = 100;
const Double_t theta_min  =   0;
const Double_t theta_max  = 180;
const Int_t    r_ndiv =  100;
const Double_t r_min  = -600;
const Double_t r_max  =  600;

//  static const Double_t  res_check_factor = 5.;
static const Double_t  res_check_factor = 10.;

// define the parametric line equation
[[maybe_unused]] void
Line(Double_t t, const Double_t *p, Double_t &x, Double_t &y, Double_t &z)
{
   x = p[0] + p[1]*t;
   y = p[2] + p[3]*t;
   z = t;
}

//_____________________________________________________________________________
// function Object to be minimized
struct Chisquared
{
  // the TGraph is a data member of the object
  TGraph2D *graph;
  std::vector<TVector3> Positions;
  std::vector<TVector3> Resolutions;

  Chisquared() : graph(new TGraph2D){
    graph->SetTitle("TPCLocalTrack::Chisquared");
    graph->GetXaxis()->SetRangeUser(-400, 400);
    graph->GetYaxis()->SetRangeUser(-400, 400);
    graph->GetZaxis()->SetRangeUser(-400, 400);
  }

  void Add(const TVector3& pos, const TVector3& res){
    Positions.push_back(pos);
    Resolutions.push_back(res);
  }

  Double_t Distance2(Int_t index, const Double_t *p){
    TVector3 xp = Positions.at(index);
    TVector3 x0(p[0], p[2], 0.);
    TVector3 x1(p[0] + p[1], p[2] + p[3], 1.);
    TVector3 u = (x1-x0).Unit();
    Double_t d2 = ((xp-x0).Cross(u)).Mag2();
    return d2;
  }

  Double_t Resolution2(Int_t index, const Double_t *p){
    TVector3 rp = Resolutions.at(index);
    return 0.3;
  }

  // implementation of the function to be minimized
  Double_t operator() (const Double_t *par){
    // assert(graph != nullptr);
    // Double_t* x = graph->GetX();
    // Double_t* y = graph->GetY();
    // Double_t* z = graph->GetZ();
    Double_t chisqr = 0;
    for(Int_t i=0, n=graph->GetN(); i<n; ++i){
      chisqr += Distance2(i, par)/Resolution2(i, par);
    }
    // std::cout << "Total Initial distance square = " << sum << std::endl;
    return chisqr;
  }
};

}

//_____________________________________________________________________________
TPCLocalTrack::TPCLocalTrack()
  : m_is_fitted(false),
    m_is_calculated(false),
    m_hit_array(),
    m_Ax(), m_Ay(), m_Au(), m_Av(),
    m_x0(), m_y0(),
    m_u0(), m_v0(),
    m_a(),  m_b(),
    m_chisqr(TMath::QuietNaN()),
    m_good_for_tracking(true),
    m_n_iteration(),
    m_de(),
    m_minuit(new TMinuit(NumOfParam))
{
  m_hit_array.reserve(ReservedNumOfHits);
  // m_cluster_array.reserve(ReservedNumOfHits);
  TROOT minexam("LinearFit","linear fit using TMinuit");
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCLocalTrack::~TPCLocalTrack()
{
  debug::ObjectCounter::decrease(ClassName());
  delete m_minuit;
}

//_____________________________________________________________________________
void
TPCLocalTrack::ClearHits()
{
  m_hit_array.clear();
  // m_cluster_array.clear();
}


//_____________________________________________________________________________
void
TPCLocalTrack::AddTPCHit(TPCLTrackHit *hit)
{
  if(hit) m_hit_array.push_back(hit);
}

//_____________________________________________________________________________
// void
// TPCLocalTrack::AddTPCCluster(TPCCluster *cluster)
// {
//   //not supported
//   if(cluster) m_cluster_array.push_back(cluster);
// }

//_____________________________________________________________________________
void
TPCLocalTrack::Calculate()
{
  if(m_is_calculated){
    hddaq::cerr << FUNC_NAME << " "
		<< "already called" << std::endl;
    return;
  }

  for(const auto& hit: m_hit_array){
    hit->SetCalX0Y0(m_x0, m_y0);
    hit->SetCalUV(m_u0, m_v0);
    hit->SetCalPosition(hit->GetLocalCalPos());
  }
  Print();
  m_is_calculated = true;
}

//_____________________________________________________________________________
static void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
                 Int_t iflag)
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
    chisqr += TVector3(d.x()/Res.x(), d.y()/Res.y(), d.z()/Res.z()).Mag();
  }
  f = chisqr/(gNumOfHits-TPCLocalTrack::NumOfParam);
}


//_____________________________________________________________________________
static void fcn2_rt(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
                    Int_t iflag)
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
    chisqr += TVector3(d.x()/Res.x(), d.y()/Res.y(), d.z()/Res.z()).Mag();
  }
  f = chisqr/(gNumOfHits-TPCLocalTrack::NumOfParam);
}

//_____________________________________________________________________________
static void fcn2_rt2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
                     Int_t iflag)
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
    chisqr += TVector3(d.x()/Res.x(), d.y()/Res.y(), d.z()/Res.z()).Mag();
  }
  f = chisqr/(gNumOfHits-TPCLocalTrack::NumOfParam);
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

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::DoFit()
{
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 800, 800);
  c1.cd();

  for(const auto& hit: m_hit_array){
    gHitPos.push_back(hit->GetLocalHitPos());
   gRes.push_back(hit->GetResolutionVect());
  }

  //using namespace ROOT::Math;

  // double xmin = 0; double ymin = 0;
  // double xmax = 10; double ymax = 10;

  Chisquared chisqr;
  // Double_t p0[NumOfParam] = { m_Ax, m_Au, m_Ay, m_Av };
  for(Int_t i=0, n=m_hit_array.size(); i<n; ++i){
    auto hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    // Line(t, p0, pos.X(), pos.Y(), pos.Z());
    chisqr.graph->SetPoint(i, pos.X(), pos.Y(), pos.Z());
    //dt->SetPointError(N,0,0,err);
  }

  ROOT::Fit::Fitter fitter;
  ROOT::Math::Functor fcn(chisqr, NumOfParam);
  // set the function and the initial parameter values
  Double_t pStart[NumOfParam] = { m_Ax, m_Au, m_Ay, m_Av };
  fitter.SetFCN(fcn, pStart);
  // set step sizes different than default ones (0.3 times parameter values)
  fitter.Config().ParSettings(0).SetStepSize(1e-3);
  fitter.Config().ParSettings(1).SetStepSize(1e-5);
  fitter.Config().ParSettings(2).SetStepSize(1e-3);
  fitter.Config().ParSettings(3).SetStepSize(1e-5);
  Bool_t ok = fitter.FitFCN();
  if(!ok){
    Error("line3Dfit","Line3D Fit failed");
    return false;
  }

  const ROOT::Fit::FitResult& result = fitter.Result();

  hddaq::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
  result.Print(hddaq::cout);

  chisqr.graph->Draw("p0");

  // get fit parameters
  // const auto parFit = result.GetParams();

  // draw the fitted line
  // int n = 1000;
  // double t0 = 0;
  // double dt = 10;
  // TPolyLine3D *l = new TPolyLine3D(n);
  // for (int i = 0; i <n;++i) {
  //   double t = t0+ dt*i/n;
  //   double x,y,z;
  //   line(t,parFit,x,y,z);
  //   l->SetPoint(i,x,y,z);
  // }
  // l->SetLineColor(kRed);
  // l->Draw("same");

  // // draw original line
  // TPolyLine3D *l0 = new TPolyLine3D(n);
  // for (int i = 0; i <n;++i) {
  //   double t = t0+ dt*i/n;
  //   double x,y,z;
  //   line(t,p0,x,y,z);
  //   l0->SetPoint(i,x,y,z);
  // }
  // l0->SetLineColor(kBlue);
  // l0->Draw("same");
#if DebugEvDisp
  c1.Modified();
  c1.Update();
  getchar();
  // app.Run();
#endif
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::DoFit(Int_t min_hits)
{
  if(m_is_fitted){
    hddaq::cerr << FUNC_NAME << " "
		<< "already called" << std::endl;
    return false;
  }
  DeleteNullHit();

  Bool_t status = DoFitLinear(min_hits);
  m_is_fitted = status;
  return status;
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::DoFitLinear(Int_t min_hits)
{
#if DebugEvDisp
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 400, 800);
  c1.Clear();
  c1.Divide(1, 2);
  TGraphErrors gxz, gyz;
  gxz.SetTitle(FUNC_NAME+" X%Z");
  gyz.SetTitle(FUNC_NAME+" Y%Z");
  // gxz.SetMarkerStyle(8);
  // gyz.SetMarkerStyle(8);
#endif

  // hddaq::cout << FUNC_NAME << std::endl;

  const Int_t n = m_hit_array.size();

  if(n<min_hits || GetNDF()<1){
    // hddaq::cerr << FUNC_NAME << " "
    //             << "Min layer should be > NDF" << std::endl;
    return false;
  }

  gNumOfHits = n;
  gHitPos.clear();
  gRes.clear();
  gHitPos.resize(n);
  gRes.resize(n);

  // r = x * TMath::Cos(theta) + y * TMath::Sin(theta)
  static TH2D hist("hist",";theta (deg.); r (mm)",
		   theta_ndiv, theta_min, theta_max, r_ndiv, r_min, r_max);
  hist.Reset();

  ///// Initial parameter using least squares method
  // Double_t lsm_x, lsm_y, lsm_u, lsm_v;
  // Double_t szx, szy, sz, sx, sy, sz2;
  // for(Int_t i=0; i<n; ++i){
  //   auto hitp = m_hit_array[i];
  //   auto pos = hitp->GetLocalHitPos();
  //   auto Res = hitp->GetResolutionVect();
  //   szx += pos.Z()*pos.X();
  //   szy += pos.Z()*pos.Y();
  //   sz  += pos.Z();
  //   sx  += pos.X();
  //   sy  += pos.Y();
  //   sz2 += pos.Z()*pos.Z();
  // }
  // lsm_u = (n*szx - sz*sx)/(n*sz2 - sz*sz);
  // lsm_x = (sz2*sx - sz*szx)/(n*sz2 - sz*sz);
  // lsm_v = (n*szy - sz*sy)/(n*sz2 - sz*sz);
  // lsm_y = (sz2*sy - sz*szy)/(n*sz2 - sz*sz);

  //hough translation for ini-param of Ay and Av
  for(Int_t i=0; i<n; ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 res = hitp->GetResolutionVect();
    gHitPos[i] = pos;
    gRes[i] = res;
    for(Int_t ti=0; ti<theta_ndiv; ti++){
      Double_t theta = theta_min+ti*(theta_max-theta_min)/theta_ndiv;
      hist.Fill(theta, TMath::Cos(theta*TMath::ACos(-1)/180.)*pos.Z()
                +TMath::Sin(theta*TMath::ACos(-1)/180.)*pos.Y());
      if(fabs(TMath::Cos(theta*TMath::ACos(-1)/180.)*pos.Z()+TMath::Sin(theta*TMath::ACos(-1)/180.)*pos.Y())>r_max)
	std::cout<<"HoughY: out of range:"<<TMath::Cos(theta*TMath::ACos(-1)/180.)*pos.Z()+TMath::Sin(theta*TMath::ACos(-1)/180.)*pos.Y()<<std::endl;
    }
    gxz.SetPoint(i, pos.Z(), pos.X());
    gyz.SetPoint(i, pos.Z(), pos.Y());
    gxz.SetPointError(i, res.Z(), res.X());
    gyz.SetPointError(i, res.Z(), res.Y());
  }
  Int_t maxbin = hist.GetMaximumBin();
  Int_t mx,my,mz;
  hist.GetBinXYZ(maxbin, mx, my, mz);
  Double_t mtheta = hist.GetXaxis()->GetBinCenter(mx)*TMath::ACos(-1)/180.;
  Double_t mr = hist.GetYaxis()->GetBinCenter(my);
  Double_t p0 = mr/TMath::Sin(mtheta);
  Double_t p1 = -TMath::Cos(mtheta)/TMath::Sin(mtheta);

  Double_t m_Ay = p0+p1*tpc::ZTarget;
  Double_t m_Av = p1;
  // std::cout<<"Hough x param:"<<m_Ax<<", y param:"<<m_Ay<<std::endl;
  // std::cout<<"Hough u param:"<<m_Au<<", v param:"<<m_Av<<std::endl;

  Double_t m_Au_atan = TMath::ATan(m_Au);
  Double_t m_Av_atan = TMath::ATan(m_Av);

  Double_t m_Ax_r = m_Ax*TMath::Cos(m_Au_atan);
  Double_t m_Ay_r = m_Ay*TMath::Cos(m_Av_atan);

  Double_t m_Au_theta = TMath::ATan(-1./m_Au);
  Double_t m_Av_theta = TMath::ATan(-1./m_Av);
  Double_t m_Ax_r2 = m_Ax*TMath::Sin(m_Au_theta);
  Double_t m_Ay_r2 = m_Ay*TMath::Sin(m_Av_theta);

  Double_t par[NumOfParam]={m_Ax, m_Au_atan, m_Ay, m_Av_atan};
  Double_t par2[NumOfParam]={m_Ax_r, m_Au_atan, m_Ay_r, m_Av_atan};
  Double_t par3[NumOfParam]={m_Ax_r2, m_Au_theta, m_Ay_r2, m_Av_theta};
  Double_t err[NumOfParam]={-999.,-999.,-999.,-999.};

  // m_minuit = new TMinuit(NumOfParam);
  // TROOT minexam("LinearFit", "Linear fit using TMinuit");

  m_minuit->SetPrintLevel(-1);
  m_minuit->SetFCN(fcn2);

  Int_t ierflg = 0;
  Double_t arglist[10];
  arglist[0] = 1;
  m_minuit->mnexcm("SET NOW", arglist,1,ierflg); // No warnings
  TString name[NumOfParam] = {"x0_r", "atan_u0", "y0_r", "atan_v0"};
  for(Int_t i = 0; i<NumOfParam; i++)
  {
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
  //m_minuit->mnprin(NumOfParam, amin);
  Int_t Err;
  Double_t bnd1, bnd2;
  for(Int_t i=0; i<NumOfParam; i++)
  {
    m_minuit->mnpout(i, name[i], par[i], err[i], bnd1, bnd2, Err);
    //std::cout<<Par[i]<<"  "<<std::endl;
  }

  Double_t chisqr1 = 100000., chisqr2=100000., chisqr3=100000.;

  m_x0=par[0];
  m_u0=tan(par[1]);
  m_y0=par[2];
  m_v0=tan(par[3]);
  CalcChi2();

  chisqr1 = m_chisqr;

  m_minuit->SetFCN (fcn2_rt2);
  for(Int_t i = 0; i<NumOfParam; i++)
  {
    m_minuit->mnparm(i, name[i], par2[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  for(Int_t i=0; i<NumOfParam; i++)
  {
    m_minuit->mnpout(i, name[i], par2[i], err[i], bnd1, bnd2, Err);
  }
  m_x0=par2[0]/TMath::Cos(par2[1]);
  m_u0=tan(par2[1]);
  m_y0=par2[2]/TMath::Cos(par2[3]);
  m_v0=tan(par2[3]);
  CalcChi2();
  chisqr2 = m_chisqr;


  m_minuit->SetFCN (fcn2_rt);
  for(Int_t i = 0; i<NumOfParam; i++)
  {
    m_minuit->mnparm(i, name[i], par3[i], FitStep[i], LowLimit[i], UpLimit[i], ierflg);
  }
  m_minuit->mnexcm("MIGRAD", arglist, 2, ierflg);
  for(Int_t i=0; i<NumOfParam; i++)
  {
    m_minuit->mnpout(i, name[i], par3[i], err[i], bnd1, bnd2, Err);
  }
  m_x0=par3[0]/TMath::Sin(par3[1]);
  m_u0=-TMath::Cos(par3[1])/TMath::Sin(par3[1]);
  m_y0=par3[2]/TMath::Sin(par3[3]);
  m_v0=-TMath::Cos(par3[3])/TMath::Sin(par3[3]);

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
    m_x0=par2[0]/TMath::Cos(par2[1]);
    m_u0=tan(par2[1]);
    m_y0=par2[2]/TMath::Cos(par2[3]);
    m_v0=tan(par2[3]);
    CalcChi2();
  }

  // m_x0 += m_u0*tpc::ZTarget;
  // m_y0 += m_v0*tpc::ZTarget;

  // std::cout<<"chisqr1="<<chisqr1<<", chisqr2="<<chisqr2<<", chisqr3="<<chisqr3<<std::endl;
  // std::cout<<"m_chisqr="<<m_chisqr<<std::endl;

  if(m_chisqr > MaxChisqr||std::isnan(m_chisqr))
    return false;

  Int_t false_layer =0;

  for(std::size_t i=0; i<m_hit_array.size(); ++i){
    TPCLTrackHit *hitp = m_hit_array[i];
    TVector3 pos = hitp->GetLocalHitPos();
    TVector3 Res = hitp->GetResolutionVect();
    if(!Residual_check(pos,Res)){
      // std::cout<<"false layer:"<<i<<", n="<<m_hit_array.size()<<std::endl;
      // std::cout<<"m_x0:"<<m_x0<<", m_y0="<<m_y0<<std::endl;
      // std::cout<<"m_u0:"<<m_u0<<", m_v0="<<m_v0<<", m_chisqr="<<m_chisqr<<std::endl;
      m_hit_array.erase(m_hit_array.begin()+i);
      ++false_layer;
      --i;
    }
    if(m_hit_array.size()<min_hits)
      return false;
  }


  if(false_layer ==0){
    //std::cout<<"return true"<<std::endl;
    // hddaq::cout << "hough init param: u=" << m_Au << " x=" << m_Ax << " v="<< m_Av << " y=" << m_Ay << std::endl;
    // hddaq::cout << "lsm   init param: u=" << lsm_u << " x=" << lsm_x << " v="<< lsm_v << " y=" << lsm_y << std::endl;
    // hddaq::cout << "fitted param    : u0=" << m_u0 << " x0=" << m_x0 << " v="<< m_v0 << " y=" << m_y0 << std::endl;
#if DebugEvDisp
    {
      TLine lxz, lyz;
      Double_t x0 = m_x0 + -400.*m_u0;
      Double_t x1 = m_x0 +  400.*m_u0;
      Double_t y0 = m_y0 + -400.*m_v0;
      Double_t y1 = m_y0 +  400.*m_v0;
      gxz.GetXaxis()->SetLimits(-400, 400);
      gxz.GetYaxis()->SetRangeUser(-400, 400);
      gyz.GetXaxis()->SetLimits(-400, 400);
      gyz.GetYaxis()->SetRangeUser(-400, 400);
      c1.cd(1);
      gxz.Draw("AP");
      lxz.DrawLine(-400, x0, 400, x1);
      c1.cd(2);
      gyz.Draw("AP");
      lyz.DrawLine(-400, y0, 400, y1);
      c1.Modified();
      c1.Update();
      c1.Print("c1.pdf");
      getchar();
    }
#endif
    return true;
  }
  else
    return DoFitLinear(min_hits);
}

//_____________________________________________________________________________
Bool_t
TPCLocalTrack::Residual_check(TVector3 pos, TVector3  Res)
{
  Bool_t status_rescheck=false;
  // TVector3 x0(m_x0, m_y0, 0.);
  // TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, 1.);
  TVector3 x0(m_x0, m_y0, tpc::ZTarget);
  TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, tpc::ZTarget+1.);
  TVector3 u = (x1-x0).Unit();
  TVector3 AP = pos-x0;
  Double_t dist_AX = u.Dot(AP);
  TVector3 AI(x0.x()+(u.x()*dist_AX),
	      x0.y()+(u.y()*dist_AX),
	      x0.z()+(u.z()*dist_AX));
  TVector3 d = pos-AI;

  if(d.Mag()<Res.Mag()*res_check_factor)
    status_rescheck = true;
  else{
    // std::cout<<"false hit"<<std::endl;
    // std::cout<<"hitpos:"<<pos<<", calcpos:"<<AI<<std::endl;
    // std::cout<<"dMag:"<<d.Mag()<<", Res:"<<Res.Mag()<<", d:"<<d<<std::endl;
  }

  return status_rescheck;
}


//_____________________________________________________________________________
void
TPCLocalTrack::CalcChi2()
{
  Double_t chisqr = 0.;
  for(const auto& hit: m_hit_array){
    const auto& pos = hit->GetLocalHitPos();
    const auto& Res = hit->GetResolutionVect();
    // TVector3 x0(m_x0, m_y0, 0.);
    // TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, 1.);
    TVector3 x0(m_x0, m_y0, tpc::ZTarget);
    TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, tpc::ZTarget+1.);
    TVector3 u = (x1-x0).Unit();
    TVector3 AP = pos-x0;
    Double_t dist_AX = u.Dot(AP);
    TVector3 AI(x0.x()+(u.x()*dist_AX),
		x0.y()+(u.y()*dist_AX),
		x0.z()+(u.z()*dist_AX));
    TVector3 d = pos-AI;
    //    TVector3 d = (pos-x0).Cross(u);
    chisqr += TVector3(d.x()/Res.x(), d.y()/Res.y(), d.z()/Res.z()).Mag();
  }
  m_chisqr = chisqr/GetNDF();
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
TPCLocalTrack::Print(const TString& arg, std::ostream& ost) const
{
  PrintHelper helper(3, std::ios::fixed, ost);

  const Int_t w = 8;
  ost << FUNC_NAME << " " << arg << std::endl
      << " X0 : " << std::setw(w) << std::left << m_x0
      << " Y0 : " << std::setw(w) << std::left << m_y0
      << " U0 : " << std::setw(w) << std::left << m_u0
      << " V0 : " << std::setw(w) << std::left << m_v0;
  ost << " Chisqr : " << std::setw(w) << m_chisqr << std::endl;
  for(Int_t i=0, n=m_hit_array.size(); i<n; ++i){
    auto hitp = m_hit_array[i];
    if(!hitp) continue;
    Int_t lnum = hitp->GetLayer();
    auto hitpos = hitp->GetLocalHitPos();
    auto calpos = hitp->GetLocalCalPos();
    auto res = hitp->GetResidualVect();
    ost << "[" << std::setw(2) << i << "]"
	<< " #"  << std::setw(2) << lnum
        << " " << hitpos << " " << calpos
	<< " -> " << res << std::endl;
  }
  ost << std::endl;
}
