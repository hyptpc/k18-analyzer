
// -*- C++ -*-
#include "PidCommon.hh"
#include "PidData.hh"

#include <TMath.h>
#include <iostream>
#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TKey.h>
#include <TROOT.h>
#include <TClass.h>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <fstream>
#include <array>
#include <std_ostream.hh>
#include <numeric>

#include "FuncName.hh"
#include "DeleteUtility.hh"

const std::array<TString, static_cast<size_t>(CorrGraph::Graph::COUNT)>
CorrGraph::GraphNames = {
  "gMeanM2", "gMeandEdx", "gSigM2", "gSigdEdx", "gRotAngle", "gYield",
  "gMeanM2Mom", "gMeandEdxMom", "gSigM2Mom", "gSigdEdxMom", "gRotAngleMom", "gYieldMom"
};

const std::array<TString, static_cast<size_t>(CorrFunc::Func::COUNT)>
CorrFunc::FuncNames = {
  "fMeanM2", "fMeandEdx", "fSigM2", "fSigdEdx", "fRotAngle", "fYield",
  "fMeanM2Mom", "fMeandEdxMom", "fSigM2Mom", "fSigdEdxMom", "fRotAngleMom", "fYieldMom"
};

Double_t pidfunc::RotGauss2D(Double_t* xy, Double_t* par){
  const Double_t x    = xy[0];
  const Double_t y    = xy[1];
  const Double_t x0   = par[0]; // m2
  const Double_t y0   = par[1]; // dEdx
  const Double_t sigx = par[2]; // sigma of m2
  const Double_t sigy = par[3]; // sigma of dEdx
  const Double_t th   = par[4];
  const Double_t tot  = par[5]; 
  // rotation
  const Double_t u = (x - x0)/sigx;
  const Double_t v = (y - y0)/sigy;
  const Double_t cost = TMath::Cos(th);
  const Double_t sint = TMath::Sin(th);
  const Double_t up =  u * cost + v * sint;    
  const Double_t vp = -u * sint + v * cost;
  // const Double_t sigxp = sigx * TMath::Cos(th) + sigy * TMath::Sin(th);    
  // const Double_t sigyp = -sigx * TMath::Sin(th) + sigy * TMath::Cos(th);    
  // exponential
  const Double_t e = TMath::Exp( -0.5*(up*up + vp*vp) );
  // normalization  1/(2pi*sigmax*sigmagy)
  const Double_t norm = 1.0 / (2.0 * TMath::Pi() * sigx * sigy);
  return tot * norm * e;
}

Double_t pidfunc::RotGauss2DFit(Double_t* xy, Double_t* par){
  const Double_t x    = xy[0];
  const Double_t y    = xy[1];
  const Double_t x0   = par[0]; // m2
  const Double_t y0   = par[1]; // dEdx
  const Double_t sigx = par[2]; // sigma of m2
  const Double_t sigy = par[3]; // sigma of dEdx
  const Double_t th   = par[4];
  const Double_t tot  = par[5];

  const double binwm2 = pidlikeli::binwm2;
  const double binwdedx = pidlikeli::binwdedx;    
    
  // rotation 
  const Double_t u = (x - x0)/sigx;
  const Double_t v = (y - y0)/sigy;
  const Double_t cost = TMath::Cos(th);
  const Double_t sint = TMath::Sin(th);
  const Double_t up =  u * cost + v * sint;
  const Double_t vp = -u * sint + v * cost;
  // const Double_t sigxp = sigx * TMath::Cos(th) + sigy * TMath::Sin(th);    
  // const Double_t sigyp = -sigx * TMath::Sin(th) + sigy * TMath::Cos(th);    
  // exponential
  const Double_t e = TMath::Exp( -0.5*(up*up + vp*vp) );
  // normalization  1/(2pi*sigmax*sigmagy)
  const Double_t norm = 1.0 / (2.0 * TMath::Pi() * sigx * sigy);
  return tot * norm * e * binwdedx * binwm2;
}  
  
Double_t pidfunc::RotGauss2DProjX(Double_t* x, Double_t* par){
  const Double_t y0   = par[1]; // sigma of 1/beta    
  const Double_t sigx = par[2]; // sigma of 1/beta
  const Double_t sigy = par[3]; // sigma of dEdx
  const Double_t th   = par[4];
  const Double_t tot  = par[5];

  const Double_t sinT = TMath::Sin(th);
  const Double_t cosT = TMath::Cos(th);
  // const Double_t sinT = TMath::Sin(0.);
  // const Double_t cosT = TMath::Cos(0.);    
  const double cv2xx = std::pow(sigx*cosT,2) + std::pow(sigy*sinT,2);
  const double norm = tot / std::sqrt(2*TMath::Pi()*cv2xx);
  const double dx   = x[0] - par[0];
  return norm * std::exp(-0.5*dx*dx / cv2xx);
}

Double_t pidfunc::RotGauss2DProjXFit(Double_t* x, Double_t* par){
  const Double_t y0   = par[1]; // sigma of 1/beta    
  const Double_t sigx = par[2]; // sigma of 1/beta
  const Double_t sigy = par[3]; // sigma of dEdx
  const Double_t th   = par[4];
  const Double_t tot  = par[5];

  const double binwm2 = pidlikeli::binwm2;

  const Double_t sinT = TMath::Sin(th);
  const Double_t cosT = TMath::Cos(th);
  // const Double_t sinT = TMath::Sin(0.);
  // const Double_t cosT = TMath::Cos(0.);    
  const double cv2xx = std::pow(sigx*cosT,2) + std::pow(sigy*sinT,2);
  const double norm = tot / std::sqrt(2*TMath::Pi()*cv2xx);
  const double dx   = x[0] - par[0];
  return norm * std::exp(-0.5*dx*dx / cv2xx) * binwm2;
  //return 400 * std::exp(-0.5*dx*dx / cv2xx);
  //double sx = 0.00242;
  //return norm * std::exp(-0.5*dx*dx/sx/sx);
}    

Double_t pidfunc::RotGauss2DProjY(Double_t* x, Double_t* par){
  const Double_t sigx = par[2]; // sigma of 1/beta
  const Double_t sigy = par[3]; // sigma of dEdx
  const Double_t th   = par[4];
  const Double_t tot  = par[5];

  const Double_t sinT = TMath::Sin(th);
  const Double_t cosT = TMath::Cos(th);
  const double cv2yy = std::pow(sigx*sinT,2) + std::pow(sigy*cosT,2);
  const double norm = tot / std::sqrt(2*TMath::Pi()*cv2yy);
  const double dx   = x[0] - par[1];
  return norm * std::exp(-0.5*dx*dx / cv2yy);
  //return 160  * std::exp(-0.5*dx*dx / cv2yy);
}

Double_t pidfunc::RotGauss2DProjYFit(Double_t* x, Double_t* par){
  const Double_t sigx = par[2]; // sigma of 1/beta
  const Double_t sigy = par[3]; // sigma of dEdx
  const Double_t th   = par[4];
  const Double_t tot  = par[5];

  const double binwdedx = pidlikeli::binwdedx;    

  const Double_t sinT = TMath::Sin(th);
  const Double_t cosT = TMath::Cos(th);
  const double cv2yy = std::pow(sigx*sinT,2) + std::pow(sigy*cosT,2);
  const double norm = tot / std::sqrt(2*TMath::Pi()*cv2yy);
  const double dx   = x[0] - par[1];
  return norm * std::exp(-0.5*dx*dx / cv2yy) * binwdedx;
  //return 160  * std::exp(-0.5*dx*dx / cv2yy);
}    
  
Double_t pidfunc::RotGauss2DBeta(Double_t* xy, Double_t* par){
  const Double_t x = xy[0];
  const Double_t y = xy[1];
  const Double_t x0   = 1/x; // 1/beta
  const Double_t y0   = Kinematics::CalcDedx(x); // dEdx
  const Double_t sigx = par[0]; // sigma of 1/beta
  const Double_t sigy = par[1]; // sigma of dEdx
  const Double_t th   = par[2];
  const Double_t tot  = par[3]; 
  // rotation
  const Double_t dx = x - x0;
  const Double_t dy = y - y0;
  const Double_t xp =  dx * TMath::Cos(th) + dy * TMath::Sin(th);
  const Double_t yp = -dx * TMath::Sin(th) + dy * TMath::Cos(th);
  const Double_t sigxp = sigx * TMath::Cos(th) + sigy * TMath::Sin(th);
  const Double_t sigyp = -sigx * TMath::Sin(th) + sigy * TMath::Cos(th);
  // exponential
  const Double_t e = TMath::Exp( -0.5*(xp*xp/(sigxp*sigxp) + yp*yp/(sigyp*sigyp)) );
  // normalization  1/(2pi*sigmax*sigmagy)
  const Double_t norm = 1.0 / (2.0 * TMath::Pi() * sigxp * sigyp);
  return tot * norm * e;
}

Double_t pidfunc::RotDoubleGauss2D(Double_t* xy, Double_t* par){
  const Double_t x    = xy[0];
  const Double_t y    = xy[1];
  const Double_t x0   = par[0]; // m2
  const Double_t y0   = par[1]; // dEdx
  const Double_t sigx = par[2]; // sigma of m2
  const Double_t sigy = par[3]; // sigma of dEdx
  const Double_t th   = par[4];
  const Double_t tot  = par[5];
  //
  const Double_t xt    = xy[2];
  const Double_t yt    = xy[3];
  const Double_t xt0   = par[6]; // m2
  const Double_t yt0   = par[7]; // dEdx
  const Double_t sigxt = par[8]; // sigma of m2
  const Double_t sigyt = par[9]; // sigma of dEdx
  const Double_t tht   = par[10];
  const Double_t tott  = par[11];
  
  // gauss 1 
  const Double_t u = (x - x0)/sigx;
  const Double_t v = (y - y0)/sigy;
  const Double_t cost = TMath::Cos(th);
  const Double_t sint = TMath::Sin(th);
  const Double_t up =  u * cost + v * sint;    
  const Double_t vp = -u * sint + v * cost;
  const Double_t e = TMath::Exp( -0.5*(up*up + vp*vp) );
  const Double_t norm = 1.0 / (2.0 * TMath::Pi() * sigx * sigy);
  const Double_t g1 = tot * norm * e;
  
  // gauss 2
  const Double_t ut = (xt - xt0)/sigxt;
  const Double_t vt = (yt - yt0)/sigyt;
  const Double_t costt = TMath::Cos(tht);
  const Double_t sintt = TMath::Sin(tht);
  const Double_t utp =  ut * costt + vt * sintt;
  const Double_t vtp = -ut * sintt + vt * costt;
  const Double_t et = TMath::Exp( -0.5*(utp*utp + vtp*vtp) );
  const Double_t normt = 1.0 / (2.0 * TMath::Pi() * sigxt * sigyt);
  const Double_t g2 = tott * normt * et;
  
  return g1 + g2;
}

Double_t pidfunc::RotFiveGauss2D(Double_t* xy, Double_t* par){
  // par[0]  - par[5]  :pion
  // par[6]  - par[11] :kaon
  // par[12] - par[17] :proton
  // par[18] - par[23] :deutron
  // par[24] - par[29] :electron  
  Double_t val = 0.;
  val += RotGauss2D(xy, &par[pidfunc::kNparamGauss*pidlikeli::kPion]); // pion
  val += RotGauss2D(xy, &par[pidfunc::kNparamGauss*pidlikeli::kKaon]); // kaon
  val += RotGauss2D(xy, &par[pidfunc::kNparamGauss*pidlikeli::kProton]); // proton
  val += RotGauss2D(xy, &par[pidfunc::kNparamGauss*pidlikeli::kDeutron]); // deutron
  val += RotGauss2D(xy, &par[pidfunc::kNparamGauss*pidlikeli::kElectron]); // electron
  return val;
}

Double_t pidfunc::SigmaDedx(Double_t* x, Double_t* par){
  const Double_t beta = x[0];
  Double_t c = par[0];
  Double_t ib = par[1]/beta/beta;
  return hypot(c,ib);
}

Double_t pidfunc::SigmaM2Pid(Double_t* x, int pid){
  return 0.;
}

double pidfunc::CalcSigM2(double mom, int pid)
{
  std::cout << "debug " << __FILE__ << " " << __LINE__ << " " << __func__ << std::endl;  
  Double_t beta = pidlikeli::MomToBetaPid(mom,pid);
  if (beta < 1e-9) { 
        return 0.1; 
  }
  std::array<double, 4> sigmaM2param;
  if(pid==pidlikeli::kPion)
    sigmaM2param = pidlikeli::sigmaM2paramPi;
  else if(pid==pidlikeli::kKaon)
    sigmaM2param = pidlikeli::sigmaM2paramK;
  else if(pid==pidlikeli::kProton)
    sigmaM2param = pidlikeli::sigmaM2paramP;
  else if(pid==pidlikeli::kDeutron)
    sigmaM2param = pidlikeli::sigmaM2paramD;
  else if(pid==pidlikeli::kElectron)
    sigmaM2param = pidlikeli::sigmaM2paramE;
  else
    sigmaM2param = {0.01,0.01,0.01,0.01};
  return ExpPlusPol1(&beta,sigmaM2param.data());  
}

double pidfunc::CalcSigdEdx(double mom, int pid)
{
  std::cout << "debug " << __FILE__ << " " << __LINE__ << " " << __func__ << std::endl;  
  Double_t beta = pidlikeli::MomToBetaPid(mom,pid);
  if (beta < 1e-9) { 
        return 0.1; 
  }  
  std::array<double, 2> sigmaDedxparam;
  if(pid==pidlikeli::kPion)
    sigmaDedxparam = pidlikeli::sigmaDedxparamPi;
  else if(pid==pidlikeli::kKaon)
    sigmaDedxparam = pidlikeli::sigmaDedxparamK;
  else if(pid==pidlikeli::kProton)
    sigmaDedxparam = pidlikeli::sigmaDedxparamP;
  else if(pid==pidlikeli::kDeutron)
    sigmaDedxparam = pidlikeli::sigmaDedxparamD;
  else if(pid==pidlikeli::kElectron)
    sigmaDedxparam = pidlikeli::sigmaDedxparamE;
  else
    sigmaDedxparam = {0.01,0.01};    
  return SigmaDedx(&beta,sigmaDedxparam.data());
}  

Double_t pidfunc::SigmaInvBeta(Double_t* x, Double_t* par){
  Double_t xx = x[0];
  return par[0]*TMath::Exp(par[1]*xx)          // p0, p1 : expo
    + par[2] + par[3]*xx;        // p2–3  : pol1
}

Double_t pidfunc::ExpPlusPol3(Double_t *x, Double_t *p) {
  Double_t xx = x[0];
  return p[0]*TMath::Exp(p[1]*xx)          // p0, p1 : expo
    + p[2] + p[3]*xx + p[4]*xx*xx       // p2–4  : pol2
    + p[5]*xx*xx*xx;                    // p5    : x³
}
Double_t pidfunc::ExpPlusPol1(Double_t *x, Double_t *p) {
  Double_t xx = x[0];
  return p[0]*TMath::Exp(p[1]*xx)          // p0, p1 : expo
    + p[2] + p[3]*xx;        // p2–4  : pol2
}

Double_t pidfunc::Pol1(Double_t *x, Double_t *p) {
  Double_t xx = x[0];
  return p[0] + p[1]*xx;        // p2–4  : pol2
}

Double_t pidfunc::RationalFunc1(Double_t *x, Double_t *p) {
  Double_t xx = x[0];
  Double_t m2 = p[0];
  Double_t nume = p[1]*xx + (1+p[2])*m2-p[1];
  Double_t deno = 1 + p[2]*xx;
  return nume/deno;
}

Double_t pidfunc::LogiFunc1(Double_t *x, Double_t *p) {
  Double_t xx = x[0];
  Double_t m2 = p[0];
  Double_t deno = 1 + TMath::Exp(-(xx-p[1])/p[2]);
  return m2/deno + p[3];
}

// struct CorrFunc
CorrFunc::CorrFunc() {
  using Func = CorrFunc::Func;
  functions[pidlikeli::scast(Func::MeanM2)]
    = new TF1("fMeanM2", pidfunc::Pol1, 0., pidlikeli::maxpoq, 2);
  functions[pidlikeli::scast(Func::MeandEdx)]
    = new TF1("fMeandEdx", pidfunc::Pol1, 0, pidlikeli::maxpoq, 2);
  functions[pidlikeli::scast(Func::SigM2)]
    = new TF1("fSigM2", pidfunc::ExpPlusPol1, 0, pidlikeli::maxpoq, 4);
  functions[pidlikeli::scast(Func::SigdEdx)]
    = new TF1("fSigdEdx", pidfunc::SigmaInvBeta, 0, pidlikeli::maxpoq, 4);
  functions[pidlikeli::scast(Func::Yield)]
    = new TF1("fYield", pidfunc::Pol1, 0, pidlikeli::maxpoq, 2);
}
CorrGraph::CorrGraph() {
  using Graph = CorrGraph::Graph;
  for (size_t i = 0; i < static_cast<size_t>(Graph::COUNT); ++i) {
    graphs[i] = new TGraphErrors();
  }    
}
