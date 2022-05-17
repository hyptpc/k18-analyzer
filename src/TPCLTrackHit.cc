// -*- C++ -*-

#include "TPCLTrackHit.hh"

#include <string>
#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <TString.h>
#include <TF1.h>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DCAnalyzer.hh"
#include "MathTools.hh"

namespace
{
const double zTgtTPC = -143.;
const double& HS_field_0 = 0.9860;
const double& HS_field_Hall_calc = ConfMan::Get<Double_t>("HSFLDCALC");
const double& HS_field_Hall = ConfMan::Get<Double_t>("HSFLDHALL");

//for Helix tracking
//[0]~[4] are the Helix parameters,
//([5],[6],[7]) = (x, y, z)
std::string s_tmp="pow([5]-([0]+([3]*cos(x))),2)+pow([6]-([1]+([3]*sin(x))),2)+pow([7]-([2]+([3]*[4]*x)),2)";
//static TF1 fint("fint",s_tmp.c_str(),-10.,10.);
static TF1 fint("fint",s_tmp.c_str(),-4.,4.);
}

//______________________________________________________________________________
TPCLTrackHit::TPCLTrackHit(TPCHit *hit)
  : m_hit(hit),
    m_x0(TMath::QuietNaN()),
    m_y0(TMath::QuietNaN()),
    m_u0(TMath::QuietNaN()),
    m_v0(TMath::QuietNaN()),
    m_cx(TMath::QuietNaN()),
    m_cy(TMath::QuietNaN()),
    m_z0(TMath::QuietNaN()),
    m_r(TMath::QuietNaN()),
    m_dz(TMath::QuietNaN())
{
  m_local_hit_pos = hit->GetPosition();
  m_cal_pos = TVector3(0.,0.,0.);
  m_res = TVector3(hit->GetResolutionX(),
		   hit->GetResolutionY(),
		   hit->GetResolutionZ());
  debug::ObjectCounter::increase(ClassName());
  m_hit->RegisterHits(this);
}

//______________________________________________________________________________
TPCLTrackHit::TPCLTrackHit(const TPCLTrackHit& right)
  : m_hit(right.m_hit),
    m_x0(right.m_x0),
    m_y0(right.m_y0),
    m_u0(right.m_u0),
    m_v0(right.m_v0),
    m_cx(right.m_cx),
    m_cy(right.m_cy),
    m_z0(right.m_z0),
    m_r(right.m_r),
    m_dz(right.m_dz)
{
  m_local_hit_pos = right.m_local_hit_pos;
  m_cal_pos = right.m_cal_pos;
  m_res = right.m_res;
  m_hit->RegisterHits(this);
  debug::ObjectCounter::increase(ClassName());
}

//______________________________________________________________________________
TPCLTrackHit::~TPCLTrackHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

// For Helix fit
//______________________________________________________________________________
TVector3
TPCLTrackHit::GetHelixPosition(double par[5], double t) const
{
  //This is the eqation of Helix
  // double  x = p[0] + p[3]*cos(t+theta0);
  // double  y = p[1] + p[3]*sin(t+theta0);
  // double  z = p[2] + (p[4]*p[3]*t);
  double  x = par[0] + par[3]*cos(t);
  double  y = par[1] + par[3]*sin(t);
  double  z = par[2] + (par[4]*par[3]*t);

  return TVector3(x, y, z);
}



//______________________________________________________________________________
TVector3
TPCLTrackHit::GetLocalCalPos() const
{
  TVector3 pos = m_local_hit_pos;
  // TVector3 x0(m_x0, m_y0, 0.);
  // TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, 1.);
  //temp

  TVector3 x0(m_x0, m_y0, zTgtTPC);
  TVector3 x1(m_x0 + m_u0, m_y0 + m_v0, zTgtTPC+1.);
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
TPCLTrackHit::GetLocalCalPosHelix() const
{
  TVector3 pos(-m_local_hit_pos.X(),
	       m_local_hit_pos.Z() - zTgtTPC,
	       m_local_hit_pos.Y());

  double par[5]={m_cx, m_cy, m_z0, m_r, m_dz};
  double fpar[8];
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
  }
  fpar[5] = pos.X();
  fpar[6] = pos.Y();
  fpar[7] = pos.Z();

  fint.SetParameters(fpar);
  double min_t = fint.GetMinimumX();
  TVector3 fittmp = GetHelixPosition(par, min_t);
  TVector3 calpos(-fittmp.X(),
		   fittmp.Z(),
		   fittmp.Y()+zTgtTPC);
  return calpos;
}


//______________________________________________________________________________
double
TPCLTrackHit::GetTcal() const
{
  TVector3 pos(-m_local_hit_pos.X(),
	       m_local_hit_pos.Z() - zTgtTPC,
	       m_local_hit_pos.Y());

  double par[5]={m_cx, m_cy, m_z0, m_r, m_dz};
  double fpar[8];
  for(int ip=0; ip<5; ++ip){
    fpar[ip] = par[ip];
  }
  fpar[5] = pos.X();
  fpar[6] = pos.Y();
  fpar[7] = pos.Z();

  fint.SetParameters(fpar);
  double min_t = fint.GetMinimumX();

  return min_t;
}


//______________________________________________________________________________
TVector3
TPCLTrackHit::GetMomentumHelix() const
{
  TVector3 pos(-m_cal_pos.X(),
   	       m_cal_pos.Z() - zTgtTPC,
   	       m_cal_pos.Y());

  const double Const = 0.299792458; // =c/10^9
  //const double dMagneticField = 1.; //T, "-1" is needed. // Should be given by field param
  const double dMagneticField = HS_field_0*(HS_field_Hall/HS_field_Hall_calc);

  double t = (pos.Z()-m_z0)/(m_r*m_dz);
  double pt = fabs(m_r)*(Const*dMagneticField); // MeV/c
  //From here!!!!
  double tmp_px = pt*(-1.*sin(t));
  double tmp_py = pt*(cos(t));
  double tmp_pz = pt*(m_dz);
  double px = -tmp_px*0.001;
  double py = tmp_pz*0.001;
  double pz = tmp_py*0.001;

  return TVector3(px,py,pz);
}



//______________________________________________________________________________
TVector3
TPCLTrackHit::GetResidualVect() const
{
  return m_cal_pos - m_local_hit_pos;
}

//______________________________________________________________________________
double
TPCLTrackHit::GetResidual() const
{
  TVector3 Res = m_cal_pos - m_local_hit_pos;
  return Res.Mag();
}

//______________________________________________________________________________
Bool_t
TPCLTrackHit::ResidualCut() const
{
  Bool_t status = false;
  TVector3 Res = m_cal_pos - m_local_hit_pos;
  double resolution = m_hit->GetResolution();
  if(Res.Mag()<resolution*5.)
    status = true;
  return status;
}

//______________________________________________________________________________
void
TPCLTrackHit::Print(const std::string& arg) const
{
  m_hit->Print(arg);
  hddaq::cout << "local_hit_pos " << m_local_hit_pos << std::endl
	      << "residual " << GetResidual() << std::endl;
}
