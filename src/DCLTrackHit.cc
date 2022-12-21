// -*- C++ -*-

#include "DCLTrackHit.hh"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <std_ostream.hh>

#include "DCAnalyzer.hh"
#include "FuncName.hh"
#include "MathTools.hh"

namespace
{
const Double_t qnan = TMath::QuietNaN();
}

//_____________________________________________________________________________
DCLTrackHit::DCLTrackHit(DCHit *hit, Double_t pos, Int_t nh)
  : m_hit(hit),
    m_nth_hit(nh),
    m_local_hit_pos(pos),
    m_cal_pos(qnan),
    m_xcal(qnan),
    m_ycal(qnan),
    m_ucal(qnan),
    m_vcal(qnan),
    m_honeycomb(false)
{
  debug::ObjectCounter::increase(ClassName());
  m_hit->RegisterHits(this);
}

//_____________________________________________________________________________
DCLTrackHit::DCLTrackHit(const DCLTrackHit& right)
  : m_hit(right.m_hit),
    m_nth_hit(right.m_nth_hit),
    m_local_hit_pos(right.m_local_hit_pos),
    m_cal_pos(right.m_cal_pos),
    m_xcal(right.m_xcal),
    m_ycal(right.m_ycal),
    m_ucal(right.m_ucal),
    m_vcal(right.m_vcal),
    m_honeycomb(right.m_honeycomb)
{
  m_hit->RegisterHits(this);
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCLTrackHit::~DCLTrackHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Double_t
DCLTrackHit::GetLocalCalPos() const
{
  Double_t angle = m_hit->GetTiltAngle()*TMath::DegToRad();
  return m_xcal*TMath::Cos(angle) + m_ycal*TMath::Sin(angle);
}

//_____________________________________________________________________________
Double_t
DCLTrackHit::GetResidual() const
{
  Double_t a    = GetTiltAngle()*TMath::DegToRad();
  Double_t dsdz = m_ucal*TMath::Cos(a)+m_vcal*TMath::Sin(a);
  Double_t coss = m_honeycomb ? TMath::Cos(TMath::ATan(dsdz)) : 1.;
  Double_t scal = GetLocalCalPos();
  Double_t wp   = GetWirePosition();
  Double_t ss   = wp+(m_local_hit_pos-wp)/coss;
  return (ss-scal)*coss;
}

//_____________________________________________________________________________
void
DCLTrackHit::Print(const TString& arg) const
{
  m_hit->Print(arg);
  hddaq::cout << "local_hit_pos " << m_local_hit_pos << std::endl
	      << "residual " << GetResidual() << std::endl
	      << "honeycomb " << m_honeycomb << std::endl;
}

//_____________________________________________________________________________
Bool_t
DCLTrackHit::ReCalc(Bool_t applyRecursively)
{
  if(applyRecursively){
    if(!m_hit->ReCalcDC(applyRecursively)) return false;
    // if(!m_hit->ReCalcDC(applyRecursively)   ||
    // 	!m_hit->ReCalcMWPC(applyRecursively)){
    //   return false;
    // }
  }

  Double_t wp = GetWirePosition();
  Double_t dl = GetDriftLength();

  if(m_local_hit_pos>wp)  m_local_hit_pos = wp+dl;
  else                    m_local_hit_pos = wp-dl;

  return true;
}
