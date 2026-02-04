// -*- C++ -*-

#include "CherenkovHit.hh"

#include <cmath>
#include <algorithm>

#include <TMath.h>

#include "DebugCounter.hh"
#include "DetectorID.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoRawHit.hh"

#include <std_ostream.hh>

namespace
{
const auto U = HodoRawHit::kUp;
const auto D = HodoRawHit::kDown;
const auto E = HodoRawHit::kExtra;
const auto& gHodo = HodoParamMan::GetInstance();
}

//_____________________________________________________________________________
CherenkovHit::CherenkovHit(HodoRawHit* rhit, Double_t max_time_diff)
  : HodoHit(rhit, max_time_diff),
    m_npe_high()
{
  m_npe_high.resize(m_n_ch);
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
CherenkovHit::~CherenkovHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
CherenkovHit::Calculate()
{
  if(m_is_calculated){
    hddaq::cerr << FUNC_NAME << " already calculated" << std::endl;
    return false;
  }

  if(!gHodo.IsReady()){
    hddaq::cerr << FUNC_NAME << " HodoParamMan must be initialized" << std::endl;
    return false;
  }

  Int_t id    = m_raw->DetectorId();
  Int_t plane = m_raw->PlaneId();
  Int_t seg   = m_raw->SegmentId();

  using data_t = std::vector<std::vector<Double_t>>;
  data_t leading(m_n_ch);
  data_t trailing(m_n_ch);
  data_t cleading(m_n_ch);   // No PHC: ctime=time. Same names as HodoHit for pipeline compatibility.
  data_t ctrailing(m_n_ch);

  for(Int_t ch=0; ch<m_n_ch; ++ch){
    // AdcHigh only (AdcLow unused, same as Hodo). GetNpe uses m_AHPContainer.
    for(const auto& adc: m_raw->GetArrayAdcHigh(ch)){
      Double_t npe = TMath::QuietNaN();
      if(gHodo.GetNpe(id, plane, seg, ch, adc, npe)){
        m_npe_high.at(ch).push_back(npe);
      }
    }
    // TDC→time. No PHC for Cherenkov, so cleading=leading.
    for(const auto& tdc: m_raw->GetArrayTdcLeading(ch)){
      Double_t time = TMath::QuietNaN();
      if(gHodo.GetTime(id, plane, seg, ch, tdc, time)){
        leading.at(ch).push_back(time);
        cleading.at(ch).push_back(time);
      }
    }
    for(const auto& tdc: m_raw->GetArrayTdcTrailing(ch)){
      Double_t time = TMath::QuietNaN();
      if(gHodo.GetTime(id, plane, seg, ch, tdc, time)){
        trailing.at(ch).push_back(time);
        ctrailing.at(ch).push_back(time);
      }
    }
    std::sort(leading.at(ch).begin(), leading.at(ch).end());
    std::sort(trailing.at(ch).begin(), trailing.at(ch).end());
    std::sort(cleading.at(ch).begin(), cleading.at(ch).end());
    std::sort(ctrailing.at(ch).begin(), ctrailing.at(ch).end());
  }

  // one-side readout
  if(m_n_ch == 1){
    m_time_leading.at(U)  = leading.at(U);
    m_ctime_leading.at(U) = cleading.at(U);
    m_is_clustered.resize(leading.at(U).size(), false);
  }
  // two-side readout
  else{
    for(Int_t ju=0; ju<(Int_t)leading.at(U).size(); ++ju){
      Double_t lu = leading.at(U).at(ju);
      Double_t clu = cleading.at(U).at(ju);
      for(Int_t jd=0; jd<(Int_t)leading.at(D).size(); ++jd){
        Double_t ld = leading.at(D).at(jd);
        if(TMath::Abs(lu - ld) < m_max_time_diff){
          Double_t cld = cleading.at(D).at(jd);
          m_time_leading.at(U).push_back(lu);
          m_time_leading.at(D).push_back(ld);
          m_ctime_leading.at(U).push_back(clu);
          m_ctime_leading.at(D).push_back(cld);
          m_is_clustered.push_back(false);
          break;
        }
      }
    }
  }

  // KVC: TDC uses kSUM only; .
  if(id == DetIdKVC){
    m_time_leading.at(E)  = leading.at(HodoRawHit::EChannelKVC::kSUM);
    m_ctime_leading.at(E) = cleading.at(HodoRawHit::EChannelKVC::kSUM);
  }

  for(Int_t ch=0; ch<m_n_ch; ++ch){
    m_time_trailing.at(ch)  = trailing.at(ch);
    m_ctime_trailing.at(ch) = ctrailing.at(ch);
  }

  m_is_calculated = true;
  // Valid hit: KVC uses E(kSUM). Others use U (one-side) or U/D (two-side).
  if(id == DetIdKVC)
    return (m_ctime_leading.at(E).size() > 0);
  return (m_ctime_leading.at(U).size() > 0);
}

//_____________________________________________________________________________
// i=EChannel (kUp=0, kDown=1, kExtra=2). Returns Npe from m_npe_high.
// KVC: i=kExtra -> kSUM(ch=4). BAC: i=kExtra -> Npe(j).
Double_t
CherenkovHit::GetDeltaEHighGain(Int_t i, Int_t j) const
{
  if(DetectorId() == DetIdBAC && i == (Int_t)HodoRawHit::kExtra)
    return Npe(j);
  Int_t ch;
  if(DetectorId() == DetIdKVC && i == (Int_t)HodoRawHit::kExtra)
    ch = HodoRawHit::EChannelKVC::kSUM;
  else if(i >= m_n_ch)
    return TMath::QuietNaN();
  else
    ch = i;
  if(ch < 0 || ch >= (Int_t)m_npe_high.size()) return TMath::QuietNaN();
  if(j < 0 || j >= (Int_t)m_npe_high.at(ch).size()) return TMath::QuietNaN();
  return m_npe_high.at(ch).at(j);
}

//_____________________________________________________________________________
// j: hit index. 1ch: that ch. 2ch(U/D): sum of U and D. KVC: sum kA+kB+kC+kD (ch0-3).
Double_t
CherenkovHit::Npe(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return m_npe_high.at(U).at(j);
    }
    if(DetectorId() == DetIdKVC){
      Double_t sum = 0;
      for(Int_t c = 0; c < 4; ++c){
        if(j < (Int_t)m_npe_high.at(c).size())
          sum += m_npe_high.at(c).at(j);
      }
      return sum;
    }
    return m_npe_high.at(U).at(j) + m_npe_high.at(D).at(j);
  } catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
// Online sum (hardware SUM). KVC: ch4. BAC: seg4 only, else NaN.
Double_t
CherenkovHit::NpeSum(Int_t j) const
{
  if(DetectorId() == DetIdKVC)
    return GetNpe(HodoRawHit::EChannelKVC::kSUM, j);
  if(DetectorId() == DetIdBAC && m_raw->SegmentId() == 4)
    return Npe(j);
  return TMath::QuietNaN();
}

//_____________________________________________________________________________
Double_t
CherenkovHit::GetNpe(Int_t ch, Int_t j) const
{
  try {
    return m_npe_high.at(ch).at(j);
  } catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
const std::vector<Double_t>&
CherenkovHit::GetArrayNpe(Int_t ch) const
{
  return m_npe_high.at(ch);
}
