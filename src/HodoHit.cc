// -*- C++ -*-

#include "HodoHit.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "PrintHelper.hh"
#include "RawData.hh"

#include <std_ostream.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

namespace
{
const auto U = HodoRawHit::kUp;
const auto D = HodoRawHit::kDown;
const auto E = HodoRawHit::kExtra;
const auto& gUnpackerConf = hddaq::unpacker::GConfig::get_instance();
const auto& gHodo = HodoParamMan::GetInstance();
const auto& gPHC = HodoPHCMan::GetInstance();
}

//_____________________________________________________________________________
HodoHit::HodoHit(const HodoRawHit *rhit, Double_t max_time_diff)
  : m_raw(rhit),
    m_is_calculated(false),
    m_max_time_diff(max_time_diff),
    m_n_ch(gUnpackerConf.get_digit_info().get_n_ch(rhit->DetectorId())),
    m_time_offset(),
    m_de_high(m_n_ch),
    m_de_low(m_n_ch),
    m_time_leading(m_n_ch),
    m_time_trailing(m_n_ch),
    m_ctime_leading(m_n_ch),
    m_ctime_trailing(m_n_ch),
    m_is_clustered()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
HodoHit::~HodoHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
HodoHit::Calculate()
{
  if(m_is_calculated){
    hddaq::cerr << FUNC_NAME << " already calculated" << std::endl;
    return false;
  }

  if(!gHodo.IsReady()){
    hddaq::cerr << FUNC_NAME << " HodoParamMan must be initialized" << std::endl;
    return false;
  }

  if(!gPHC.IsReady()){
    hddaq::cout << FUNC_NAME << " HodoPHCMan must be initialized" << std::endl;
    return false;
  }

  Int_t id    = m_raw->DetectorId();
  Int_t plane = m_raw->PlaneId();
  Int_t seg   = m_raw->SegmentId();

  data_t leading(m_n_ch);
  data_t trailing(m_n_ch);
  data_t cleading(m_n_ch);
  data_t ctrailing(m_n_ch);

  for(Int_t ch=0; ch<m_n_ch; ++ch){
    // adc
    for(const auto& adc: m_raw->GetArrayAdcHigh(ch)){
      Double_t de = TMath::QuietNaN();
      if(gHodo.GetDeHighGain(id, plane, seg, ch, adc, de)){
        m_de_high.at(ch).push_back(de);
      }
    }
    for(const auto& adc: m_raw->GetArrayAdcLow(ch)){
      Double_t de = TMath::QuietNaN();
      if(gHodo.GetDeLowGain(id, plane, seg, ch, adc, de)){
        m_de_low.at(ch).push_back(de);
      }
    }
    // tdc
    Double_t de =
      m_de_high.at(ch).size() > 0 ?
      m_de_high.at(ch).at(0) : TMath::QuietNaN();
    for(const auto& tdc: m_raw->GetArrayTdcLeading(ch)){
      Double_t time = TMath::QuietNaN();
      if(gHodo.GetTime(id, plane, seg, ch, tdc, time)){
        leading.at(ch).push_back(time);
        Double_t ctime = TMath::QuietNaN();
        gPHC.DoCorrection(id, plane, seg, ch, time, de, ctime);
        cleading.at(ch).push_back(ctime);
      }
    }
    for(const auto& tdc: m_raw->GetArrayTdcTrailing(ch)){
      Double_t time = TMath::QuietNaN();
      if(gHodo.GetTime(id, plane, seg, ch, tdc, time)){
        trailing.at(ch).push_back(time);
        Double_t ctime = TMath::QuietNaN();
        gPHC.DoCorrection(id, plane, seg, ch, time, de, ctime);
        ctrailing.at(ch).push_back(ctime);
      }
    }
    std::sort(leading.at(ch).begin(), leading.at(ch).end());
    std::sort(trailing.at(ch).begin(), trailing.at(ch).end());
    std::sort(cleading.at(ch).begin(), cleading.at(ch).end());
    std::sort(ctrailing.at(ch).begin(), ctrailing.at(ch).end());
  }

  // Double_t offset_vtof = 0.;
  // if(m_name.EqualTo("TOF")){
  //   gHodo.GetTime(id, plane, seg, 2, 0., offset_vtof);
  // }

  // one-side readout
  if(m_n_ch == 1){
    m_time_leading.at(U) = leading.at(U);
    m_ctime_leading.at(U) = cleading.at(U);
  }
  // two-side readout
  else{
    for(Int_t ju=0; ju<leading.at(U).size(); ++ju){
      Double_t lu = leading.at(U).at(ju);
      Double_t clu = cleading.at(U).at(ju);
      for(Int_t jd=0; jd<leading.at(D).size(); ++jd){
        Double_t ld = leading.at(D).at(jd);
        if(TMath::Abs(lu - ld) < m_max_time_diff){
          Double_t cld = cleading.at(D).at(jd);
          m_time_leading.at(U).push_back(lu);
          m_time_leading.at(D).push_back(ld);
          m_ctime_leading.at(U).push_back(clu);
          m_ctime_leading.at(D).push_back(cld);
          m_is_clustered.push_back(false);
          break;
        }else{
          ;
        }
      }
    }
  }

  // extra channel remains
  if(m_n_ch == HodoRawHit::kNChannel){
    m_time_leading.at(E) = leading.at(E);
    m_ctime_leading.at(E) = cleading.at(E);
  }

  for(Int_t ch=0; ch<m_n_ch; ++ch){
    m_time_trailing.at(ch) = trailing.at(ch);
    m_ctime_trailing.at(ch) = trailing.at(ch);
  }

  /*
    HodoHit considers only the coinsidence of Up/Down LEADINGs
    and does not care about the presence/absence of TRAILINGs or the counts.
  */

  m_is_calculated = true;
  return (m_ctime_leading.at(U).size() > 0);
}

//_____________________________________________________________________________
Double_t
HodoHit::DeltaEHighGain(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return m_de_high.at(U).at(j);
    }else{
      return TMath::Sqrt(
        TMath::Abs(m_de_high.at(U).at(j) *
                   m_de_high.at(D).at(j)));
    }
  }catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoHit::DeltaELowGain(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return m_de_low.at(U).at(j);
    }else{
      return TMath::Sqrt(
        TMath::Abs(m_de_low.at(U).at(j) *
                   m_de_low.at(D).at(j)));
    }
  }catch(const std::out_of_range& e){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoHit::MeanTime(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return m_time_leading.at(U).at(j);
    }else{
      return 0.5*(m_time_leading.at(U).at(j) +
                  m_time_leading.at(D).at(j));
    }
  }catch(const std::out_of_range& e){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoHit::CMeanTime(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return m_ctime_leading.at(U).at(j);
    }else{
      return 0.5*(m_ctime_leading.at(U).at(j) +
                  m_ctime_leading.at(D).at(j));
    }
  }catch(const std::out_of_range& e){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoHit::TimeDiff(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return TMath::QuietNaN();
    }else{
      return (m_time_leading.at(D).at(j) -
              m_time_leading.at(U).at(j));
    }
  }catch(const std::out_of_range& e){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoHit::CTimeDiff(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return 0.;
    }else{
      return (m_ctime_leading.at(D).at(j) -
              m_ctime_leading.at(U).at(j));
    }
  }catch(const std::out_of_range& e){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Bool_t
HodoHit::IsClusteredAll()
{
  Bool_t ret = true;
  for(const auto& f: m_is_clustered){
    ret &= f;
  }
  return ret;
}

//_____________________________________________________________________________
void
HodoHit::Print(Option_t* arg) const
{
  PrintHelper helper(3, std::ios::fixed);
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
	      << "detector_name = " << m_raw->DetectorName() << std::endl
	      << "detector_id   = " << m_raw->DetectorId() << std::endl
	      << "plane_id      = " << m_raw->PlaneId()    << std::endl
	      << "segment_id    = " << m_raw->SegmentId()  << std::endl
              << "n_ch          = " << m_n_ch              << std::endl
              << "de            = " << DeltaE() << std::endl
              << "time offset   = " << m_time_offset << std::endl
              << "mt/cmt        = " << MeanTime()
              << " / " << CMeanTime() << std::endl
              << "tdiff/ctdiff  = " << TimeDiff()
              << " / " << CTimeDiff() << std::endl;
  for(const auto data_map: std::map<TString, data_t>
        {{"de-hi  ", m_de_high},      {"de-lo  ", m_de_low},
         {"time-l ", m_time_leading}, {"time-t ", m_time_trailing},
         {"ctime-l", m_ctime_leading}, {"ctime-t", m_ctime_trailing}
        }){
    for(const auto& cont: data_map.second){
      hddaq::cout << " " << data_map.first << ":" << cont.size()
                  << " ";
      std::copy(cont.begin(), cont.end(),
                std::ostream_iterator<Double_t>(hddaq::cout," "));
    }
    hddaq::cout << std::endl;
  }
}
