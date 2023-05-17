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
const auto& gConf = hddaq::unpacker::GConfig::get_instance();
const auto& gHodo = HodoParamMan::GetInstance();
const auto& gPHC = HodoPHCMan::GetInstance();
}

//_____________________________________________________________________________
HodoHit::HodoHit(HodoRawHit *rhit, Double_t max_time_diff)
  : m_raw(rhit),
    m_is_calculated(false),
    m_max_time_diff(max_time_diff),
    m_n_ch(gConf.get_digit_info().get_n_ch(rhit->DetectorId())),
    m_de_high(HodoRawHit::kNChannel),
    m_de_low(HodoRawHit::kNChannel),
    m_time_leading(HodoRawHit::kNChannel),
    m_time_trailing(HodoRawHit::kNChannel),
    m_ctime_leading(HodoRawHit::kNChannel),
    m_ctime_trailing(HodoRawHit::kNChannel),
    m_flag_join()
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

  data_t leading(HodoRawHit::kNChannel);
  data_t trailing(HodoRawHit::kNChannel);
  data_t cleading(HodoRawHit::kNChannel);
  data_t ctrailing(HodoRawHit::kNChannel);

  for(Int_t ch=0; ch<HodoRawHit::kNChannel; ++ch){
    // adc
    for(const auto& adc: m_raw->GetArrayAdcHigh(ch)){
      Double_t de = TMath::QuietNaN();
      if(gHodo.GetDe(id, plane, seg, ch, adc, de)){
        m_de_high.at(ch).push_back(de);
      }
    }
    for(const auto& adc: m_raw->GetArrayAdcLow(ch)){
      Double_t de = TMath::QuietNaN();
      if(gHodo.GetDe(id, plane, seg, ch, adc, de)){
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

  auto& lu_cont = leading.at(HodoRawHit::kUp);
  auto& ld_cont = leading.at(HodoRawHit::kDown);
  auto& clu_cont = cleading.at(HodoRawHit::kUp);
  auto& cld_cont = cleading.at(HodoRawHit::kDown);

  auto& tu_cont = trailing.at(HodoRawHit::kUp);
  auto& td_cont = trailing.at(HodoRawHit::kDown);
  auto& ctu_cont = ctrailing.at(HodoRawHit::kUp);
  auto& ctd_cont = ctrailing.at(HodoRawHit::kDown);

  for(Int_t ju=0; ju<lu_cont.size(); ++ju){
    Double_t lu = lu_cont.at(ju);
    Double_t clu = clu_cont.at(ju);
    // one-side readout
    if(m_n_ch == 1){
      m_time_leading.at(HodoRawHit::kUp).push_back(lu);
      m_ctime_leading.at(HodoRawHit::kUp).push_back(clu);
    }
    // two-side readout
    else{
      for(Int_t jd=0; jd<ld_cont.size(); ++jd){
        Double_t ld = ld_cont.at(jd);
        if(TMath::Abs(lu - ld) < m_max_time_diff){
          m_time_leading.at(HodoRawHit::kUp).push_back(lu);
          m_time_leading.at(HodoRawHit::kDown).push_back(ld);
          Double_t cld = cld_cont.at(jd);
          m_ctime_leading.at(HodoRawHit::kUp).push_back(clu);
          m_ctime_leading.at(HodoRawHit::kDown).push_back(cld);
          break;
        }else{
          ;
        }
      }
    }
  }

  m_time_leading.at(HodoRawHit::kExtra) = leading.at(HodoRawHit::kExtra);
  m_ctime_leading.at(HodoRawHit::kExtra) = cleading.at(HodoRawHit::kExtra);

  m_is_calculated = true;
  return (m_ctime_leading.at(HodoRawHit::kUp).size() > 0);
}

//_____________________________________________________________________________
Double_t
HodoHit::DeltaE(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return m_de_high.at(HodoRawHit::kUp).at(j);
    }else{
      return TMath::Sqrt(
        TMath::Abs(m_de_high.at(HodoRawHit::kUp).at(j) *
                   m_de_high.at(HodoRawHit::kDown).at(j)));
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
      return m_time_leading.at(HodoRawHit::kUp).at(j);
    }else{
      return 0.5*(m_time_leading.at(HodoRawHit::kUp).at(j) +
                  m_time_leading.at(HodoRawHit::kDown).at(j));
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
      return m_ctime_leading.at(HodoRawHit::kUp).at(j);
    }else{
      return 0.5*(m_ctime_leading.at(HodoRawHit::kUp).at(j) +
                  m_ctime_leading.at(HodoRawHit::kDown).at(j));
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
      return (m_time_leading.at(HodoRawHit::kDown).at(j) -
              m_time_leading.at(HodoRawHit::kUp).at(j));
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
      return (m_ctime_leading.at(HodoRawHit::kDown).at(j) -
              m_ctime_leading.at(HodoRawHit::kUp).at(j));
    }
  }catch(const std::out_of_range& e){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Bool_t
HodoHit::JoinedAllMhit()
{
  Bool_t ret = true;
  for(Int_t i=0, n=m_flag_join.size(); i<n; ++i){
    ret = ret & m_flag_join[i];
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
