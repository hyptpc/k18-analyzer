// -*- C++ -*-

#include "HodoRawHit.hh"

#include <iterator>

#include <TMath.h>

#include <UnpackerManager.hh>

#include "DebugCounter.hh"
#include "FuncName.hh"

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();
}

//_____________________________________________________________________________
HodoRawHit::HodoRawHit(const TString& detector_name,
                       Int_t plane_id, Int_t segment_id)
  : m_detector_name(detector_name),
    m_detector_id(gUnpacker.get_device_id(detector_name)),
    m_plane_id(plane_id),
    m_segment_id(segment_id),
    m_adc_high(kNChannel),
    m_adc_low(kNChannel),
    m_tdc_leading(kNChannel),
    m_tdc_trailing(kNChannel),
    m_tdc_is_overflow(false)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
HodoRawHit::~HodoRawHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
HodoRawHit::Clear()
{
  m_adc_high.clear();
  m_adc_low.clear();
  m_tdc_leading.clear();
  m_tdc_trailing.clear();
  m_adc_high.resize(kNChannel);
  m_adc_low.resize(kNChannel);
  m_tdc_leading.resize(kNChannel);
  m_tdc_trailing.resize(kNChannel);
}

//_____________________________________________________________________________
Double_t
HodoRawHit::GetAdcHigh(Int_t i, Int_t j) const
{
  try{
    return m_adc_high.at(i).at(j);
  }catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoRawHit::GetAdcLow(Int_t i, Int_t j) const
{
  try{
    return m_adc_low.at(i).at(j);
  }catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoRawHit::GetTdcLeading(Int_t i, Int_t j) const
{
  try{
    return m_tdc_leading.at(i).at(j);
  }catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoRawHit::GetTdcTrailing(Int_t i, Int_t j) const
{
  try{
    return m_tdc_trailing.at(i).at(j);
  }catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
void
HodoRawHit::Print(const TString& arg) const
{
  hddaq::cerr << FUNC_NAME << " " << arg << std::endl
	      << "detector_name = " << m_detector_name << std::endl
	      << "detector_id   = " << m_detector_id   << std::endl
	      << "plane_id      = " << m_plane_id      << std::endl
	      << "segment_id    = " << m_segment_id    << std::endl;
  for(const auto data_map: std::map<TString, data_t>
        {{"adc-hi", m_adc_high}, {"adc-lo", m_adc_low},
         {"tdc-l ", m_tdc_leading}, {"tdc-t ", m_tdc_trailing}}){
    for(const auto& cont: data_map.second){
      hddaq::cout << " " << data_map.first << ":" << cont.size()
                  << " ";
      std::copy(cont.begin(), cont.end(),
                std::ostream_iterator<UInt_t>(hddaq::cout," "));
    }
    hddaq::cout << std::endl;
  }
}
