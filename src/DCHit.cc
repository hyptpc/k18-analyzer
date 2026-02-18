// -*- C++ -*-
#include "DCHit.hh"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iomanip>

#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <TString.h>

#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCLTrackHit.hh"
#include "DCParameters.hh"
#include "DCRawHit.hh"
#include "DCTdcCalibMan.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "PrintHelper.hh"
#include "RootHelper.hh"
#include "std_ostream.hh"
#include "UserParamMan.hh"

#include <spdlog/spdlog.h>

namespace
{
const auto& gDrift = DCDriftParamMan::GetInstance();
const auto& gGeom  = DCGeomMan::GetInstance();
const auto& gTdc   = DCTdcCalibMan::GetInstance();
const auto& gUser  = UserParamMan::GetInstance();
const Double_t qnan = TMath::QuietNaN();
const Bool_t SelectTDC1st = false;
}

//_____________________________________________________________________________
DCHit::DCHit(const DCRawHit* rhit)
  : m_raw_hit(rhit),
    m_plane(rhit->PlaneId()),
    m_layer(rhit->DCGeomLayerId()),
    m_wire(rhit->WireId()),
    m_tdc(),
    m_adc(),
    m_trailing(),
    m_drift_time(),
    m_drift_length(),
    m_tot(),
    m_belong_to_track(),
    m_is_good(),
    m_wpos(qnan),
    m_angle(0.),
    m_z(qnan),
    m_register_container()

{
  for(const auto& t: rhit->GetTdcArray())
    m_tdc.push_back(t);
  for(const auto& t: rhit->GetTrailingArray())
    m_trailing.push_back(t);
  std::sort(m_tdc.begin(), m_tdc.end(), std::greater<Int_t>());
  std::sort(m_trailing.begin(), m_trailing.end(), std::greater<Int_t>());
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCHit::DCHit(Int_t layer)
  : m_raw_hit(nullptr),
    m_plane(-1),
    m_layer(layer),
    m_wire(-1),
    m_tdc(),
    m_adc(),
    m_trailing(),
    m_drift_time(),
    m_drift_length(),
    m_tot(),
    m_belong_to_track(),
    m_is_good(),
    m_wpos(qnan),
    m_angle(0.),
    m_z(qnan),
    m_register_container()

{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCHit::DCHit(Int_t layer, Double_t wire)
  : m_raw_hit(nullptr),
    m_plane(-1),
    m_layer(layer),
    m_wire(wire),
    m_tdc(),
    m_adc(),
    m_trailing(),
    m_drift_time(),
    m_drift_length(),
    m_tot(),
    m_belong_to_track(),
    m_is_good(),
    m_wpos(qnan),
    m_angle(0.),
    m_z(qnan),
    m_register_container()

{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCHit::~DCHit()
{
  ClearRegisteredHits();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
DCHit::ClearDCData()
{
  m_drift_time.clear();
  m_drift_length.clear();
  m_tot.clear();
  m_belong_to_track.clear();
  m_is_good.clear();
}

//_____________________________________________________________________________
void
DCHit::EraseDCData(Int_t i)
{
  m_tdc.erase(m_tdc.begin() + i);
  m_trailing.erase(m_trailing.begin() + i);
  m_drift_time.erase(m_drift_time.begin() + i);
  m_drift_length.erase(m_drift_length.begin() + i);
  m_tot.erase(m_tot.begin() + i);
  m_belong_to_track.erase(m_belong_to_track.begin() + i);
  m_is_good.erase(m_is_good.begin() + i);
}

//_____________________________________________________________________________
void
DCHit::SetDCData(Double_t dt, Double_t dl, Double_t tot,
                 Bool_t belong_to_track, Bool_t is_good)
{
  m_drift_time.push_back(dt);
  m_drift_length.push_back(dl);
  m_tot.push_back(tot);
  m_belong_to_track.push_back(belong_to_track);
  m_is_good.push_back(is_good);
}

//_____________________________________________________________________________
void
DCHit::ClearRegisteredHits()
{
  del::ClearContainer(m_register_container);
}

//_____________________________________________________________________________
Bool_t
DCHit::CalcDCObservables()
{
  if(false
     || !gGeom.IsReady()
     || !gTdc.IsReady()
     || !gDrift.IsReady()){
    return false;
  }

  m_wpos  = gGeom.CalcWirePosition(m_layer, m_wire);
  m_angle = gGeom.GetTiltAngle(m_layer);
  m_z     = gGeom.GetLocalZ(m_layer);

  std::sort(m_tdc.begin(), m_tdc.end(), std::greater<Double_t>());
  std::sort(m_trailing.begin(), m_trailing.end(), std::greater<Double_t>());

  const auto detector_id = m_raw_hit->DetectorId();

  data_t leading, trailing;
  for(Int_t il=0, nl=m_tdc.size(); il<nl; ++il){
    Double_t l = m_tdc[il];
    Double_t l_next = (il+1) != nl ? m_tdc[il+1] : DBL_MIN;
    Double_t buf = DBL_MAX; // TMath::QuietNaN();
    for(const auto& t: m_trailing){
      if(l_next<t && t<l){
        buf = t;
        break;
      }
    }
    leading.push_back(l);
    trailing.push_back(buf);
    Double_t ctime = qnan;
    gTdc.GetTime(detector_id, m_plane, m_wire, l, ctime);
    Double_t dt = ctime;
    Double_t dl = qnan;
    gDrift.CalcDrift(m_raw_hit->DetectorName(),
                     m_plane, m_wire, ctime, dt, dl);
    Double_t ctime_trailing = qnan;
    gTdc.GetTime(detector_id, m_plane, m_wire, buf, ctime_trailing);
    // Double_t tot = ctime - ctime_trailing;
    Double_t tot = l - buf;
    Bool_t dl_is_good = true;
    if(!SelectTDC1st){
      SetDCData(dt, dl, tot, false, dl_is_good);
    }else if(dl_is_good){
      SetDCData(dt, dl, tot, false, dl_is_good);
      break;
    }
  }
  m_tdc = leading;
  m_trailing = trailing;
  return true;
}

//_____________________________________________________________________________
Int_t
DCHit::GetTdc1st() const
{
  if(m_tdc.empty())
    return TMath::QuietNaN();
  else
    return m_tdc.front();
}

//_____________________________________________________________________________
Double_t
DCHit::GetResolution() const
{
  return gGeom.GetResolution(m_layer);
}

//_____________________________________________________________________________
void
DCHit::DriftTimeCut(Double_t min, Double_t max, Bool_t select_1st)
{
  for(Int_t i=GetEntries()-1; i>=0; --i){
    if(m_drift_time[i] < min || max < m_drift_time[i]){
      EraseDCData(i);
    }
  }
  if(select_1st){
    for(Int_t i=GetEntries()-1; i>0; --i){
      EraseDCData(i);
    }
  }
}

//_____________________________________________________________________________
void
DCHit::TotCut(Double_t min, Bool_t keep_nan)
{
  for(Int_t i=GetEntries()-1; i>=0; --i){
    if(keep_nan && TMath::IsNaN(m_tot[i]))
      continue;
    if(m_tot[i] < min){
      EraseDCData(i);
    }
  }
}

//_____________________________________________________________________________
void
DCHit::Print(Option_t* arg) const
{
  PrintHelper helper(3, std::ios::fixed);

  const Int_t w = 16;
  TString dname = m_raw_hit ? m_raw_hit->DetectorName() : "Unknown";
  Double_t z = gGeom.GetLocalZ(m_layer);

  hddaq::cout << FUNC_NAME << " " << arg << std::endl
              << std::setw(w) << std::left << "detector" << dname << std::endl
              << std::setw(w) << std::left << "plane"    << m_plane << std::endl
              << std::setw(w) << std::left << "layer"    << m_layer << std::endl
              << std::setw(w) << std::left << "wire"     << m_wire  << std::endl
              << std::setw(w) << std::left << "wpos"     << m_wpos  << std::endl
              << std::setw(w) << std::left << "angle"    << m_angle << std::endl
              << std::setw(w) << std::left << "z"        << z       << std::endl;

  hddaq::cout << std::setw(w) << std::left << "tdc" << (Int_t)m_tdc.size() << " : ";
  for(const auto& t: m_tdc) hddaq::cout << t << " ";
  hddaq::cout << std::endl;

  hddaq::cout << std::setw(w) << std::left << "trailing" << (Int_t)m_trailing.size() << " : ";
  for(const auto& t: m_trailing) hddaq::cout << t << " ";
  hddaq::cout << std::endl;

  for(const auto& data_map: std::map<TString, data_t>
        {{"dt", m_drift_time},
         {"dl", m_drift_length},
         {"tot", m_tot},
        }){
    const auto& cont = data_map.second;
    hddaq::cout << std::setw(w) << std::left << data_map.first << (Int_t)cont.size()
                << " : ";
    for(const auto& val: cont) hddaq::cout << val << " ";
    hddaq::cout << std::endl;
  }
  for(const auto& data_map: std::map<TString, flag_t>
        {{"belong_to_track", m_belong_to_track},
         {"is_good", m_is_good},
        }){
    const auto& cont = data_map.second;
    hddaq::cout << std::setw(w) << std::left << data_map.first << (Int_t)cont.size()
                << " : ";
    for(const auto& val: cont) hddaq::cout << val << " ";
    hddaq::cout << std::endl;
  }
}
