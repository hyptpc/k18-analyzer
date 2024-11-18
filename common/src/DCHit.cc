// -*- C++ -*-

#include "DCHit.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCParameters.hh"
#include "DCRawHit.hh"
#include "DCTdcCalibMan.hh"
#include "DCLTrackHit.hh"
#include "UserParamMan.hh"
#include "DebugCounter.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "RootHelper.hh"

namespace
{
const Double_t qnan = TMath::QuietNaN();
const auto& gGeom  = DCGeomMan::GetInstance();
const auto& gTdc   = DCTdcCalibMan::GetInstance();
const auto& gDrift = DCDriftParamMan::GetInstance();
const auto& gUser  = UserParamMan::GetInstance();
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
    m_cluster_size(0.),
    m_mwpc_flag(false)
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
  : m_layer(layer),
    m_wire(-1),
    m_wpos(qnan),
    m_angle(0.),
    m_cluster_size(0.),
    m_mwpc_flag(false)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCHit::DCHit(Int_t layer, Double_t wire)
  : m_layer(layer),
    m_wire(wire),
    m_wpos(qnan),
    m_angle(0.),
    m_cluster_size(0.),
    m_mwpc_flag(false)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCHit::DCHit(Int_t plane, Int_t layer, Int_t wire) // for Geant4
  : m_plane(plane),
    m_layer(layer),
    m_wire(wire),
    m_lpos(),
    m_de(),
    m_tot(),
    m_belong_to_track(),
    m_is_good(),
    m_wpos(qnan),
    m_angle(0.),
    m_cluster_size(0.),
    m_mwpc_flag(false)
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
DCHit::SetDCDataGeant4(TVector3 lpos, Double_t de)
{
  m_lpos.push_back(lpos);
  m_de.push_back(de);
}

//_____________________________________________________________________________
void
DCHit::ClearRegisteredHits()
{
  Int_t n = m_register_container.size();
  for(Int_t i=0; i<n; ++i){
    delete m_register_container[i];
  }
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

  std::sort(m_tdc.begin(), m_tdc.end(), std::greater<Int_t>());
  std::sort(m_trailing.begin(), m_trailing.end(), std::greater<Int_t>());

  data_t leading, trailing;
  for(Int_t il=0, nl=m_tdc.size(); il<nl; ++il){
    Double_t l = m_tdc[il];
    Double_t l_next = (il+1) != nl ? m_tdc[il+1] : DBL_MIN;
    Double_t buf = TMath::QuietNaN();
    for(const auto& t: m_trailing){
      if(l_next<t && t<l){
        buf = t;
        break;
      }
    }
    leading.push_back(l);
    trailing.push_back(buf);
    Double_t ctime = TMath::QuietNaN();
    gTdc.GetTime(m_layer, m_wire, l, ctime);
    Double_t dt = TMath::QuietNaN();
    Double_t dl = TMath::QuietNaN();
    // gDrift.CalcDrift(m_layer, m_wire, ctime, dt, dl);
    Double_t ctime_trailing = TMath::QuietNaN();
    gTdc.GetTime(m_layer, m_wire, buf, ctime_trailing);
    // Double_t tot = ctime - ctime_trailing;
    Double_t tot = l - buf;
    Bool_t dl_is_good = false;
    switch(m_layer){
      // BC3,4
    case 113: case 114: case 115: case 116: case 117: case 118:
    case 119: case 120: case 121: case 122: case 123: case 124:
      if(MinDLBc[m_layer-100] < dl && dl < MaxDLBc[m_layer-100]){
	dl_is_good = true;
      }
      break;
      // SDC1,2,3,4,5
    case 1: case 2: case 3: case 4: case 5: case 6:
    case 7: case 8: case 9: case 10:
    case 31: case 32: case 33: case 34:
    case 35: case 36: case 37: case 38:
    case 39: case 40: case 41: case 42:
      if(MinDLSdc[m_layer] < dl && dl < MaxDLSdc[m_layer]){
      	dl_is_good = true;
      }
      break;
    default:
      hddaq::cout << FUNC_NAME << " "
		  << "invalid layer id : " << m_layer << std::endl;
      return false;
    }

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
Bool_t
DCHit::CalcDCObservablesGeant4()
{

  if(false
     || !gGeom.IsReady()){
    return false;
  }

  m_wpos  = gGeom.CalcWirePosition(m_layer, m_wire);
  m_angle = gGeom.GetTiltAngle(m_layer);
  m_z     = gGeom.GetLocalZ(m_layer);

  // gRandom->SetSeed(TDatime().Convert());
  // Double_t res = gUser.GetParameter(Form("ResolutionLayer%d", m_layer));
  for(Int_t i=0, n=m_lpos.size(); i<n; ++i){
    Double_t dt = TMath::QuietNaN();
    Double_t a  = m_angle*TMath::DegToRad();
    Double_t s  = m_lpos[i].x()*TMath::Cos(a) + m_lpos[i].y()*TMath::Sin(a);
    Double_t dl = TMath::Abs(s-m_wpos);
    // dl = gRandom->Gaus(dl, res);
    Double_t tot = m_de[i];
    Bool_t dl_is_good = false;
    switch(m_layer){
      // BC3,4
    case 113: case 114: case 115: case 116: case 117: case 118:
    case 119: case 120: case 121: case 122: case 123: case 124:
      if(MinDLBc[m_layer-100] < dl && dl < MaxDLBc[m_layer-100]){
	dl_is_good = true;
      }
      break;
      // SDC1,2,3,4,5
    case 1: case 2: case 3: case 4: case 5: case 6:
    case 7: case 8: case 9: case 10:
    case 31: case 32: case 33: case 34:
    case 35: case 36: case 37: case 38:
    case 39: case 40: case 41: case 42:
      if(MinDLSdc[m_layer] < dl && dl < MaxDLSdc[m_layer]){
	dl_is_good = true;
      }
      break;
    default:
      hddaq::cout << FUNC_NAME << " "
		  << "invalid layer id : " << m_layer << std::endl;
      return false;
    }

    if(!SelectTDC1st){
      SetDCData(dt, dl, tot, false, dl_is_good);
    }else if(dl_is_good){
      SetDCData(dt, dl, tot, false, dl_is_good);
      break;
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
DCHit::CalcFiberObservables()
{
  if(!gGeom.IsReady())
    return false;
  m_angle = gGeom.GetTiltAngle(m_layer);
  m_z     = gGeom.GetLocalZ(m_layer);
  for(const auto& tdc: m_tdc){
    m_drift_time.push_back(tdc);
    m_drift_length.push_back(0.);
    m_tot.push_back(qnan);
    m_belong_to_track.push_back(false);
    m_is_good.push_back(true);
  }
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
  const Int_t w = 16;
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
              << std::setw(w) << std::left << "detector"
              << m_raw_hit->DetectorName() << std::endl
              << std::setw(w) << std::left << "plane" << m_plane << std::endl
              << std::setw(w) << std::left << "layer" << m_layer << std::endl
              << std::setw(w) << std::left << "wire"  << m_wire  << std::endl
              << std::setw(w) << std::left << "wpos"  << m_wpos  << std::endl
              << std::setw(w) << std::left << "angle" << m_angle << std::endl
              << std::setw(w) << std::left << "z"     << m_z     << std::endl;

  hddaq::cout << std::setw(w) << std::left << "tdc" << m_tdc.size() << " : ";
  std::copy(m_tdc.begin(), m_tdc.end(),
            std::ostream_iterator<Int_t>(hddaq::cout, " "));
  hddaq::cout << std::endl;
  hddaq::cout << std::setw(w) << std::left << "trailing" << m_trailing.size() << " : ";
  std::copy(m_trailing.begin(), m_trailing.end(),
            std::ostream_iterator<Int_t>(hddaq::cout, " "));
  hddaq::cout << std::endl;
  for(const auto& data_map: std::map<TString, data_t>
        {{"dt", m_drift_time},
         {"dl", m_drift_length},
         {"tot", m_tot},
        }){
    const auto& cont = data_map.second;
    hddaq::cout << std::setw(w) << std::left << data_map.first << cont.size()
                << " : ";
    std::copy(cont.begin(), cont.end(),
              std::ostream_iterator<Double_t>(hddaq::cout, " "));
    hddaq::cout << std::endl;
  }
  for(const auto& data_map: std::map<TString, flag_t>
        {{"belong_to_track", m_belong_to_track},
         {"is_good", m_is_good},
        }){
    const auto& cont = data_map.second;
    hddaq::cout << std::setw(w) << std::left << data_map.first << cont.size()
                << " : ";
    std::copy(cont.begin(), cont.end(),
              std::ostream_iterator<Bool_t>(hddaq::cout, " "));
    hddaq::cout << std::endl;
  }

  if(m_mwpc_flag){
    hddaq::cout << std::endl
                << std::setw(w) << std::left << "clsize" << m_cluster_size << std::endl
                << std::setw(w) << std::left << "mean wire" << m_mwpc_wire << std::endl
                << std::setw(w) << std::left << "mean pos"  << m_mwpc_wpos << std::endl;
  }
}
