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
    m_data(),
    m_wpos(qnan),
    m_angle(0.),
    m_cluster_size(0.),
    m_mwpc_flag(false),
    m_ofs_dt(0.)
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
    m_mwpc_flag(false),
    m_ofs_dt(0.)
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
    m_mwpc_flag(false),
    m_ofs_dt(0.),
    m_hitnum(0)
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
DCHit::SetTdcCFT(Int_t tdc)
{
  m_tdc.push_back(tdc);
  //m_belong_track.push_back(false);
}

//_____________________________________________________________________________
void
DCHit::SetDummyData()
{
  data_t a_data = { 0., 0., qnan, false, true };
  m_data.push_back(a_data);
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
DCHit::Calculate()
{
  return CalcDCObservables();
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

  Bool_t status = true;

  DVec_t leading, trailing;
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
    ctime += m_ofs_dt;
    Double_t dt = TMath::QuietNaN();
    Double_t dl = TMath::QuietNaN();
    gDrift.CalcDrift(m_layer, m_wire, ctime, dt, dl);
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
      status = false;
      break;
    }
    data_t a_data;
    a_data.drift_time = dt;
    a_data.drift_length = dl;
    a_data.tot = tot;
    a_data.belong_to_track = false;
    a_data.dl_is_good = dl_is_good;
    m_data.push_back(a_data);
    if(SelectTDC1st && dl_is_good){
      m_data.clear();
      m_data.push_back(a_data);
      break;
    }
  }

  return status;
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
    data_t a_data = { static_cast<Double_t>(tdc), 0., qnan, false, true };
    m_data.push_back(a_data);
  }
  return true;
}

#if 0
//_____________________________________________________________________________
Bool_t
DCHit::CalcMWPCObservables()
{
  if(!gGeom.IsReady() || !gTdc.IsReady())
    return false;

  m_angle = gGeom.GetTiltAngle(m_layer);
  m_wpos  = gGeom.CalcWirePosition(m_layer, m_wire);
  m_z     = gGeom.GetLocalZ(m_layer);

  Bool_t status = true;
  Int_t  nh_tdc      = m_tdc.size();
  Int_t  nh_trailing = m_trailing.size();

  IVec_t leading_cont, trailing_cont;

  // Prepare
  {
    for(Int_t m = 0; m < nh_tdc; ++m) {
      leading_cont.push_back(m_tdc.at(m));
    }
    for(Int_t m = 0; m < nh_trailing; ++m) {
      trailing_cont.push_back(m_trailing.at(m));
    }

    std::sort(leading_cont.begin(),  leading_cont.end(),  std::greater<Int_t>());
    std::sort(trailing_cont.begin(), trailing_cont.end(), std::greater<Int_t>());

    Int_t i_t = 0;
    for(Int_t i = 0; i<nh_tdc; ++i){
      data_t a_data = {0., 0., qnan, qnan, -1, false, false};

      Int_t leading  = leading_cont.at(i);
      while(i_t < nh_trailing){
	Int_t trailing = trailing_cont.at(i_t);

	if(leading > trailing){
	  a_data.index_t = i_t;
	  m_data.push_back(a_data);
	  break;
	}else{
	  ++i_t;
	}// Goto next trailing
      }

      if(i_t == nh_trailing){
	a_data.index_t = -1;
	m_data.push_back(a_data);
	continue;
      }// no more trailing data
    }// for(i)
  }

  // Delete duplication index_t
  for(Int_t i = 0; i<nh_tdc-1; ++i){
    if(true
       && m_data.at(i).index_t != -1
       && m_data.at(i).index_t == m_data.at(i+1).index_t)
    {
      m_data.at(i).index_t = -1;
    }
  }// for(i)


  for(Int_t i=0; i<nh_tdc; ++i) {
    Double_t ctime;
    if(!gTdc.GetTime(m_layer, m_wire, leading_cont.at(i), ctime)){
      return false;
    }

    Double_t dtime, dlength;
    if(!gDrift.CalcDrift(m_layer, m_wire, ctime, dtime, dlength)){
      status = false;
    }

    m_data.at(i).drift_time   = dtime;
    m_data.at(i).drift_length = dlength;

    if(m_data.at(i).index_t != -1){
      Double_t trailing_ctime;
      gTdc.GetTime(m_layer, m_wire, trailing_cont.at(i), trailing_ctime);
      m_data.at(i).trailing_time = trailing_ctime;
      m_data.at(i).tot           = ctime - trailing_ctime;
    }else{
      m_data.at(i).trailing_time = qnan;
      m_data.at(i).tot           = qnan;
    }

    if(m_data.at(i).drift_time > MinDLBc[m_layer-100] && m_data.at(i).drift_time < MaxDLBc[m_layer-100]){
      m_data.at(i).dl_range = true;
    }else{
      status = false;
    }
  }

  return status;
}

//_____________________________________________________________________________
Bool_t
DCHit::CalcCFTObservables()
{
  if(!gGeom.IsReady()) return false;

  //m_angle = gGeom.GetTiltAngle(m_layer);
  //m_z     = gGeom.GetLocalZ(m_layer);

  Bool_t status = true;

  std::size_t nh_tdc = m_tdc.size();
  for(std::size_t i=0; i<nh_tdc; i++){
    data_t a_data = {(Double_t)m_tdc[i], 0., qnan, qnan, //-1, false, true};
      -1, false, false};
    m_data.push_back(a_data);
  }


  return status;
}
#endif

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
DCHit::TotCut(Double_t min_tot, Bool_t adopt_nan)
{
  auto itr_new_end =
    std::remove_if(m_data.begin(), m_data.end(),
		   [min_tot, adopt_nan](data_t a_data)->Bool_t
                     {return(isnan(a_data.tot) && adopt_nan) ? false : !(a_data.tot > min_tot);}
      );
  m_data.erase(itr_new_end, m_data.end());
}

//_____________________________________________________________________________
void
DCHit::GateDriftTime(Double_t min, Double_t max, Bool_t select_1st)
{
  auto itr_new_end =
    std::remove_if(m_data.begin(), m_data.end(),
		   [min, max](data_t a_data)->Bool_t
                     {return !(min < a_data.drift_time && a_data.drift_time < max);}
      );
  m_data.erase(itr_new_end, m_data.end());
  if(0 == m_data.size()) return;

  if(select_1st) m_data.erase(m_data.begin()+1, m_data.end());
}

//_____________________________________________________________________________
void
DCHit::Print(const TString& arg, std::ostream& ost) const
{
  const Int_t w = 16;
  ost << FUNC_NAME << " " << arg << std::endl
      << std::setw(w) << std::left << "layer" << m_layer << std::endl
      << std::setw(w) << std::left << "wire"  << m_wire  << std::endl
      << std::setw(w) << std::left << "wpos"  << m_wpos  << std::endl
      << std::setw(w) << std::left << "angle" << m_angle << std::endl
      << std::setw(w) << std::left << "z"     << m_z     << std::endl;

  ost << std::setw(w) << std::left << "tdc" << m_tdc.size() << " : ";
  std::copy(m_tdc.begin(), m_tdc.end(),
            std::ostream_iterator<Int_t>(ost, " "));
  ost << std::endl;
  ost << std::setw(w) << std::left << "trailing" << m_trailing.size() << " : ";
  std::copy(m_trailing.begin(), m_trailing.end(),
            std::ostream_iterator<Int_t>(ost, " "));
  ost << std::endl;
  hddaq::cout << std::setw(w) << std::left << "data.size()"
              << m_data.size() << std::endl;
  for(const auto& d: m_data){
    hddaq::cout << "\tdt:" << d.drift_time
                << "\tdl:" << d.drift_length
                << "\ttot:" << d.tot
                << "\tbelong_to_track:" << d.belong_to_track
                << "\tdl_is_good:" << d.dl_is_good
                << std::endl;
  }

  if(m_mwpc_flag){
    ost << std::endl
	<< std::setw(w) << std::left << "clsize" << m_cluster_size << std::endl
	<< std::setw(w) << std::left << "mean wire" << m_mwpc_wire << std::endl
	<< std::setw(w) << std::left << "mean pos"  << m_mwpc_wpos << std::endl;
  }
}
