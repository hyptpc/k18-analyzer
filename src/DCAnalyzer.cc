// -*- C++ -*-
#include "DCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <TString.h>

#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "DCTrackSearch.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

#include <spdlog/spdlog.h>
#include <UnpackerConfig.hh>
#include <UnpackerXMLReadDigit.hh>

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

// Tracking routine selection __________________________________________________
/* BcInTracking */
// #define UseBcIn    0 // not supported
/* BcOutTracking */
#define BcOut_XUV  0 // XUV Tracking (slow but accerate)
#define BcOut_Pair 1 // Pair plane Tracking (fast but bad for large angle track)

namespace
{
const auto& gGeom   = DCGeomMan::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
DCAnalyzer::DCAnalyzer(const RawData& raw_data)
  : m_raw_data(&raw_data),
    m_dc_hit_collection(),
    m_TempBcInHC(NumOfLayersBcIn+1),

    m_BcInHC(NumOfLayersBcIn+1),
    m_BcOutHC(),
    m_BcInTC(),
    m_BcOutTC()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCAnalyzer::DCAnalyzer()
  : m_raw_data(nullptr),
    m_dc_hit_collection(),
    m_TempBcInHC(NumOfLayersBcIn+1),

    m_BcInHC(NumOfLayersBcIn+1),
    m_BcOutHC(),
    m_BcInTC(),
    m_BcOutTC()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCAnalyzer::~DCAnalyzer()
{
  for(auto& elem: m_dc_hit_collection)
    del::ClearContainer(elem.second);

  ClearTracksBcIn();
  ClearTracksBcOut();
  ClearDCHits();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeBcInHits()
{
  static const auto& digit_info =
    hddaq::unpacker::GConfig::get_instance().get_digit_info();
  m_BcInHC.clear();
  Int_t plane_offset = 0;
  for(const auto& name: DCNameList.at("BcIn")){
    Int_t id = digit_info.get_device_id(name.Data());
    Int_t n_plane = digit_info.get_n_plane(id);
    m_BcInHC.resize(n_plane + m_BcInHC.size());
    DecodeHits(name);
    for(const auto& hit: m_dc_hit_collection.at(name)){
      m_BcInHC[hit->PlaneId() + plane_offset].push_back(hit);
    }
    plane_offset += n_plane;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeBcOutHits()
{
  static const auto& digit_info =
    hddaq::unpacker::GConfig::get_instance().get_digit_info();
  m_BcOutHC.clear();
  Int_t plane_offset = 0;
  for(const auto& name: DCNameList.at("BcOut")){
    Int_t id = digit_info.get_device_id(name.Data());
    Int_t n_plane = digit_info.get_n_plane(id);
    m_BcOutHC.resize(n_plane + m_BcOutHC.size());
    DecodeHits(name);
    for(const auto& hit: m_dc_hit_collection.at(name)){
      m_BcOutHC[hit->PlaneId() + plane_offset].push_back(hit);
    }
    plane_offset += n_plane;
  }
  return true;
}

//_____________________________________________________________________________
void
DCAnalyzer::DecodeHits(const TString& name)
{
  if(m_dc_hit_collection.find(name) != m_dc_hit_collection.end()){
    spdlog::warn("{} {} is already decoded", FUNC_NAME.Data(), name.Data());
    return;
  }
  auto& HitCont = m_dc_hit_collection[name];
  del::ClearContainer(HitCont);
  for(const auto& rhit: m_raw_data->GetDCRawHitContainer(name)){
    auto hit = new DCHit(rhit);
    if(hit && hit->CalcDCObservables()){
      HitCont.push_back(hit);
    }else{
      delete hit;
    }
  }
  std::sort(HitCont.begin(), HitCont.end(), DCHit::Compare);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeRawHits()
{
  ClearDCHits();
  DecodeBcInHits();
  DecodeBcOutHits();
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchBcIn(Int_t T0Seg)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerBcIn");

  Int_t ntrack = track::LocalTrackSearch(m_BcInHC, PPInfoBcIn, NPPInfoBcIn,
                                         m_BcInTC, MinLayer, T0Seg);

  return ntrack == -1 ? false : true;
}


//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchBcOut(Int_t T0Seg)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerBcOut");

#if BcOut_Pair //Pair Plane Tracking Routine for BcOut
  Int_t ntrack = track::LocalTrackSearch(m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
                                         m_BcOutTC, MinLayer, T0Seg);
  return ntrack == -1 ? false : true;
#endif

#if BcOut_XUV  //XUV Tracking Routine for BcOut
  Int_t ntrack = track::LocalTrackSearchVUX(m_BcOutHC, PPInfoBcOut, NPPInfoBcOut,
                                            m_BcOutTC, MinLayer, T0Seg);
  return ntrack == -1 ? false : true;
#endif

  return false;
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearDCHits()
{
#if UseBcIn
  ClearBcInHits();
#endif
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearBcInHits()
{
  del::ClearContainerAll(m_TempBcInHC);
  del::ClearContainerAll(m_BcInHC);
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearBcOutHits()
{
  del::ClearContainerAll(m_BcOutHC);
}


//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksBcIn()
{
  del::ClearContainer(m_BcInTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOut()
{
  del::ClearContainer(m_BcOutTC);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcDCHits(std::vector<DCHC>& cont,
                         Bool_t applyRecursively)
{
  const Int_t n = (Int_t)cont.size();
  for(Int_t l=0; l<n; ++l){
    const Int_t m = (Int_t)cont[l].size();
    for(Int_t i=0; i<m; ++i){
      DCHit *hit = (cont[l])[i];
      if(!hit) continue;
      hit->ReCalcDC(applyRecursively);
    }
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcDCHits(Bool_t applyRecursively)
{

  ReCalcDCHits(m_BcOutHC, applyRecursively);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrack(DCLocalTC& cont,
                        Bool_t applyRecursively)
{
  const Int_t n = (Int_t)cont.size();
  for(Int_t i=0; i<n; ++i){
    DCLocalTrack *track = cont[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::ReCalcTrackBcIn(Bool_t applyRecursively)
{
  return ReCalcTrack(m_BcInTC, applyRecursively);
}
#endif

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrackBcOut(Bool_t applyRecursively)
{
  return ReCalcTrack(m_BcOutTC, applyRecursively);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcAll()
{
  ReCalcDCHits();
#if UseBcIn
  ReCalcTrackBcIn();
#endif
  ReCalcTrackBcOut();

  return true;
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutBcOut(Double_t chisqr)
{
  ChiSqrCut(m_BcOutTC, chisqr);
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCut(DCLocalTC& TrackCont,
                      Double_t chisqr)
{
  DCLocalTC ValidCand;
  DCLocalTC DeleteCand;
  for(auto& tempTrack : TrackCont){
    if(tempTrack->GetChiSquare() > chisqr){
      DeleteCand.push_back(tempTrack);
    }else{
      ValidCand.push_back(tempTrack);
    }
  }

  del::ClearContainer(DeleteCand);

  TrackCont.clear();
  TrackCont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), TrackCont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
void
DCAnalyzer::EraseEmptyHits(std::vector<DCHC>& HitCont)
{
  for(auto& hc: HitCont){
    for(auto it = hc.begin(); it != hc.end(); ){
      if((*it)->IsEmpty()) it = hc.erase(it);
      else ++it;
    }
  }
}

//_____________________________________________________________________________
void
DCAnalyzer::EraseEmptyHits(const TString& name)
{
  for(const auto& dcname: DCNameList){
    if(std::find(dcname.second.begin(), dcname.second.end(), name)
       != dcname.second.end()){
      if(dcname.first == "BcOut") EraseEmptyHits(m_BcOutHC);
    }
  }
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutBCOut(Double_t min_tot)
{
  for(const auto& name: DCNameList.at("BcOut")){
    TotCut(name, min_tot, true);
  }
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCut(const TString& name)
{
  Double_t min_tot = gUser.Get(Form("%s_TOT", name.Data()));
  TotCut(name, min_tot);
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCut(const TString& name, Double_t min_tot, Bool_t keep_nan)
{
  try {
    auto& HitCont = m_dc_hit_collection.at(name);
    for(auto& hit: HitCont){
      hit->TotCut(min_tot, keep_nan);
    }
    EraseEmptyHits(name);
  } catch (const std::out_of_range&) {
    spdlog::error("{} {} is not listed in dc_hit_collection",
                  FUNC_NAME.Data(), name.Data());
  }
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCut(const TString& name)
{
  Double_t min_dt = gUser.Get(Form("%s_DT", name.Data()), 0);
  Double_t max_dt = gUser.Get(Form("%s_DT", name.Data()), 1);
  DriftTimeCut(name, min_dt, max_dt, false);
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCut(const TString& name,
                         Double_t min_dt, Double_t max_dt, Bool_t select_1st)
{
  try {
    auto& HitCont = m_dc_hit_collection.at(name);
    for(auto& hit: HitCont){
      hit->DriftTimeCut(min_dt, max_dt, select_1st);
    }
    EraseEmptyHits(name);
  } catch (const std::out_of_range&) {
    spdlog::error("{} {} is not listed in dc_hit_collection",
                  FUNC_NAME.Data(), name.Data());
  }
}
