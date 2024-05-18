// -*- C++ -*-

#include "DCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <TH2D.h>

#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

#include "ConfMan.hh"
#include "DCDriftParamMan.hh"
#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCLocalTrack.hh"
#include "DCRawHit.hh"
#include "DCTrackSearch.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "FiberCluster.hh"
#include "FuncName.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "HodoCluster.hh"
#include "K18Parameters.hh"
//#include "K18TrackU2D.hh"
#include "K18TrackD2U.hh"
#include "S2sTrack.hh"
#include "MathTools.hh"
#include "MWPCCluster.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"

#define DefStatic
#include "DCParameters.hh"
#undef DefStatic

// Tracking routine selection __________________________________________________
/* BcInTracking */
#define UseBcIn    0 // not supported
/* BcOutTracking */
#define BcOut_XUV  0 // XUV Tracking (slow but accerate)
#define BcOut_Pair 1 // Pair plane Tracking (fast but bad for large angle track)
/* SdcInTracking */
#define SdcIn_XUV         0 // XUV Tracking (not used in S2S)
#define SdcIn_Pair        1 // Pair plane Tracking (fast but bad for large angle track)
#define SdcIn_Deletion    1 // Deletion method for too many combinations

namespace
{
const auto& gUnpacker = hddaq::unpacker::GUnpacker::get_instance();

using namespace K18Parameter;
const auto& gConf   = ConfMan::GetInstance();
const auto& gGeom   = DCGeomMan::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();

//_____________________________________________________________________________
const auto& valueNMR  = ConfMan::Get<Double_t>("FLDNMR");
const Double_t& pK18 = ConfMan::Get<Double_t>("PK18");
const Int_t& IdTOFUX = gGeom.DetectorId("TOF-UX");
const Int_t& IdTOFUY = gGeom.DetectorId("TOF-UY");
const Int_t& IdTOFDX = gGeom.DetectorId("TOF-DX");
const Int_t& IdTOFDY = gGeom.DetectorId("TOF-DY");

const Double_t TimeDiffToYTOF = 77.3511; // [mm/ns]

const Double_t MaxChiSqrS2sTrack = 10000.;
const Double_t MaxTimeDifMWPC       =   100.;

const Double_t kMWPCClusteringWireExtension =  1.0; // [mm]
const Double_t kMWPCClusteringTimeExtension = 10.0; // [nsec]

//_____________________________________________________________________________
inline Bool_t /* for MWPCCluster */
isConnectable(Double_t wire1, Double_t leading1, Double_t trailing1,
              Double_t wire2, Double_t leading2, Double_t trailing2,
              Double_t wExt,  Double_t tExt)
{
  Double_t w1Min = wire1 - wExt;
  Double_t w1Max = wire1 + wExt;
  Double_t t1Min = leading1  - tExt;
  Double_t t1Max = trailing1 + tExt;
  Double_t w2Min = wire2 - wExt;
  Double_t w2Max = wire2 + wExt;
  Double_t t2Min = leading2  - tExt;
  Double_t t2Max = trailing2 + tExt;
  Bool_t isWireOk = !(w1Min>w2Max || w1Max<w2Min);
  Bool_t isTimeOk = !(t1Min>t2Max || t1Max<t2Min);
#if 0
  hddaq::cout << __func__ << std::endl
              << " w1 = " << wire1
              << " le1 = " << leading1
              << " tr1 = " << trailing1 << "\n"
              << " w2 = " << wire2
              << " le2 = " << leading2
              << " tr2 = " << trailing2 << "\n"
              << " w1(" << w1Min << " -- " << w1Max << "), t1("
              << t1Min << " -- " << t1Max << ")\n"
              << " w2(" << w2Min << " -- " << w2Max << "), t2("
              << t2Min << " -- " << t2Max << ")\n"
              << " wire : " << isWireOk
              << ", time : " << isTimeOk
              << std::endl;
#endif
  return (isWireOk && isTimeOk);
}

//_____________________________________________________________________________
inline void
printConnectionFlag(const std::vector<std::deque<Bool_t> >& flag)
{
  for(Int_t i=0, n=flag.size(); i<n; ++i){
    hddaq::cout << std::endl;
    for(Int_t j=0, m=flag[i].size(); j<m; ++j){
      hddaq::cout << " " << flag[i][j];
    }
  }
  hddaq::cout << std::endl;
}
}

//_____________________________________________________________________________
DCAnalyzer::DCAnalyzer(const RawData& raw_data)
  : m_raw_data(&raw_data),
    m_dc_hit_collection(),
    m_max_v0diff(90.),
    m_TempBcInHC(NumOfLayersBcIn+1),
    m_BcInHC(NumOfLayersBcIn+1),
    m_BcOutHC(),
    m_SdcInHC(),
    m_SdcOutHC(),
    m_SdcInExTC(NumOfLayersSdcIn),
    m_SdcOutExTC(NumOfLayersSdcOut+1),
    m_MWPCClCont(NumOfLayersBcIn+1)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCAnalyzer::DCAnalyzer()
  : m_raw_data(),
    m_dc_hit_collection(),
    m_max_v0diff(90.),
    m_TempBcInHC(NumOfLayersBcIn+1),
    m_BcInHC(NumOfLayersBcIn+1),
    m_BcOutHC(),
    m_SdcInHC(),
    m_SdcOutHC(),
    m_SdcInExTC(NumOfLayersSdcIn),
    m_SdcOutExTC(NumOfLayersSdcOut+1),
    m_MWPCClCont(NumOfLayersBcIn+1)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
DCAnalyzer::~DCAnalyzer()
{
  for(auto& elem: m_dc_hit_collection)
    del::ClearContainer(elem.second);

  ClearS2sTracks();
#if UseBcIn
  ClearK18TracksU2D();
  ClearTracksBcIn();
#endif
  ClearK18TracksD2U();
  ClearTracksSdcOut();
  ClearTracksSdcIn();
  ClearTracksBcOut();
  ClearTracksBcOutSdcIn();
  ClearTracksSdcInSdcOut();
  ClearDCHits();
  ClearVtxHits();
  debug::ObjectCounter::decrease(ClassName());
}


//_____________________________________________________________________________
void
DCAnalyzer::PrintS2s(const TString& arg) const
{
  Int_t nn = m_S2sTC.size();
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
              << "   S2sTC.size : " << nn << std::endl;
  for(const auto& track: m_S2sTC){
    hddaq::cout << " Niter=" << std::setw(3) << track->Niteration()
                << " ChiSqr=" << track->ChiSquare()
                << " P=" << track->PrimaryMomentum().Mag()
                << " PL(TOF)=" << track->PathLengthToTOF()
                << std::endl;
  }
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::DecodeBcInHits()
{
  ClearBcInHits();

  for(Int_t layer=1; layer<=NumOfLayersBcIn; ++layer){
    const DCRHitContainer &RHitCont = m_raw_data->GetBcInRawHC(layer);
    Int_t nh = RHitCont.size();
    for(Int_t i=0; i<nh; ++i){
      DCRawHit *rhit  = RHitCont[i];
      DCHit    *thit  = new DCHit(rhit->PlaneId()+PlOffsBc, rhit->WireId());
      Int_t       nhtdc = rhit->GetTdcSize();
      if(!thit) continue;
      for(Int_t j=0; j<nhtdc; ++j){
        thit->SetTdcVal(rhit->GetTdc(j));
        thit->SetTdcTrailing(rhit->GetTrailing(j));
      }

      if(thit->CalcMWPCObservables())
        m_TempBcInHC[layer].push_back(thit);
      else
        delete thit;
    }

    // hddaq::cout<<"*************************************"<<std::endl;
    Int_t ncl = clusterizeMWPCHit(m_TempBcInHC[layer], m_MWPCClCont[layer]);
    // hddaq::cout<<"numCl="<< ncl << std::endl;
    for(Int_t i=0; i<ncl; ++i){
      MWPCCluster *p = m_MWPCClCont[layer][i];
      if(!p) continue;

      const MWPCCluster::Statistics& mean  = p->GetMean();
      const MWPCCluster::Statistics& first = p->GetFirst();
      Double_t mwire    = mean.m_wire;
      Double_t mwirepos = mean.m_wpos;
      Double_t mtime    = mean.m_leading;
      Double_t mtrail   = mean.m_trailing;

      DCHit *hit = new DCHit(layer+PlOffsBc, mwire);
      if(!hit) continue;
      hit->SetClusterSize(p->GetClusterSize());
      hit->SetMWPCFlag(true);
      hit->SetWire(mwire);
      hit->SetMeanWire(mwire);
      hit->SetMeanWirePosition(mwirepos);
      hit->SetTrailing(mtrail);
      hit->SetNomalizedData();
      hit->SetTdcVal(0);
      hit->SetTdcTrailing(0);

      if(hit->CalcMWPCObservables())
        m_BcInHC[layer].push_back(hit);
      else
        delete hit;
    }
    // hddaq::cout << "nh="<< m_BcInHC[layer].size() <<std::endl;
  }

  return true;
}
#endif

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
Bool_t
DCAnalyzer::DecodeSdcInHits()
{
  static const auto& digit_info =
    hddaq::unpacker::GConfig::get_instance().get_digit_info();
  m_SdcInHC.clear();
  Int_t plane_offset = 0;
  for(const auto& name: DCNameList.at("SdcIn")){
    Int_t id = digit_info.get_device_id(name.Data());
    Int_t n_plane = digit_info.get_n_plane(id);
    m_SdcInHC.resize(n_plane + m_SdcInHC.size());
    DecodeHits(name);
    for(const auto& hit: m_dc_hit_collection.at(name)){
      m_SdcInHC[hit->PlaneId() + plane_offset].push_back(hit);
    }
    plane_offset += n_plane;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeSdcOutHits(Double_t ofs_dt/*=0.*/)
{
  static const auto& digit_info =
    hddaq::unpacker::GConfig::get_instance().get_digit_info();
  m_SdcOutHC.clear();
  Int_t plane_offset = 0;
  for(const auto& name: DCNameList.at("SdcOut")){
    Int_t id = digit_info.get_device_id(name.Data());
    Int_t n_plane = digit_info.get_n_plane(id);
    m_SdcOutHC.resize(n_plane + m_SdcOutHC.size());
    DecodeHits(name);
    for(const auto& hit: m_dc_hit_collection.at(name)){
      m_SdcOutHC[hit->PlaneId() + plane_offset].push_back(hit);
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
    hddaq::cerr << FUNC_NAME << std::endl
                << " " << name << " is already decoded." << std::endl;
    return;
  }
  auto& HitCont = m_dc_hit_collection[name];
  del::ClearContainer(HitCont);
  HitCont.clear();
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
DCAnalyzer::DecodeSdcInHitsGeant4(const std::vector<Int_t>& nhit,
				  const std::vector<std::vector<TVector3>>& pos,
				  const std::vector<std::vector<Double_t>>& de)
{
  static const auto& digit_info =
    hddaq::unpacker::GConfig::get_instance().get_digit_info();
  m_SdcInHC.clear();
  Int_t plane_offset = 0;
  for(const auto& name: DCNameList.at("SdcIn")){
    Int_t id = digit_info.get_device_id(name.Data());
    Int_t n_plane = digit_info.get_n_plane(id);
    std::vector<Int_t> plane, layer;
    std::vector<TVector3> gpos;
    std::vector<Double_t> dE;
    for( Int_t i = 0; i < n_plane; i++ ){
      Int_t planeId = plane_offset + i;
      Int_t layerId = PlMinSdcIn + planeId;
      for( Int_t ih = 0; ih < nhit[planeId]; ih++ ){
	plane.push_back(i);
	layer.push_back(layerId);
	gpos.push_back(TVector3( pos[planeId][ih].x(),
				 pos[planeId][ih].y(),
				 pos[planeId][ih].z() ));
	dE.push_back(de[planeId][ih]);
      }
    }
    DecodeHitsGeant4(name, plane, layer, gpos, dE);
    m_SdcInHC.resize(n_plane + m_SdcInHC.size());
    for(const auto& hit: m_dc_hit_collection.at(name)){
      m_SdcInHC[hit->PlaneId() + plane_offset].push_back(hit);
    }
    plane_offset += n_plane;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeSdcOutHitsGeant4(const std::vector<Int_t>& nhit,
				   const std::vector<std::vector<TVector3>>& pos,
				   const std::vector<std::vector<Double_t>>& de)
{
  static const auto& digit_info =
    hddaq::unpacker::GConfig::get_instance().get_digit_info();
  m_SdcOutHC.clear();
  Int_t plane_offset = 0;
  for(const auto& name: DCNameList.at("SdcOut")){
    Int_t id = digit_info.get_device_id(name.Data());
    Int_t n_plane = digit_info.get_n_plane(id);
    std::vector<Int_t> plane, layer;
    std::vector<TVector3> gpos;
    std::vector<Double_t> dE;
    for( Int_t i = 0; i < n_plane; i++ ){
      Int_t planeId = plane_offset + i;
      Int_t layerId = PlMinSdcOut + planeId;
      for( Int_t ih = 0; ih < nhit[planeId]; ih++ ){
	plane.push_back(i);
	layer.push_back(layerId);
	gpos.push_back(TVector3( pos[planeId][ih].x(),
				 pos[planeId][ih].y(),
				 pos[planeId][ih].z() ));
	dE.push_back(de[planeId][ih]);
      }
    }
    DecodeHitsGeant4(name, plane, layer, gpos, dE);
    m_SdcOutHC.resize(n_plane + m_SdcOutHC.size());
    for(const auto& hit: m_dc_hit_collection.at(name)){
      m_SdcOutHC[hit->PlaneId() + plane_offset].push_back(hit);
    }
    plane_offset += n_plane;
  }
  return true;
}

//_____________________________________________________________________________
void
DCAnalyzer::DecodeHitsGeant4(const TString& name,
			     const std::vector<Int_t>& planeId, const std::vector<Int_t>& layerId,
			     const std::vector<TVector3>& gpos, const std::vector<Double_t>& de)
{
  if(m_dc_hit_collection.find(name) != m_dc_hit_collection.end()){
    hddaq::cerr << FUNC_NAME << std::endl
                << " " << name << " is already decoded." << std::endl;
    return;
  }
  auto& HitCont = m_dc_hit_collection[name];
  del::ClearContainer(HitCont);
  HitCont.clear();
  for( Int_t ih = 0, nh = planeId.size(); ih < nh; ih++ ){
    Int_t plane = planeId[ih];
    Int_t layer = layerId[ih];
    TVector3 pos  = gGeom.Global2LocalPos(layer, gpos[ih]);
    Int_t       wire = gGeom.CalcWireNumber(layer, pos.x());
    Double_t    wpos = gGeom.CalcWirePosition(layer, wire);
    Double_t    dl   = TMath::Abs(pos.x()-wpos);
    DCHit *p = nullptr;
    for(Int_t i=0, n=HitCont.size(); i<n; ++i){
      DCHit* q = HitCont[i];
      if(true
	 && q->PlaneId() == plane
	 && q->LayerId() == layer
	 && q->WireId()  == wire){
	p=q; break;
      }
    }
    if(!p){
      p = new DCHit(plane, layer, wire, wpos);
      HitCont.push_back(p);
    }
    p->SetDCObservablesGeant4(dl, de[ih]);
  }
  std::sort(HitCont.begin(), HitCont.end(), DCHit::Compare);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeRawHits()
{
  ClearDCHits();
#if UseBcIn
  DecodeBcInHits();
#endif
  DecodeBcOutHits();
  DecodeSdcInHits();
  DecodeSdcOutHits();
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeTOFHits(const HodoHC& HitCont)
{
  ClearTOFHits();

  // for the tilting plane case
  static const Double_t RA2 = gGeom.GetRotAngle2("TOF");
  static const TVector3 TOFPos[2] = { gGeom.GetGlobalPosition("TOF-UX"),
				      gGeom.GetGlobalPosition("TOF-DX") };

  for(const auto& hodo_hit: HitCont){
    const Double_t seg = hodo_hit->SegmentId()+1;
    const Double_t dt  = hodo_hit->TimeDiff();
    Int_t layer_x = -1;
    Int_t layer_y = -1;
    if((Int_t)seg%2==0){
      layer_x  = IdTOFUX;
      layer_y  = IdTOFUY;
    }
    if((Int_t)seg%2==1){
      layer_x  = IdTOFDX;
      layer_y  = IdTOFDY;
    }
    Double_t wpos = gGeom.CalcWirePosition(layer_x, seg-1);
    TVector3 w(wpos, 0., 0.);
    w.RotateY(RA2*TMath::DegToRad()); // for the tilting plane case
    const TVector3 hit_pos = TOFPos[(Int_t)seg%2] + w
      + TVector3(0., dt*TimeDiffToYTOF, 0.);
    // X
    DCHit *dc_hit_x = new DCHit(layer_x, seg-1);
    dc_hit_x->SetWirePosition(hit_pos.x());
    dc_hit_x->SetZ(hit_pos.z());
    dc_hit_x->SetTiltAngle(0.);
    dc_hit_x->SetDCData();
    m_TOFHC.push_back(dc_hit_x);
    // Y
    DCHit *dc_hit_y = new DCHit(layer_y, seg-1);
    dc_hit_y->SetWirePosition(hit_pos.y());
    dc_hit_y->SetZ(hit_pos.z());
    dc_hit_y->SetTiltAngle(90.);
    dc_hit_y->SetDCData();
    m_TOFHC.push_back(dc_hit_y);
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::DecodeTOFHits(const HodoClusterContainer& ClCont)
{
  ClearTOFHits();

  static const Double_t RA2 = gGeom.GetRotAngle2("TOF");
  static const TVector3 TOFPos[2] = {
    gGeom.GetGlobalPosition("TOF-UX"),
    gGeom.GetGlobalPosition("TOF-DX") };

  for(const auto& hodo_cluster: ClCont){
    const Double_t seg = hodo_cluster->MeanSeg()+1;
    const Double_t dt  = hodo_cluster->TimeDiff();
    Int_t layer_x = -1;
    Int_t layer_y = -1;
    if((Int_t)seg%2==0){
      layer_x  = IdTOFUX;
      layer_y  = IdTOFUY;
    }
    if((Int_t)seg%2==1){
      layer_x  = IdTOFDX;
      layer_y  = IdTOFDY;
    }
    Double_t wpos = gGeom.CalcWirePosition(layer_x, seg-1);
    TVector3 w(wpos, 0., 0.);
    w.RotateY(RA2*TMath::DegToRad());
    const TVector3& hit_pos = TOFPos[(Int_t)seg%2] + w
      + TVector3(0., dt*TimeDiffToYTOF, 0.);
    // X
    DCHit *dc_hit_x = new DCHit(layer_x, seg-1);
    dc_hit_x->SetWirePosition(hit_pos.x());
    dc_hit_x->SetZ(hit_pos.z());
    dc_hit_x->SetTiltAngle(0.);
    dc_hit_x->SetDCData();
    m_TOFHC.push_back(dc_hit_x);
    // Y
    DCHit *dc_hit_y = new DCHit(layer_y, seg-1);
    dc_hit_y->SetWirePosition(hit_pos.y());
    dc_hit_y->SetZ(hit_pos.z());
    dc_hit_y->SetTiltAngle(90.);
    dc_hit_y->SetDCData();
    m_TOFHC.push_back(dc_hit_y);
  }

  return true;
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::TrackSearchBcIn()
{
  track::MWPCLocalTrackSearch(&(m_BcInHC[1]), m_BcInTC);
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchBcIn(const std::vector<std::vector<DCHC> >& hc)
{
  track::MWPCLocalTrackSearch(hc, m_BcInTC);
  return true;
}
#endif

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
                                            m_BcOutTC, MinLayer);
  return ntrack == -1 ? false : true;
#endif

  return false;
}

//_____________________________________________________________________________
// Use with BH2Filter
Bool_t
DCAnalyzer::TrackSearchBcOut(const std::vector<std::vector<DCHC> >& hc, Int_t T0Seg)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerBcOut");

#if BcOut_Pair //Pair Plane Tracking Routine for BcOut
  Int_t ntrack = track::LocalTrackSearch(hc, PPInfoBcOut, NPPInfoBcOut, m_BcOutTC, MinLayer, T0Seg);
  return ntrack == -1 ? false : true;
#endif

#if BcOut_XUV  //XUV Tracking Routine for BcOut
  Int_t ntrack = track::LocalTrackSearchVUX(hc, PPInfoBcOut, NPPInfoBcOut, m_BcOutTC, MinLayer);
  return ntrack == -1 ? false : true;
#endif

  return false;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchSdcIn()
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerSdcIn");
  track::LocalTrackSearch(m_SdcInHC, PPInfoSdcIn, NPPInfoSdcIn, m_SdcInTC, MinLayer);
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchSdcOut()
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerSdcOut");

  track::LocalTrackSearchSdcOut(m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut,
                                m_SdcOutTC, MinLayer);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchSdcOut(const HodoHC& TOFCont)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerSdcOut");

  if(!DecodeTOFHits(TOFCont)) return false;

#if 0
  for(std::size_t i=0, n=m_TOFHC.size(); i<n; ++i){
    m_TOFHC[i]->Print();
  }
#endif

  track::LocalTrackSearchSdcOut(m_TOFHC, m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut+2,
                                m_SdcOutTC, MinLayer);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchSdcOut(const HodoClusterContainer& TOFCont)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerSdcOut");

  if(!DecodeTOFHits(TOFCont)) return false;

  track::LocalTrackSearchSdcOut(m_TOFHC, m_SdcOutHC, PPInfoSdcOut, NPPInfoSdcOut+2,
                                m_SdcOutTC, MinLayer);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchBcOutSdcIn()
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerBcOutSdcIn");

  track::LocalTrackSearchBcOutSdcIn(m_BcOutHC, PPInfoBcOut,
                                    m_SdcInHC, PPInfoSdcIn,
                                    NPPInfoBcOut, NPPInfoSdcIn,
                                    m_BcOutSdcInTC, MinLayer);

  return true;
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::TrackSearchK18U2D()
{
  ClearK18TracksU2D();

  Int_t nIn  = m_BcInTC.size();
  Int_t nOut = m_BcOutTC.size();

# if 0
  hddaq::cout << "**************************************" << std::endl;
  hddaq::cout << FUNC_NAME << ": #TracksIn=" << std::setw(3) << nIn
              << " #TracksOut=" << std::setw(3) << nOut << std::endl;
# endif

  if(nIn==0 || nOut==0) return true;

  for(Int_t iIn=0; iIn<nIn; ++iIn){
    DCLocalTrack *trIn = m_BcInTC[iIn];
# if 0
    hddaq::cout << "TrackIn  :" << std::setw(2) << iIn
                << " X0=" << trIn->GetX0() << " Y0=" << trIn->GetY0()
                << " U0=" << trIn->GetU0() << " V0=" << trIn->GetV0()
                << std::endl;
# endif
    if(!trIn->GoodForTracking() ||
       trIn->GetX0()<MinK18InX || trIn->GetX0()>MaxK18InX ||
       trIn->GetY0()<MinK18InY || trIn->GetY0()>MaxK18InY ||
       trIn->GetU0()<MinK18InU || trIn->GetU0()>MaxK18InU ||
       trIn->GetV0()<MinK18InV || trIn->GetV0()>MaxK18InV) continue;
    for(Int_t iOut=0; iOut<nOut; ++iOut){
      DCLocalTrack *trOut=m_BcOutTC[iOut];
# if 0
      hddaq::cout << "TrackOut :" << std::setw(2) << iOut
                  << " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
                  << " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
                  << std::endl;
# endif
      if(!trOut->GoodForTracking() ||
         trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
         trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
         trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
         trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV) continue;

# if 0
      hddaq::cout << FUNC_NAME << ": In -> " << trIn->GetChiSquare()
                  << " (" << std::setw(2) << trIn->GetNHit() << ") "
                  << "Out -> " << trOut->GetChiSquare()
                  << " (" << std::setw(2) << trOut->GetNHit() << ") "
                  << std::endl;
# endif

      K18TrackU2D *track=new K18TrackU2D(trIn, trOut, TMath::Abs(pK18));
      if(track && track->DoFit())
        m_K18U2DTC.push_back(track);
      else
        delete track;
    }
  }

# if 0
  hddaq::cout << "********************" << std::endl;
  {
    Int_t nn = m_K18U2DTC.size();
    hddaq::cout << FUNC_NAME << ": Before sorting. #Track="
                << nn << std::endl;
    for(Int_t i=0; i<nn; ++i){
      K18TrackU2D *tp = m_K18U2DTC[i];
      hddaq::cout << std::setw(3) << i
                  << " ChiSqr=" << tp->chisquare()
                  << " Delta=" << tp->Delta()
                  << " P=" << tp->P() << "\n";
      //      hddaq::cout<<"********************"<<std::endl;
      //      hddaq::cout << "In :"
      //         << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
      //         << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
      //         << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
      //         << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
      //         << "\n";
      //      hddaq::cout << "Out:"
      //         << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
      //         << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
      //         << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
      //         << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
      //         << std::endl;
    }
  }
# endif

  std::sort(m_K18U2DTC.begin(), m_K18U2DTC.end(), K18TrackU2DComp());

  return true;
}
#endif // UseBcIn

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchK18D2U(const std::vector<Double_t>& XinCont)
{
  ClearK18TracksD2U();

  std::size_t nIn  = XinCont.size();
  std::size_t nOut = m_BcOutTC.size();

  if(nIn==0 || nOut==0)
    return true;

  for(std::size_t iIn=0; iIn<nIn; ++iIn){
    Double_t LocalX = XinCont.at(iIn);
    for(std::size_t iOut=0; iOut<nOut; ++iOut){
      DCLocalTrack *trOut = m_BcOutTC[iOut];
#if 0
      hddaq::cout << "TrackOut :" << std::setw(2) << iOut
                  << " X0=" << trOut->GetX0() << " Y0=" << trOut->GetY0()
                  << " U0=" << trOut->GetU0() << " V0=" << trOut->GetV0()
                  << std::endl;
#endif
      if(!trOut->GoodForTracking() ||
         trOut->GetX0()<MinK18OutX || trOut->GetX0()>MaxK18OutX ||
         trOut->GetY0()<MinK18OutY || trOut->GetY0()>MaxK18OutY ||
         trOut->GetU0()<MinK18OutU || trOut->GetU0()>MaxK18OutU ||
         trOut->GetV0()<MinK18OutV || trOut->GetV0()>MaxK18OutV) continue;

#if 0
      hddaq::cout << FUNC_NAME
                  << "Out -> " << trOut->GetChiSquare()
                  << " (" << std::setw(2) << trOut->GetNHit() << ") "
                  << std::endl;
#endif

      K18TrackD2U *track = new K18TrackD2U(LocalX, trOut, TMath::Abs(pK18));
      if(track && track->CalcMomentumD2U())
	m_K18D2UTC.push_back(track);
      else
	delete track;
    }
  }

#if 0
  hddaq::cout<<"********************"<<std::endl;
  {
    Int_t nn = m_K18D2UTC.size();
    hddaq::cout << FUNC_NAME << ": Before sorting. #Track="
                << nn << std::endl;
    for(Int_t i=0; i<nn; ++i){
      K18TrackD2U *tp = m_K18D2UTC[i];
      hddaq::cout << std::setw(3) << i
                  << " ChiSqr=" << tp->chisquare()
                  << " Delta=" << tp->Delta()
                  << " P=" << tp->P() << "\n";
      //      hddaq::cout<<"********************"<<std::endl;
      //      hddaq::cout << "In :"
      //         << " X " << tp->Xin() << "(" << tp->TrackIn()->GetX0() << ")"
      //         << " Y " << tp->Yin() << "(" << tp->TrackIn()->GetY0() << ")"
      //         << " U " << tp->Uin() << "(" << tp->TrackIn()->GetU0() << ")"
      //         << " V " << tp->Vin() << "(" << tp->TrackIn()->GetV0() << ")"
      //         << "\n";
      //      hddaq::cout << "Out:"
      //         << " X " << tp->Xout() << "(" << tp->TrackOut()->GetX0() << ")"
      //         << " Y " << tp->Yout() << "(" << tp->TrackOut()->GetY0() << ")"
      //         << " U " << tp->Uout() << "(" << tp->TrackOut()->GetU0() << ")"
      //         << " V " << tp->Vout() << "(" << tp->TrackOut()->GetV0() << ")"
      //         << std::endl;
    }
  }
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchS2s()
{
  ClearS2sTracks();

  auto nIn = m_SdcInTC.size();
  auto nOut = m_SdcOutTC.size();
  if(nIn==0 || nOut==0) return true;
  for(Int_t iIn=0; iIn<nIn; ++iIn){
    const auto& trIn = GetTrackSdcIn(iIn);
    if(!trIn || !trIn->GoodForTracking()) continue;
    for(Int_t iOut=0; iOut<nOut; ++iOut){
      const auto& trOut = GetTrackSdcOut(iOut);
      if(!trOut || !trOut->GoodForTracking()) continue;
      auto trS2s = new S2sTrack(trIn, trOut);
      if(!trS2s) continue;
      Double_t u0In    = trIn->GetU0();
      Double_t u0Out   = trOut->GetU0();
      // Double_t v0In    = trIn->GetV0();
      // Double_t v0Out   = trOut->GetV0();
      Double_t bending = u0Out - u0In;
      // Double_t p[3] = { 0.08493, 0.2227, 0.01572 };
      // Double_t initial_momentum = p[0] + p[1]/(bending-p[2]);
      Double_t s = valueNMR/TMath::Abs(valueNMR);
      // Double_t initial_momentum = s*pK18;
      Double_t initial_momentum = 1.4;
      if(false
         && bending>0. && initial_momentum>0.){
        trS2s->SetInitialMomentum(initial_momentum);
      } else {
        trS2s->SetInitialMomentum(initial_momentum);
      }
      if(true
         && trS2s->DoFit()
         && trS2s->ChiSquare()<MaxChiSqrS2sTrack){
        // trS2s->Print("in "+FUNC_NAME);
        m_S2sTC.push_back(trS2s);
      }
      else{
        // trS2s->Print("in "+FUNC_NAME);
        delete trS2s;
      }
    }
  }

  std::sort(m_S2sTC.begin(), m_S2sTC.end(), S2sTrackComp());

#if 0
  PrintS2s("Before Deleting");
#endif

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::TrackSearchS2s(Double_t initial_momentum)
{
  ClearS2sTracks();

  Int_t nIn  = GetNtracksSdcIn();
  Int_t nOut = GetNtracksSdcOut();

  if(nIn==0 || nOut==0) return true;

  for(Int_t iIn=0; iIn<nIn; ++iIn){
    const auto& trIn = GetTrackSdcIn(iIn);
    if(!trIn->GoodForTracking()) continue;
    for(Int_t iOut=0; iOut<nOut; ++iOut){
      const auto& trOut = GetTrackSdcOut(iOut);
      if(!trOut->GoodForTracking()) continue;
      S2sTrack *trS2s = new S2sTrack(trIn, trOut);
      if(!trS2s) continue;
      trS2s->SetInitialMomentum(initial_momentum);
      if(trS2s->DoFit() && trS2s->ChiSquare()<MaxChiSqrS2sTrack){
        m_S2sTC.push_back(trS2s);
      }
      else{
        // trS2s->Print(" in "+FUNC_NAME);
        delete trS2s;
      }
    }// for(iOut)
  }// for(iIn)

  std::sort(m_S2sTC.begin(), m_S2sTC.end(), S2sTrackComp());

#if 0
  PrintS2s("Before Deleting");
#endif

  return true;
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearDCHits()
{
#if UseBcIn
  ClearBcInHits();
#endif
  ClearTOFHits();
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearBcInHits()
{
  del::ClearContainerAll(m_TempBcInHC);
  del::ClearContainerAll(m_BcInHC);
  del::ClearContainerAll(m_MWPCClCont);
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
DCAnalyzer::ClearSdcInHits()
{
  del::ClearContainerAll(m_SdcInHC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearSdcOutHits()
{
  del::ClearContainerAll(m_SdcOutHC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearVtxHits()
{
  del::ClearContainer(m_VtxPoint);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTOFHits()
{
  del::ClearContainer(m_TOFHC);
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearTracksBcIn()
{
  del::ClearContainer(m_BcInTC);
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOut()
{
  del::ClearContainer(m_BcOutTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcIn()
{
  del::ClearContainer(m_SdcInTC);
  del::ClearContainerAll(m_SdcInExTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcOut()
{
  del::ClearContainer(m_SdcOutTC);
  del::ClearContainerAll(m_SdcOutExTC);
}

//_____________________________________________________________________________
#if UseBcIn
void
DCAnalyzer::ClearK18TracksU2D()
{
  del::ClearContainer(m_K18U2DTC);
}
#endif

//_____________________________________________________________________________
void
DCAnalyzer::ClearK18TracksD2U()
{
  del::ClearContainer(m_K18D2UTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearS2sTracks()
{
  del::ClearContainer(m_S2sTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksBcOutSdcIn()
{
  del::ClearContainer(m_BcOutSdcInTC);
}

//_____________________________________________________________________________
void
DCAnalyzer::ClearTracksSdcInSdcOut()
{
  del::ClearContainer(m_SdcInSdcOutTC);
}

#if UseBcIn
//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcMWPCHits(std::vector<DCHC>& cont,
                           Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t l=0; l<n; ++l){
    const std::size_t m = cont[l].size();
    for(std::size_t i=0; i<m; ++i){
      DCHit *hit = (cont[l])[i];
      if(!hit) continue;
      hit->ReCalcMWPC(applyRecursively);
    }
  }
  return true;
}
#endif

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcDCHits(std::vector<DCHC>& cont,
                         Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t l=0; l<n; ++l){
    const std::size_t m = cont[l].size();
    for(std::size_t i=0; i<m; ++i){
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
#if UseBcIn
  ReCalcMWPCHits(m_TempBcInHC, applyRecursively);
  ReCalcMWPCHits(m_BcInHC, applyRecursively);
#endif

  ReCalcDCHits(m_BcOutHC, applyRecursively);
  ReCalcDCHits(m_SdcInHC, applyRecursively);
  ReCalcDCHits(m_SdcOutHC, applyRecursively);

  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrack(DCLocalTC& cont,
                        Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t i=0; i<n; ++i){
    DCLocalTrack *track = cont[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrack(K18TC& cont, Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t i=0; i<n; ++i){
    K18TrackD2U *track = cont[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrack(S2sTC& cont,
                        Bool_t applyRecursively)
{
  const std::size_t n = cont.size();
  for(std::size_t i=0; i<n; ++i){
    S2sTrack *track = cont[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::ReCalcTrackBcIn(Bool_t applyRecursively)
{
  return ReCalcTrack(m_BcInTC);
}
#endif

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrackBcOut(Bool_t applyRecursively)
{
  return ReCalcTrack(m_BcOutTC);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrackSdcIn(Bool_t applyRecursively)
{
  return ReCalcTrack(m_SdcInTC);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcTrackSdcOut(Bool_t applyRecursively)
{
  return ReCalcTrack(m_SdcOutTC);
}

//_____________________________________________________________________________
#if UseBcIn
Bool_t
DCAnalyzer::ReCalcK18TrackU2D(Bool_t applyRecursively)
{
  Int_t n = m_K18U2DTC.size();
  for(Int_t i=0; i<n; ++i){
    K18TrackU2D *track = m_K18U2DTC[i];
    if(track) track->ReCalc(applyRecursively);
  }
  return true;
}
#endif

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcK18TrackD2U(Bool_t applyRecursively)
{
  return ReCalcTrack(m_K18D2UTC, applyRecursively);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcS2sTrack(Bool_t applyRecursively)
{
  return ReCalcTrack(m_S2sTC, applyRecursively);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::ReCalcAll()
{
  ReCalcDCHits();
#if UseBcIn
  ReCalcTrackBcIn();
  ReCalcK18TrackU2D();
#endif
  ReCalcTrackBcOut();
  ReCalcTrackSdcIn();
  ReCalcTrackSdcOut();

  //ReCalcK18TrackD2U();
  ReCalcS2sTrack();

  return true;
}

//_____________________________________________________________________________
// Int_t
// clusterizeMWPCHit(const DCHC& hits,
//     MWPCClusterContainer& clusters)
// {
//   if (!clusters.empty()){
//       std::for_each(clusters.begin(), clusters.end(), DeleteObject());
//       clusters.clear();
//   }

//   const Int_t nhits = hits.size();
//   //   hddaq::cout << __func__ << " " << nhits << std::endl;
//   if (nhits==0)
//     return 0;

//   Int_t n = 0;
//   for (Int_t i=0; i<nhits; ++i){
//     const DCHit* h = hits[i];
//     if (!h)
//       continue;
//     n += h->GetTdcSize();
//   }

//   DCHC singleHits;
//   singleHits.reserve(n);
//   for (Int_t i=0; i<nhits; ++i){
//     const DCHit* h = hits[i];
//     if (!h)
//       continue;
//     Int_t nn = h->GetTdcSize();
//     for (Int_t ii=0; ii<nn; ++ii){
//       DCHit* htmp = new DCHit(h->GetLayer(), h->GetWire());
//       htmp->SetTdcVal(h->GetTdcVal());
//       htmp->SetTdcTrailing(h->GetTdcTrailing());
//       htmp->SetTrailingTime(h->GetTrailingTime());
//       htmp->SetDriftTime(h->GetDriftTime());
//       htmp->SetDriftLength(h->GetDriftLength());
//       htmp->SetTiltAngle(h->GetTiltAngle());
//       htmp->SetWirePosition(h->GetWirePosition());
//       htmp->setRangeCheckStatus(h->rangecheck(), 0);
//       singleHits.push_back(htmp);
//     }
//   }

//   std::vector<std::deque<Bool_t> > flag(n, std::deque<Bool_t>(n, false));
//   n = singleHits.size();
//   for (Int_t i=0;  i<n; ++i){
//     flag[i][i] = true;
//     const DCHit* h1 = singleHits[i];
//     //       h1->print("h1");
//     for (Int_t j=i+1; j<n; ++j){
//       const DCHit* h2 = singleHits[j];
//       //    h2->print("h2");
//       //    hddaq::cout << " (i,j) = (" << i << ", " << j << ")" << std::endl;
//       Bool_t val
//  = isConnectable(h1->GetWirePosition(),
//    h1->GetDriftTime(),
//    h1->GetTrailingTime(),
//    h2->GetWirePosition(),
//    h2->GetDriftTime(),
//    h2->GetTrailingTime(),
//    kMWPCClusteringWireExtension,
//    kMWPCClusteringTimeExtension);
//       //    hddaq::cout << "#D val = " << val << std::endl;
//       flag[i][j] = val;
//       flag[j][i] = val;
//     }
//   }

//   //   hddaq::cout << __func__ << "  before " << std::endl;
//   //   printConnectionFlag(flag);

//   const Int_t maxLoop = static_cast<Int_t>(std::log(x)/std::log(2.))+1;
//   for (Int_t loop=0; loop<maxLoop; ++loop){
//     std::vector<std::deque<Bool_t> > tmp(n, std::deque<Bool_t>(n, false));
//     for (Int_t i=0; i<n; ++i){
//       for (Int_t j=i; j<n; ++j){
//  for (Int_t k=0; k<n; ++k){
//    tmp[i][j] |= (flag[i][k] && flag[k][j]);
//    tmp[j][i] = tmp[i][j];
//  }
//       }
//     }
//     flag = tmp;
//     //       hddaq::cout << " n iteration = " << loop << std::endl;
//     //       printConnectionFlag(flag);
//   }

//   //   hddaq::cout << __func__ << "  after " << std::endl;
//   //   printConnectionFlag(flag);

//   std::set<Int_t> checked;
//   for (Int_t i=0; i<n; ++i){
//     if (checked.find(i)!=checked.end())
//       continue;
//     MWPCCluster* c = 0;
//     for (Int_t j=i; j<n; ++j){
//       if (flag[i][j]){
//  checked.insert(j);
//  if (!c) {
//    c = new MWPCCluster;
//    //     hddaq::cout << " new cluster " << std::endl;
//  }
//  //        hddaq::cout << " " << i << "---" << j << std::endl;
//  c->Add(singleHits[j]);
//       }
//     }

//     if (c){
//       c->Calculate();
//       clusters.push_back(c);
//     }
//   }

//   //   hddaq::cout << " end of " << __func__
//   //      << " : n = " << n << ", " << checked.size()
//   //      << std::endl;

//   //   hddaq::cout << __func__ << " n clusters = " << clusters.size() << std::endl;

//   return clusters.size();
// }

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutBcOut(Double_t chisqr)
{
  ChiSqrCut(m_BcOutTC, chisqr);
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSdcIn(Double_t chisqr)
{
  ChiSqrCut(m_SdcInTC, chisqr);
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCutSdcOut(Double_t chisqr)
{
  ChiSqrCut(m_SdcOutTC, chisqr);
}

//_____________________________________________________________________________
void
DCAnalyzer::ChiSqrCut(DCLocalTC& TrackCont,
                      Double_t chisqr)
{
  DCLocalTC DeleteCand;
  DCLocalTC ValidCand;
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
    auto i = hc.begin();
    while(i != hc.end()){
      if((*i)->IsEmpty()) i = hc.erase(i);
      else ++i;
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
      if(dcname.first == "SdcIn") EraseEmptyHits(m_SdcInHC);
      if(dcname.first == "SdcOut") EraseEmptyHits(m_SdcOutHC);
    }
  }
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutBCOut(Double_t min_tot)
{
  for(const auto& name: DCNameList.at("BcOut")){
    TotCut(name, min_tot, true);
  }// for(i)
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC1(Double_t min_tot)
{
  TotCut("SDC1", min_tot, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC2(Double_t min_tot)
{
  TotCut("SDC2", min_tot, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC3(Double_t min_tot)
{
  TotCut("SDC3", min_tot, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC4(Double_t min_tot)
{
  TotCut("SDC4", min_tot, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCutSDC5(Double_t min_tot)
{
  TotCut("SDC5", min_tot, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::TotCut(const TString& name, Double_t min_tot, Bool_t keep_nan)
{
  DCHC& HitCont = m_dc_hit_collection.at(name);
  for(auto& hit: HitCont){
    hit->TotCut(min_tot, keep_nan);
  }
  EraseEmptyHits(name);
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutBC34(Double_t min_dt, Double_t max_dt)
{
  for(const auto& name: DCNameList.at("BcOut")){
    DriftTimeCut(name, min_dt, max_dt, true);
  }
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC1(Double_t min_dt, Double_t max_dt)
{
  DriftTimeCut("SDC1", min_dt, max_dt, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC2(Double_t min_dt, Double_t max_dt)
{
  DriftTimeCut("SDC2", min_dt, max_dt, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC3(Double_t min_dt, Double_t max_dt)
{
  DriftTimeCut("SDC3", min_dt, max_dt, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC4(Double_t min_dt, Double_t max_dt)
{
  DriftTimeCut("SDC4", min_dt, max_dt, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCutSDC5(Double_t min_dt, Double_t max_dt)
{
  DriftTimeCut("SDC5", min_dt, max_dt, true);
}

//_____________________________________________________________________________
void
DCAnalyzer::DriftTimeCut(const TString& name,
                         Double_t min_dt, Double_t max_dt, Bool_t select_1st)
{
  DCHC& HitCont = m_dc_hit_collection.at(name);
  for(auto& hit: HitCont){
    hit->DriftTimeCut(min_dt, max_dt, select_1st);
  }
  EraseEmptyHits(name);
}

//_____________________________________________________________________________
Bool_t
DCAnalyzer::MakeBH2DCHit(Int_t t0seg)
{
  static const Double_t centerbh2[] = {
    -41.8, -19.3, -10.7, -3.6, 3.6, 10.7, 19.3, 41.8
  };

  Bool_t status = true;

  Double_t bh2pos = centerbh2[t0seg];
  DCHit *dchit = new DCHit(125, t0seg);
  dchit->SetTdcVal(0.);
  if(dchit->CalcFiberObservables()){
    dchit->SetWirePosition(bh2pos);
    m_BcOutHC[13].push_back(dchit);
  }else{
    delete dchit;
    status = false;
  }

  return status;
}
