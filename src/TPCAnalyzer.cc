// -*- C++ -*-

#include "TPCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "DCGeomMan.hh"
#include "DCGeomRecord.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "TPCRawData.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "TPCPadHelper.hh"
#include "TPCRawHit.hh"
#include "TPCHit.hh"
#include "TRandom3.h"
#include "TPCCluster.hh"
#include "TPCTrackSearch.hh"
#include "TPCLocalTrack.hh"

/* TPCTracking */
#define UseTpcCluster 1 // 1 : Common clustering method, 0 : Cluster size=1 no clustering

namespace
{
TRandom3 RandGen;

const auto& gUser   = UserParamMan::GetInstance();
const auto& gGeom   = DCGeomMan::GetInstance();

const double ztgt = tpc::Z_TARGET;
}
//_____________________________________________________________________________
TPCAnalyzer::TPCAnalyzer()
  : m_is_decoded(n_type),
    m_TPCHitCont(NumOfLayersTPC+1),
    m_TPCClCont(NumOfLayersTPC)
{
  for(Int_t i=0; i<n_type; ++i){
    m_is_decoded[i] = false;
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCAnalyzer::~TPCAnalyzer()
{
  ClearTPCHits();
  ClearTPCClusters();
  ClearTPCTracks();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::DecodeTPCHits(TPCRawData &TPCrawData,
			   const std::vector<Double_t> clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearTPCHits();

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    for(const auto& rhit: TPCrawData.GetTPCCorHitContainer(layer)){
      auto hit = new TPCHit(rhit);
      Int_t row  = hit -> GetRow();
      Int_t cobo_id = tpc::GetCoBoId(layer, row);
      if(cobo_id >= clock.size()){
        delete hit;
        return false; //No cobo input
      }

      if(hit->DoFit() && hit->Calculate(clock[cobo_id])){
        m_TPCHitCont[layer].push_back(hit);
      }else{
        delete hit;
      }
    }
  }

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::MakeUpTPCClusters(const TPCHitContainer& HitCont,
			       TPCClusterContainer& ClCont,
			       Double_t maxdy)
{
  static const Double_t MinClusterDe   = gUser.GetParameter("MinClusterDeTPC");
  static const Int_t    MinClusterSize = gUser.GetParameter("MinClusterSizeTPC");
  static const Double_t MinClusterYPos = gUser.GetParameter("MinClusterYPosTPC");
  static const Double_t MaxClusterYPos = gUser.GetParameter("MaxClusterYPosTPC");

  const auto nh = HitCont.size();
  if(nh==0) return false;

  std::vector<Int_t> joined(nh, 0);
  for(Int_t i=0; i<nh; ++i){
    if(joined[i] > 0) continue;
    TPCHitContainer CandCont;
    TPCHit* hit = HitCont[i];
    if(!hit || !hit->IsGood()) continue;
    Int_t layer = hit->GetLayer();
    CandCont.push_back(hit);
    joined[i]++;
    Double_t padlength = hit -> GetPadLength();
    TVector3 dist2tgt = hit -> GetPosition() - TVector3(0., 0., tpc::Z_TARGET);
    Double_t verticalpathlength_forpad =
      padlength*dist2tgt.y()/TMath::Hypot(dist2tgt.x(), dist2tgt.z());
    maxdy = TMath::Max(maxdy, verticalpathlength_forpad);
#if UseTpcCluster
    for(Int_t j=0; j<nh; ++j){
      if(i==j || joined[j]>0) continue;
      TPCHit* thit = HitCont[j];
      if(!thit || !thit->IsGood()) continue;
      Int_t rowID = thit->GetRow();
      for(const auto& c_hit: CandCont){
        Int_t c_rowID = c_hit->GetRow();
        if(tpc::IsClusterable(layer, rowID, c_rowID)
           && TMath::Abs(thit->GetY() - c_hit->GetY()) < maxdy){
          CandCont.push_back(thit);
          joined[j]++;
          break;
        }
      }
    }
#endif

    TPCCluster* cluster = new TPCCluster(layer, CandCont);
    if(!cluster) continue;
    if(cluster->Calculate()
       && cluster->GetDe()>=MinClusterDe && cluster->GetClusterSize()>=MinClusterSize
       && cluster->GetY()>=MinClusterYPos && cluster->GetY()<=MaxClusterYPos){
      ClCont.push_back(cluster);
    }else{
      delete cluster;
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::ReCalcTPCHits(const Int_t nhits,
			   const std::vector<Int_t>& pad,
			   const std::vector<Double_t>& time,
			   const std::vector<Double_t>& de,
			   const std::vector<Double_t>& clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

  static const Double_t MinCDe = gUser.GetParameter("MinCDeTPC");

  ClearTPCHits();
  ClearTPCClusters();

  if(nhits != pad.size() || nhits != time.size() || nhits != de.size()){
    hddaq::cerr << FUNC_NAME << " vector size mismatch" << std::endl;
    return false;
  }

  for(Int_t ih=0; ih<nhits; ih++){
    const Int_t layer = tpc::getLayerID(pad[ih]);
    const Double_t row = tpc::getRowID(pad[ih]);
    auto hit = new TPCHit(layer, row);
    hit->AddHit(de[ih], time[ih]);
    Int_t cobo_id = tpc::GetCoBoId(layer, row);
    if(cobo_id >= clock.size()){
      delete hit;
      return false; //No cobo input
    }

    if(hit->Calculate(clock[cobo_id]) && hit->GetCDe()>=MinCDe && hit->IsGood()){
      m_TPCHitCont[layer].push_back(hit);
    }else{
      delete hit;
    }
  }

#if 1
  static const Double_t MaxYDif = gUser.GetParameter("MaxYDifClusterTPC");
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    MakeUpTPCClusters(m_TPCHitCont[layer], m_TPCClCont[layer], MaxYDif);
  }
#endif

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
//HS-OFF: Track searching
Bool_t
TPCAnalyzer::TrackSearchTPC(Bool_t exclusive)
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::LocalTrackSearch(m_TPCClCont, m_TPCTC, m_TPCTCFailed, exclusive, MinLayer);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
//HS-On: Track Searching
Bool_t
TPCAnalyzer::TrackSearchTPCHelix(Bool_t exclusive)
{
  if(m_is_decoded[kTPCTracking]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");
  
  tpc::LocalTrackSearchHelix(m_TPCClCont, m_TPCTCHelix, m_TPCTCHelixInverted, m_TPCTCHelixFailed, m_TPCVC, m_TPCVCClustered, exclusive, MinLayer);

  m_is_decoded[kTPCTracking] = true;
  return true;
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCHits()
{
  del::ClearContainerAll(m_TPCHitCont);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCClusters()
{
  del::ClearContainerAll(m_TPCClCont);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCTracks()
{
  del::ClearContainer(m_TPCTC);
  del::ClearContainer(m_TPCTCFailed);
}
