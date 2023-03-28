// -*- C++ -*-

#include "TPCAnalyzer.hh"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

//#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DebugTimer.hh"
#include "FuncName.hh"
#include "MathTools.hh"
#include "RawData.hh"
#include "UserParamMan.hh"
#include "DeleteUtility.hh"
#include "TPCPadHelper.hh"
//#include "TPCParamMan.hh"
//#include "TPCPositionCorrector.hh"
#include "TPCRawHit.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCTrackSearch.hh"

/* TPCTracking */
#define UseTpcCluster 1 // 1 : Common clustering method, 0 : Cluster size=1 no clustering

namespace
{

//const auto& gConf   = ConfMan::GetInstance();
//const auto& gGeom   = DCGeomMan::GetInstance();
//const auto& gTPC  = TPCParamMan::GetInstance();
//const auto& gTPCPos = TPCPositionCorrector::GetInstance();
const auto& gUser   = UserParamMan::GetInstance();

const Int_t MaxRowDifTPC = 2; // for cluster

}

//_____________________________________________________________________________
TPCAnalyzer::TPCAnalyzer()
  : m_is_decoded(n_type),
    m_TPCHitCont(NumOfLayersTPC+1),
    m_TempTPCHitCont(NumOfLayersTPC+1),
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
  ClearTracksTPC();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::MakeUpTPCClusters(const TPCHitContainer& HitCont,
                              TPCClusterContainer& ClCont,
                              Double_t maxdy)
{
  static const Double_t MinClusterDe = gUser.GetParameter("MinClusterDeTPC");
  static const Int_t MinClusterSize = gUser.GetParameter("MinClusterSizeTPC");
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
#if UseTpcCluster
    for(Int_t j=0; j<nh; ++j){
      if(i==j || joined[j]>0) continue;
      TPCHit* thit = HitCont[j];
      if(!thit || !thit->IsGood()) continue;
      Int_t rowID = thit->GetRow();
      for(const auto& c_hit: CandCont){
        Int_t c_rowID = c_hit->GetRow();
        if(tpc::IsClusterable(layer, rowID, c_rowID, MaxRowDifTPC)
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
    if(cluster->Calculate() && cluster->GetDe()>=MinClusterDe && cluster->GetClusterSize()>=MinClusterSize){
      ClCont.push_back(cluster);
    }else{
      delete cluster;
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::DecodeTPCHits(RawData *rawData, Double_t clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearTPCHits();
  ClearTPCClusters();

  for(Int_t layer=0; layer<=NumOfLayersTPC; ++layer){
    for(const auto& rhit: rawData->GetTPCCorHC(layer)){
      auto hit = new TPCHit(rhit);
      if(hit->DoFit() && hit->Calculate(clock)){
        m_TPCHitCont[layer].push_back(hit);
      }else{
        delete hit;
      }
    }
  }

#if 0 // Cluster analysis will be done by RecalcTPCHits() in Dst.
  static const Double_t MaxYDif = gUser.GetParameter("MaxYDifClusterTPC");
  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    MakeUpTPCClusters(m_TPCHitCont[layer], m_TPCClCont[layer], MaxYDif);
  }
#endif

  m_is_decoded[kTPC] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::ReCalcTPCHits(const Int_t nhits,
                          const std::vector<Int_t>& pad,
                          const std::vector<Double_t>& time,
                          const std::vector<Double_t>& de,
                          Double_t clock)
{
  if(m_is_decoded[kTPC]){
    hddaq::cerr << FUNC_NAME << " already decoded" << std::endl;
    return false;
  }

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
    if(hit->Calculate(clock)){
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
Bool_t
TPCAnalyzer::DecodeTPCHitsGeant4(const Int_t nhits,
                                const Double_t *x, const Double_t *y,
                                const Double_t *z, const Double_t *de)
{
  if(m_is_decoded[kTPC]){
    hddaq::cout << FUNC_NAME << " "
                << "already decoded" << std::endl;
    return true;
  }

  ClearTPCHits();
  ClearTPCClusters();

  for(Int_t i=0; i<nhits; i++){
    Int_t pad = tpc::findPadID(z[i], x[i]);
    Int_t layer = tpc::getLayerID(pad);
    Int_t row = tpc::getRowID(pad);
    auto hit = new TPCHit(layer, row);
    hit->AddHit(TMath::QuietNaN(), TMath::QuietNaN()); // allocate hit
    // tentative treatment
    if(de[i] == 0. || de[i] == TMath::QuietNaN()){
      hit->SetDe(1.e-3);
    }else{
      hit->SetDe(de[i]);
    }
    // end of tentative treatment
    hit->SetPosition(TVector3(x[i], y[i], z[i]));
    // hit->Print();
    m_TPCHitCont[layer].push_back(hit);
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
Bool_t
TPCAnalyzer::TrackSearchTPC(Bool_t exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::LocalTrackSearch(m_TPCClCont, m_TPCTC, m_TPCTCFailed, exclusive, MinLayer);
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TrackSearchTPCHelix(Bool_t exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::LocalTrackSearchHelix(m_TPCClCont, m_TPCTCHelix, m_TPCTCHelixFailed, exclusive, MinLayer);
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TrackSearchTPCHelix(std::vector<std::vector<TVector3>> K18VPs,
				std::vector<std::vector<TVector3>> KuramaVPs,
				Bool_t exclusive)
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::LocalTrackSearchHelix(K18VPs, KuramaVPs, m_TPCClCont, m_TPCTCHelix, m_TPCTCHelixFailed, exclusive, MinLayer);
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TestHoughTransform()
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::HoughTransformTest(m_TPCClCont, m_TPCTC, MinLayer);
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCAnalyzer::TestHoughTransformHelix()
{
  static const Int_t MinLayer = gUser.GetParameter("MinLayerTPC");

  tpc::HoughTransformTestHelix(m_TPCClCont, m_TPCTCHelix, MinLayer);
  return true;
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCHits()
{
  del::ClearContainerAll(m_TPCHitCont);
  del::ClearContainerAll(m_TempTPCHitCont);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTPCClusters()
{
  del::ClearContainerAll(m_TPCClCont);
}

//_____________________________________________________________________________
void
TPCAnalyzer::ClearTracksTPC()
{
  del::ClearContainer(m_TPCTC);
  del::ClearContainer(m_TPCTCFailed);
  del::ClearContainer(m_TPCTCHelix);
  del::ClearContainer(m_TPCTCHelixFailed);
}
