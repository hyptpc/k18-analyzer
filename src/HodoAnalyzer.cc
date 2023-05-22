// -*- C++ -*-

#include "HodoAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "BH2Cluster.hh"
#include "BH2Hit.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
// #include "FLHit.hh"
#include "FuncName.hh"
#include "HodoHit.hh"
#include "HodoCluster.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

#define REQDE   0

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const Double_t MaxTimeDifBH1   =  2.0;
const Double_t MaxTimeDifBH2   =  2.0;
const Double_t MaxTimeDifBAC   = -1.0;
const Double_t MaxTimeDifTOF   = -1.0;
const Double_t MaxTimeDifLAC   = -1.0;
const Double_t MaxTimeDifWC    = -1.0;
const Double_t MaxTimeDifWCSUM = -1.0;
const Double_t MaxTimeDifBFT   =  8.0;
const Int_t    MaxSizeCl       = 8;
}

//_____________________________________________________________________________
HodoAnalyzer::HodoAnalyzer(RawData* raw_data)
  : m_raw_data(raw_data),
    m_hodo_hit_collection(),
    m_hodo_cluster_collection()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
HodoAnalyzer::~HodoAnalyzer()
{
  for(auto& elem: m_hodo_hit_collection)
    del::ClearContainer(elem.second);
  for(auto& elem: m_hodo_cluster_collection)
    del::ClearContainer(elem.second);
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
// Bool_t
// HodoAnalyzer::DecodeBFTHits()
// {
//   ClearBFTHits();
//   for(auto& hit: m_raw_data->GetHodoRawHitContainer("BFT")){
//     if(!hit) continue;
//     auto hp = new FiberHit(hit);
//     if(hp && hp->Calculate()){
//       delete hp;
//       //m_BFTCont.at(p).push_back(hp);
//     }else{
//       delete hp;
//     }
//   }
//     // std::sort(m_BFTCont.at(p).begin(), m_BFTCont.at(p).end(),
//     //           FiberHit::CompFiberHit);

// #if 0 // Cluster
//   FiberHitContainer cont_merge(m_BFTCont.at(0));
//   cont_merge.reserve(m_BFTCont.at(0).size() + m_BFTCont.at(1).size());
//   cont_merge.insert(cont_merge.end(), m_BFTCont.at(1).begin(),
//                     m_BFTCont.at(1).end());
//   std::sort(cont_merge.begin(), cont_merge.end(), FiberHit::CompFiberHit);
//   MakeUpClusters(cont_merge, m_BFTClCont, MaxTimeDifBFT, 3);
// #endif

//   return true;
// }

//_____________________________________________________________________________
// for HodoHit
Int_t
HodoAnalyzer::MakeUpClusters(const HodoHitContainer& HitCont,
                             HodoClusterContainer& ClusterCont,
                             Int_t MaxClusterSize,
                             Double_t MaxTimeDiff)
{
  del::ClearContainer(ClusterCont);

  if(HitCont.empty())
    return 0;

  for(Int_t iA=0, n=HitCont.size(); iA<n; ++iA){
    const auto& hitA = HitCont[iA];
    if(hitA->IsClusteredAll())
      continue;
    Double_t segA = hitA->SegmentId();
    Double_t seg = segA;
    for(Int_t jA=0, mA=hitA->GetEntries(); jA<mA; ++jA){
      if(hitA->IsClustered(jA))
        continue;
      HodoHitContainer CandCont;
      std::vector<Int_t> index;
      hitA->JoinCluster(jA);
      CandCont.push_back(hitA);
      index.push_back(jA);
      Double_t cmtA = hitA->CMeanTime(jA);
      for(Int_t iB=iA+1; iB<n; ++iB){
        const auto& hitB = HitCont[iB];
        Double_t segB = hitB->SegmentId();
        for(Int_t jB=0, mB=hitB->GetEntries(); jB<mB; ++jB){
          if(hitB->IsClustered(jB))
            continue;
          Double_t cmtB = hitB->CMeanTime(jB);
          if(true
             && TMath::Abs(segB - seg) <= 1
             && TMath::Abs(cmtA - cmtB) < MaxTimeDiff){
            seg = segB;
            hitB->JoinCluster(jB);
            CandCont.push_back(hitB);
            index.push_back(jB);
            break;
          }
        }
      }
      auto cluster = new HodoCluster(CandCont, index);
      if(cluster && cluster->IsGood()){
        ClusterCont.push_back(cluster);
      }else{
        delete cluster;
      }
    }
  }
  return ClusterCont.size();
}

//_____________________________________________________________________________
Int_t
HodoAnalyzer::MakeUpClusters(const BH2HitContainer& HitCont,
                             BH2ClusterContainer& ClusterCont,
                             Double_t maxTimeDif)
{
  del::ClearContainer(ClusterCont);

  // for(Int_t i=0, nh=HitCont.size(); i<nh; ++i){
  //   BH2Hit* hitA = HitCont[i];
  //   Int_t   segA = hitA->SegmentId();
  //   BH2Hit* hitB = 0;
  //   Int_t     iB = -1;
  //   Int_t     mB = -1;
  //   Int_t     mC = -1;
  //   Int_t   segB = -1;
  //   if(hitA->JoinedAllMhit()) continue;
  //   for(Int_t ma=0, n_mhitA=hitA->GetEntries(); ma<n_mhitA; ++ma){
  //     if(hitA->Joined(ma)) continue;
  //     Double_t cmtA = hitA->CMeanTime(ma);
  //     Double_t cmtB = -9999;
  //     for(int j=i+1; j<nh; ++j){
  //       BH2Hit *hit = HitCont[j];
  //       Int_t   seg = hit->SegmentId();
  //       for(Int_t mb=0, n_mhitB=hit->GetEntries(); mb<n_mhitB; ++mb){
  //         if(hit->Joined(mb)) continue;
  //         Double_t cmt = hit->CMeanTime(mb);
  //         if(std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif){
  //           hitB = hit; iB = j; mB = mb; segB = seg; cmtB = cmt; break;
  //         }
  //       }// for(mb:hitB)
  //     }// for(j:hitB)
  //     if(hitB){
  //       BH2Hit *hitC = nullptr;
  //       for(Int_t j=i+1; j<nh; ++j){
  //         if(j==iB) continue;
  //         BH2Hit *hit = HitCont[j];
  //         Int_t   seg = hit->SegmentId();
  //         for(Int_t mc=0, n_mhitC=hit->GetEntries(); mc<n_mhitC; ++mc){
  //           if(hit->Joined(mc)) continue;
  //           Double_t cmt = hit->CMeanTime(mc);
  //           if((std::abs(seg-segA)==1 && std::abs(cmt-cmtA)<maxTimeDif) ||
  //              (std::abs(seg-segB)==1 && std::abs(cmt-cmtB)<maxTimeDif)){
  //             hitC = hit; mC=mc; break;
  //           }
  //         }
  //       }
  //       if(hitC){
  //         hitA->SetJoined(ma);
  //         hitB->SetJoined(mB);
  //         hitC->SetJoined(mC);
  //         BH2Cluster *cluster = new BH2Cluster(hitA, hitB, hitC);
  //         cluster->SetIndex(ma, mB, mC);
  //         cluster->Calculate();
  //         if(cluster) ClusterCont.push_back(cluster);
  //       } else {
  //         hitA->SetJoined(ma);
  //         hitB->SetJoined(mB);
  //         BH2Cluster *cluster = new BH2Cluster(hitA, hitB);
  //         cluster->SetIndex(ma, mB);
  //         cluster->Calculate();
  //         if(cluster) ClusterCont.push_back(cluster);
  //       }
  //     } else {
  //       hitA->SetJoined(ma);
  //       BH2Cluster *cluster = new BH2Cluster(hitA);
  //       cluster->SetIndex(ma);
  //       cluster->Calculate();
  //       if(cluster) ClusterCont.push_back(cluster);
  //     }
  //   }// for(ma:hitA)
  // }// for(i:hitA)

  return ClusterCont.size();
}

//_____________________________________________________________________________
Int_t
HodoAnalyzer::MakeUpClusters(const FiberHitContainer& cont,
                             FiberClusterContainer& ClusterCont,
                             Double_t maxTimeDif,
                             Int_t DifPairId)
{
  del::ClearContainer(ClusterCont);

  for(Int_t seg=0, NofSeg=cont.size(); seg<NofSeg; ++seg){
    FiberHit* HitA = cont.at(seg);
    Bool_t fl_ClCandA = false;
    if(seg != (NofSeg -1)){
      if(DifPairId > (cont.at(seg+1)->PairId() - HitA->PairId())){
	fl_ClCandA = true;
      }
    }

    for(Int_t mhitA=0, NofHitA=HitA->GetEntries(); mhitA<NofHitA; ++mhitA){
      if(HitA->IsClustered(mhitA)) continue;
      FiberCluster *cluster = new FiberCluster;
      // cluster->push_back(new FLHit(HitA, mhitA));
      if(!fl_ClCandA){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	} else {
	  delete cluster;
	  cluster = nullptr;
	}
	continue;
      }
      // Start Search HitB
      Double_t cmtA    = (Double_t)HitA->GetCTime(mhitA);
      Int_t    NofHitB = cont.at(seg+1)->GetEntries();
      Bool_t   fl_HitB = false;
      Double_t cmtB    = -1;
      Int_t    CurrentPair = HitA->PairId();
      for(Int_t mhitB = 0; mhitB<NofHitB; ++mhitB){
	if(cont.at(seg+1)->IsClustered(mhitB)) continue;
	FiberHit* HitB = cont.at(seg+1);
	cmtB = (Double_t)HitB->GetCTime(mhitB);
	if(std::abs(cmtB-cmtA)<maxTimeDif){
	  // cluster->push_back(new FLHit(HitB, mhitB));
	  CurrentPair = HitB->PairId();
	  fl_HitB = true;
	  break;
	}
      }
      Bool_t fl_ClCandB  = false;
      if((seg+1) != (NofSeg -1)){
	if(DifPairId > (cont.at(seg+2)->PairId() - CurrentPair)){
	  fl_ClCandB = true;
	}
      }
      if(!fl_ClCandB){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	} else {
	  delete cluster;
	  cluster = nullptr;
	}
	continue;
      }
      // Start Search HitC
      Bool_t   fl_HitC = false;
      Double_t cmtC    = -1;
      for(Int_t mhitC=0, NofHitC=cont.at(seg+2)->GetEntries();
          mhitC<NofHitC; ++mhitC){
	if(cont.at(seg+2)->IsClustered(mhitC)) continue;
	FiberHit* HitC = cont.at(seg+2);
	cmtC = (Double_t)HitC->GetCTime(mhitC);
	if(true
           && std::abs(cmtC-cmtA)<maxTimeDif
           && !(fl_HitB && (std::abs(cmtC-cmtB)>maxTimeDif))
          ){
	  // cluster->push_back(new FLHit(HitC, mhitC));
	  CurrentPair = HitC->PairId();
	  fl_HitC = true;
	  break;
	}
      }
      Bool_t fl_ClCandC = false;
      if((seg+2) != (NofSeg -1)){
	if(DifPairId > (cont.at(seg+3)->PairId() - CurrentPair)){
	  fl_ClCandC = true;
	}
      }
      if(!fl_ClCandC){
	// there is no more candidates
	if(cluster->Calculate()){
	  ClusterCont.push_back(cluster);
	} else {
	  delete cluster;
	  cluster = nullptr;
	}
	continue;
      }
      // Start Search HitD
      Double_t cmtD = -1;
      for(Int_t mhitD=0, NofHitD=cont.at(seg+3)->GetEntries();
          mhitD<NofHitD; ++mhitD){
	if(cont.at(seg+3)->IsClustered(mhitD)) continue;
	FiberHit* HitD = cont.at(seg+3);
	cmtD = (Double_t)HitD->GetCTime(mhitD);
	if(true
           && std::abs(cmtD-cmtA)<maxTimeDif
           && !(fl_HitB && (std::abs(cmtD-cmtB)>maxTimeDif))
           && !(fl_HitC && (std::abs(cmtD-cmtC)>maxTimeDif))
          ){
	  // cluster->push_back(new FLHit(HitD, mhitD));
	  break;
	}
      }
      // Finish
      if(cluster->Calculate()){
	ClusterCont.push_back(cluster);
      } else {
	delete cluster;
	cluster = nullptr;
      }
    }
  }
  return ClusterCont.size();
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcHit(const TString& name, Bool_t applyRecursively)
{
  for(auto& cl: GetHitContainer(name)){
    cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcCluster(const TString& name, Bool_t applyRecursively)
{
  for(auto& cl: GetClusterContainer(name)){
    cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcAll()
{
  return true;
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCut(const TString& name, Double_t min, Double_t max)
{
  TimeCut(m_hodo_cluster_collection[name], min, max);
}

//_____________________________________________________________________________
// Implementation of Time cut for the cluster container
template <typename T>
void
HodoAnalyzer::TimeCut(std::vector<T>& cont,
                      Double_t min, Double_t max)
{
  std::vector<T> DeleteCand;
  std::vector<T> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t ctime = cont.at(i)->CMeanTime();
    if(min < ctime && ctime < max){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
void
HodoAnalyzer::TotCut(const TString& name, Double_t min, Double_t max,
                     Bool_t adopt_nan)
{
  TotCut(m_hodo_cluster_collection[name], min, max, adopt_nan);
}

//_____________________________________________________________________________
// Implementation of tot cut for the cluster container
template <typename T>
void
HodoAnalyzer::TotCut(std::vector<T>& cont,
                     Double_t min, Double_t max, Bool_t adopt_nan)
{
  std::vector<T> DeleteCand;
  std::vector<T> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t tot = cont.at(i)->TOT();
    if(TMath::IsNaN(tot) && adopt_nan){
      ValidCand.push_back(cont.at(i));
    }else if(min < tot && tot < max){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
// Implementation of tot cut for the cluster container
// template <typename T>
// void
// HodoAnalyzer::WidthCutR(std::vector<T>& cont,
//                         Double_t min, Double_t max, Bool_t adopt_nan)
// {
//   std::vector<T> DeleteCand;
//   std::vector<T> ValidCand;
//   for(Int_t i=0, n=cont.size(); i<n; ++i){
//     Double_t width = cont.at(i)->Width();
//     //Double_t width = -cont.at(i)->minWidth();//reverse

//     if(TMath::IsNaN(width) && adopt_nan){
//       ValidCand.push_back(cont.at(i));
//     }else if(min < width && width < max){
//       ValidCand.push_back(cont.at(i));
//     }else{
//       DeleteCand.push_back(cont.at(i));
//     }
//   }

//   del::ClearContainer(DeleteCand);

//   cont.clear();
//   cont.resize(ValidCand.size());
//   std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
//   ValidCand.clear();
// }

//_____________________________________________________________________________
// Implementation of ADC cut for the cluster container
template <typename T>
void
HodoAnalyzer::AdcCut(std::vector<T>& cont,
                     Double_t amin, Double_t amax)
{
  std::vector<T> DeleteCand;
  std::vector<T> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t adc = cont.at(i)->MaxAdcLow();
    if(amin < adc && adc < amax){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
BH2Cluster*
HodoAnalyzer::GetTime0BH2Cluster()
{
  static const Double_t MinMt = gUser.GetParameter("MtBH2", 0);
  static const Double_t MaxMt = gUser.GetParameter("MtBH2", 1);
#if REQDE
  static const Double_t MinDe = gUser.GetParameter("DeBH2", 0);
  static const Double_t MaxDe = gUser.GetParameter("DeBH2", 1);
#endif

  BH2Cluster* time0_cluster = nullptr;
  Double_t min_mt = -9999;
  for(const auto& cluster : GetClusterContainer("BH2")){
    Double_t mt = cluster->MeanTime();
    if(true
       && std::abs(mt) < std::abs(min_mt)
       && MinMt < mt && mt < MaxMt
#if REQDE
       && (MinDe < cluster->DeltaE() && cluster->DeltaE() < MaxDe)
#endif
      ){
      min_mt        = mt;
      time0_cluster = dynamic_cast<BH2Cluster*>(cluster);
    }// T0 selection
  }// for

  return time0_cluster;
}

//_____________________________________________________________________________
HodoCluster*
HodoAnalyzer::GetBtof0BH1Cluster(Double_t time0)
{
  static const Double_t MinBtof = gUser.GetParameter("BTOF",  0);
  static const Double_t MaxBtof = gUser.GetParameter("BTOF",  1);
#if REQDE
  static const Double_t MinDe   = gUser.GetParameter("DeBH1", 0);
  static const Double_t MaxDe   = gUser.GetParameter("DeBH1", 1);
#endif

  HodoCluster* time0_cluster = nullptr;
  Double_t min_btof            = -9999;
  for(const auto& cluster : GetClusterContainer("BH1")){
    Double_t cmt  = cluster->CMeanTime();
    Double_t btof = cmt - time0;
    if(true
       && std::abs(btof) < std::abs(min_btof)
       && MinBtof < btof && btof < MaxBtof
#if REQDE
       && (MinDe < cluster->DeltaE() && cluster->DeltaE() < MaxDe)
#endif
      ){
      min_btof      = btof;
      time0_cluster = cluster;
    }// T0 selection
  }// for

  return time0_cluster;
}


//_____________________________________________________________________________
const HodoHitContainer&
HodoAnalyzer::GetHitContainer(const TString& name) const
{
  auto itr = m_hodo_hit_collection.find(name);
  if(itr == m_hodo_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoHitContainer null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
const HodoClusterContainer&
HodoAnalyzer::GetClusterContainer(const TString& name) const
{
  auto itr = m_hodo_cluster_collection.find(name);
  if(itr == m_hodo_cluster_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoClusterContainer null_container;
    return null_container;
  }else{
    return itr->second;
  }
}
