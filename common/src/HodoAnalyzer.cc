// -*- C++ -*-

#include "HodoAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

// #include "BH2Cluster.hh"
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
}

//_____________________________________________________________________________
HodoAnalyzer::HodoAnalyzer(const RawData& raw_data)
  : m_raw_data(&raw_data),
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
void
HodoAnalyzer::DeCut(const TString& name, Double_t min, Double_t max)
{
  DeCut(m_hodo_cluster_collection[name], min, max);
}

//_____________________________________________________________________________
// Implementation of Time cut for the cluster container
template <typename T>
void
HodoAnalyzer::DeCut(std::vector<T>& cont,
                      Double_t min, Double_t max)
{
  std::vector<T> DeleteCand;
  std::vector<T> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t de = cont.at(i)->DeltaE();
    if(min < de && de < max){
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
const HodoCluster*
HodoAnalyzer::GetTime0BH2Cluster() const
{
  static const Double_t MinMt = gUser.GetParameter("MtBH2", 0);
  static const Double_t MaxMt = gUser.GetParameter("MtBH2", 1);
#if REQDE
  static const Double_t MinDe = gUser.GetParameter("DeBH2", 0);
  static const Double_t MaxDe = gUser.GetParameter("DeBH2", 1);
#endif

  HodoCluster* time0_cluster = nullptr;
  Double_t min_mt = DBL_MAX;
  for(const auto& cluster : GetClusterContainer("BH2")){
    Double_t mt = cluster->MeanTime();
    if(true
       && TMath::Abs(mt) < TMath::Abs(min_mt)
       && MinMt < mt && mt < MaxMt
#if REQDE
       && (MinDe < cluster->DeltaE() && cluster->DeltaE() < MaxDe)
#endif
      ){
      min_mt        = mt;
      time0_cluster = cluster;
    }
  }

  return time0_cluster;
}

//_____________________________________________________________________________
const HodoCluster*
HodoAnalyzer::GetBtof0BH1Cluster() const
{
  static const Double_t MinBtof = gUser.GetParameter("BTOF",  0);
  static const Double_t MaxBtof = gUser.GetParameter("BTOF",  1);
#if REQDE
  static const Double_t MinDe   = gUser.GetParameter("DeBH1", 0);
  static const Double_t MaxDe   = gUser.GetParameter("DeBH1", 1);
#endif

  HodoCluster* btof0_cluster = nullptr;

  Double_t time0 = Time0();
  if(TMath::IsNaN(time0)) return nullptr;
  Double_t min_btof = DBL_MAX;
  for(const auto& cluster : GetClusterContainer("BH1")){
    Double_t cmt  = cluster->CMeanTime();
    Double_t btof = cmt - time0;
    if(true
       && TMath::Abs(btof) < TMath::Abs(min_btof)
       && MinBtof < btof && btof < MaxBtof
#if REQDE
       && (MinDe < cluster->DeltaE() && cluster->DeltaE() < MaxDe)
#endif
      ){
      min_btof      = btof;
      btof0_cluster = cluster;
    }// T0 selection
  }// for

  return btof0_cluster;
}

//_____________________________________________________________________________
const HodoHC&
HodoAnalyzer::GetHitContainer(const TString& name) const
{
  auto itr = m_hodo_hit_collection.find(name);
  if(itr == m_hodo_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoHC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
const HodoCC&
HodoAnalyzer::GetClusterContainer(const TString& name) const
{
  auto itr = m_hodo_cluster_collection.find(name);
  if(itr == m_hodo_cluster_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoCC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Time0() const
{
  const auto& cl = GetTime0BH2Cluster();
  if(cl) return cl->CTime0();
  else   return TMath::QuietNaN();
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Btof0() const
{
  const auto& cl = GetBtof0BH1Cluster();
  if(cl) return cl->CMeanTime() - Time0();
  else   return TMath::QuietNaN();
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Time0Seg() const
{
  const auto& cl = GetTime0BH2Cluster();
  if(cl) return cl->MeanSeg();
  else   return TMath::QuietNaN();
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Btof0Seg() const
{
  const auto& cl = GetBtof0BH1Cluster();
  if(cl) return cl->MeanSeg();
  else   return TMath::QuietNaN();
}
