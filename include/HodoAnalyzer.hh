// -*- C++ -*-

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

#include <TString.h>

#include "DetectorID.hh"
#include "HodoHit.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

class RawData;
class BH2Hit;
class FiberHit;
class HodoCluster;
class BH2Cluster;
class FiberCluster;

using HodoHitContainer = std::vector<HodoHit*>;
using FiberHitContainer = std::vector<FiberHit*>;
typedef std::vector<BH2Hit*>   BH2HitContainer;
typedef std::vector< std::vector<FiberHit*> >
MultiPlaneFiberHitContainer;

typedef std::vector<HodoCluster*>  HodoClusterContainer;
typedef std::vector<BH2Cluster*>   BH2ClusterContainer;
typedef std::vector<FiberCluster*> FiberClusterContainer;
typedef std::vector< std::vector<FiberCluster*> >
MultiPlaneFiberClusterContainer;

//_____________________________________________________________________________
class HodoAnalyzer
{
public:
  static TString ClassName();
  explicit HodoAnalyzer(RawData* raw_data);
  ~HodoAnalyzer();

private:
  HodoAnalyzer(const HodoAnalyzer&);
  HodoAnalyzer& operator =(const HodoAnalyzer&);

private:
  template <typename T> using map_t = std::map<TString, T>;
  RawData*                    m_raw_data;
  map_t<HodoHitContainer>     m_hodo_hit_collection;
  map_t<HodoClusterContainer> m_hodo_cluster_collection;

public:
  template <typename T=HodoHit>
  Bool_t DecodeHits(const TString& name, Double_t max_time_diff=10.);
  const HodoHitContainer& GetHitContainer(const TString& name) const;
  const HodoClusterContainer& GetClusterContainer(const TString& name) const;
  Int_t  GetNHits(const TString& name) const
    { return GetHitContainer(name).size(); };
  Int_t  GetNClusters(const TString& name) const
    { return GetClusterContainer(name).size(); };

  template <typename T=HodoHit>
  T* GetHit(const TString& name, Int_t i=0) const
    { return dynamic_cast<T*>(GetHitContainer(name).at(i)); }

  template <typename T=HodoCluster>
  T* GetCluster(const TString& name, Int_t i=0) const
    { return dynamic_cast<T*>(GetClusterContainer(name).at(i)); }

  Bool_t ReCalcHit(const TString& name, Bool_t applyRecursively=false);
  Bool_t ReCalcCluster(const TString& name, Bool_t applyRecursively=false);
  Bool_t ReCalcAll();

  void TimeCut(const TString& name, Double_t min, Double_t max);
  void TotCut(const TString& name, Double_t min, Double_t max,
              Bool_t adopt_nan=true);

  BH2Cluster*  GetTime0BH2Cluster();
  HodoCluster* GetBtof0BH1Cluster(Double_t time0);

private:
  void ClearBH1Hits();
  void ClearBH2Hits();
  void ClearBACHits();
  void ClearTOFHits();
  void ClearLACHits();
  void ClearWCHits();
  void ClearWCSUMHits();
  void ClearBFTHits();

  template <typename T>
  void TimeCut(std::vector<T>& cont, Double_t min, Double_t max);
  template <typename T>
  void TotCut(std::vector<T>& cont,
              Double_t min, Double_t max, Bool_t adopt_nan);
  template <typename T>
  void WidthCutR(std::vector<T>& cont,
                 Double_t min, Double_t max, Bool_t adopt_nan);
  template <typename T>
  void AdcCut(std::vector<T>& cont, Double_t min, Double_t max);

  static Bool_t PairPlane(TString planeA, TString planeB); // pair or same plane
  static Int_t MakeUpClusters(const HodoHitContainer& HitCont,
                              HodoClusterContainer& ClusterCont,
                              Int_t MaxClusterSize,
                              Double_t MaxTimeDif);
  static Int_t MakeUpClusters(const BH2HitContainer& HitCont,
                              BH2ClusterContainer& ClusterCont,
                              Double_t maxTimeDif);
  static Int_t MakeUpClusters(const FiberHitContainer& cont,
                              FiberClusterContainer& ClusterCont,
                              Double_t maxTimeDif,
                              Int_t DifPairId);
};

//_____________________________________________________________________________
inline TString
HodoAnalyzer::ClassName()
{
  static TString s_name("HodoAnalyzer");
  return s_name;
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
HodoAnalyzer::DecodeHits(const TString& name, Double_t max_time_diff)
{
  auto& cont = m_hodo_hit_collection[name];
  for(auto& hit: cont){
    delete hit;
  }
  cont.clear();

  for(auto& rhit: m_raw_data->GetHodoRawHitContainer(name)){
    if(!rhit) continue;
    auto hit = new T(rhit);
    if(hit && hit->Calculate()){
      cont.push_back(hit);
    }else{
      delete hit;
    }
  }

  std::sort(cont.begin(), cont.end(), T::Compare);

#if 1
  // static const auto& gUser = UserParamMan::GetInstance();
  const auto MaxClusterSize = 3;
  const auto MaxTimeDiff    = 100;
  MakeUpClusters(cont, m_hodo_cluster_collection[name],
                 MaxClusterSize, MaxTimeDiff);
#endif

  return true;
}

//_____________________________________________________________________________
inline Bool_t
HodoAnalyzer::PairPlane(TString planeA, TString planeB)
{
  if(planeA.EndsWith("P", TString::kIgnoreCase))
    planeA.Chop();
  if(planeB.EndsWith("P", TString::kIgnoreCase))
    planeB.Chop();
  return planeA.EqualTo(planeB, TString::kIgnoreCase);
}

#endif
