// -*- C++ -*-

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

#include <TString.h>

#include "DetectorID.hh"
#include "RawData.hh"

class RawData;
class HodoHit;
class BH2Hit;
class FiberHit;
class FLHit;
class HodoCluster;
class BH2Cluster;
class FiberCluster;

using HodoHitContainer = std::vector<HodoHit*>;
typedef std::vector<BH2Hit*>   BH2HitContainer;
typedef std::vector<FiberHit*> FiberHitContainer;
typedef std::vector<FLHit*>    FLHitContainer;
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
  MultiPlaneFiberHitContainer m_BFTCont;
  FiberClusterContainer       m_BFTClCont;

public:
  Bool_t DecodeHits(const TString& name, Double_t max_time_diff=10.);
  Bool_t DecodeBFTHits();
  const HodoHitContainer& GetHitContainer(const TString& name) const;
  const HodoClusterContainer& GetClusterContainer(const TString& name) const;
  Int_t  GetNHits(const TString& name) const
    { return GetHitContainer(name).size(); };
  Int_t  GetNClusters(const TString& name) const
    { return GetClusterContainer(name).size(); };
  Int_t  GetNHitsBFT(Int_t plane) const
  { return m_BFTCont.at(plane).size(); };

  template <typename T=HodoHit>
  T* GetHit(const TString& name, Int_t i=0) const
    { return dynamic_cast<T*>(GetHitContainer(name).at(i)); }

  template <typename T=HodoCluster>
  T* GetCluster(const TString& name, Int_t i=0) const
    { return dynamic_cast<T*>(GetClusterContainer(name).at(i)); }

  FiberHit* GetHitBFT(Int_t plane, UInt_t seg) const { return m_BFTCont.at(plane).at(seg); }

  Int_t GetNClustersBFT() const { return m_BFTClCont.size(); };
  FiberCluster* GetClusterBFT(UInt_t i) const { return m_BFTClCont.at(i); }
  const FiberClusterContainer& GetClustersBFT() const { return m_BFTClCont; }

  Bool_t ReCalcHit(const TString& name, Bool_t applyRecursively=false);
  Bool_t ReCalcCluster(const TString& name, Bool_t applyRecursively=false);
  Bool_t ReCalcAll();

  void TimeCut(const TString& name, Double_t tmin, Double_t tmax);
  void TimeCutBFT(Double_t tmin, Double_t tmax);

  void WidthCutBFT(Double_t min_width, Double_t max_width);
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

  template<typename TypeCluster>
  void TimeCut(std::vector<TypeCluster>& cont, Double_t tmin, Double_t tmax);
  template<typename TypeCluster>
  void WidthCut(std::vector<TypeCluster>& cont,
		 Double_t min_width, Double_t max_width, Bool_t adopt_nan);
  template<typename TypeCluster>
  void WidthCutR(std::vector<TypeCluster>& cont,
		  Double_t min_width, Double_t max_width, Bool_t adopt_nan);
  template<typename TypeCluster>
  void AdcCut(std::vector<TypeCluster>& cont, Double_t amin, Double_t amax);
  static Int_t MakeUpClusters(const HodoHitContainer& HitCont,
                              HodoClusterContainer& ClusterCont,
                              Double_t maxTimeDif);
  static Int_t MakeUpClusters(const BH2HitContainer& HitCont,
                              BH2ClusterContainer& ClusterCont,
                              Double_t maxTimeDif);
  static Int_t MakeUpClusters(const FiberHitContainer& cont,
                              FiberClusterContainer& ClusterCont,
                              Double_t maxTimeDif,
                              Int_t DifPairId);
  static Int_t MakeUpClusters(const FLHitContainer& cont,
                              FiberClusterContainer& ClusterCont,
                              Double_t maxTimeDif,
                              Int_t DifPairId);
  static Int_t MakeUpCoincidence(const FiberHitContainer& cont,
                                 FLHitContainer& CoinCont,
                                 Double_t maxTimeDif);
};

//_____________________________________________________________________________
inline TString
HodoAnalyzer::ClassName()
{
  static TString s_name("HodoAnalyzer");
  return s_name;
}

#endif
