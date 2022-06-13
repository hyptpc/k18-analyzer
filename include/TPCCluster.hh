// -*- C++ -*-

#ifndef TPC_CLUSTER_HH
#define TPC_CLUSTER_HH

#include <cstddef>
#include <string>
#include <vector>

#include <TString.h>
#include <TVector3.h>

#include "TPCHit.hh"

typedef std::vector<TPCHit*>       TPCHitContainer;

//_____________________________________________________________________________
class TPCCluster
{
public:
  static const TString& ClassName();
  TPCCluster(Int_t layer, const TPCHitContainer& HitCont);
  TPCCluster(Double_t x, Double_t y, Double_t z, Double_t de); // for MC data
  ~TPCCluster();

private:
  Bool_t          m_is_good;
  Int_t           m_layer;
  Int_t           m_pad_id;
  Double_t        m_cluster_de;
  TVector3        m_cluster_position;
  TVector3        m_pos_center;
  TPCHitContainer m_hit_array;
  Bool_t          m_pos_calculated;
  Double_t        m_mrow;
  Int_t           m_mrow_int;
  Double_t        m_cluster_de_center;
  TPCHit*         m_mean_hit; // representative hit for tracking

public:
  void     AddTPCHit(TPCHit* hit);
  Bool_t   Calculate();
  Bool_t   CalculateWeightedMean(); // unused
  Bool_t   CalculateWeightedMeanTheta(); // unused
  void     ClearTPCHits();
  Int_t    GetClusterSize() const { return m_hit_array.size(); }
  Double_t GetDe() const { return m_cluster_de; }
  TPCHit*  GetMeanHit() const { return m_mean_hit; }
  Bool_t   IsGood() const { return m_is_good; }
  Int_t    MeanPadId() const { return m_pad_id; }
  Double_t MeanRow() const { return m_mrow; }
  Double_t GetRow() const { return MeanRow(); }
  Double_t GetMRow() const { return MeanRow(); }
  // Double_t GetDe_center() const { return m_cluster_de_center; }
  void     Print(Option_t* opt="") const;
  const TVector3& GetPosition() const { return m_cluster_position; }
  // TVector3 GetPos_center() const { return Position_CLcenter(); }
  // TVector3 Position_CLcenter() const;
  Double_t GetX() const { return m_cluster_position.X(); }
  Double_t GetY() const { return m_cluster_position.Y(); }
  Double_t GetZ() const { return m_cluster_position.Z(); }
  Double_t ResX() const;
  Double_t ResY() const;
  Double_t ResZ() const ;
  const TPCHitContainer& GetHitContainer() const { return m_hit_array; }
};

//_____________________________________________________________________________
inline const TString&
TPCCluster::ClassName()
{
  static TString s_name("TPCCluster");
  return s_name;
}

#endif
