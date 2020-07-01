// -*- C++ -*-

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

#include <TString.h>

#include "DetectorID.hh"
#include "RawData.hh"

class RawData;
class Hodo1Hit;
class Hodo2Hit;
class BH2Hit;
class FiberHit;
class FLHit;
class HodoCluster;
class BH2Cluster;
class FiberCluster;

typedef std::vector<Hodo1Hit*> Hodo1HitContainer;
typedef std::vector<Hodo2Hit*> Hodo2HitContainer;
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
  static TString       ClassName( void );
  static HodoAnalyzer& GetInstance( void );
  HodoAnalyzer( void );
  ~HodoAnalyzer( void );

private:
  HodoAnalyzer( const HodoAnalyzer& );
  HodoAnalyzer& operator =( const HodoAnalyzer& );

private:
  Hodo2HitContainer           m_BH1Cont;
  BH2HitContainer             m_BH2Cont;
  Hodo1HitContainer           m_BACCont;
  BH2HitContainer             m_E42BH2Cont;
  Hodo1HitContainer           m_SACCont;
  Hodo2HitContainer           m_TOFCont;
  Hodo1HitContainer           m_HtTOFCont;
  Hodo1HitContainer           m_LACCont;
  Hodo1HitContainer           m_LCCont;
  MultiPlaneFiberHitContainer m_BFTCont;
  FiberHitContainer           m_SCHCont;
  MultiPlaneFiberHitContainer m_SFTCont;
  MultiPlaneFiberHitContainer m_FBT1UCont;
  MultiPlaneFiberHitContainer m_FBT1DCont;
  MultiPlaneFiberHitContainer m_FBT2UCont;
  MultiPlaneFiberHitContainer m_FBT2DCont;

  HodoClusterContainer            m_BH1ClCont;
  BH2ClusterContainer             m_BH2ClCont;
  BH2ClusterContainer             m_E42BH2ClCont;
  HodoClusterContainer            m_BACClCont;
  HodoClusterContainer            m_SACClCont;
  HodoClusterContainer            m_TOFClCont;
  HodoClusterContainer            m_HtTOFClCont;
  HodoClusterContainer            m_LACClCont;
  HodoClusterContainer            m_LCClCont;
  FiberClusterContainer           m_BFTClCont;
  FiberClusterContainer           m_SCHClCont;
  MultiPlaneFiberClusterContainer m_SFTClCont;
  MultiPlaneFiberClusterContainer m_FBT1UClCont;
  MultiPlaneFiberClusterContainer m_FBT1DClCont;
  MultiPlaneFiberClusterContainer m_FBT2UClCont;
  MultiPlaneFiberClusterContainer m_FBT2DClCont;

public:
  Bool_t DecodeRawHits( RawData* rawData );
  Bool_t DecodeBH1Hits( RawData* rawData );
  Bool_t DecodeBH2Hits( RawData* rawData );
  Bool_t DecodeBACHits( RawData* rawData );
  Bool_t DecodeE42BH2Hits( RawData* rawData );
  Bool_t DecodeSACHits( RawData* rawData );
  Bool_t DecodeTOFHits( RawData* rawData );
  Bool_t DecodeHtTOFHits( RawData* rawData );
  Bool_t DecodeLACHits( RawData* rawData );
  Bool_t DecodeLCHits( RawData* rawData );
  Bool_t DecodeBFTHits( RawData* rawData );
  Bool_t DecodeSFTHits( RawData* rawData );
  Bool_t DecodeSCHHits( RawData* rawData );
  Bool_t DecodeFBT1Hits( RawData* rawData );
  Bool_t DecodeFBT2Hits( RawData* rawData );
  Int_t  GetNHitsBH1( void ) const { return m_BH1Cont.size(); };
  Int_t  GetNHitsBH2( void ) const { return m_BH2Cont.size(); };
  Int_t  GetNHitsBAC( void ) const { return m_BACCont.size(); };
  Int_t  GetNHitsE42BH2( void ) const { return m_E42BH2Cont.size(); };
  Int_t  GetNHitsSAC( void ) const { return m_SACCont.size(); };
  Int_t  GetNHitsTOF( void ) const { return m_TOFCont.size(); };
  Int_t  GetNHitsHtTOF( void ) const { return m_HtTOFCont.size(); };
  Int_t  GetNHitsLAC( void ) const { return m_LACCont.size(); };
  Int_t  GetNHitsLC( void ) const { return m_LCCont.size(); };
  Int_t  GetNHitsBFT( Int_t plane) const
  { return m_BFTCont.at( plane ).size(); };
  Int_t  GetNHitsSFT( Int_t plane ) const
  { return m_SFTCont.at( plane ).size(); };
  Int_t  GetNHitsSCH( void )  const { return m_SCHCont.size(); };
  Int_t  GetNHitsFBT1( Int_t layer, Int_t UorD) const
  { return UorD==0? m_FBT1UCont.size() : m_FBT1DCont.size(); }
  Int_t  GetNHitsFBT2( Int_t layer, Int_t UorD) const
  { return UorD==0? m_FBT2UCont.size() : m_FBT2DCont.size(); }

  inline Hodo2Hit* GetHitBH1( UInt_t i ) const;
  inline BH2Hit*   GetHitBH2( UInt_t i ) const;
  inline Hodo1Hit* GetHitBAC( UInt_t i ) const;
  inline BH2Hit*   GetHitE42BH2( UInt_t i ) const;
  inline Hodo1Hit* GetHitSAC( UInt_t i ) const;
  inline Hodo2Hit* GetHitTOF( UInt_t i ) const;
  inline Hodo1Hit* GetHitHtTOF( UInt_t i ) const;
  inline Hodo1Hit* GetHitLAC( UInt_t i ) const;
  inline Hodo1Hit* GetHitLC( UInt_t i ) const;
  inline FiberHit* GetHitBFT( Int_t plane, UInt_t seg ) const;
  inline FiberHit* GetHitSFT( Int_t plane, UInt_t seg ) const;
  inline FiberHit* GetHitSCH( UInt_t seg ) const;
  inline FiberHit* GetHitFBT1( Int_t layer, Int_t UorD, UInt_t seg ) const;
  inline FiberHit* GetHitFBT2( Int_t layer, Int_t UorD, UInt_t seg ) const;

  Int_t GetNClustersBH1( void ) const { return m_BH1ClCont.size(); }
  Int_t GetNClustersBH2( void ) const { return m_BH2ClCont.size(); }
  Int_t GetNClustersBAC( void ) const { return m_BACClCont.size(); }
  Int_t GetNClustersE42BH2( void ) const { return m_E42BH2ClCont.size(); }
  Int_t GetNClustersSAC( void ) const { return m_SACClCont.size(); }
  Int_t GetNClustersTOF( void ) const { return m_TOFClCont.size(); }
  Int_t GetNClustersHtTOF( void ) const{ return m_HtTOFClCont.size(); }
  Int_t GetNClustersLAC( void ) const { return m_LACClCont.size(); }
  Int_t GetNClustersLC( void ) const { return m_LCClCont.size(); }
  Int_t GetNClustersBFT( void ) const { return m_BFTClCont.size(); };
  Int_t GetNClustersSFT( Int_t layer ) const
  { return m_SFTClCont.at( layer ).size(); };
  Int_t GetNClustersSCH( void ) const { return m_SCHClCont.size(); };
  Int_t  GetNClustersFBT1( Int_t layer, Int_t UorD)  const
  { return UorD==0 ? m_FBT1UClCont.at(layer).size()
      : m_FBT1DClCont.at(layer).size(); }
  Int_t  GetNClustersFBT2( Int_t layer, Int_t UorD)  const
  { return UorD==0 ? m_FBT2UClCont.at(layer).size()
      : m_FBT2DClCont.at(layer).size(); }
  inline HodoCluster*  GetClusterBH1( UInt_t i ) const;
  inline BH2Cluster*   GetClusterBH2( UInt_t i ) const;
  inline HodoCluster*  GetClusterBAC( UInt_t i ) const;
  inline BH2Cluster*   GetClusterE42BH2( UInt_t i ) const;
  inline HodoCluster*  GetClusterSAC( UInt_t i ) const;
  inline HodoCluster*  GetClusterTOF( UInt_t i ) const;
  inline HodoCluster*  GetClusterHtTOF( UInt_t i ) const;
  inline HodoCluster*  GetClusterLAC( UInt_t i ) const;
  inline HodoCluster*  GetClusterLC( UInt_t i ) const;
  inline FiberCluster* GetClusterBFT( UInt_t i ) const;
  inline FiberCluster* GetClusterSFT( Int_t layer, UInt_t i ) const;
  inline FiberCluster* GetClusterSCH( UInt_t i ) const;
  inline FiberCluster* GetClusterFBT1( Int_t layer, Int_t UorD, UInt_t i ) const;
  inline FiberCluster* GetClusterFBT2( Int_t layer, Int_t UorD, UInt_t i ) const;

  Bool_t ReCalcBH1Hits( Bool_t applyRecursively=false );
  Bool_t ReCalcBH2Hits( Bool_t applyRecursively=false );
  Bool_t ReCalcBACHits( Bool_t applyRecursively=false );
  Bool_t ReCalcE42BH2Hits( Bool_t applyRecursively=false );
  Bool_t ReCalcSACHits( Bool_t applyRecursively=false );
  Bool_t ReCalcTOFHits( Bool_t applyRecursively=false );
  Bool_t ReCalcHtTOFHits( Bool_t applyRecursively=false );
  Bool_t ReCalcLACHits( Bool_t applyRecursively=false );
  Bool_t ReCalcLCHits( Bool_t applyRecursively=false );
  Bool_t ReCalcBH1Clusters( Bool_t applyRecursively=false );
  Bool_t ReCalcBH2Clusters( Bool_t applyRecursively=false );
  Bool_t ReCalcBACClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcE42BH2Clusters( Bool_t applyRecursively=false );
  Bool_t ReCalcSACClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcTOFClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcHtTOFClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcLACClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcLCClusters( Bool_t applyRecursively=false );
  Bool_t ReCalcFBT1Clusters( Bool_t applyRecursively=false );
  Bool_t ReCalcFBT2Clusters( Bool_t applyRecursively=false );
  Bool_t ReCalcAll( void );

  void TimeCutBH1( Double_t tmin, Double_t tmax );
  void TimeCutBH2( Double_t tmin, Double_t tmax );
  void TimeCutTOF( Double_t tmin, Double_t tmax );
  void TimeCutE42BH2( Double_t tmin, Double_t tmax );
  void TimeCutBFT( Double_t tmin, Double_t tmax );
  void TimeCutSFT( Int_t layer, Double_t tmin, Double_t tmax );
  void TimeCutSCH( Double_t tmin, Double_t tmax );
  void TimeCutFBT1( Int_t layer, Int_t UorD, Double_t tmin, Double_t tmax );
  void TimeCutFBT2( Int_t layer, Int_t UorD, Double_t tmin, Double_t tmax );

  void WidthCutBFT( Double_t min_width, Double_t max_width );
  void WidthCutSFT( Int_t layer, Double_t min_width, Double_t max_width );
  void WidthCutSCH( Double_t min_width, Double_t max_width );
  void WidthCutFBT1( Int_t layer, Int_t UorD,
		     Double_t min_width, Double_t max_width );
  void WidthCutFBT2( Int_t layer, Int_t UorD,
		     Double_t min_width, Double_t max_width );
  BH2Cluster*  GetTime0BH2Cluster( void );
  BH2Cluster*  GetTime0E42BH2Cluster( void );
  HodoCluster* GetBtof0BH1Cluster( Double_t time0 );

private:
  void ClearBH1Hits( void );
  void ClearBH2Hits( void );
  void ClearBACHits( void );
  void ClearE42BH2Hits( void );
  void ClearSACHits( void );
  void ClearTOFHits( void );
  void ClearHtTOFHits( void );
  void ClearLACHits( void );
  void ClearLCHits( void );
  void ClearBFTHits( void );
  void ClearSFTHits( void );
  void ClearSCHHits( void );
  void ClearFBT1Hits( void );
  void ClearFBT2Hits( void );

  template<typename TypeCluster>
  void TimeCut( std::vector<TypeCluster>& cont, Double_t tmin, Double_t tmax );
  template<typename TypeCluster>
  void WidthCut( std::vector<TypeCluster>& cont,
		 Double_t min_width, Double_t max_width, Bool_t adopt_nan );
  template<typename TypeCluster>
  void WidthCutR( std::vector<TypeCluster>& cont,
		  Double_t min_width, Double_t max_width, Bool_t adopt_nan );
  template<typename TypeCluster>
  void AdcCut( std::vector<TypeCluster>& cont, Double_t amin, Double_t amax );
  static Int_t MakeUpClusters( const Hodo1HitContainer& HitCont,
			       HodoClusterContainer& ClusterCont,
			       Double_t maxTimeDif );
  static Int_t MakeUpClusters( const Hodo2HitContainer& HitCont,
			       HodoClusterContainer& ClusterCont,
			       Double_t maxTimeDif );
  static Int_t MakeUpClusters( const BH2HitContainer& HitCont,
			       BH2ClusterContainer& ClusterCont,
			       Double_t maxTimeDif );
  static Int_t MakeUpClusters( const FiberHitContainer& cont,
			       FiberClusterContainer& ClusterCont,
			       Double_t maxTimeDif,
			       Int_t DifPairId );
  static Int_t MakeUpClusters( const FLHitContainer& cont,
			       FiberClusterContainer& ClusterCont,
			       Double_t maxTimeDif,
			       Int_t DifPairId );
  static Int_t MakeUpCoincidence( const FiberHitContainer& cont,
				  FLHitContainer& CoinCont,
				  Double_t maxTimeDif);
};

//_____________________________________________________________________________
inline TString
HodoAnalyzer::ClassName( void )
{
  static TString s_name("HodoAnalyzer");
  return s_name;
}

//_____________________________________________________________________________
inline HodoAnalyzer&
HodoAnalyzer::GetInstance( void )
{
  static HodoAnalyzer g_instance;
  return g_instance;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterBH1( UInt_t i ) const
{
  if( i<m_BH1ClCont.size() )
    return m_BH1ClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline BH2Cluster*
HodoAnalyzer::GetClusterBH2( UInt_t i ) const
{
  if( i<m_BH2ClCont.size() )
    return m_BH2ClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterBAC( UInt_t i ) const
{
  if( i<m_BACClCont.size() )
    return m_BACClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline BH2Cluster*
HodoAnalyzer::GetClusterE42BH2( UInt_t i ) const
{
  if( i<m_E42BH2ClCont.size() )
    return m_E42BH2ClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterSAC( UInt_t i ) const
{
  if( i<m_SACClCont.size() )
    return m_SACClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterTOF( UInt_t i ) const
{
  if( i<m_TOFClCont.size() )
    return m_TOFClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterHtTOF( UInt_t i ) const
{
  if( i<m_HtTOFClCont.size() )
    return m_HtTOFClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterLAC( UInt_t i ) const
{
  if( i<m_LACClCont.size() )
    return m_LACClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterLC( UInt_t i ) const
{
  if( i<m_LCClCont.size() )
    return m_LCClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterBFT( UInt_t i ) const
{
  if( i<m_BFTClCont.size() )
    return m_BFTClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterSFT( Int_t layer, UInt_t i ) const
{
  if( i<m_SFTClCont.at( layer ).size() )
    return m_SFTClCont.at( layer ).at( i );
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterSCH( UInt_t i ) const
{
  if( i<m_SCHClCont.size() )
    return m_SCHClCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterFBT1( Int_t layer, Int_t UorD, UInt_t seg ) const
{
  if(!(0 <= UorD && UorD <=1))
    return nullptr;
  if(UorD==0){
    if( seg<m_FBT1UClCont.at( layer ).size() )
      return m_FBT1UClCont.at( layer ).at( seg );
    else
      return nullptr;
  }else{
    if( seg<m_FBT1DClCont.at( layer ).size() )
      return m_FBT1DClCont.at( layer ).at( seg );
    else
      return nullptr;
  }
}

//_____________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterFBT2( Int_t layer, Int_t UorD, UInt_t seg ) const
{
  if(!(0 <= UorD && UorD <=1))
    return nullptr;
  if(UorD==0){
    if( seg<m_FBT2UClCont.at( layer ).size() )
      return m_FBT2UClCont.at( layer ).at( seg );
    else
      return nullptr;
  }else{
    if( seg<m_FBT2DClCont.at( layer ).size() )
      return m_FBT2DClCont.at( layer ).at( seg );
    else
      return nullptr;
  }
}

//_____________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHitBH1( UInt_t i ) const
{
  if( i<m_BH1Cont.size() )
    return m_BH1Cont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline BH2Hit*
HodoAnalyzer::GetHitBH2( UInt_t i ) const
{
  if( i<m_BH2Cont.size() )
    return m_BH2Cont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline BH2Hit*
HodoAnalyzer::GetHitE42BH2( UInt_t i ) const
{
  if( i<m_E42BH2Cont.size() )
    return m_E42BH2Cont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitBAC( UInt_t i ) const
{
  if( i<m_BACCont.size() )
    return m_BACCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitSAC( UInt_t i ) const
{
  if( i<m_SACCont.size() )
    return m_SACCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHitTOF( UInt_t i ) const
{
  if( i<m_TOFCont.size() )
    return m_TOFCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitHtTOF( UInt_t i ) const
{
  if( i<m_HtTOFCont.size() )
    return m_HtTOFCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitLAC( UInt_t i ) const
{
  if( i<m_LACCont.size() )
    return m_LACCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitLC( UInt_t i ) const
{
  if( i<m_LCCont.size() )
    return m_LCCont[i];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitBFT( Int_t plane, UInt_t seg ) const
{
  if( seg<m_BFTCont.at(plane).size() )
    return m_BFTCont.at(plane).at(seg);
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitSFT( Int_t plane, UInt_t seg ) const
{
  if( seg<m_SFTCont.at( plane ).size() )
    return m_SFTCont.at( plane ).at( seg );
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitSCH( UInt_t seg ) const
{
  if( seg<m_SCHCont.size() )
    return m_SCHCont[seg];
  else
    return nullptr;
}

//_____________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitFBT1( Int_t layer, Int_t UorD, UInt_t seg ) const
{
  if(!(0 <= UorD && UorD <=1))
    return nullptr;

  if(UorD==0){
    if( seg<m_FBT1UCont.at( layer ).size() )
      return m_FBT1UCont.at( layer ).at( seg );
    else
      return nullptr;
  }else{
    if( seg<m_FBT1DCont.at( layer ).size() )
      return m_FBT1DCont.at( layer ).at( seg );
    else
      return nullptr;
  }
}

//_____________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitFBT2( Int_t layer, Int_t UorD, UInt_t seg ) const
{
  if(!(0 <= UorD && UorD <=1))
    return nullptr;

  if(UorD==0){
    if( seg<m_FBT2UCont.at( layer ).size() )
      return m_FBT2UCont.at( layer ).at( seg );
    else
      return nullptr;
  }else{
    if( seg<m_FBT2DCont.at( layer ).size() )
      return m_FBT2DCont.at( layer ).at( seg );
    else
      return nullptr;
  }
}

#endif
