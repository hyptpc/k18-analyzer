/**
 *  file: HodoAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO_ANALYZER_HH
#define HODO_ANALYZER_HH

#include <vector>

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

typedef std::vector <Hodo1Hit*> Hodo1HitContainer;
typedef std::vector <Hodo2Hit*> Hodo2HitContainer;
typedef std::vector <BH2Hit*>   BH2HitContainer;
typedef std::vector <FiberHit*> FiberHitContainer;
typedef std::vector < std::vector< FiberHit*> > MultiPlaneFiberHitContainer;
typedef std::vector <FLHit*>    FLHitContainer;

typedef std::vector <HodoCluster*>  HodoClusterContainer;
typedef std::vector <BH2Cluster*>   BH2ClusterContainer;
typedef std::vector <FiberCluster*> FiberClusterContainer;
typedef std::vector < std::vector< FiberCluster*> > MultiPlaneFiberClusterContainer;

//______________________________________________________________________________
class HodoAnalyzer
{
public:
  HodoAnalyzer( void );
  ~HodoAnalyzer( void );

  static HodoAnalyzer& GetInstance( void );

private:
  HodoAnalyzer( const HodoAnalyzer& );
  HodoAnalyzer& operator =( const HodoAnalyzer& );

private:
  Hodo2HitContainer     m_BH1Cont;
  BH2HitContainer       m_BH2Cont;
  Hodo1HitContainer     m_BACCont;
  Hodo1HitContainer     m_PVACCont;
  Hodo1HitContainer     m_SACCont;
  Hodo1HitContainer     m_FACCont;
  Hodo1HitContainer     m_SSDTCont;
  Hodo2HitContainer     m_TOFCont;
  Hodo1HitContainer     m_LCCont;
  MultiPlaneFiberHitContainer m_BFTCont;
  FiberHitContainer     m_SCHCont;
  FiberHitContainer     m_FBHCont;
  MultiPlaneFiberHitContainer m_SFTCont;
  FLHitContainer        m_FBHCoinCont;
  HodoClusterContainer  m_BH1ClCont;
  BH2ClusterContainer   m_BH2ClCont;
  HodoClusterContainer  m_TOFClCont;
  HodoClusterContainer  m_LCClCont;
  FiberClusterContainer m_BFTClCont;
  FiberClusterContainer m_SCHClCont;
  FiberClusterContainer m_FBHClCont;
  MultiPlaneFiberClusterContainer m_SFTClCont;

public:
  bool DecodeRawHits( RawData* rawData );
  bool DecodeBH1Hits( RawData* rawData );
  bool DecodeBH2Hits( RawData* rawData );
  bool DecodeBACHits( RawData* rawData );
  bool DecodePVACHits( RawData* rawData );
  bool DecodeSACHits( RawData* rawData );
  bool DecodeFACHits( RawData* rawData );
  bool DecodeSSDTHits( RawData* rawData );
  bool DecodeTOFHits( RawData* rawData );
  bool DecodeLCHits( RawData* rawData );
  bool DecodeBFTHits( RawData* rawData );
  bool DecodeSFTHits( RawData* rawData );
  bool DecodeSCHHits( RawData* rawData );
  bool DecodeFBHHits( RawData* rawData );
  int  GetNHitsBH1( void )  const { return m_BH1Cont.size();  };
  int  GetNHitsBH2( void )  const { return m_BH2Cont.size();  };
  int  GetNHitsBAC( void )  const { return m_BACCont.size();  };
  int  GetNHitsPVAC( void ) const { return m_PVACCont.size(); };
  int  GetNHitsSAC( void )  const { return m_SACCont.size();  };
  int  GetNHitsFAC( void )  const { return m_FACCont.size();  };
  int  GetNHitsSSDT( void ) const { return m_SSDTCont.size(); };
  int  GetNHitsTOF( void )  const { return m_TOFCont.size();  };
  int  GetNHitsLC( void )   const { return m_LCCont.size();  };
  int  GetNHitsBFT( int plane)  const { return m_BFTCont.at( plane ).size(); };
  int  GetNHitsSFT( int plane ) const { return m_SFTCont.at( plane ).size(); };
  int  GetNHitsSCH( void )  const { return m_SCHCont.size();  };
  int  GetNHitsFBH( void )  const { return m_FBHCont.size();  };
  int  GetNHitsFBHCoin( void ) const { return m_FBHCoinCont.size();  };

  inline Hodo2Hit * GetHitBH1( std::size_t i )  const;
  inline BH2Hit   * GetHitBH2( std::size_t i )  const;
  inline Hodo1Hit * GetHitBAC( std::size_t i )  const;
  inline Hodo1Hit * GetHitPVAC( std::size_t i ) const;
  inline Hodo1Hit * GetHitFAC( std::size_t i )  const;
  inline Hodo1Hit * GetHitSAC( std::size_t i )  const;
  inline Hodo1Hit * GetHitSSDT( std::size_t i )  const;
  inline Hodo2Hit * GetHitTOF( std::size_t i )  const;
  inline Hodo1Hit * GetHitLC( std::size_t i )  const;
  inline FiberHit * GetHitBFT( int plane, std::size_t seg ) const;
  inline FiberHit * GetHitSFT( int plane, std::size_t seg ) const;
  inline FiberHit * GetHitSCH( std::size_t seg ) const;
  inline FiberHit * GetHitFBH( std::size_t seg ) const;
  inline FLHit    * GetHitFBHCoin( std::size_t seg ) const;

  int GetNClustersBH1( void ) const { return m_BH1ClCont.size(); };
  int GetNClustersBH2( void ) const { return m_BH2ClCont.size(); };
  int GetNClustersTOF( void ) const { return m_TOFClCont.size(); }
  int GetNClustersLC( void )  const { return m_LCClCont.size(); }
  int GetNClustersSAC( void ) const { return m_TOFClCont.size(); }
  int GetNClustersBFT( void ) const { return m_BFTClCont.size(); };
  int GetNClustersSFT( int layer ) const { return m_SFTClCont.at( layer ).size(); };
  int GetNClustersSCH( void ) const { return m_SCHClCont.size(); };
  int GetNClustersFBH( void ) const { return m_FBHClCont.size(); };

  inline HodoCluster  * GetClusterBH1( std::size_t i ) const;
  inline BH2Cluster   * GetClusterBH2( std::size_t i ) const;
  inline HodoCluster  * GetClusterTOF( std::size_t i ) const;
  inline HodoCluster  * GetClusterLC( std::size_t i ) const;
  inline HodoCluster  * GetClusterSAC( std::size_t i ) const;
  inline FiberCluster * GetClusterBFT( std::size_t i ) const;
  inline FiberCluster * GetClusterSFT( int layer, std::size_t i ) const;
  inline FiberCluster * GetClusterSCH( std::size_t i ) const;
  inline FiberCluster * GetClusterFBH( std::size_t i ) const;

  bool ReCalcBH1Hits( bool applyRecursively=false );
  bool ReCalcBH2Hits( bool applyRecursively=false );
  bool ReCalcBACHits( bool applyRecursively=false );
  bool ReCalcPVACHits( bool applyRecursively=false );
  bool ReCalcFACHits( bool applyRecursively=false );
  bool ReCalcSACHits( bool applyRecursively=false );
  bool ReCalcSSDTHits( bool applyRecursively=false );
  bool ReCalcTOFHits( bool applyRecursively=false );
  bool ReCalcLCHits( bool applyRecursively=false );
  bool ReCalcBH1Clusters( bool applyRecursively=false );
  bool ReCalcBH2Clusters( bool applyRecursively=false );
  bool ReCalcTOFClusters( bool applyRecursively=false );
  bool ReCalcLCClusters( bool applyRecursively=false );
  bool ReCalcAll( void );

  void TimeCutBH1(double tmin, double tmax);
  void TimeCutBH2(double tmin, double tmax);
  void TimeCutBFT(double tmin, double tmax);
  void TimeCutSFT( int layer, double tmin, double tmax );
  void TimeCutSCH(double tmin, double tmax);
  void TimeCutFBH(double tmin, double tmax);

  void WidthCutBFT(double min_width, double max_width);
  void WidthCutSFT( int layer, double min_width, double max_width);
  void WidthCutSCH(double min_width, double max_width);

private:
  void ClearBH1Hits( void );
  void ClearBH2Hits( void );
  void ClearBACHits( void );
  void ClearPVACHits( void );
  void ClearFACHits( void );
  void ClearSACHits( void );
  void ClearSSDTHits( void );
  void ClearTOFHits( void );
  void ClearLCHits( void );
  void ClearBFTHits();
  void ClearSFTHits();
  void ClearSCHHits();
  void ClearFBHHits();

  template<typename TypeCluster>
  void TimeCut(std::vector<TypeCluster>& cont, double tmin, double tmax);

  template<typename TypeCluster>
  void WidthCut(std::vector<TypeCluster>& cont, double min_width, double max_width, bool adopt_nan);

  static int MakeUpClusters( const Hodo2HitContainer& HitCont,
			     HodoClusterContainer& ClusterCont,
			     double maxTimeDif );

  static int MakeUpClusters( const BH2HitContainer& HitCont,
			     BH2ClusterContainer& ClusterCont,
			     double maxTimeDif );

  static int MakeUpClusters( const FiberHitContainer& cont,
			     FiberClusterContainer& ClusterCont,
			     double maxTimeDif,
			     int DifPairId);

  static int MakeUpClusters( const FLHitContainer& cont,
			     FiberClusterContainer& ClusterCont,
			     double maxTimeDif,
			     int DifPairId);

  static int MakeUpCoincidence( const FiberHitContainer& cont,
				FLHitContainer& CoinCont,
				double maxTimeDif);
};

//______________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterBH1( std::size_t i ) const
{
  if( i<m_BH1ClCont.size() )
    return m_BH1ClCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline BH2Cluster*
HodoAnalyzer::GetClusterBH2( std::size_t i ) const
{
  if( i<m_BH2ClCont.size() )
    return m_BH2ClCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterTOF( std::size_t i ) const
{
  if( i<m_TOFClCont.size() )
    return m_TOFClCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline HodoCluster*
HodoAnalyzer::GetClusterLC( std::size_t i ) const
{
  if( i<m_LCClCont.size() )
    return m_LCClCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterBFT( std::size_t i ) const
{
  if( i<m_BFTClCont.size() )
    return m_BFTClCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterSFT( int layer, std::size_t i ) const
{
  if( i<m_SFTClCont.at( layer ).size() )
    return m_SFTClCont.at( layer ).at( i );
  else
    return 0;
}

//______________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterSCH( std::size_t i ) const
{
  if( i<m_SCHClCont.size() )
    return m_SCHClCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline FiberCluster*
HodoAnalyzer::GetClusterFBH( std::size_t i ) const
{
  if( i<m_FBHClCont.size() )
    return m_FBHClCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHitBH1( std::size_t i ) const
{
  if( i<m_BH1Cont.size() )
    return m_BH1Cont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline BH2Hit*
HodoAnalyzer::GetHitBH2( std::size_t i ) const
{
  if( i<m_BH2Cont.size() )
    return m_BH2Cont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitBAC( std::size_t i ) const
{
  if( i<m_BACCont.size() )
    return m_BACCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitPVAC( std::size_t i ) const
{
  if( i<m_PVACCont.size() )
    return m_PVACCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitSAC( std::size_t i ) const
{
  if( i<m_SACCont.size() )
    return m_SACCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitFAC( std::size_t i ) const
{
  if( i<m_FACCont.size() )
    return m_FACCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitSSDT( std::size_t i ) const
{
  if( i<m_SSDTCont.size() )
    return m_SSDTCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline Hodo2Hit*
HodoAnalyzer::GetHitTOF( std::size_t i ) const
{
  if( i<m_TOFCont.size() )
    return m_TOFCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline Hodo1Hit*
HodoAnalyzer::GetHitLC( std::size_t i ) const
{
  if( i<m_LCCont.size() )
    return m_LCCont[i];
  else
    return 0;
}

//______________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitBFT( int plane, std::size_t seg ) const
{
  if( seg<m_BFTCont.at(plane).size() )
    return m_BFTCont.at(plane).at(seg);
  else
    return NULL;
}

//______________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitSFT( int plane, std::size_t seg ) const
{
  if( seg<m_SFTCont.at( plane ).size() )
    return m_SFTCont.at( plane ).at( seg );
  else
    return NULL;
}

//______________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitSCH( std::size_t seg ) const
{
  if( seg<m_SCHCont.size() )
    return m_SCHCont[seg];
  else
    return NULL;
}

//______________________________________________________________________________
inline FiberHit*
HodoAnalyzer::GetHitFBH( std::size_t seg ) const
{
  if( seg<m_FBHCont.size() )
    return m_FBHCont[seg];
  else
    return NULL;
}

//______________________________________________________________________________
inline FLHit*
HodoAnalyzer::GetHitFBHCoin( std::size_t seg ) const
{
  if( seg<m_FBHCoinCont.size() )
    return m_FBHCoinCont[seg];
  else
    return NULL;
}

inline HodoAnalyzer&
HodoAnalyzer::GetInstance( void )
{
    static HodoAnalyzer g_instance;
    return g_instance;
}
#endif
