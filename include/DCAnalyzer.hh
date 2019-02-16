/**
 *  file: DCAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef DC_ANALYZER_HH
#define DC_ANALYZER_HH

#include "DetectorID.hh"
#include "ThreeVector.hh"
#include <vector>

class DCHit;
class DCLocalTrack;
class K18TrackU2D;
class K18TrackD2U;
class KuramaTrack;
class RawData;
class MWPCCluster;
class FiberCluster;
class SsdCluster;
class HodoCluster;

class Hodo1Hit;
class Hodo2Hit;
class HodoAnalyzer;

typedef std::vector<DCHit*>        DCHitContainer;
typedef std::vector<MWPCCluster*>  MWPCClusterContainer;
typedef std::vector<SsdCluster*>   SsdClusterContainer;
typedef std::vector<DCLocalTrack*> DCLocalTrackContainer;
typedef std::vector<K18TrackU2D*>  K18TrackU2DContainer;
typedef std::vector<K18TrackD2U*>  K18TrackD2UContainer;
typedef std::vector<KuramaTrack*>  KuramaTrackContainer;

typedef std::vector<Hodo1Hit*> Hodo1HitContainer;
typedef std::vector<Hodo2Hit*> Hodo2HitContainer;
typedef std::vector<HodoCluster*> HodoClusterContainer;

//______________________________________________________________________________
class DCAnalyzer
{
public:
  DCAnalyzer( void );
  ~DCAnalyzer( void );

private:
  DCAnalyzer( const DCAnalyzer& );
  DCAnalyzer& operator =( const DCAnalyzer& );

private:
  enum e_type { k_BcIn,  k_BcOut,
		k_SdcIn, k_SdcOut,
		k_SsdIn, k_SsdOut,
		// k_SftIn, k_SftOut,
		k_TOF, n_type };
  std::vector<bool>     m_is_decoded;
  std::vector<int>      m_much_combi;
  std::vector<MWPCClusterContainer> m_MWPCClCont;
  std::vector<DCHitContainer>       m_TempBcInHC;
  std::vector<DCHitContainer>       m_BcInHC;
  std::vector<DCHitContainer>       m_BcOutHC;
  std::vector<DCHitContainer>       m_SdcInHC;
  std::vector<DCHitContainer>       m_SdcOutHC;
  std::vector<DCHitContainer>       m_SsdInHC;
  std::vector<DCHitContainer>       m_SsdOutHC;
  DCHitContainer        m_TOFHC;
  std::vector<SsdClusterContainer>  m_SsdInClCont;
  std::vector<SsdClusterContainer>  m_SsdOutClCont;
  DCHitContainer        m_VtxPoint;
  DCLocalTrackContainer m_BcInTC;
  DCLocalTrackContainer m_BcOutTC;
  DCLocalTrackContainer m_SdcInTC;
  DCLocalTrackContainer m_SdcOutTC;
  DCLocalTrackContainer m_SsdInTC;
  DCLocalTrackContainer m_SsdOutTC;
  DCLocalTrackContainer m_SsdXTC;   // for PreTracking
  DCLocalTrackContainer m_SsdYTC;   // for PreTracking
  DCLocalTrackContainer m_SsdInXTC; // for XiTracking
  DCLocalTrackContainer m_SsdInYTC; // for XiTracking
  K18TrackU2DContainer  m_K18U2DTC;
  K18TrackD2UContainer  m_K18D2UTC;
  KuramaTrackContainer  m_KuramaTC;
  DCLocalTrackContainer m_BcOutSdcInTC;
  DCLocalTrackContainer m_BcOutSsdInTC;
  DCLocalTrackContainer m_SsdOutSdcInTC;
  DCLocalTrackContainer m_SdcInSdcOutTC;
  DCLocalTrackContainer m_SsdInSsdOutTC;
  // Exclusive Tracks
  std::vector<DCLocalTrackContainer> m_SdcInExTC;
  std::vector<DCLocalTrackContainer> m_SdcOutExTC;

public:
  int  MuchCombinationSdcIn( void ) const { return m_much_combi[k_SdcIn]; }
  bool DecodeRawHits( RawData* rawData );
  // bool DecodeFiberHits( FiberCluster* FiberCl, int layer );
  bool DecodeFiberHits( RawData* rawData );
  bool DecodeBcInHits( RawData* rawData );
  bool DecodeBcOutHits( RawData* rawData );
  bool DecodeSdcInHits( RawData* rawData );
  bool DecodeSdcOutHits( RawData* rawData );
  bool DecodeSsdInHits( RawData* rawData );
  bool DecodeSsdOutHits( RawData* rawData );
  bool DecodeTOFHits( const Hodo2HitContainer& HitCont );
  bool DecodeTOFHits( const HodoClusterContainer& ClCont );
  //bool DecodeSimuHits( SimuData *simuData );
  int  ClusterizeMWPCHit( const DCHitContainer& hits,
			  MWPCClusterContainer& clusters );
  bool ClusterizeSsd( void );
  bool ClusterizeSsd( const DCHitContainer& HitCont,
		      SsdClusterContainer& ClCont,
		      const double MaxTimeDiff=25. );

  inline const DCHitContainer& GetTempBcInHC( int layer ) const;
  inline const DCHitContainer& GetBcInHC( int layer ) const;
  inline const DCHitContainer& GetBcOutHC( int layer ) const;
  inline const DCHitContainer& GetSdcInHC( int layer ) const;
  inline const DCHitContainer& GetSdcOutHC( int layer ) const;
  inline const DCHitContainer& GetSsdInHC( int layer ) const;
  inline const DCHitContainer& GetSsdOutHC( int layer ) const;
  inline const DCHitContainer& GetTOFHC( void ) const;

  bool TrackSearchBcIn( void );
  bool TrackSearchBcIn( const std::vector< std::vector<DCHitContainer> >& hc );
  bool TrackSearchBcOut( void );
  bool TrackSearchBcOut( const std::vector< std::vector<DCHitContainer> >& hc );
  bool TrackSearchSdcIn( void );
  bool TrackSearchSdcInFiber( void );
  bool TrackSearchSdcOut( void );
  bool TrackSearchSdcOut( const Hodo2HitContainer& HitCont );
  bool TrackSearchSdcOut( const HodoClusterContainer& ClCont );
  bool TrackSearchSsdIn( void );
  bool TrackSearchSsdOut( void );
  bool TrackSearchSsdInXY( void );

  int GetNtracksBcIn( void )   const { return m_BcInTC.size(); }
  int GetNtracksBcOut( void )  const { return m_BcOutTC.size(); }
  int GetNtracksSdcIn( void )  const { return m_SdcInTC.size(); }
  // int GetNtracksSdcIn( void )  const { return m_SsdOutSdcInTC.size(); }
  int GetNtracksSdcOut( void ) const { return m_SdcOutTC.size(); }
  int GetNtracksSsdIn( void )  const { return m_SsdInTC.size(); }
  int GetNtracksSsdOut( void ) const { return m_SsdOutTC.size(); }
  int GetNtracksSsdInX( void ) const { return m_SsdInXTC.size(); }
  int GetNtracksSsdInY( void ) const { return m_SsdInYTC.size(); }
  // for PreTracking
  int GetNtracksSsdX( void ) const { return m_SsdXTC.size(); }
  int GetNtracksSsdY( void ) const { return m_SsdYTC.size(); }
  // Exclusive Tracks
  int GetNtracksSdcInEx( int layer ) const { return m_SdcInExTC[layer].size(); }
  int GetNtracksSdcOutEx( int layer ) const { return m_SdcOutExTC[layer].size(); }

  inline DCLocalTrack* GetTrackBcIn( int i ) const;
  inline DCLocalTrack* GetTrackBcOut( int i ) const;
  inline DCLocalTrack* GetTrackSdcIn( int i ) const;
  inline DCLocalTrack* GetTrackSdcOut( int i ) const;
  inline DCLocalTrack* GetTrackSsdIn( int i ) const;
  inline DCLocalTrack* GetTrackSsdOut( int i ) const;
  inline DCLocalTrack* GetTrackSsdInX( int i ) const;
  inline DCLocalTrack* GetTrackSsdInY( int i ) const;
  // for PreTracking
  inline DCLocalTrack* GetTrackSsdX( int i ) const;
  inline DCLocalTrack* GetTrackSsdY( int i ) const;
  // Exclusive Tracks
  inline DCLocalTrack* GetTrackSdcInEx( int layer, int i ) const;
  inline DCLocalTrack* GetTrackSdcOutEx( int layer, int i ) const;

  bool TrackSearchK18U2D( void );
  bool TrackSearchK18D2U( const std::vector<double>& XinCont );
  bool TrackSearchKurama( double initial_momentum );
  bool TrackSearchKurama( void );

  ////////// Filters for SSD
  void DoTimeCorrectionSsd( const std::vector<double>& t0 );
  void DoTimeCorrectionSsdIn( const std::vector<double>& t0 );
  void DoTimeCorrectionSsdOut( const std::vector<double>& t0 );
  void ResetStatusSsd( void );
  void ResetStatusSsdIn( void );
  void ResetStatusSsdOut( void );
  void SlopeFilterSsd( void );
  void SlopeFilterSsdIn( void );
  void SlopeFilterSsdOut( void );
  void AdcPeakHeightFilterSsd( double min, double max );
  void AdcPeakHeightFilterSsdIn( double min, double max );
  void AdcPeakHeightFilterSsdOut( double min, double max );
  void AmplitudeFilterSsd( double min, double max, bool cluster_flag=false );
  void AmplitudeFilterSsdIn(double min, double max, bool cluster_flag=false );
  void AmplitudeFilterSsdOut( double min, double max, bool cluster_flag=false );
  void DeltaEFilterSsd( double min, double max, bool cluster_flag=false );
  void DeltaEFilterSsdIn( double min, double max, bool cluster_flag=false );
  void DeltaEFilterSsdOut( double min, double max, bool cluster_flag=false );
  void RmsFilterSsd( double min, double max );
  void RmsFilterSsdIn( double min, double max );
  void RmsFilterSsdOut( double min, double max );
  void DeviationFilterSsd( double min, double max );
  void DeviationFilterSsdIn( double min, double max );
  void DeviationFilterSsdOut( double min, double max );
  void TimeFilterSsd( double min, double max, bool cluster_flag=false );
  void TimeFilterSsdIn( double min, double max, bool cluster_flag=false );
  void TimeFilterSsdOut( double min, double max, bool cluster_flag=false );
  void ChisqrFilterSsd( double maxchi2 );
  void ChisqrFilterSsdIn( double maxchi2 );
  void ChisqrFilterSsdOut( double maxchi2 );

  void ChiSqrCutBcOut( double chisqr );
  void ChiSqrCutSdcIn( double chisqr );
  void ChiSqrCutSdcOut( double chisqr );
  void ChiSqrCutSsdIn( double chisqr );
  void ChiSqrCutSsdOut( double chisqr );

  void TotCutSDC2( double min_tot );
  void TotCutSDC3( double min_tot );

  int GetNTracksK18U2D( void ) const { return m_K18U2DTC.size(); }
  int GetNTracksK18D2U( void ) const { return m_K18D2UTC.size(); }
  int GetNTracksKurama( void ) const { return m_KuramaTC.size(); }

  inline K18TrackU2D  * GetK18TrackU2D( int i ) const;
  inline K18TrackD2U  * GetK18TrackD2U( int i ) const;
  inline KuramaTrack  * GetKuramaTrack( int i )    const;

  int GetNClustersMWPC( int layer ) const { return m_MWPCClCont[layer].size(); };
  int GetNClustersSsdIn( int layer ) const { return m_SsdInClCont[layer].size(); };
  int GetNClustersSsdOut( int layer ) const { return m_SsdOutClCont[layer].size(); };

  inline const MWPCClusterContainer & GetClusterMWPC( int layer ) const;
  inline const SsdClusterContainer  & GetClusterSsdIn( int layer ) const;
  inline const SsdClusterContainer  & GetClusterSsdOut( int layer ) const;

  void PrintKurama( const std::string& arg="" ) const;

  bool ReCalcMWPCHits( std::vector<DCHitContainer>& cont,
		       bool applyRecursively=false );
  bool ReCalcDCHits( std::vector<DCHitContainer>& cont,
		     bool applyRecursively=false );
  bool ReCalcDCHits( bool applyRecursively=false );

  bool ReCalcTrack( DCLocalTrackContainer& cont, bool applyRecursively=false );
  bool ReCalcTrack( K18TrackD2UContainer& cont, bool applyRecursively=false );
  bool ReCalcTrack( KuramaTrackContainer& cont, bool applyRecursively=false );

  bool ReCalcTrackBcIn( bool applyRecursively=false );
  bool ReCalcTrackBcOut( bool applyRecursively=false );
  bool ReCalcTrackSdcIn( bool applyRecursively=false );
  bool ReCalcTrackSdcOut( bool applyRecursively=false );
  bool ReCalcTrackSsdIn( bool applyRecursively=false );
  bool ReCalcTrackSsdOut( bool applyRecursively=false );

  bool ReCalcK18TrackD2U( bool applyRecursively=false );
  // bool ReCalcK18TrackU2D( bool applyRecursively=false );
  bool ReCalcKuramaTrack( bool applyRecursively=false );

  bool ReCalcAll( void );

  bool TrackSearchBcOutSdcIn( void );
  bool TrackSearchBcOutSsdIn( void );
  bool TrackSearchSsdOutSdcIn( void );
  bool TrackSearchSdcInSdcOut( void );
  bool TrackSearchSsdInSsdOut( void );
  int GetNtracksBcOutSdcIn( void ) const { return m_BcOutSdcInTC.size(); }
  int GetNtracksBcOutSsdIn( void ) const { return m_BcOutSsdInTC.size(); }
  int GetNtracksSsdOutSdcIn( void ) const { return m_SsdOutSdcInTC.size(); }
  int GetNtracksSdcInSdcOut( void ) const { return m_SdcInSdcOutTC.size(); }
  int GetNtracksSsdInSsdOut( void ) const { return m_SsdInSsdOutTC.size(); }
  inline DCLocalTrack * GetTrackBcOutSdcIn( int i ) const;
  inline DCLocalTrack * GetTrackBcOutSsdIn( int i ) const;
  inline DCLocalTrack * GetTrackSsdOutSdcIn( int i ) const;
  inline DCLocalTrack * GetTrackSdcInSdcOut( int i ) const;
  inline DCLocalTrack * GetTrackSsdInSsdOut( int i ) const;

protected:
  void ClearDCHits( void );
  void ClearBcInHits( void );
  void ClearBcOutHits( void );
  void ClearSdcInHits( void );
  void ClearSdcOutHits( void );
  void ClearSsdInHits( void );
  void ClearSsdOutHits( void );
  // void ClearSftHits( void );
  void ClearTOFHits( void );
  void ClearVtxHits( void );
  void ClearTracksBcIn( void );
  void ClearTracksBcOut( void );
  void ClearTracksSdcIn( void );
  void ClearTracksSdcOut( void );
  void ClearTracksSsdIn( void );
  void ClearTracksSsdOut( void );
  void ClearTracksSsdXY( void );
  void ClearTracksBcOutSdcIn( void );
  void ClearTracksBcOutSsdIn( void );
  void ClearTracksSsdOutSdcIn( void );
  void ClearTracksSdcInSdcOut( void );
  void ClearTracksSsdInSsdOut( void );
  void ClearK18TracksU2D( void );
  void ClearK18TracksD2U( void );
  void ClearKuramaTracks( void );
  void ChiSqrCut( DCLocalTrackContainer& cont, double chisqr );
  void TotCut( DCHitContainer& cont, double min_tot, bool adopt_nan );
  static int MakeUpMWPCClusters( const DCHitContainer& HitCont,
  				 MWPCClusterContainer& ClusterCont,
  				 double maxTimeDif );
public:
  void ResetTracksBcIn( void )        { ClearTracksBcIn();        }
  void ResetTracksBcOut( void )       { ClearTracksBcOut();       }
  void ResetTracksSdcIn( void )       { ClearTracksSdcIn();       }
  void ResetTracksSdcOut( void )      { ClearTracksSdcOut();      }
  void ResetTracksSsdIn( void )       { ClearTracksSsdIn();       }
  void ResetTracksSsdXY( void )       { ClearTracksSsdXY();       }
  void ResetTracksSsdOut( void )      { ClearTracksSsdOut();      }
  void ResetTracksBcOutSdcIn( void )  { ClearTracksBcOutSdcIn();  }
  void ResetTracksBcOutSsdIn( void )  { ClearTracksBcOutSsdIn();  }
  void ResetTracksSsdOutSdcIn( void ) { ClearTracksSsdOutSdcIn(); }
  void ResetTracksSdcInSdcOut( void )  { ClearTracksSdcInSdcOut();  }
  void ApplyBh1SegmentCut(const std::vector<double>& validBh1Cluster);
  void ApplyBh2SegmentCut(const double Time0_Cluster);

};

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetTempBcInHC( int layer ) const
{
  if( layer>NumOfLayersBcIn ) layer=0;
  return m_TempBcInHC[layer];
}

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetBcInHC( int layer ) const
{
  if( layer>NumOfLayersBcIn ) layer=0;
  return m_BcInHC[layer];
}

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetBcOutHC( int layer ) const
{
  if( layer>NumOfLayersBcOut ) layer=0;
  return m_BcOutHC[layer];
}

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetSdcInHC( int layer ) const
{
  if( layer>NumOfLayersSdcIn ) layer=0;
  return m_SdcInHC[layer];
}

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetSdcOutHC( int layer ) const
{
  if( layer>NumOfLayersSdcOut ) layer=0;
  return m_SdcOutHC[layer];
}

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetSsdInHC( int layer ) const
{
  if( layer>NumOfLayersSsdIn ) layer=0;
  return m_SsdInHC[layer];
}

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetSsdOutHC( int layer ) const
{
  if( layer>NumOfLayersSsdOut ) layer=0;
  return m_SsdOutHC[layer];
}

//______________________________________________________________________________
inline const DCHitContainer&
DCAnalyzer::GetTOFHC( void ) const
{
  return m_TOFHC;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackBcIn( int i ) const
{
  if( i<m_BcInTC.size() )
    return m_BcInTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackBcOut( int i ) const
{
  if( i<m_BcOutTC.size() )
    return m_BcOutTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcIn( int i ) const
{
  if( i<m_SdcInTC.size() )
    return m_SdcInTC[i];
  // if( i<m_SsdOutSdcInTC.size() )
  //   return m_SsdOutSdcInTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcOut( int i ) const
{
  if( i<m_SdcOutTC.size() )
    return m_SdcOutTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSsdIn( int i ) const
{
  if( i<m_SsdInTC.size() )
    return m_SsdInTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSsdOut( int i ) const
{
  if( i<m_SsdOutTC.size() )
    return m_SsdOutTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSsdInX( int i ) const
{
  if( i<m_SsdInXTC.size() )
    return m_SsdInXTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSsdInY( int i ) const
{
  if( i<m_SsdInYTC.size() )
    return m_SsdInYTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSsdX( int i ) const
{
  if( i<m_SsdXTC.size() )
    return m_SsdXTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSsdY( int i ) const
{
  if( i<m_SsdYTC.size() )
    return m_SsdYTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcInEx( int layer, int i ) const
{
  if( i<m_SdcInExTC[layer].size() )
    return m_SdcInExTC[layer][i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcOutEx( int layer, int i ) const
{
  if( i<m_SdcOutExTC[layer].size() )
    return m_SdcOutExTC[layer][i];
  else
    return 0;
}

//______________________________________________________________________________
inline K18TrackU2D*
DCAnalyzer::GetK18TrackU2D( int i ) const
{
  if( i<m_K18U2DTC.size() )
    return m_K18U2DTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline K18TrackD2U*
DCAnalyzer::GetK18TrackD2U( int i ) const
{
  if( i<m_K18D2UTC.size() )
    return m_K18D2UTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline KuramaTrack*
DCAnalyzer::GetKuramaTrack( int i ) const
{
  if( i<m_KuramaTC.size() )
    return m_KuramaTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackBcOutSdcIn( int i ) const
{
  if( i<m_BcOutSdcInTC.size() )
    return m_BcOutSdcInTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackBcOutSsdIn( int i ) const
{
  if( i<m_BcOutSsdInTC.size() )
    return m_BcOutSsdInTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSsdOutSdcIn( int i ) const
{
  if( i<m_SsdOutSdcInTC.size() )
    return m_SsdOutSdcInTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSdcInSdcOut( int i ) const
{
  if( i<m_SdcInSdcOutTC.size() )
    return m_SdcInSdcOutTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline DCLocalTrack*
DCAnalyzer::GetTrackSsdInSsdOut( int i ) const
{
  if( i<m_SsdInSsdOutTC.size() )
    return m_SsdInSsdOutTC[i];
  else
    return 0;
}

//______________________________________________________________________________
inline const MWPCClusterContainer&
DCAnalyzer::GetClusterMWPC( int layer ) const
{
  if( layer>NumOfLayersBcIn ) layer=0;
  return m_MWPCClCont[layer];
}

//______________________________________________________________________________
inline const SsdClusterContainer&
DCAnalyzer::GetClusterSsdIn( int layer ) const
{
  if( layer>NumOfLayersSsdIn ) layer=0;
  return m_SsdInClCont[layer];
}

//______________________________________________________________________________
inline const SsdClusterContainer&
DCAnalyzer::GetClusterSsdOut( int layer ) const
{
  if( layer>NumOfLayersSsdOut ) layer=0;
  return m_SsdOutClCont[layer];
}

#endif
