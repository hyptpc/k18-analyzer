/**
 *  file: DCTrackSearch.hh
 *  date: 2017.04.10
 *
 */

#ifndef DC_TRACK_SEARCH_HH
#define DC_TRACK_SEARCH_HH

#include "DCAnalyzer.hh"

#include <vector>

struct DCPairPlaneInfo;
class  DCPairHitCluster;
class  DCLocalTrack;
class  DCLTrackHit;
class  MWPCCluster;

typedef std::vector<DCPairHitCluster*> ClusterList;
typedef std::vector<int>               IndexList;

namespace track
{
  //______________________________________________________________________________
  int LocalTrackSearch( const std::vector<DCHitContainer>& HC,
			const DCPairPlaneInfo *PpInfo,
			int npp, std::vector<DCLocalTrack*>& TrackCont,
			int MinNumOfHits=6 );

  //______________________________________________________________________________
  int LocalTrackSearch( const std::vector< std::vector<DCHitContainer> >& hcAssemble,
			const DCPairPlaneInfo *PpInfo,
			int npp, std::vector<DCLocalTrack*>& TrackCont,
			int MinNumOfHits=6 );

  //______________________________________________________________________________
  int LocalTrackSearchSdcInFiber( const std::vector<DCHitContainer>& HC,
				  const DCPairPlaneInfo *PpInfo,
				  int npp, std::vector<DCLocalTrack*>& trackCont,
				  int MinNumOfHits=6 );

  //______________________________________________________________________________
  int LocalTrackSearchVUX( const std::vector<DCHitContainer>& HC,
			   const DCPairPlaneInfo *PpInfo,
			   int npp, std::vector<DCLocalTrack*>& TrackCont,
			   int MinNumOfHits=6 );

  //______________________________________________________________________________
  int LocalTrackSearchSdcOut( const std::vector<DCHitContainer>& SdcOutHC,
			      const DCPairPlaneInfo *PpInfo, int npp,
			      std::vector<DCLocalTrack*>& TrackCont,
			      int MinNumOfHits=6 );

  //______________________________________________________________________________
  int LocalTrackSearchSdcOut( const DCHitContainer& TOFHC,
			      const std::vector<DCHitContainer>& SdcOutHC,
			      const DCPairPlaneInfo *PpInfo,
			      int npp,
			      std::vector<DCLocalTrack*>& TrackCont,
			      int MinNumOfHits=6 );

  //______________________________________________________________________________
  int LocalTrackSearchSsdIn( const std::vector<DCHitContainer>& HC,
			     std::vector<DCLocalTrack*>& TrackCont,
			     int MinNumOfHits=4 );

  //______________________________________________________________________________
  int LocalTrackSearchSsdIn( const std::vector<SsdClusterContainer>& SsdInClCont,
			     std::vector<DCLocalTrack*>& TrackCont,
			     int MinNumOfHits=4 );

  //______________________________________________________________________________
  int LocalTrackSearchSsdOut( const std::vector<DCHitContainer>& HC,
			      std::vector<DCLocalTrack*>& TrackCont,
			      int MinNumOfHits=4 );

  //______________________________________________________________________________
  int LocalTrackSearchSsdOut( const std::vector<SsdClusterContainer>& SsdOutClCont,
			      std::vector<DCLocalTrack*>& TrackCont,
			      int MinNumOfHits=4 );

  //______________________________________________________________________________
  int LocalTrackSearchSsdInXY( const std::vector<DCHitContainer>& HC,
			       std::vector<DCLocalTrack*>& TrackContX,
			       std::vector<DCLocalTrack*>& TrackContY );

  //______________________________________________________________________________
  int LocalTrackSearchSsdInXY( const std::vector<SsdClusterContainer>& SsdInClCont,
			       std::vector<DCLocalTrack*>& TrackContX,
			       std::vector<DCLocalTrack*>& TrackContY );

  //______________________________________________________________________________
  int PreTrackSearchSsdXY( const std::vector<SsdClusterContainer>& SsdInClCont,
			   const std::vector<SsdClusterContainer>& SsdOutClCont,
			   std::vector<DCLocalTrack*>& TrackContX,
			   std::vector<DCLocalTrack*>& TrackContY );

  //______________________________________________________________________________
  int LocalTrackSearchBcOutSdcIn( const std::vector<DCHitContainer>& BcHC,
				  const DCPairPlaneInfo *BcPpInfo,
				  const std::vector<DCHitContainer>& SdcHC,
				  const DCPairPlaneInfo *SdcPpInfo,
				  int BcNpp, int SdcNpp,
				  std::vector<DCLocalTrack*>& TrackCont,
				  int MinNumOfHits=18 );

  //______________________________________________________________________________
  int LocalTrackSearchBcOutSsdIn( const std::vector<DCHitContainer>& BcHC,
				  const DCPairPlaneInfo *BcPpInfo,
				  int BcNpp,
				  const std::vector<DCHitContainer>& SsdHC,
				  std::vector<DCLocalTrack*>& TrackCont,
				  int MinNumOfHits=14 );

  //______________________________________________________________________________
  int LocalTrackSearchSsdInSsdOut( const std::vector<DCHitContainer>& SsdInHC,
				   const std::vector<DCHitContainer>& SsdOutHC,
				   std::vector<DCLocalTrack*>& TrackCont,
				   int MinNumOfHits=6 );

  //______________________________________________________________________________
  int LocalTrackSearchSsdOutSdcIn( const std::vector<DCHitContainer>& SsdInHC,
				   const std::vector<DCHitContainer>& SsdOutHC,
				   const std::vector<DCHitContainer>& SdcHC,
				   const DCPairPlaneInfo *SdcPpInfo,
				   int SdcNpp,
				   std::vector<DCLocalTrack*>& TrackCont,
				   int MinNumOfHits=10 );

  //______________________________________________________________________________
  int LocalTrackSearchSsdOutSdcIn( const std::vector<SsdClusterContainer>& SsdInClCont,
				   const std::vector<SsdClusterContainer>& SsdOutClCont,
				   const std::vector<DCHitContainer>& SdcHC,
				   const DCPairPlaneInfo *SdcPpInfo,
				   int SdcNpp,
				   std::vector<DCLocalTrack*>& TrackCont,
				   int MinNumOfHits=10,
				   // Delete layers having many hits
				   bool DeleteFlag=false );

  //______________________________________________________________________________
  // w/SsdPreTracking
  int LocalTrackSearchSsdOutSdcIn( const std::vector<DCHitContainer>& SdcInHC,
				   const DCPairPlaneInfo *PpInfo,
				   int SdcInNpp,
				   std::vector<DCLocalTrack*>& SsdXTC,
				   std::vector<DCLocalTrack*>& SsdYTC,
				   std::vector<DCLocalTrack*>& TrackCont,
				   int MinNumOfHits=10 );

  //______________________________________________________________________________
  int LocalTrackSearchSdcInSdcOut( const std::vector<DCHitContainer>& SdcInHC,
				   const DCPairPlaneInfo *SdcInPpInfo,
				   const std::vector<DCHitContainer>& SdcOutHC,
				   const DCPairPlaneInfo *SdcOutPpInfo,
				   int SdcInNpp, int SdcOutNpp,
				   std::vector<DCLocalTrack*>& TrackCont,
				   int MinNumOfHits=12 );
  //______________________________________________________________________________
  int MWPCLocalTrackSearch( const std::vector<DCHitContainer>& HC,
			    std::vector<DCLocalTrack*>& TrackCont );

  //______________________________________________________________________________
  int MWPCLocalTrackSearch( const std::vector< std::vector<DCHitContainer> >& hcList,
			    std::vector<DCLocalTrack*>& trackCont );

}

#endif
