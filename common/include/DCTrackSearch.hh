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
  int LocalTrackSearchVUX( const std::vector<DCHitContainer>& HC,
			   const DCPairPlaneInfo *PpInfo,
			   int npp, std::vector<DCLocalTrack*>& TrackCont,
			   int MinNumOfHits=6 );
  //______________________________________________________________________________
  int MWPCLocalTrackSearch( const std::vector<DCHitContainer>& HC,
			    std::vector<DCLocalTrack*>& TrackCont );

  //______________________________________________________________________________
  int MWPCLocalTrackSearch( const std::vector< std::vector<DCHitContainer> >& hcList,
			    std::vector<DCLocalTrack*>& trackCont );

}

#endif
