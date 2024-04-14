/**
 *  file: DCTrackTools.hh
 *  date: 2017.04.10
 *
 */

#ifndef DC_TRACK_TOOLS_HH
#define DC_TRACK_TOOLS_HH

#include "DCAnalyzer.hh"

#include <vector>

class  DCPairHitCluster;
class  LocalTrack;

typedef std::vector<DCPairHitCluster*> ClusterList;
typedef std::vector<int>               IndexList;

namespace dctrack
{
  //______________________________________________________________________________
  int LocalTrackSearch( const std::vector<DCHitContainer>& HC,
			std::vector<LocalTrack*>& TrackCont,
			int MinNumOfHits=6 );
}

#endif
