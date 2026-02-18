// -*- C++ -*-

#ifndef DC_TRACK_SEARCH_HH
#define DC_TRACK_SEARCH_HH

#include <vector>

#include <TString.h>

#include "DCAnalyzer.hh"

struct DCPairPlaneInfo;
class DCLocalTrack;
class DCLTrackHit;
class DCPairHitCluster;

using ClusterList = std::vector<DCPairHitCluster*>;
using IndexList = std::vector<Int_t>;

namespace track
{
inline const TString& ClassName()
{
  static TString s_name("DCTrackSearch");
  return s_name;
}

std::vector<IndexList> MakeIndex(Int_t ndim, const Int_t *index1, Bool_t& status);
std::vector<IndexList> MakeIndex(Int_t ndim, const IndexList& index1, Bool_t& status);
std::vector<IndexList> MakeIndex_VXU(Int_t ndim, Int_t maximumHit, const Int_t *index1);
std::vector<IndexList> MakeIndex_VXU(Int_t ndim, Int_t maximumHit, const IndexList& index1);
DCLocalTrack*          MakeTrack(const std::vector<ClusterList>& CandCont,
                                 const IndexList& combination);

Int_t LocalTrackSearch(const std::vector<DCHC>& HC,
                       const DCPairPlaneInfo *PpInfo,
                       Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
                       Int_t MinNumOfHits=6, Int_t T0Seg=-1);
Int_t LocalTrackSearchVUX(const std::vector<DCHC>& HC,
                          const DCPairPlaneInfo *PpInfo,
                          Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
                          Int_t MinNumOfHits=6, Int_t T0Seg=-1);
} // namespace track

#endif
