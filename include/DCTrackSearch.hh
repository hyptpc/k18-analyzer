// -*- C++ -*-

#ifndef DC_TRACK_SEARCH_HH
#define DC_TRACK_SEARCH_HH

#include "DCAnalyzer.hh"

#include <vector>

#include <TString.h>

struct DCPairPlaneInfo;
class DCPairHitCluster;
class DCLocalTrack;
class DCLTrackHit;
class MWPCCluster;

using ClusterList = std::vector<DCPairHitCluster*>;
using IndexList = std::vector<Int_t>;

namespace track
{
inline const TString& ClassName()
{
  static TString s_name("DCTrackSearch");
  return s_name;
}

std::vector<IndexList> MakeIndex(Int_t ndim, const Int_t *index1, bool& status);
std::vector<IndexList> MakeIndex(Int_t ndim, const IndexList& index1, bool& status);
std::vector<IndexList> MakeIndex_VXU(Int_t ndim, Int_t maximumHit, const Int_t *index1);
std::vector<IndexList> MakeIndex_VXU(Int_t ndim, Int_t maximumHit, const IndexList& index1);
DCLocalTrack*          MakeTrack(const std::vector<ClusterList>& CandCont,
                                 const IndexList& combination);

Int_t LocalTrackSearch(const std::vector<DCHC>& HC,
                       const DCPairPlaneInfo *PpInfo,
                       Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
                       Int_t MinNumOfHits=6, Int_t T0Seg = -1);
Int_t LocalTrackSearch(const std::vector< std::vector<DCHC> >& hcAssemble,
                       const DCPairPlaneInfo *PpInfo,
                       Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
                       Int_t MinNumOfHits=6, Int_t T0Seg = -1);
Int_t LocalTrackSearchSdcInFiber(const std::vector<DCHC>& HC,
                                 const DCPairPlaneInfo *PpInfo,
                                 Int_t npp, std::vector<DCLocalTrack*>& trackCont,
                                 Int_t MinNumOfHits=6);
Int_t LocalTrackSearchVUX(const std::vector<DCHC>& HC,
                          const DCPairPlaneInfo *PpInfo,
                          Int_t npp, std::vector<DCLocalTrack*>& TrackCont,
                          Int_t MinNumOfHits=6);
Int_t LocalTrackSearchSdcOut(const std::vector<DCHC>& SdcOutHC,
                             const DCPairPlaneInfo *PpInfo, Int_t npp,
                             std::vector<DCLocalTrack*>& TrackCont,
                             Int_t MinNumOfHits=6);
Int_t LocalTrackSearchSdcOut(const DCHC& TOFHC,
                             const std::vector<DCHC>& SdcOutHC,
                             const DCPairPlaneInfo *PpInfo,
                             Int_t npp,
                             std::vector<DCLocalTrack*>& TrackCont,
                             Int_t MinNumOfHits=6);
Int_t MakeLocalTrackGeant4(const std::vector<DCHC>& HC,
			   std::vector<DCLocalTrack*>& TrackCont,
			   Int_t MinNumOfHits=6);
Int_t LocalTrackSearchBcOutSdcIn(const std::vector<DCHC>& BcHC,
                                 const DCPairPlaneInfo *BcPpInfo,
                                 const std::vector<DCHC>& SdcHC,
                                 const DCPairPlaneInfo *SdcPpInfo,
                                 Int_t BcNpp, Int_t SdcNpp,
                                 std::vector<DCLocalTrack*>& TrackCont,
                                 Int_t MinNumOfHits=18);
Int_t LocalTrackSearchSdcInSdcOut(const std::vector<DCHC>& SdcInHC,
                                  const DCPairPlaneInfo *SdcInPpInfo,
                                  const std::vector<DCHC>& SdcOutHC,
                                  const DCPairPlaneInfo *SdcOutPpInfo,
                                  Int_t SdcInNpp, Int_t SdcOutNpp,
                                  std::vector<DCLocalTrack*>& TrackCont,
                                  Int_t MinNumOfHits=12);
Int_t MWPCLocalTrackSearch(const std::vector<DCHC>& HC,
                           std::vector<DCLocalTrack*>& TrackCont);
Int_t MWPCLocalTrackSearch(const std::vector< std::vector<DCHC> >& hcList,
                           std::vector<DCLocalTrack*>& trackCont);
Int_t LocalTrackSearchCFT(const std::vector<DCHC>& HC,
                          const DCPairPlaneInfo *PpInfo,
                          Int_t npp, std::vector<DCLocalTrack*>& trackCont,
                          Int_t MinNumOfHits=3);
Int_t LocalTrackSearchCFTppPhi(const std::vector<DCHC>& HC,
                               const DCPairPlaneInfo *PpInfo,
                               Int_t npp, std::vector<DCLocalTrack*>& trackCont,
                               Int_t MinNumOfHits=3);
}

#endif
