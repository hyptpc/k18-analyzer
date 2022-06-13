// -*- C++ -*-

#ifndef TPC_TRACK_SEARCH_HH
#define TPC_TRACK_SEARCH_HH

#include <vector>

#include <TString.h>

#include "DCAnalyzer.hh"

// struct DCPairPlaneInfo;
// class DCPairHitCluster;
// class DCLocalTrack;
// class DCLTrackHit;
// class MWPCCluster;
class TPCCluster;
class TPCLocalTrack;
class TPCLocalTrackHelix;

namespace tpc
{
inline const TString& ClassName()
{
  static TString s_name("TPCTrackSearch");
  return s_name;
}

Int_t LocalTrackSearch(const std::vector<TPCHitContainer>& HitCont,
                       std::vector<TPCLocalTrack*>& TrackCont,
                       Int_t MinNumOfHits=8);
Int_t LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
                       std::vector<TPCLocalTrack*>& TrackCont,
                       Int_t MinNumOfHits=8);
Int_t LocalTrackSearchHelix(const std::vector<TPCHitContainer>& HitCont,
                            std::vector<TPCLocalTrackHelix*>& TrackCont,
                            Int_t MinNumOfHits=8);
Int_t LocalTrackSearchHelix(const std::vector<TPCClusterContainer>& ClCont,
                            std::vector<TPCLocalTrackHelix*>& TrackCont,
                            Int_t MinNumOfHits=8);
}

#endif
