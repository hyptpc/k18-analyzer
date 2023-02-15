// -*- C++ -*-

#ifndef TPC_TRACK_SEARCH_HH
#define TPC_TRACK_SEARCH_HH

#include <vector>

#include <TString.h>

#include "DCAnalyzer.hh"
#include <TVector3.h>

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

//functions for the HS-OFF data
//Track searching
Int_t LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<TPCLocalTrack*>& TrackCont,
		       std::vector<TPCLocalTrack*>& TrackContFailed,
		       Int_t MinNumOfHits=8);

//Hough distance calculation
Bool_t LinearHoughDistCheck(TPCLocalTrack *track,
			    const std::vector<TPCClusterContainer>& ClCont,
			    Double_t *LinearPar,
			    Double_t MaxHoughWindowY);

//Track fitting
void FitLinearTrack(TPCLocalTrack *Track,
		    const std::vector<TPCClusterContainer>& ClCont,
		    std::vector<TPCLocalTrack*>& TrackCont,
		    std::vector<TPCLocalTrack*>& TrackContFailed,
		    Int_t MinNumOfHits);
void ExclusiveTracking(std::vector<TPCLocalTrack*>& TrackCont);

//functions for the HS-ON data
//Track searching
Int_t LocalTrackSearchHelix(std::vector<std::vector<TVector3>> K18VPs,
			    std::vector<std::vector<TVector3>> KuramaVPs,
			    const std::vector<TPCClusterContainer>& ClCont,
			    std::vector<TPCLocalTrackHelix*>& TrackCont,
			    std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			    Int_t MinNumOfHits=8);
//commom
void HelixTrackSearch(Int_t Beamflag, Int_t Houghflag,
		      const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      Int_t MinNumOfHits=8);

//Kurama tracks
void KuramaTrackSearch(std::vector<std::vector<TVector3>> VPs,
		       const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<TPCLocalTrackHelix*>& TrackCont,
		       std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		       Int_t MinNumOfHits=8);
//K1.8 tracks
void K18TrackSearch(std::vector<std::vector<TVector3>> VPs,
		    const std::vector<TPCClusterContainer>& ClCont,
		    std::vector<TPCLocalTrackHelix*>& TrackCont,
		    std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		    Int_t MinNumOfHits=8);
//Accidental beam tracks
void AccidentalBeamSearch(const std::vector<TPCClusterContainer>& ClCont,
			  std::vector<TPCLocalTrackHelix*>& TrackCont,
			  std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			  Int_t MinNumOfHits=10);

//Hough distance calculation
Bool_t HelixHoughDistCheck(TPCLocalTrackHelix *track,
			   const std::vector<TPCClusterContainer>& ClCont,
			   Double_t *HelixPar,
			   Double_t MaxHoughWindow,
			   Double_t MaxHoughWindowY);

//Track fitting
void FitHelixTrack(TPCLocalTrackHelix *Track,
		   Int_t Houghflag,
		   const std::vector<TPCClusterContainer>& ClCont,
		   std::vector<TPCLocalTrackHelix*>& TrackCont,
		   std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		   Int_t MinNumOfHits);
void ExclusiveTrackingHelix(std::vector<TPCLocalTrackHelix*>& TrackCont);

//undev
void TestRemainingHits(const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<TPCLocalTrackHelix*>& TrackCont,
		       Int_t MinNumOfHits);

}
#endif
