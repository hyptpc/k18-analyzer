// -*- C++ -*-

#ifndef TPC_TRACK_SEARCH_HH
#define TPC_TRACK_SEARCH_HH

#include <vector>

#include <TString.h>

#include "DCAnalyzer.hh"
#include <TVector3.h>
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"

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

//Track fitting
template <typename T>
void FitTrack(T *Track, Int_t Houghflag,
	      const std::vector<TPCClusterContainer>& ClCont,
	      std::vector<T*>& TrackCont,
	      std::vector<T*>& TrackContFailed,
	      Int_t MinNumOfHits);

//functions for the HS-OFF data
//for DCAnalyzer
Int_t LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<TPCLocalTrack*>& TrackCont,
		       std::vector<TPCLocalTrack*>& TrackContFailed,
		       bool Exclusive,
		       Int_t MinNumOfHits=8);
//Hough distance calculation
Bool_t MakeLinearTrack(TPCLocalTrack *Track, Bool_t &VtxFlag,
		       const std::vector<TPCClusterContainer>& ClCont,
		       Double_t *LinearPar,
		       Double_t MaxHoughWindowY);

//functions for the HS-ON data
//Track searching
//for DCAnalyzer
Int_t LocalTrackSearchHelix(const std::vector<TPCClusterContainer>& ClCont,
			    std::vector<TPCLocalTrackHelix*>& TrackCont,
			    std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			    bool Exclusive,
			    Int_t MinNumOfHits=8);
Int_t LocalTrackSearchHelix(std::vector<std::vector<TVector3>> K18VPs,
			    std::vector<std::vector<TVector3>> KuramaVPs,
			    const std::vector<TPCClusterContainer>& ClCont,
			    std::vector<TPCLocalTrackHelix*>& TrackCont,
			    std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			    bool Exclusive,
			    Int_t MinNumOfHits=8);
//Commom helix track
void HelixTrackSearch(Int_t Trackflag, Int_t Houghflag,
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
void UseBeamRemover(const std::vector<TPCClusterContainer>& ClCont,
		     std::vector<TPCLocalTrackHelix*>& TrackCont,
		     std::vector<TPCLocalTrackHelix*>& TrackContFailed);


void AccidentalBeamSearchTemp(const std::vector<TPCClusterContainer>& ClCont,
			      std::vector<TPCLocalTrackHelix*>& TrackCont,
			      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			      Int_t MinNumOfHits=10);
//Hough distance calculation
Bool_t MakeHelixTrack(TPCLocalTrackHelix *Track, Bool_t &VtxFlag,
		      const std::vector<TPCClusterContainer>& ClCont,
		      Double_t *HelixPar,
		      Double_t MaxHoughWindow,
		      Double_t MaxHoughWindowY);

//functions for tracking study
void HoughTransformTest(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrack*>& TrackCont,
			Int_t MinNumOfHits /*=8*/);
void HoughTransformTestHelix(const std::vector<TPCClusterContainer>& ClCont,
			     std::vector<TPCLocalTrackHelix*>& TrackCont,
			     Int_t MinNumOfHits /*=8*/);

//undev
void TestRemainingHits(const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<TPCLocalTrackHelix*>& TrackCont,
		       Int_t MinNumOfHits);

}
#endif
