// -*- C++ -*-

#ifndef TPC_TRACK_SEARCH_HH
#define TPC_TRACK_SEARCH_HH

#include <vector>

#include <TString.h>

#include "DCAnalyzer.hh"
#include <TVector3.h>
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "TPCVertexHelix.hh"

class TPCCluster;
class TPCLocalTrack;
class TPCLocalTrackHelix;
class TPCVertexHelix;

namespace tpc
{

inline const TString& ClassName()
{
  static TString s_name("TPCTrackSearch");
  return s_name;
}

//Functions for the TPCAnalyzer.
//HS-OFF data
Int_t LocalTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<TPCLocalTrack*>& TrackCont,
		       std::vector<TPCLocalTrack*>& TrackContFailed,
		       Bool_t Exclusive,
		       Int_t MinNumOfHits);

//HS-ON data
//w/o K1.8, Kurama tracking information
Int_t LocalTrackSearchHelix(const std::vector<TPCClusterContainer>& ClCont,
			    std::vector<TPCLocalTrackHelix*>& TrackCont,
			    std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			    std::vector<TPCVertexHelix*>& VertexCont,
			    Bool_t Exclusive,
			    Int_t MinNumOfHits);

//for common runs(w/ K1.8, Kurama tracking information)
Int_t LocalTrackSearchHelix(std::vector<std::vector<TVector3>> K18VPs,
			    std::vector<std::vector<TVector3>> KuramaVPs,
			    const std::vector<TPCClusterContainer>& ClCont,
			    std::vector<TPCLocalTrackHelix*>& TrackCont,
			    std::vector<TPCLocalTrackHelix*>& TrackContVP,
			    std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			    std::vector<TPCVertexHelix*>& VertexCont,
			    Bool_t Exclusive,
			    Int_t MinNumOfHits);

//Track fitting
template <typename T>
Bool_t FitStep(T* Track,
	       const std::vector<TPCClusterContainer>& ClCont,
	       std::vector<T*>& TrackContFailed,
	       Int_t MinNumOfHits);

template <typename T>
void FitTrack(T *Track, Int_t Houghflag,
	      const std::vector<TPCClusterContainer>& ClCont,
	      std::vector<T*>& TrackCont,
	      std::vector<T*>& TrackContFailed,
	      Int_t MinNumOfHits);

//HS-OFF data
//Make a track after Hough-transform
Bool_t MakeLinearTrack(TPCLocalTrack *Track, Bool_t &VtxFlag,
		       const std::vector<TPCClusterContainer>& ClCont,
		       Double_t *LinearPar,
		       Double_t MaxHoughWindowY);

//HS-ON data
//Make a track after Hough-transform
Bool_t MakeHelixTrack(TPCLocalTrackHelix *Track, Bool_t &VtxFlag,
		      const std::vector<TPCClusterContainer>& ClCont,
		      Double_t *HelixPar,
		      Double_t MaxHoughWindow,
		      Double_t MaxHoughWindowY);
//Vertex finding
void VertexSearch(std::vector<TPCLocalTrackHelix*>& TrackCont,
		  std::vector<TPCVertexHelix*>& VertexCont);

//Merge fragmented tracks
template <typename T> void
RestoreFragmentedTracks(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<T*>& TrackCont,
			std::vector<T*>& TrackContFailed,
			std::vector<TPCVertexHelix*>& VertexCont,
			Bool_t Exclusive,
			Int_t MinNumOfHits);

//Commom helix track searching
void HelixTrackSearch(Int_t Trackflag, Int_t Houghflag,
		      const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      Int_t MinNumOfHits);

//High momentum helix track searching with linear tracking
void HighMomHelixTrackSearch(const std::vector<TPCClusterContainer>& ClCont,
			     std::vector<TPCLocalTrackHelix*>& TrackCont,
			     std::vector<TPCLocalTrackHelix*>& TrackContFailed,
			     Int_t MinNumOfHits);

//With constraints for K1.8 & Kurama tracks
//Track fitting with constraints
/* Legacy : frequently failed to fit.
  template <typename T>
void FitTrack(T *Track, Int_t Houghflag,
	      const std::vector<TPCClusterContainer>& ClCont,
	      std::vector<T*>& TrackCont,
	      std::vector<T*>& TrackContFailed,
	      Int_t ChargeConstraint, Double_t RKHelixParam[5],
	      Int_t MinNumOfHits);

void HelixTrackSearch(Int_t Trackflag, Int_t Houghflag,
		      const std::vector<TPCClusterContainer>& ClCont,
		      std::vector<TPCLocalTrackHelix*>& TrackCont,
		      std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		      Int_t ChargeConstraint, Double_t RKHelixParam[5],
		      Int_t MinNumOfHits);
*/

//Kurama scattered track finding
void KuramaTrackSearch(std::vector<std::vector<TVector3>> VPs,
		       const std::vector<TPCClusterContainer>& ClCont,
		       std::vector<TPCLocalTrackHelix*>& TrackCont,
		       std::vector<TPCLocalTrackHelix*>& TrackContFailed,
		       std::vector<TPCLocalTrackHelix*>& TrackContVP,
		       std::vector<TPCVertexHelix*>& VertexCont,
		       Bool_t Exclusive,
		       Int_t MinNumOfHits);
//K1.8 beam track finding
void K18TrackSearch(std::vector<std::vector<TVector3>> VPs,
		    const std::vector<TPCClusterContainer>& ClCont,
		    std::vector<TPCLocalTrackHelix*>& TrackCont,
		    std::vector<TPCLocalTrackHelix*>& TrackContVP,
		    Int_t MinNumOfHits=2);

//functions for R&D of tracking
void HoughTransformTest(const std::vector<TPCClusterContainer>& ClCont,
			std::vector<TPCLocalTrack*>& TrackCont,
			Int_t MinNumOfHits);
void HoughTransformTestHelix(const std::vector<TPCClusterContainer>& ClCont,
			     std::vector<TPCLocalTrackHelix*>& TrackCont,
			     Int_t MinNumOfHits);
}
#endif
