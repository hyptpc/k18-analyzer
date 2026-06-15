// -*- C++ -*-

#ifndef TPC_VERTEX_HH
#define TPC_VERTEX_HH

#include <vector>

#include <TMath.h>
#include <TString.h>
#include <TVector3.h>

#include "TPCReconstructor.hh"

class TPCLocalTrack;
class TPCLocalTrackHelix;

//_____________________________________________________________________________
class TPCVertex
{
public:
  static const TString& ClassName();
  TPCVertex(Int_t id1, Int_t id2);
  TPCVertex(TVector3 vertex, std::vector<Int_t> trackid);
  ~TPCVertex();

private:

  Bool_t m_is_calculated; // flag of Calculate() m_vertex;
  Int_t m_is_accidental;
  TVector3 m_vertex;
  Double_t m_angle;
  Double_t m_distance;
  std::vector<Int_t>  m_track_id;
  std::vector<Int_t>  m_track_pid;
  std::vector<Int_t>  m_track_charge;
  UInt_t m_scatter_track_flags;
  std::vector<Double_t> m_track_chisqr;
  std::vector<Int_t> m_track_nhit;
  std::vector<Int_t> m_track_fit_flag;
  std::vector<TVector3> m_track_pos;
  std::vector<TVector3> m_track_mom;
  std::vector<Double_t> m_track_theta;
  std::vector<TPCRecoCandidate> m_reco_candidates;
  Bool_t m_has_lambda_candidate{false};

public:
  // Pair-level reconstruction flags packed in one word.
  // Bits 3-7 are intentionally reserved for future per-pair extensions.
  enum e_scatter_flag : UInt_t {
    kTrack1IsK18        = 1u << 0,
    kTrack1IsBeam       = 1u << 1,
    kTrack1IsAccidental = 1u << 2,
    kTrack2IsK18        = 1u << 8,
    kTrack2IsBeam       = 1u << 9,
    kTrack2IsAccidental = 1u << 10,
    kScatterTrackFlagsUnknown  = 1u << 31
  };

  void Calculate(TPCLocalTrackHelix* track1,
		 TPCLocalTrackHelix* track2);
  void Calculate(TPCLocalTrack* track1,
		 TPCLocalTrack* track2);
  Bool_t IsCalculated() const { return m_is_calculated; }
  Int_t GetIsAccidental() const { return m_is_accidental; }
  void SetIsAccidental(Int_t flag=1) { m_is_accidental=flag; }

  TVector3 GetVertex() const { return m_vertex; }
  Double_t GetOpeningAngle() const { return m_angle; }
  Double_t GetClosestDist() const { return m_distance; }

  Int_t GetTrackId(Int_t i) const { return m_track_id.at(i); }
  Int_t GetTrackCharge(Int_t i) const { return m_track_charge.at(i); }
  Int_t GetTrackPid(Int_t i) const { return m_track_pid.at(i); }
  UInt_t GetScatterTrackFlags() const { return m_scatter_track_flags; }
  Double_t GetTrackChiSquare(Int_t i) const { return m_track_chisqr.at(i); }
  Int_t GetTrackNHit(Int_t i) const { return m_track_nhit.at(i); }
  Int_t GetTrackFitFlag(Int_t i) const { return m_track_fit_flag.at(i); }
  TVector3 GetTrackPos(Int_t i) const { return m_track_pos.at(i); }
  TVector3 GetTrackMom(Int_t i) const { return m_track_mom.at(i); }
  Double_t GetTrackTheta(Int_t i) const { return m_track_theta.at(i); }

  //For clustered tracks
  Int_t GetNTracks() const { return m_track_id.size(); }

  // For decay candidate reconstruction.
  void ClearRecoCandidates()
  {
    m_reco_candidates.clear();
    m_has_lambda_candidate = false;
  }
  void AddRecoCandidate(const TPCRecoCandidate& cand) { m_reco_candidates.push_back(cand); }
  Bool_t GetHasLambdaCandidate() const { return m_has_lambda_candidate; }
  void SetHasLambdaCandidate(Bool_t v) { m_has_lambda_candidate = v; }
  Int_t GetNRecoCandidates() const { return static_cast<Int_t>(m_reco_candidates.size()); }
  const TPCRecoCandidate& GetRecoCandidate(Int_t i) const { return m_reco_candidates.at(i); }

  void Print(const TString& arg="", Bool_t print_all=false) const;
};

//_____________________________________________________________________________
inline const TString&
TPCVertex::ClassName()
{
  static TString s_name("TPCVertex");
  return s_name;
}

#endif
