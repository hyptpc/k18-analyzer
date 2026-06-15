// -*- C++ -*-

#ifndef TPC_RECONSTRUCTOR_HH
#define TPC_RECONSTRUCTOR_HH

#include <vector>

#include <Rtypes.h>
#include <TMath.h>
#include <TString.h>
#include <TVector3.h>

class TPCVertex;
class TPCLocalTrackHelix;

//_____________________________________________________________________________
/// One reconstructed decay candidate (Lambda, K0s, etc.) for Dst / histograms.
class TPCRecoCandidate
{
public:
  TPCRecoCandidate() = default;
  TPCRecoCandidate(Int_t mother_pdg, Double_t mass,
                   const TVector3& vertex, const TVector3& momentum,
                   Double_t closest_dist,
                   Int_t track1_id, Int_t track2_id,
                   Int_t track1_pid, Int_t track2_pid,
                   Int_t track1_charge, Int_t track2_charge);

  Int_t GetMotherPdg() const { return m_mother_pdg; }
  Double_t GetMass() const { return m_mass; }
  TVector3 GetVertex() const { return m_vertex; }
  TVector3 GetMomentum() const { return m_momentum; }
  Double_t GetClosestDist() const { return m_closest_dist; }
  Int_t GetTrackId1() const { return m_track1_id; }
  Int_t GetTrackId2() const { return m_track2_id; }
  Int_t GetTrack1Pid() const { return m_track1_pid; }
  Int_t GetTrack2Pid() const { return m_track2_pid; }
  Int_t GetTrack1Charge() const { return m_track1_charge; }
  Int_t GetTrack2Charge() const { return m_track2_charge; }
  void Print(const TString& label="") const;

private:
  Int_t m_mother_pdg{0};
  Double_t m_mass{TMath::QuietNaN()};
  TVector3 m_vertex;
  TVector3 m_momentum;
  Double_t m_closest_dist{TMath::QuietNaN()};
  Int_t m_track1_id{-1};
  Int_t m_track2_id{-1};
  Int_t m_track1_pid{0};
  Int_t m_track2_pid{0};
  Int_t m_track1_charge{0};
  Int_t m_track2_charge{0};
};

//_____________________________________________________________________________
/// Decay reconstruction from vertex geometry + PID (writes candidates into TPCVertex).
class TPCReconstructor
{
public:
  enum RecoMode : UInt_t {
    kRecoNone    = 0u,
    kRecoLambda  = 1u << 0,
    kRecoK0Short = 1u << 1,
    // kRecoParticle = 1u << n, <- you can add like this
    kRecoAll     = kRecoLambda | kRecoK0Short
  };

  static const TString& ClassName();

  TPCReconstructor();
  ~TPCReconstructor();

  /// Kinematic Lambda check without storing candidates (e.g. temporary vertex in track search).
  /// After ReconstructLambda, use TPCVertex::GetHasLambdaCandidate() instead.
  static Bool_t HasLambdaCandidate(const TPCVertex* vertex);
  static Bool_t HasK0ShortCandidate(const TPCVertex* vertex);

  void Clear();
  void Calculate(const std::vector<TPCVertex*>& vertices,
                 const std::vector<TPCLocalTrackHelix*>* helix_tracks = nullptr,
                 UInt_t reco_mode = kRecoAll);

private:
  Bool_t RefitLambdaWithVertex(const TPCVertex* vertex,
                     const std::vector<TPCLocalTrackHelix*>* helix_tracks,
                     TVector3& out_vertex,
                     TVector3& out_mom1,
                     TVector3& out_mom2,
                     Double_t& out_dist) const;
  void ReconstructLambda(TPCVertex* vertex,
                          const std::vector<TPCLocalTrackHelix*>* helix_tracks);
  void ReconstructK0Short(TPCVertex* vertex);
};

//_____________________________________________________________________________
inline const TString&
TPCReconstructor::ClassName()
{
  static TString s_name("TPCReconstructor");
  return s_name;
}

#endif
