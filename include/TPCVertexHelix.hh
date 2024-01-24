// -*- C++ -*-

#ifndef TPC_VERTEX_HELIX_HH
#define TPC_VERTEX_HELIX_HH

#include <vector>

#include "TPCLocalTrackHelix.hh"
#include <TString.h>

class TPCLocalTrackHelix;

//_____________________________________________________________________________
class TPCVertexHelix
{
public:
  static const TString& ClassName();
  TPCVertexHelix(Int_t id1, Int_t id2);
  ~TPCVertexHelix();

private:
  TVector3 m_vertex;
  Int_t  m_track1_id;
  Int_t  m_track2_id;
  TVector3 m_track1_pos;
  TVector3 m_track2_pos;
  TVector3 m_track1_mom;
  TVector3 m_track2_mom;
  Double_t m_track1_theta;
  Double_t m_track2_theta;
  Double_t m_angle;
  Double_t m_distance;

public:

  void Calculate(TPCLocalTrackHelix* track1,
		 TPCLocalTrackHelix* track2);

  void SetTrack1Id(Int_t id){ m_track1_id = id; }
  void SetTrack2Id(Int_t id){ m_track2_id = id; }
  void SetTrack1Pos(TVector3 vec) { m_track1_pos = vec; }
  void SetTrack2Pos(TVector3 vec) { m_track2_pos = vec; }
  void SetTrack1Mom(TVector3 vec) { m_track1_mom = vec; }
  void SetTrack2Mom(TVector3 vec) { m_track2_mom = vec; }
  void SetTrack1Theta(Double_t theta) { m_track1_theta = theta; }
  void SetTrack2Theta(Double_t theta) { m_track2_theta = theta; }
  void SetOpeningAngle(Double_t angle) { m_angle = angle; }
  void SetClosestDist(Double_t dist) { m_distance = dist; }

  Int_t GetTrack1Id() const { return m_track1_id; }
  Int_t GetTrack2Id() const { return m_track2_id; }
  TVector3 GetTrack1Pos() const { return m_track1_pos; }
  TVector3 GetTrack2Pos() const { return m_track2_pos; }
  TVector3 GetTrack1Mom() const { return m_track1_mom; }
  TVector3 GetTrack2Mom() const { return m_track2_mom; }
  Double_t GetTrack1Theta() const { return m_track1_theta; }
  Double_t GetTrack2Theta() const { return m_track2_theta; }

  TVector3 GetVertex() const { return m_vertex; }
  Double_t GetClosestDist() const { return m_distance; }
  Double_t GetOpeningAngle() const { return m_angle; }

  void   Print(const TString& arg="", Bool_t print_all=false) const;
};

//_____________________________________________________________________________
inline const TString&
TPCVertexHelix::ClassName()
{
  static TString s_name("TPCVertexHelix");
  return s_name;
}

#endif
