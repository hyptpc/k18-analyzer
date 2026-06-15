// -*- C++ -*-
/*
  //comment by wooseung

  This class has two purposes.
  The first purpose is to find the vertex between two tracks, which is used in Lambda reconstruction and the tracking process.
  The second purpose is for clustering in multitrack for accidental concidence evnet veto.
  These two cases are distinguished by using different constructors.
 */

#include "TPCVertex.hh"

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "Kinematics.hh"
#include "ThreeVector.hh"
#include "TPCLocalTrack.hh"
#include "TPCLocalTrackHelix.hh"
#include "UserParamMan.hh"

#include <std_ostream.hh>

namespace
{
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
TPCVertex::TPCVertex(Int_t id1, Int_t id2)
  : m_is_calculated(false), m_is_accidental(),
    m_vertex(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_angle(TMath::QuietNaN()), m_distance(1.e+10),
    m_track_id(), m_track_pid(), m_track_charge(),
    m_scatter_track_flags(kScatterTrackFlagsUnknown),
    m_track_chisqr(), m_track_nhit(), m_track_fit_flag(),
    m_track_pos(), m_track_mom(), m_track_theta()

{
  debug::ObjectCounter::increase(ClassName());
  m_track_id.push_back(id1);
  m_track_id.push_back(id2);
}

//_____________________________________________________________________________
TPCVertex::TPCVertex(TVector3 vertex, std::vector<Int_t> trackid)
  : m_is_calculated(false), m_is_accidental(),
    m_vertex(vertex.x(), vertex.y(), vertex.z()),
    m_angle(TMath::QuietNaN()), m_distance(1.e+10),
    m_track_id(), m_track_pid(), m_track_charge(),
    m_scatter_track_flags(kScatterTrackFlagsUnknown),
    m_track_chisqr(), m_track_nhit(), m_track_fit_flag(),
    m_track_pos(), m_track_mom(), m_track_theta()
{
  debug::ObjectCounter::increase(ClassName());
  m_track_id.reserve(trackid.size());
  for (const auto& id : trackid) {
    m_track_id.push_back(id);
  }
  // This ctor is used for clustered accidental vertices where a representative vertex position is already known.
  m_is_calculated = true;
  m_distance = 0.;
}

//_____________________________________________________________________________
TPCVertex::~TPCVertex()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCVertex::Calculate(TPCLocalTrackHelix* track1, TPCLocalTrackHelix* track2)
{

  // For vertex finding (Closest point betweepn helix tracks)
  // Scanning range (-VertexScanRange, Tracklength + VertexScanRange) for each helix track
  static const Double_t VertexScanRange = gUser.GetParameter("VertexScanRange"); //mm

  Double_t helix_par1[5];
  track1->GetParam(helix_par1);
  Double_t scantheta1 = VertexScanRange/helix_par1[kHelixR]; // mm -> rad.
  Double_t range_theta1[2] = {track1->GetMint() - scantheta1,
			                        track1->GetMaxt() + scantheta1};
  Double_t helix_par2[5];
  track2->GetParam(helix_par2);
  Double_t scantheta2 = VertexScanRange/helix_par2[kHelixR]; // mm -> rad.
  Double_t range_theta2[2] = {track2->GetMint() - scantheta2,
			                        track2->GetMaxt() + scantheta2};

  Double_t theta1, theta2, dist;
  TVector3 vertex = Kinematics::VertexPointHelix(
    helix_par1, helix_par2,
    range_theta1[0], range_theta1[1],
    range_theta2[0], range_theta2[1],
    theta1, theta2, dist
  );

  // check within the TPC volume
  if(!TMath::IsNaN(dist) &&
     TMath::Abs(vertex.x()) < 250. &&
     TMath::Abs(vertex.z()) < 250. &&
     TMath::Abs(vertex.y()) < 250.) m_is_calculated = true;

  if(m_is_calculated){
    m_scatter_track_flags = 0;
    m_vertex = vertex;
    m_track_charge.push_back(track1->GetCharge());
    m_track_charge.push_back(track2->GetCharge());
    m_track_pid.push_back(track1->GetPid());
    m_track_pid.push_back(track2->GetPid());
    if (track1->GetIsK18() != 0) m_scatter_track_flags |= kTrack1IsK18;
    if (track1->GetIsBeam() != 0) m_scatter_track_flags |= kTrack1IsBeam;
    if (track1->GetIsAccidental() != 0) m_scatter_track_flags |= kTrack1IsAccidental;
    if (track2->GetIsK18() != 0) m_scatter_track_flags |= kTrack2IsK18;
    if (track2->GetIsBeam() != 0) m_scatter_track_flags |= kTrack2IsBeam;
    if (track2->GetIsAccidental() != 0) m_scatter_track_flags |= kTrack2IsAccidental;
    m_track_chisqr.push_back(track1->GetChiSquare());
    m_track_chisqr.push_back(track2->GetChiSquare());
    m_track_nhit.push_back(track1->GetNHit());
    m_track_nhit.push_back(track2->GetNHit());
    m_track_fit_flag.push_back(track1->GetFitFlag());
    m_track_fit_flag.push_back(track2->GetFitFlag());
    m_track_pos.push_back(track1->GetPosition(helix_par1, theta1));
    m_track_pos.push_back(track2->GetPosition(helix_par2, theta2));
    TVector3 mom1 = track1->CalcHelixMom(helix_par1, theta1);
    TVector3 mom2 = track2->CalcHelixMom(helix_par2, theta2);
    m_track_mom.push_back(mom1);
    m_track_mom.push_back(mom2);
    m_track_theta.push_back(theta1);
    m_track_theta.push_back(theta2);

    m_angle = m_track_mom[0].Angle(m_track_mom[1]);
    m_distance = dist;
  }

}

//not supported (E42 is not using this)
//_____________________________________________________________________________
void
TPCVertex::Calculate(TPCLocalTrack* track1, TPCLocalTrack* track2)
{

  /*
  // For vertex finding (Closest point betweepn helix tracks)
  // Scanning range (-VertexScanRange, Tracklength + VertexScanRange) for each helix track
  static const Double_t VertexScanRange = gUser.GetParameter("VertexScanRange"); //mm

  Double_t par1[5];
  track1 -> GetParam(par1);
  Double_t par2[5];
  track2 -> GetParam(par2);

  Double_t dist;
  TVector3 vertex = Kinematics::VertexPoint();

  if(!TMath::IsNaN(dist) &&
     TMath::Abs(vertex.x()) < 250. &&
     TMath::Abs(vertex.z()) < 250. &&
     TMath::Abs(vertex.y()) < 250.) m_is_calculated = true;

  if(m_is_calculated){
    m_vertex = vertex;
    m_track_pos.push_back(track1 -> GetPosition(par1, theta1));
    m_track_pos.push_back(track2 -> GetPosition(par2, theta2));
    m_angle = ;
    m_distance = dist;
  }
  */
}

//_____________________________________________________________________________
void
TPCVertex::Print(const TString& arg, Bool_t print_all) const
{
  hddaq::cerr << arg << std::endl
	      << "Vertex point = " << m_vertex << std::endl
	      << "Closest distance = " << m_distance << " mm, "
	      << "Opening angle = " << m_angle << std::endl;
  if(print_all){
    hddaq::cerr << " Track1 #id : "<< m_track_id[0] << std::endl
		<< " closest point " << m_track_pos[0]
		<< " mm, theta "<<m_track_theta[0] << std::endl
		<< " mom "<<m_track_mom[0] << " GeV/c"<<std::endl;
    hddaq::cerr << " Track2 #id : "<< m_track_id[1] << std::endl
		<< " closest point " << m_track_pos[1]
		<< " mm, theta "<<m_track_theta[1] << std::endl
		<< " mom "<<m_track_mom[1] << " GeV/c"<<std::endl;
  }

}
