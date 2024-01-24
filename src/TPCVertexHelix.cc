// -*- C++ -*-

#include "TPCVertexHelix.hh"

#include <iostream>
#include <iterator>

#include "std_ostream.hh"

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "Kinematics.hh"
#include "UserParamMan.hh"

namespace
{
const auto& gUser = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
TPCVertexHelix::TPCVertexHelix(Int_t id1, Int_t id2)
  : m_vertex(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_track1_id(id1),
    m_track2_id(id2),
    m_track1_pos(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_track2_pos(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_track1_mom(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_track2_mom(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_track1_theta(TMath::QuietNaN()),
    m_track2_theta(TMath::QuietNaN()),
    m_angle(TMath::QuietNaN()),
    m_distance(TMath::QuietNaN())
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCVertexHelix::~TPCVertexHelix()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCVertexHelix::Calculate(TPCLocalTrackHelix* track1, TPCLocalTrackHelix* track2)
{

  // For vertex finding (Closest point between helix tracks)
  // Scanning range (-VertexScanRange, Tracklength + VertexScanRange) for each helix track
  static const Double_t VertexScanRange = gUser.GetParameter("VertexScanRange"); //mm

  Double_t helix_par1[5];
  track1 -> GetParam(helix_par1);
  Double_t scantheta1 = VertexScanRange/helix_par1[3]; //mm -> rad.
  Double_t range_theta1[2] = {track1 -> GetMint() - scantheta1,
			      track1 -> GetMaxt() + scantheta1};
  Double_t helix_par2[5];
  track2 -> GetParam(helix_par2);
  Double_t scantheta2 = VertexScanRange/helix_par2[3]; //mm -> rad.
  Double_t range_theta2[2] = {track2 -> GetMint() - scantheta2,
			      track2 -> GetMaxt() + scantheta2};

  Double_t theta1, theta2, dist;
  m_vertex = Kinematics::VertexPointHelix(helix_par1, helix_par2,
					  range_theta1[0], range_theta1[1],
					  range_theta2[0], range_theta2[1],
					  theta1, theta2, dist);

  m_track1_pos = track1 -> GetPosition(helix_par1, theta1);
  m_track1_mom = track1 -> CalcHelixMom(helix_par1, theta1);
  m_track1_theta = theta1;

  m_track2_pos = track2 -> GetPosition(helix_par2, theta2);
  m_track2_mom = track2 -> CalcHelixMom(helix_par2, theta2);
  m_track2_theta = theta2;

  m_angle = m_track1_mom.Angle(m_track2_mom);
  m_distance = dist;

}

//_____________________________________________________________________________
void
TPCVertexHelix::Print(const TString& arg, Bool_t print_all) const
{
  hddaq::cerr << arg << std::endl
	      << "Vertex point = " << m_vertex << std::endl
	      << "Closest distance = " << m_distance << " mm, "
	      << "Opening angle = " << m_angle << std::endl;
  if(print_all){
    hddaq::cerr << " Track1 #id : "<< m_track1_id << std::endl
		<< " closest point " << m_track1_pos
		<< " mm, theta "<<m_track1_theta << std::endl
		<< " mom "<<m_track1_mom << " GeV/c"<<std::endl;
    hddaq::cerr << " Track2 #id : "<< m_track2_id << std::endl
		<< " closest point " << m_track2_pos
		<< " mm, theta "<<m_track2_theta << std::endl
		<< " mom "<<m_track2_mom << " GeV/c"<<std::endl;
  }

}
