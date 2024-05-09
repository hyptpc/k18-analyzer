// -*- C++ -*-

#include "TPCVertexHelix.hh"

#include <iostream>
#include <iterator>

#include "std_ostream.hh"

#include <TLorentzVector.h>
#include "DebugCounter.hh"
#include "FuncName.hh"
#include "Kinematics.hh"
#include "DatabasePDG.hh"
#include "UserParamMan.hh"

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const Double_t ppi_distcut = 10.; //Closest distance for p, pi at the vertex point
}

//_____________________________________________________________________________
TPCVertexHelix::TPCVertexHelix(Int_t id1, Int_t id2)
  : m_is_calculated(false),
    m_vertex(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_angle(TMath::QuietNaN()), m_distance(1.e+10),
    m_track_id(), m_track_pid(), m_track_charge(),
    m_track_pos(), m_track_mom(), m_track_theta(),
    m_is_lambda(false),
    m_lambda_vertex(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN()),
    m_lambda_mass(), m_lambda_mom(),
    m_proton_id(), m_proton_mom(),
    m_pion_id(), m_pion_mom()

{
  debug::ObjectCounter::increase(ClassName());
  m_track_id.push_back(id1);
  m_track_id.push_back(id2);


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

  // For vertex finding (Closest point betweepn helix tracks)
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
  TVector3 vertex = Kinematics::VertexPointHelix(helix_par1, helix_par2,
					  range_theta1[0], range_theta1[1],
					  range_theta2[0], range_theta2[1],
					  theta1, theta2, dist);
  if(!TMath::IsNaN(dist) &&
     TMath::Abs(vertex.x()) < 250. &&
     TMath::Abs(vertex.z()) < 250. &&
     TMath::Abs(vertex.y()) < 250.) m_is_calculated = true;

  if(m_is_calculated){
    m_vertex = vertex;
    m_track_charge.push_back(track1 -> GetCharge());
    m_track_charge.push_back(track2 -> GetCharge());
    m_track_pid.push_back(track1 -> GetPid());
    m_track_pid.push_back(track2 -> GetPid());
    m_track_pos.push_back(track1 -> GetPosition(helix_par1, theta1));
    m_track_pos.push_back(track2 -> GetPosition(helix_par2, theta2));
    m_track_mom.push_back(track1 -> CalcHelixMom(helix_par1, theta1));
    m_track_mom.push_back(track2 -> CalcHelixMom(helix_par2, theta2));
    m_track_theta.push_back(theta1);
    m_track_theta.push_back(theta2);

    m_angle = m_track_mom[0].Angle(m_track_mom[1]);
    m_distance = dist;

    if(track1 -> GetIsK18()==0 && track2 -> GetIsK18()==0 &&
       track1 -> GetIsBeam()==0 && track2 -> GetIsBeam()==0 &&
       track1 -> GetIsKurama()==0 && track2 -> GetIsKurama()==0 &&
       track1 -> GetIsAccidental()==0 && track2 -> GetIsAccidental()==0) IsLambda();
  }
}

//_____________________________________________________________________________
Bool_t
TPCVertexHelix::ReconstructLambda(TPCLocalTrackHelix* track1, TPCLocalTrackHelix* track2)
{

  if(!IsCalculated()) return false;
  if(!IsLambda()) return false;

  return true;
}

//_____________________________________________________________________________
Bool_t
TPCVertexHelix::IsLambda()
{

  Double_t lambda_masscut = 0.1; //ref
  static const auto PionMass    = pdg::PionMass();
  //static const auto KaonMass    = pdg::KaonMass();
  static const auto ProtonMass  = pdg::ProtonMass();
  static const auto LambdaMass  = pdg::LambdaMass();

  if(!IsCalculated()) return false;
  if(m_distance > ppi_distcut) return false;

  m_lambda_vertex = m_vertex;
  for(Int_t i=0;i<2;i++){ //p, pi- or pi- p combination
    Int_t p_id = m_track_id[i];
    Int_t pi_id = m_track_id[1-i];
    Int_t p_pid = m_track_pid[i];
    Int_t pi_pid = m_track_pid[1-i];
    Int_t p_charge = m_track_charge[i];
    Int_t pi_charge = m_track_charge[1-i];
    TVector3 p_mom = m_track_mom[i];
    TVector3 pi_mom = m_track_mom[1-i];

    if((p_pid&4)!=4 || p_charge!=1) continue;
    if((pi_pid&1)!=1 || pi_charge!=-1) continue;

    TLorentzVector Lpi(pi_mom, TMath::Hypot(pi_mom.Mag(), PionMass));
    TLorentzVector Lp(p_mom, TMath::Hypot(p_mom.Mag(), ProtonMass));
    TLorentzVector Llambda = Lp + Lpi;
    TVector3 lambda_mom = pi_mom + p_mom;

    if(TMath::Abs(Llambda.M() - LambdaMass) > lambda_masscut) continue;
    m_is_lambda = true;
    m_lambda_mass.push_back(Llambda.M());
    m_lambda_mom.push_back(lambda_mom);
    m_proton_mom.push_back(p_mom);
    m_pion_mom.push_back(pi_mom);
    m_proton_id.push_back(p_id);
    m_pion_id.push_back(pi_id);
  }

  return m_is_lambda;
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
