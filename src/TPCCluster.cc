// -*- C++ -*-

#include <iostream>
#include <iterator>

#include <escape_sequence.hh>
#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "PrintHelper.hh"
#include "TPCCluster.hh"
#include "TPCPadHelper.hh"
#include "TPCPositionCorrector.hh"
#include "ThreeVector.hh"

#define WeightedMean 0
#define WeightedMeanTheta 1

namespace
{
const auto& gTPCPos = TPCPositionCorrector::GetInstance();
}

//_____________________________________________________________________________
TPCCluster::TPCCluster(Int_t layer, const TPCHitContainer& HitCont)
  : m_is_good(false),
    m_layer(layer),
    m_cluster_de(),
    m_cluster_position(),
    m_hit_array(HitCont), // shallow copy
    m_mean_row(),
    m_mean_phi(),
    // m_pos_center(),
    // m_cluster_de_center(),
    m_mean_hit(new TPCHit(layer, TMath::QuietNaN()))
{
  auto itr = m_hit_array.begin();
  while(itr != m_hit_array.end()){
    if(!*itr || !(*itr)->IsGood()){
      itr = m_hit_array.erase(itr);
    }else{
      itr++;
    }
  }
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCCluster::~TPCCluster()
{
  ClearTPCHits();
  delete m_mean_hit;
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
TPCCluster::ClearTPCHits()
{
  m_hit_array.clear();
}

//_____________________________________________________________________________
void
TPCCluster::AddTPCHit(TPCHit* hit)
{
  if(hit) m_hit_array.push_back(hit);
  m_is_good = false;
}

//_____________________________________________________________________________
Bool_t
TPCCluster::Calculate()
{
  static const TVector2 target_center(0., tpc::ZTarget); // (X, Z)
  const Double_t R = tpc::GetRadius(m_layer);
 	int max_row= tpc::padParameter[m_layer][1];
  m_cluster_de = 0.;
  m_cluster_position.SetXYZ(0., 0., 0.);
  m_mean_row = 0.;
  m_mean_phi = 0.;
  Double_t mean_y = 0.;
  Double_t buf_phi = TMath::QuietNaN();
  Double_t buf_row = TMath::QuietNaN();
  bool bf1= false;//branching flag
  bool bf2= false;//branching flag
	double branch=0;
	double branch_row=0;
	for(const auto& hit: m_hit_array){
    const auto& pos = hit->GetPosition();
    TVector2 xz_vector(pos.X(), pos.Z());
    xz_vector -= target_center;
#if 0
    if(TMath::Abs(xz_vector.Mod() - R) > 1e-10){
      hit->Print();
      throw Exception(FUNC_NAME + Form(" found invalid radius %lf/%lf",
                                       R, xz_vector.Mod()));
    }
#endif
    const Double_t de = hit->GetCDe();
    Double_t phi = xz_vector.Phi();
    Double_t row = hit->GetRow();
		if(phi<TMath::Pi()/2){
			bf1=true;//any small enough angle would be OK.
			branch+=2*TMath::Pi()*de;
			branch_row+=max_row*de;
		}
		if(phi>3*TMath::Pi()/2) bf2=true;//any small enough angle would be OK.
    m_mean_row += row * de;
    m_mean_phi += phi * de;
    mean_y += pos.Y() * de;
    m_cluster_de += de;
    buf_phi = phi;
    buf_row = row;
  }
	if(bf1&&bf2){
		m_mean_phi = fmod((m_mean_phi+branch)/m_cluster_de,2*TMath::Pi());
		m_mean_row = fmod((m_mean_row+branch_row)/m_cluster_de,max_row);
	}
	else{
		m_mean_phi *= 1./m_cluster_de;
		m_mean_row *= 1./m_cluster_de;
	}
	mean_y *= 1./m_cluster_de;
  TVector2 xz_vector;
  xz_vector.SetMagPhi(R, m_mean_phi);
  xz_vector += target_center;
  m_cluster_position.SetXYZ(xz_vector.X(), mean_y, xz_vector.Y());
//  if(m_mean_row<0) m_mean_row += tpc::padParameter[m_layer][1];
  m_mean_hit->SetMRow(m_mean_row);
  m_mean_hit->AddHit(0., 0.);
  m_mean_hit->SetDe(m_cluster_de);
  m_mean_hit->SetPosition(m_cluster_position);
  m_mean_hit->SetParentCluster(this);
  m_is_good = true;
  return true;
}

//_____________________________________________________________________________
TPCHit*
TPCCluster::GetCenterHit() const
{

  double rowdiff = 10; int id = -1;
  for(Int_t i=0; i<m_hit_array.size(); ++i){
    if(!m_hit_array[i]) continue;
    double row = (double) m_hit_array[i] -> GetRow();
    if(TMath::Abs(m_mean_row-row)<rowdiff){
      rowdiff = TMath::Abs(m_mean_row-row);
      id = i;
    }
  }
  return m_hit_array[id];
}

//_____________________________________________________________________________
void
TPCCluster::Print(Option_t* opt) const
{
  PrintHelper helper(1, std::ios::fixed);
  const Double_t R = tpc::GetRadius(m_layer);
  hddaq::cout << FUNC_NAME << " " << std::endl
              << "is good = "  << m_is_good  << std::endl
              << "de = "  << m_cluster_de << std::endl
              << "position = " << m_cluster_position << std::endl
              << "Radius = " << R << std::endl
              << "mean row = " << m_mean_row << std::endl
              << "mean phi = " << m_mean_phi*TMath::RadToDeg() << " (in XZ plane)" << std::endl;
  hddaq::cout << "L" << std::setw(2) << m_layer << " size="
              << std::setw(3) << m_hit_array.size() << "  ";
  for(const auto& hit: m_hit_array){
    hddaq::cout << hit->GetRow()
                << "(de=" << hddaq::unpacker::esc::k_purple
                << hit->GetCDe() << hddaq::unpacker::esc::k_default_color
                << ", y=" << hddaq::unpacker::esc::k_cyan
                << hit->GetDriftLength() << hddaq::unpacker::esc::k_default_color
                << ")" << " ";
    const auto& pos = hit->GetPosition();
    TVector2 xz_vector(pos.X(), pos.Z());
    TVector2 target_position(0., tpc::ZTarget);
    auto phi = (xz_vector - target_position).Phi()*TMath::RadToDeg();
    auto residual = tpc::ArcLength(m_layer, hit->GetRow(), m_mean_row);
    hddaq::cout << " pos=" << pos << ", phi=" << phi << ", res=" << residual << std::endl;
  }
  hddaq::cout << std::endl;
}
