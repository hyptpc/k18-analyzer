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
  : m_layer(layer),
    m_pad_id(TMath::QuietNaN()),
    m_cluster_de(),
    m_cluster_position(),
    m_pos_center(),
    m_hit_array(HitCont), // shallow copy
    m_pos_calculated(false),
    m_mrow(),
    m_mrow_int(),
    m_cluster_de_center(),
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
TPCCluster::TPCCluster(Double_t x, Double_t y, Double_t z, Double_t de)
  : m_layer(tpc::getLayerID(tpc::findPadID(z, x))),
    m_pad_id(tpc::findPadID(z, x)),
    m_cluster_de(de),
    m_cluster_position(x, y, z),
    m_pos_center(x, y, z),
    m_hit_array(),
    m_pos_calculated(true), // for MC data
    m_mrow(),
    m_mrow_int(tpc::getRowID(tpc::findPadID(z, x))),
    m_cluster_de_center(),
    m_mean_hit(new TPCHit(m_layer, TMath::QuietNaN()))
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
TPCCluster::~TPCCluster()
{
  ClearTPCHits();
  debug::ObjectCounter::decrease(ClassName());
  delete m_mean_hit;
}

//_____________________________________________________________________________
void
TPCCluster::ClearTPCHits()
{
  m_hit_array.clear();
}

//_____________________________________________________________________________
// void
// TPCCluster::AddTPCHit(TPCHit* hit)
// {
//   if(hit){
//     m_hit_array.push_back(hit);
//     m_cluster_de += hit->GetCDe();
//     m_pos_calculated = false;
//   }
// }

//_____________________________________________________________________________
Bool_t
TPCCluster::Calculate()
{
  static const TVector2 target_center(0., tpc::ZTarget); // (X, Z)
  m_cluster_de = 0.;
  m_cluster_position.SetXYZ(0., 0., 0.);
  m_mrow = 0.;
  Double_t mean_phi = 0.;
  Double_t mean_y = 0.;
  const Double_t R = tpc::GetRadius(m_layer);
  Double_t max_phi=-1.;
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
    if(max_phi > 0 && (phi - max_phi) > TMath::Pi()/2.){
      phi -= 2.*TMath::Pi();
      row -= tpc::padParameter[m_layer][1];
    }
    m_mrow += row * de;
    mean_phi += phi * de;
    mean_y += pos.Y() * de;
    m_cluster_de += de;
    max_phi = TMath::Max(max_phi, phi);
  }
  m_mrow *= 1./m_cluster_de;
  mean_phi *= 1./m_cluster_de;
  mean_y *= 1./m_cluster_de;
  TVector2 xz_vector;
  xz_vector.SetMagPhi(R, mean_phi);
  xz_vector += target_center;
  m_cluster_position.SetXYZ(xz_vector.X(), mean_y, xz_vector.Y());

  m_mean_hit->SetMRow(m_mrow);
  m_mean_hit->AddHit(0., 0.);
  m_mean_hit->SetDe(m_cluster_de);
  m_mean_hit->SetPosition(m_cluster_position);
  m_mean_hit->SetParentCluster(this);

  m_is_good = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCCluster::CalculateWeightedMean()
{
  if(m_pos_calculated) return false;
  Double_t x=0, y=0, z=0, de=0, dummy_padid=0, mrow=0;
  for(Int_t hiti=0; hiti<m_hit_array.size(); hiti++) {
    x+=m_hit_array[hiti]->GetX()*m_hit_array[hiti]->GetCDe();
    y+=m_hit_array[hiti]->GetY()*m_hit_array[hiti]->GetCDe();
    z+=m_hit_array[hiti]->GetZ()*m_hit_array[hiti]->GetCDe();

    dummy_padid+=(Double_t)(m_hit_array[hiti]->GetPad())*m_hit_array[hiti]->GetCDe();
    mrow+=(Double_t)(m_hit_array[hiti]->GetRow())*m_hit_array[hiti]->GetCDe();
    de+=m_hit_array[hiti]->GetCDe();
  }
  if(de){
    m_cluster_position.SetXYZ(x/de, y/de, z/de);
    m_pad_id = (Int_t)(dummy_padid/de);
    m_mrow = mrow/de;
    m_cluster_de = de;

    if(m_mrow-(Int_t)m_mrow<0.5)
      m_mrow_int = (Int_t)m_mrow;
    else
      m_mrow_int = 1+(Int_t)m_mrow;
    for(Int_t hiti=0; hiti<m_hit_array.size(); hiti++) {
      if(m_mrow_int==m_hit_array[hiti]->GetRow()){
	m_cluster_de_center = m_hit_array[hiti]->GetCDe();
	m_pos_center = m_hit_array[hiti]->GetPosition();
      }
    }
  }
  else {
    m_cluster_position.SetXYZ(0, 0, 0);
    m_pad_id = 0;
    m_mrow = 0;
  }
  m_pos_calculated = true;
  m_is_good = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCCluster::CalculateWeightedMeanTheta()
{
  if(m_pos_calculated) return false;
  m_cluster_de = 0.;
  [[maybe_unused]] Double_t x=0, y=0, z=0, dummy_padid=0, mrow=0;
  for(const auto& hit: m_hit_array){
    const Double_t de = hit->GetCDe();
    y += hit->GetY()*de;

    TVector3 pos_check = tpc::getPosition(hit->GetLayer(), (Double_t)hit->GetRow());
    // std::cout<<"(x,z)=("<<hit->GetX()
    //    	     <<","<<hit->GetZ()<<")"<<std::endl;
    // std::cout<<"theta, (x,z)=("<<pos_check.x()
    //  	     <<","<<pos_check.z()<<")"<<std::endl;

    dummy_padid+=(Double_t)(hit->GetPad())*hit->GetCDe();
    mrow+=(Double_t)(hit->GetRow())*hit->GetCDe();
    m_cluster_de += de;
  }
  if(m_cluster_de > 0.){
    //m_cluster_position.SetXYZ(x/de, y/de, z/de);
    //TVector3 m_pos_dummy(x/de, y/de, z/de);
    m_pad_id = (Int_t)(dummy_padid/m_cluster_de);
    m_mrow = mrow/m_cluster_de;
    TVector3 pos_xz = tpc::getPosition(tpc::getLayerID(m_pad_id), m_mrow);
    m_cluster_position.SetXYZ(pos_xz.x(), y/m_cluster_de, pos_xz.z());
    // std::cout<<"Weight (x,z) = ("<<m_pos_dummy.x()<<", "
    // 	     <<m_pos_dummy.z()<<")"<<std::endl;
    // std::cout<<"Theta (x,z) = ("<<m_cluster_position.x()<<", "
    // 	     <<m_cluster_position.z()<<")"<<std::endl;
    if(m_mrow-(Int_t)m_mrow<0.5)
      m_mrow_int = (Int_t)m_mrow;
    else
      m_mrow_int = 1+(Int_t)m_mrow;
    for(const auto& hit: m_hit_array){
      if(m_mrow_int==hit->GetRow()){
	m_cluster_de_center = hit->GetCDe();
	m_pos_center = hit->GetPosition();
      }
    }
  }
  else {
    m_cluster_position.SetXYZ(0, 0, 0);
    m_pad_id = 0;
    m_mrow = 0;
  }
  m_pos_calculated = true;
  return true;
}

//_____________________________________________________________________________
Double_t
TPCCluster::ResX() const
{
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  Double_t y_pos= m_cluster_position.Y();
  Double_t s0 = 0.204;// mm HIMAC result //To do parameter
  Double_t Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  Double_t L_D = 30.+(y_pos*0.1);//cm
  Double_t N_eff = 42.8;
  Double_t A = 0.0582*0.01;//m-1 -> cm-1
  Double_t e_ALD = exp(-1.*A*L_D);
  Double_t sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  Double_t sT = sqrt(sT2);

  return sT;
}

//_____________________________________________________________________________
Double_t
TPCCluster::ResY() const
{
  //Temp value (need to change)
  Double_t y_res = 0.5;
  return y_res;
}

//_____________________________________________________________________________
Double_t
TPCCluster::ResZ() const
{
  //calculated by using NIM paper
  //To do:change the resolution by checking cluster size
  Double_t y_pos= m_cluster_position.Y();
  Double_t s0 = 0.204;// mm HIMAC result //To do parameter
  Double_t Dt = 0.18;//mm/sqrt(cm) at 1T //To do parameter
  Double_t L_D = 30.+(y_pos*0.1);//cm
  Double_t N_eff = 42.8;
  Double_t A = 0.0582*0.01;//m-1 -> cm-1
  Double_t e_ALD = exp(-1.*A*L_D);
  Double_t sT2 = s0*s0 + (Dt*Dt*L_D/(N_eff*e_ALD));
  Double_t sT = sqrt(sT2);

  return sT;
}


//_____________________________________________________________________________
void
TPCCluster::Print(Option_t* opt) const
{
  PrintHelper helper(1, std::ios::fixed);
  hddaq::cout << FUNC_NAME << " " << std::endl
              // << "layerID= "  << m_layer_id  << std::endl
              << "de = "  << m_cluster_de << std::endl
              << "position = " << m_cluster_position << std::endl;
  hddaq::cout << "L" << std::setw(2) << m_layer << " size="
              << std::setw(3) << m_hit_array.size() << "  ";
  for(const auto& hit: m_hit_array){
    hddaq::cout << hit->GetRow()
                << "(" << hddaq::unpacker::esc::k_purple
                << hit->GetCDe() << hddaq::unpacker::esc::k_default_color
                << "," << hddaq::unpacker::esc::k_cyan
                << hit->GetDriftLength() << hddaq::unpacker::esc::k_default_color
                << ")" << " ";
    // const auto& pos = hit->GetPosition();
    // TVector2 xz_vector(pos.X(), pos.Z());
    // TVector2 target_position(0., tpc::ZTarget);
    // hddaq::cout << " R=" << (xz_vector - target_position).Mod() << std::endl;
  }
  hddaq::cout << std::endl;
}
