//  Authors: Wooseung Jung

//GenKEK
#include "HypTPCHit.hh"

//k18-analyzer
#include "TPCLTrackHit.hh"

ClassImp(HypTPCHit)

HypTPCHit::HypTPCHit(TPCLTrackHit& right)
  : m_hit(right)
{}

HypTPCHit::~HypTPCHit(){}

TPCLTrackHit &HypTPCHit::GetHit() const{
  return m_hit;
}
