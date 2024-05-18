// -*- C++ -*-

#include "FiberCluster.hh"

#include <cmath>
#include <string>
#include <limits>

#include <TMath.h>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "FiberHit.hh"
#include "FuncName.hh"

namespace
{
const Bool_t reject_nan = false;
}

//_____________________________________________________________________________
FiberCluster::FiberCluster(const HodoHC& cont,
                           const index_t& index)
  : HodoCluster(cont, index)
{
  Calculate();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
FiberCluster::~FiberCluster()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
FiberCluster::Calculate()
{
  if(!m_is_good)
    return;

  m_tot           = 0.;
  m_mean_position = 0.;
  for(Int_t i=0; i<m_cluster_size; ++i){
    const auto& hit = dynamic_cast<FiberHit*>(m_hit_container[i]);
    const auto& index = m_index[i];
    m_tot = TMath::Max(m_tot, hit->MeanTOT(index));
    // m_tot += hit->MeanTOT(index);
    m_mean_position += hit->Position();
  }
  m_mean_position /= Double_t(m_cluster_size);
}
