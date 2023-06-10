// -*- C++ -*-

#include "HodoCluster.hh"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "HodoHit.hh"
#include "HodoAnalyzer.hh"
#include "PrintHelper.hh"

//_____________________________________________________________________________
HodoCluster::HodoCluster(const HodoHC& cont,
                         const index_t& index)
  : m_is_good(false),
    m_hit_container(cont.size()),
    m_index(index.size()),
    m_cluster_size(cont.size()),
    m_mean_time(),
    m_ctime(),
    m_time_diff(),
    m_de(),
    m_tot(),
    m_mean_position(),
    m_segment(),
    m_1st_seg(TMath::QuietNaN()),
    m_1st_time(DBL_MAX)
{
  if(cont.size() != index.size()){
    hddaq::cerr << FUNC_NAME << " cont size mismatch" << std::endl;
    return;
  }
  std::copy(cont.begin(), cont.end(), m_hit_container.begin());
  Calculate();
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
HodoCluster::~HodoCluster()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
HodoCluster::Calculate()
{
  m_is_good = false;

  if(m_hit_container.empty())
    return;

  m_mean_time = 0.;
  m_ctime     = 0.;
  m_time_diff = 0.;
  m_de        = 0.;
  m_segment   = 0.;
  m_1st_seg   = TMath::QuietNaN();
  m_1st_time  = DBL_MAX;

  for(Int_t i=0; i<m_cluster_size; ++i){
    const auto& hit = m_hit_container[i];
    const auto& index = m_index[i];
    m_segment   += hit->SegmentId();
    m_mean_time += hit->MeanTime(index);
    m_ctime     += hit->CMeanTime(index);
    m_time_diff += hit->TimeDiff(index);
    m_de        += hit->DeltaE();
    if(hit->CMeanTime(index) < m_1st_time){
      m_1st_seg  = hit->SegmentId();
      m_1st_time = hit->CMeanTime(index);
    }
  }

  m_segment   /= Double_t(m_cluster_size);
  m_mean_time /= Double_t(m_cluster_size);
  m_ctime     /= Double_t(m_cluster_size);
  m_time_diff /= Double_t(m_cluster_size);

  m_is_good = true;
}

//_____________________________________________________________________________
HodoHit*
HodoCluster::GetHit(Int_t i) const
{
  // try {
    return m_hit_container.at(i);
  // }catch(const std::out_of_range&){
  //   return nullptr;
  // }
}

//_____________________________________________________________________________
void
HodoCluster::Print(Option_t*) const
{
  PrintHelper helper(3, std::ios::fixed);
  hddaq::cout << FUNC_NAME << std::endl
              << " detector name : " << m_hit_container.at(0)->DetectorName() << std::endl
              << " mean time     : " << m_mean_time << std::endl
              << " ctime         : " << m_ctime << std::endl
              << " de            : " << m_de << std::endl
              << " tot           : " << m_tot << std::endl
              << " cluster size  : " << m_cluster_size << std::endl;
  for(Int_t i=0; i<m_cluster_size; ++i){
    const auto& hit = m_hit_container[i];
    const auto& j = m_index[i];
    hddaq::cout << "  " << hit->PlaneName();
    hddaq::cout << " " << hit->SegmentId();
    hddaq::cout << " "  << hit->CMeanTime(j);
  }
  hddaq::cout << std::endl;
}

//_____________________________________________________________________________
Bool_t
HodoCluster::ReCalc(Bool_t applyRecursively)
{
  if(applyRecursively){
    for(auto& hit: m_hit_container){
      hit->ReCalc(applyRecursively);
    }
  }
  Calculate();
  return m_is_good;
}
