// -*- C++ -*-

#include "HodoCluster.hh"

#include <algorithm>
#include <iostream>

#include <TMath.h>

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "HodoHit.hh"
#include "PrintHelper.hh"

#include <std_ostream.hh>

#define USE_WEIGHTED_AVERAGE 0

//_____________________________________________________________________________
HodoCluster::HodoCluster(const HodoHC& cont,
                         const index_t& index)
  : m_is_good(false),
    m_hit_container(cont.size()),
    m_index(index),
    m_cluster_size(cont.size()),
    m_mean_time(),
    m_ctime(),
    m_time_diff(),
    m_de(),
    m_tot(),
    m_mean_position(),
    m_segment(),
    m_1st_seg(TMath::QuietNaN()),
    m_1st_time(TMath::QuietNaN()),
    m_time0(),
    m_ctime0()
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

  m_segment   = 0.;
  m_mean_time = 0.;
  m_ctime     = 0.;
  m_time0     = 0.;
  m_ctime0    = 0.;
  m_time_diff = 0.;
  m_de        = 0.;
  m_1st_seg   = TMath::QuietNaN();
  m_1st_time  = (DetectorName() == "BHT") ? DBL_MIN : DBL_MAX;
  
  Double_t sum_w = 0.;

  for(Int_t i=0; i<m_cluster_size; ++i){
    const auto& hit = m_hit_container[i];
    const auto& index = m_index[i];
    Double_t seg      = hit->SegmentId();
    Double_t mt       = hit->MeanTime(index);
    Double_t ctime    = hit->CMeanTime(index);
    Double_t time0    = hit->Time0(index);
    Double_t ctime0   = hit->CTime0(index);
    Double_t time_diff= hit->TimeDiff(index);
    Double_t de       = hit->DeltaE();
    
    // Add weights and accumulate
#ifdef USE_WEIGHTED_AVERAGE
    Double_t w = (de > 0.0) ? de : 0.0;
#else
    Double_t w = 1.0;
#endif
    sum_w += w; // Accumulate total weight (for normalization)

    m_segment   += seg * w;
    m_mean_time += mt * w;
    m_ctime     += ctime * w;
    m_time0     += time0 * w;
    m_ctime0    += ctime0 * w;
    m_time_diff += time_diff * w;

    // Find 1st Hit (Time) and calculate dE
    if(hit->GetName() == "BHT"){
      m_de = TMath::Max(m_de, de);
      if(m_1st_time < ctime){ 
        m_1st_seg  = hit->SegmentId();
        m_1st_time = ctime;
      }
    }else{
      m_de += de;
      if(ctime < m_1st_time){
        m_1st_seg  = hit->SegmentId();
        m_1st_time = ctime;
      }
    }
  }

  // Normalize by total weight
  if (sum_w <= 0.) {
#if 0
    hddaq::cerr << FUNC_NAME << " sum_w is non-positive : " << sum_w << std::endl
                << "  " << DetectorName() << " Size:" << m_cluster_size << std::endl;
    for(Int_t i=0; i<m_cluster_size; ++i){
       hddaq::cerr << "   Seg:" << m_hit_container[i]->SegmentId() 
                   << " dE:" << m_hit_container[i]->DeltaE() << std::endl;
    }
#endif
    m_is_good = false;
    return;
  }

  m_segment   /= sum_w;
  m_mean_time /= sum_w;
  m_ctime     /= sum_w;
  m_time0     /= sum_w;
  m_ctime0    /= sum_w;
  m_time_diff /= sum_w;

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
              << " cluster size  : " << m_cluster_size << std::endl
              << " mean position : " << m_mean_position << std::endl
              << " segment       : " << m_segment << std::endl
              << " 1st segment   : " << m_1st_seg << std::endl
              << " 1st time      : " << m_1st_time << std::endl
              << " time0         : " << m_time0 << std::endl
              << " ctime0        : " << m_ctime0 << std::endl;
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
