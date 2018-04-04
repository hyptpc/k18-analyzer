/**
 *  file: BH2Cluster.cc
 *  date: 2017.04.10
 *
 */

#include "BH2Cluster.hh"

#include <cmath>
#include <cstring>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "BH2Hit.hh"
#include "DebugCounter.hh"
#include "HodoAnalyzer.hh"

namespace
{
  const std::string& class_name("BH2Cluster");
}

//______________________________________________________________________________
BH2Cluster::BH2Cluster( BH2Hit *hitA, BH2Hit *hitB, BH2Hit *hitC ) :
  m_hitA(hitA),
  m_hitB(hitB),
  m_hitC(hitC),
  m_cluster_size(0),
  m_good_for_analysis(true)
{
  if(hitA) ++m_cluster_size;
  if(hitB) ++m_cluster_size;
  if(hitC) ++m_cluster_size;
  Calculate();
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
BH2Cluster::~BH2Cluster( void )
{
  debug::ObjectCounter::decrease(class_name);
};

//______________________________________________________________________________
BH2Hit*
BH2Cluster::GetHit( int i ) const
{
  if(i==0) return m_hitA;
  if(i==1) return m_hitB;
  if(i==2) return m_hitC;
  return 0;
}

//______________________________________________________________________________
void
BH2Cluster::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  double ms=0., mt=0., de=0., t0=0., dt=0;
  if( m_hitA ){
    ms += m_hitA->SegmentId();
    mt += m_hitA->CMeanTime();
    de += m_hitA->DeltaE();
    t0 += m_hitA->CTime0();
    dt += ( m_hitA->GetTDown() - m_hitA->GetTUp() );
  }
  if( m_hitB ){
    ms += m_hitB->SegmentId();
    mt += m_hitB->CMeanTime();
    de += m_hitB->DeltaE();
    t0 += m_hitB->CTime0();
    dt += ( m_hitB->GetTDown() - m_hitB->GetTUp() );
  }
  if( m_hitC ){
    ms += m_hitC->SegmentId();
    mt += m_hitC->CMeanTime();
    de += m_hitC->DeltaE();
    t0 += m_hitA->CTime0();
    dt += ( m_hitC->GetTDown() - m_hitB->GetTUp() );
  }

  ms /= double(m_cluster_size);
  mt /= double(m_cluster_size);
  t0 /= double(m_cluster_size);
  dt /= double(m_cluster_size);

  m_mean_seg  = ms;
  m_mean_time = mt;
  m_de        = de;
  m_time0     = t0;
  m_time_diff = dt;

}

//______________________________________________________________________________
bool
BH2Cluster::GoodForAnalysis( bool status )
{
  bool pre_status = m_good_for_analysis;
  m_good_for_analysis = status;
  return pre_status;
}

//______________________________________________________________________________
bool
BH2Cluster::ReCalc( bool applyRecursively )
{
  if( applyRecursively ){
    if(m_hitA) m_hitA->ReCalc(applyRecursively);
    if(m_hitB) m_hitB->ReCalc(applyRecursively);
    if(m_hitC) m_hitC->ReCalc(applyRecursively);
  }
  Calculate();
  return true;
}
