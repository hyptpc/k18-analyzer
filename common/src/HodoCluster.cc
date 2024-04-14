/**
 *  file: HodoCluster.cc
 *  date: 2017.04.10
 *
 */

#include "HodoCluster.hh"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

#include "DebugCounter.hh"
#include "Hodo2Hit.hh"
#include "HodoAnalyzer.hh"

namespace
{
  const std::string& class_name("HodoCluster");
}

//______________________________________________________________________________
HodoCluster::HodoCluster()
  : m_1st_index(0), m_good_for_analysis(true)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
HodoCluster::~HodoCluster( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
int
HodoCluster::AddHit( Hodo2Hit* hit, int i )
{
  if(!hit) return 0;
  m_hit.push_back(hit);
  m_nthhit.push_back(i);
  Calculate();
  return ClusterSize();
}

//______________________________________________________________________________
Hodo2Hit*
HodoCluster::GetHit( int i ) const
{
  if(i<m_hit.size()) return m_hit.at(i);
  return 0;
}

//______________________________________________________________________________
void
HodoCluster::Calculate( void )
{
  double ms = 0.;
  double mt = 0.;
  double de = 0.;
  double dt = 0.;
  double s1 = 0.;
  double t1 = 1e9;
  double de1 = 0.;
  
  for(int i=0;i<m_hit.size();i++){
    Hodo2Hit* tmphit=m_hit.at(i);
    int j=m_nthhit.at(i);
    ms += tmphit->SegmentId();
    mt += tmphit->CMeanTime(j);
    double tmpde=tmphit->DeltaE();
    if(tmphit->DetectorId()==DetIdBHT) tmpde=tmphit->DeltaE(j);
    de += tmpde;   
    dt += tmphit->CTimeDiff(j);
    if(tmphit->CMeanTime(j)<t1){
      s1  = tmphit->SegmentId();
      t1  = tmphit->CMeanTime(j);
      de1  = tmphit->DeltaCE(j);
      m_1st_index=i;
    }
  }
  ms /= double(ClusterSize());
  mt /= double(ClusterSize());
  dt /= double(ClusterSize());
  
  m_mean_time = mt;
  m_de        = de;
  m_mean_seg  = ms;
  m_time_diff = dt;
  m_1st_seg   = s1;
  m_1st_time  = t1;
  m_1st_de    = de1;
}

//______________________________________________________________________________
bool
HodoCluster::ReCalc( bool applyRecursively )
{
  if( applyRecursively ){
    for(int i=0;i<m_hit.size();i++){
      m_hit.at(i)->ReCalc(applyRecursively);
    }
  }
  Calculate();
  return true;
}

//______________________________________________________________________________
void
HodoCluster::Print( )
{
  for(int i=0;i<m_hit.size();i++){
    m_hit.at(i)->Print(m_nthhit.at(i));
  }
  std::cout<<"+++++"<<std::endl;
  std::cout<<"CMeanTime = "<<CMeanTime()<<std::endl;;
  std::cout<<"DeltaE    = "<<DeltaE()   <<std::endl;;
  std::cout<<"C1stSeg   = "<<C1stSeg()  <<std::endl;;
  std::cout<<"C1stTime  = "<<C1stTime() <<std::endl;;
  std::cout<<"+++++"<<std::endl;
}
