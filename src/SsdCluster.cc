/**
 *  file: SsdCluster.cc
 *  date: 2017.04.10
 *
 */

#include "SsdCluster.hh"

#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>

#include <std_ostream.hh>

#include "DCGeomMan.hh"
#include "DCHit.hh"
#include "DCAnalyzer.hh"
#include "DeleteUtility.hh"

namespace
{
  const std::string& class_name("SsdCluster");
  const DCGeomMan& gGeom = DCGeomMan::GetInstance();
}


//______________________________________________________________________________
SsdCluster::SsdCluster( const DCHitContainer& HitCont )
  : m_rep_hit(0),
    m_good_for_analysis( true )
{
  m_hit_array.resize( HitCont.size() );
  std::copy( HitCont.begin(), HitCont.end(), m_hit_array.begin() );
  m_cluster_size = m_hit_array.size();

  Calculate();
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
SsdCluster::~SsdCluster( void )
{
  del::DeleteObject( m_rep_hit );
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
SsdCluster::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_cluster_size>0 )
    m_layer = m_hit_array[0]->GetLayer();

  double mean_time = 0.;
  double sum_de    = 0.;
  double sum_amp   = 0.;
  double mean_seg  = 0.;
  double mean_pos  = 0.;

  for( int i=0; i<m_cluster_size; ++i ){
    if( !m_hit_array[i] ) continue;
    double de   = m_hit_array[i]->GetDe();
    double amp  = m_hit_array[i]->GetAmplitude();
    double time = m_hit_array[i]->GetPeakTime();
    mean_seg  += m_hit_array[i]->GetWire() * de;
    mean_time += time * de;
    sum_de    += de;
    sum_amp   += amp;
    mean_pos  += m_hit_array[i]->GetWirePosition() * de;
  }

  mean_seg  /= sum_de;
  mean_time /= sum_de;
  mean_pos  /= sum_de;

  m_mean_time = mean_time;
  m_de        = sum_de;
  m_amplitude = sum_amp;
  m_mean_seg  = mean_seg;
  m_mean_pos  = mean_pos;

  double angle  = gGeom.GetTiltAngle( m_layer );
  double localz = gGeom.GetLocalZ( m_layer );

  // repesentative hit for tracking
  m_rep_hit = new DCHit( m_layer, m_mean_seg );
  m_rep_hit->SetTiltAngle( angle );
  m_rep_hit->SetWire( m_mean_seg );
  m_rep_hit->SetWirePosition( m_mean_pos );
  m_rep_hit->SetZ( localz );
  m_rep_hit->SetDummyPair();
  m_rep_hit->SetSsdFlag();
  m_rep_hit->SetAmplitude( m_amplitude );
  m_rep_hit->SetDe( m_de );
  m_rep_hit->SetPeakTime( m_mean_time );
}

//______________________________________________________________________________
bool
SsdCluster::GoodForAnalysis( bool status )
{
  bool pre_status = m_good_for_analysis;
  m_good_for_analysis = status;
  return pre_status;
}

//______________________________________________________________________________
void
SsdCluster::Print( const std::string& arg )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  hddaq::cout << "#D " << func_name << " " << arg << std::endl
	      << "ClusterSize     : " << m_cluster_size << std::endl
	      << "Layer           : " << m_layer        << std::endl
	      << "MeanTime        : " << m_mean_time    << std::endl
	      << "DeltaE          : " << m_de           << std::endl
	      << "Amplitude       : " << m_amplitude    << std::endl
	      << "MeanSeg         : " << m_mean_seg     << std::endl
	      << "MeanPosition    : " << m_mean_pos     << std::endl
	      << "GoodForAnalysis : " << m_good_for_analysis << std::endl;
  m_rep_hit->Print();
#if 0
  for( DCHit *hit : m_hit_array ) hit->Print();
#endif
}
