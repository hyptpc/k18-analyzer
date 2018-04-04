/**
 *  file: SsdCluster.hh
 *  date: 2017.04.10
 *
 */

#ifndef SSD_CLUSTER_HH
#define SSD_CLUSTER_HH

#include <cstddef>
#include <string>

#include "DCAnalyzer.hh"
#include "DCHit.hh"

class DCAnalyzer;

//______________________________________________________________________________
class SsdCluster
{
public:
  SsdCluster( const DCHitContainer& HitCont );
  ~SsdCluster( void );

private:
  SsdCluster( const SsdCluster & );
  SsdCluster & operator =( const SsdCluster & );

private:
  // if kMaxClusterSize = 0, cluster size is unlimited.
  static const std::size_t kMaxClusterSize = 0;
  DCHitContainer  m_hit_array;
  DCHit          *m_rep_hit;
  int             m_cluster_size;
  int             m_layer;
  double          m_mean_time;
  double          m_de;
  double          m_amplitude;
  double          m_mean_seg;
  double          m_mean_pos;
  bool            m_good_for_analysis;

public:
  inline static const
  std::size_t MaxClusterSize( void ) { return kMaxClusterSize; }
  DCHit* GetHit( int i )         const { return m_hit_array[i]; }
  DCHit* GetHit( void )          const { return m_rep_hit;      }
  int    ClusterSize( void )     const { return m_cluster_size; }
  int    LayerId( void )         const { return m_layer;        }
  double TiltAngle( void )       const { return m_rep_hit->GetTiltAngle(); }
  double Time( void )            const { return m_mean_time;    }
  double DeltaE( void )          const { return m_de;           }
  double Amplitude( void )       const { return m_amplitude;    }
  double MeanSeg( void )         const { return m_mean_seg;     }
  double Position( void )        const { return m_mean_pos;     }
  bool   BelongToTrack( void )   const { return m_rep_hit->BelongToTrack(); }
  bool   GoodForAnalysis( void ) const { return m_good_for_analysis; }
  bool   GoodForAnalysis( bool status );
  void   Print( const std::string& arg="" );
  void   JoinKaonTrack( void ) { m_rep_hit->JoinKaonTrack(); }
  void   QuitKaonTrack( void ) { m_rep_hit->QuitKaonTrack(); }
  bool   BelongToKaonTrack( void ){ return m_rep_hit->BelongToKaonTrack(); }
private:
  void   Calculate( void );

};

#endif
