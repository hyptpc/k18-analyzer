/**
 *  file: HodoCluster.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO_CLUSTER_HH
#define HODO_CLUSTER_HH

#include <cstddef>

class Hodo2Hit;
class HodoAnalyzer;

//______________________________________________________________________________
class HodoCluster
{
public:
  HodoCluster( Hodo2Hit *hitA, Hodo2Hit *hitB=0, Hodo2Hit *hitC=0 );
  virtual ~HodoCluster( void );

private:
  HodoCluster( const HodoCluster & );
  HodoCluster & operator = ( const HodoCluster & );

private:
  Hodo2Hit *m_hitA;
  Hodo2Hit *m_hitB;
  Hodo2Hit *m_hitC;
  int       m_cluster_size;
  double    m_mean_time;
  double    m_de;
  double    m_mean_seg;
  double    m_time_diff;
  double    m_1st_seg;
  double    m_1st_time;
  bool      m_good_for_analysis;

public:
  Hodo2Hit* GetHit( int i )         const;
  int       ClusterSize( void )     const { return m_cluster_size; }
  double    CMeanTime( void )       const { return m_mean_time;    }
  double    DeltaE( void )          const { return m_de;           }
  double    MeanSeg( void )         const { return m_mean_seg;     }
  double    TimeDif( void )         const { return m_time_diff;    }
  double    C1stSeg( void )         const { return m_1st_seg;      }
  double    C1stTime( void )        const { return m_1st_time;     }
  bool      GoodForAnalysis( void ) const { return m_good_for_analysis; }
  bool      GoodForAnalysis( bool status )
  {
    bool pre_status = m_good_for_analysis;
    m_good_for_analysis = status;
    return pre_status;
  }

  bool ReCalc( bool applyRecusively=false );

private:
  void Calculate( void );
};
#endif
