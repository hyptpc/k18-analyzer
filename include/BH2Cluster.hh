/**
 *  file: BH2Cluster.hh
 *  date: 2017.04.10
 *
 */

#ifndef BH2_CLUSTER_HH
#define BH2_CLUSTER_HH

class BH2Hit;
class HodoAnalyzer;

#include <cstddef>

//______________________________________________________________________________
class BH2Cluster
{
public:
  BH2Cluster( BH2Hit *hitA, BH2Hit *hitB=0, BH2Hit *hitC=0 );
  ~BH2Cluster( void );

private:
  BH2Cluster( const BH2Cluster& );
  BH2Cluster& operator =( const BH2Cluster& );

private:
  BH2Hit* m_hitA;
  BH2Hit* m_hitB;
  BH2Hit* m_hitC;
  int     m_cluster_size;
  double  m_mean_time;
  double  m_de;
  double  m_mean_seg;
  double  m_time_diff;
  double  m_time0;
  bool    m_good_for_analysis;

public:
  int     ClusterSize( void )     const { return m_cluster_size; }
  double  CMeanTime( void )       const { return m_mean_time;    }
  double  DeltaE( void )          const { return m_de;           }
  double  MeanSeg( void )         const { return m_mean_seg;     }
  double  CTime0( void )          const { return m_time0;        }
  double  TimeDif( void )         const { return m_time_diff;    }
  BH2Hit* GetHit( int i )         const;
  bool    GoodForAnalysis( void ) const { return m_good_for_analysis; }
  bool    GoodForAnalysis( bool status );
  bool    ReCalc( bool applyRecusively=false );

protected:
  void    Calculate( void );

};
#endif
