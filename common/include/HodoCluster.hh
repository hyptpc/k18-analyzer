/**
 *  file: HodoCluster.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO_CLUSTER_HH
#define HODO_CLUSTER_HH

#include <cstddef>
#include <vector>

class Hodo2Hit;
class HodoAnalyzer;

//______________________________________________________________________________
class HodoCluster
{
public:
  HodoCluster(  );
  //  HodoCluster( Hodo2Hit *hit, int nth );
  virtual ~HodoCluster( void );

private:
  HodoCluster( const HodoCluster & );
  HodoCluster & operator = ( const HodoCluster & );

private:
  std::vector<Hodo2Hit*> m_hit;
  std::vector<int> m_nthhit;

  double    m_mean_time;
  double    m_de;
  double    m_mean_seg;
  double    m_time_diff;
  int       m_1st_index;
  double    m_1st_seg;
  double    m_1st_time;
  double    m_1st_de;
  bool      m_good_for_analysis;

public:
  int       AddHit( Hodo2Hit *hit, int nth );
  Hodo2Hit* GetHit( int i )         const;
  int       GetNthHit( int i )      const { return i<ClusterSize() ? m_nthhit.at(i) : -1; }
  Hodo2Hit* Get1stHit( )            const { return m_hit.at(m_1st_index); }
  int       Get1stNthHit( )         const { return m_nthhit.at(m_1st_index); }
  int       ClusterSize( void )     const { return (int)m_hit.size(); }
  double    CMeanTime( void )       const { return m_mean_time;    }
  double    DeltaE( void )          const { return m_de;           }
  double    MeanSeg( void )         const { return m_mean_seg;     }
  double    TimeDif( void )         const { return m_time_diff;    }
  double    C1stSeg( void )         const { return m_1st_seg;      }
  double    C1stTime( void )        const { return m_1st_time;     }
  double    C1stDe( void )          const { return m_1st_de;     }
  bool      GoodForAnalysis( void ) const { return m_good_for_analysis; }
  bool      GoodForAnalysis( bool status )
  {
    bool pre_status = m_good_for_analysis;
    m_good_for_analysis = status;
    return pre_status;
  }

  bool ReCalc( bool applyRecusively=false );
  void Print();

private:
  void Calculate( void );
};
#endif
