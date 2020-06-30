// -*- C++ -*-

#ifndef FIBER_CLUSTER_HH
#define FIBER_CLUSTER_HH

#include <vector>

class FLHit;

//______________________________________________________________________________
class FiberCluster
{
public:
  FiberCluster( void );
  virtual ~FiberCluster( void );

private:
  FiberCluster( const FiberCluster& object );
  FiberCluster& operator =( const FiberCluster& object );

protected:
  typedef std::vector<FLHit*> HitContainer;
  enum FlagsFiber { Initialized, gfastatus, sizeFlagsFiber };
  HitContainer m_hit_container;
  int          m_cluster_size;
  int          m_cluster_id;
  int          m_max_cluster_id;
  double       m_mean_time;
  double       m_max_time;
  double       m_real_mean_time;
  // real mean (not a closest value of CTime)
  double       m_max_width;
  double       m_min_width;
  double       m_mean_seg;
  double       m_mean_pos;
  double       m_sum_adc_lg;
  double       m_sum_mip_lg;
  double       m_sum_de_lg;
  double       m_max_adc_hg;
  double       m_max_adc_lg;
  double       m_max_mip_lg;
  double       m_max_de_lg;
  double       m_max_seg;
  bool         m_flag[sizeFlagsFiber];

public:
  bool   Calculate( void );
  void   push_back( FLHit* hit ) { m_hit_container.push_back(hit); };
  int    VectorSize( void )      const { return m_hit_container.size(); }
  int    ClusterId( void )       const { return m_cluster_id;      }
  int    ClusterSize( void )     const { return m_cluster_size;    }
  int    GetMaxClusterId( void ) const { return m_max_cluster_id;  }
  double CMeanTime( void )       const { return m_mean_time;       }
  double CMaxTime( void )        const { return m_max_time;        }
  double RCMeanTime( void )      const { return m_real_mean_time;  }
  double Width( void )           const { return m_max_width;       }
  double minWidth( void )        const { return m_min_width;       }
  double Tot( void )             const { return Width();           }
  double MeanPosition( void )    const { return m_mean_pos;        }
  double SumAdcLG( void )        const { return m_sum_adc_lg;      }
  double SumMipLG( void )        const { return m_sum_mip_lg;      }
  double SumDeLG( void )         const { return m_sum_de_lg;       }
  double MeanSeg( void )         const { return m_mean_seg;        }
  double MaxAdcHi( void )        const { return m_max_adc_hg;      }
  double MaxAdcLG( void )        const { return m_max_adc_lg;      }
  double MaxMipLG( void )        const { return m_max_mip_lg;      }
  double MaxDeLG( void )         const { return m_max_de_lg;       }
  double MaxSeg( void )          const { return m_max_seg;         }
  FLHit* GetHit( int i ) const;
  bool   GoodForAnalysis( void ) const { return m_flag[gfastatus]; }
  bool   GoodForAnalysis( bool status );
  bool   ReCalc( bool applyRecusively=false );

private:
  void Debug( void );

};

#endif
