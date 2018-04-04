/**
 *  file: DCHit.hh
 *  date: 2017.04.10
 *
 */

#ifndef DC_HIT_HH
#define DC_HIT_HH

#include <cmath>
#include <vector>
#include <deque>
#include <string>
#include <numeric>

#include <std_ostream.hh>

#include "DebugCounter.hh"

//typedef std::vector<bool>   BoolVec;
typedef std::deque<bool>    BoolVec;
typedef std::vector<int>    IntVec;
typedef std::vector<double> DoubleVec;

class DCLTrackHit;

//______________________________________________________________________________
class DCHit
{
public:
  DCHit( void );
  DCHit( int layer );
  DCHit( int layer, double wire );
  ~DCHit( void );

private:
  DCHit( const DCHit& );
  DCHit& operator =( const DCHit& );

protected:
  int       m_layer;
  double    m_wire;
  IntVec    m_tdc;
  IntVec    m_adc;
  IntVec    m_trailing;

  //  DoubleVec m_dt;
  //  DoubleVec m_dl;
  //  DoubleVec m_trailing_time;

  // For DC with HUL MH-TDC
  struct data_pair{
    double drift_time;
    double drift_length;
    double trailing_time;
    double tot;
    int    index_t;
    bool   belong_track;
    bool   dl_range;
  };

  std::vector<data_pair> m_pair_cont;

  double  m_wpos;
  double  m_angle;
  //  BoolVec m_belong_track;
  //  BoolVec m_dl_range;

  ///// for MWPC
  int    m_cluster_size;
  bool   m_mwpc_flag;
  double m_mwpc_wire;
  double m_mwpc_wpos;

  ///// for SSD
  bool      m_is_ssd;
  bool      m_zero_suppressed;
  bool      m_time_corrected;
  bool      m_good_waveform;
  int       m_pedestal;
  int       m_peak_height;
  int       m_peak_position;
  double    m_deviation;
  double    m_amplitude;
  double    m_peak_time;
  double    m_adc_sum;
  double    m_de;
  double    m_rms;
  double    m_chisqr;
  DoubleVec m_time;
  DoubleVec m_waveform;
  bool      m_belong_kaon;
  ///// for TOF
  double    m_z;

  mutable std::vector <DCLTrackHit *> m_register_container;

public:
  bool CalcDCObservables( void );
  bool CalcMWPCObservables( void );
  bool CalcFiberObservables( void );
  bool CalcSsdObservables( void );
  //  bool CalcObservablesSimulation( double dlength);

  void SetLayer( int layer )              { m_layer = layer;                    }
  void SetWire( double wire )             { m_wire  = wire;                     }
  void SetTdcVal( int tdc )               { m_tdc.push_back(tdc);               }
  void SetAdcVal( int adc )               { m_adc.push_back(adc);               }
  void SetTdcTrailing(int tdc)            { m_trailing.push_back(tdc);          }
  void SetDummyPair();
  void SetDriftLength( int ith, double dl ){m_pair_cont.at(ith).drift_length = dl;}
  void SetDriftTime( int ith, double dt ) { m_pair_cont.at(ith).drift_time = dt;  }
  void SetTiltAngle( double angleDegree ) { m_angle = angleDegree;              }

  void SetClusterSize( int size )          { m_cluster_size   = size;           }
  void SetMWPCFlag( bool flag )            { m_mwpc_flag = flag;                }
  void SetMeanWire( double mwire )         { m_mwpc_wire = mwire;               }
  void SetMeanWirePosition( double mwpos ) { m_mwpc_wpos = mwpos;               }
  void SetWirePosition( double wpos )      { m_wpos      = wpos;                }

  ///// for SSD
  void SetSsdFlag( bool flag=true )          { m_is_ssd        = flag;          }
  void SetGoodWaveForm( bool good=true )     { m_good_waveform = good;          }
  void SetPedestal( int pedestal )           { m_pedestal      = pedestal;      }
  void SetRms( double rms )                  { m_rms           = rms;           }
  void SetAdcSum( double sum )               { m_adc_sum       = sum;           }
  void SetDe( double de )                    { m_de            = de;            }
  void SetDeviation( double deviation )      { m_deviation     = deviation;     }
  void SetTime( DoubleVec time )             { m_time          = time;          }
  void SetWaveform( DoubleVec waveform )     { m_waveform      = waveform;      }
  void SetAmplitude( double amplitude )      { m_amplitude     = amplitude;     }
  void SetPeakTime( double peaktime )        { m_peak_time     = peaktime;      }
  void SetPeakHeight( int height )           { m_peak_height   = height;        }
  void SetPeakPosition( int position )       { m_peak_position = position;      }
  void SetChisquare( double chisqr )         { m_chisqr        = chisqr;        }

  ///// for TOF
  void SetZ( double z ) { m_z = z; }

  int GetLayer( void ) const { return m_layer; }
  double GetWire( void )  const {
    if( m_mwpc_flag ) return m_mwpc_wire;
    else return int(m_wire);
  }

  int GetTdcSize( void )             const { return m_tdc.size(); }
  int GetAdcSize( void )             const { return m_adc.size(); }
  int GetDriftTimeSize( void )       const { return m_pair_cont.size(); }
  int GetDriftLengthSize( void )     const { return m_pair_cont.size(); }
  int GetTdcVal( int nh=0 )          const { return m_tdc[nh]; }
  int GetAdcVal( int nh=0 )          const { return m_adc[nh]; }
  int GetTdcTrailing( int nh=0 )     const { return m_trailing[nh]; }
  int GetTdcTrailingSize( void )     const { return m_trailing.size(); }

  double GetResolution( void )       const;

  double GetDriftTime( int nh=0 )    const { return m_pair_cont.at(nh).drift_time; }
  double GetDriftLength( int nh=0 )  const { return m_pair_cont.at(nh).drift_length; }
  double GetTrailingTime( int nh=0 ) const { return m_pair_cont.at(nh).trailing_time; }
  double GetTot( int nh=0 )          const { return m_pair_cont.at(nh).tot; }

  double GetTiltAngle( void ) const { return m_angle; }
  double GetWirePosition( void ) const {
    if( m_mwpc_flag ) return m_mwpc_wpos;
    else return m_wpos;
  }

  int GetClusterSize( void ) const { return m_cluster_size; }
  double GetMeamWire( void ) const { return m_mwpc_wire; }
  double GetMeamWirePosition( void ) const { return m_mwpc_wpos; }

  ///// for SSD
  bool      IsSsd( void )              const { return m_is_ssd;          }
  bool      IsTimeCorrected( void )    const { return m_time_corrected;  }
  bool      IsGoodWaveForm( void )     const { return m_good_waveform;   }
  bool      IsZeroSuppressed( void )   const { return m_zero_suppressed; }
  int       GetPedestal( void )        const { return m_pedestal;        }
  DoubleVec GetTime( void )            const { return m_time;            }
  DoubleVec GetWaveform( void )        const { return m_waveform;        }
  double    GetAmplitude( void )       const { return m_amplitude;       }
  double    GetDeviation( void )       const { return m_deviation;       }
  double    GetAdcSum( void )          const { return m_adc_sum;         }
  double    GetDe( void )              const { return m_de;              }
  double    GetPeakTime( void )        const { return m_peak_time;       }
  double    GetRms( void )             const { return m_rms;             }
  double    GetAdcPeakHeight( void )   const { return m_peak_height;     }
  double    GetAdcPeakPosition( void ) const { return m_peak_position;   }
  double    GetChisquare( void )       const { return m_chisqr;          }
  bool      DoTimeCorrection( double offset );
  void      JoinKaonTrack( void ) { m_belong_kaon = true; }
  void      QuitKaonTrack( void ) { m_belong_kaon = false; }
  bool      BelongToKaonTrack( void ) const { return m_belong_kaon; }

  ///// for TOF
  double GetZ( void ) const { return m_z; }

  void JoinTrack( int nh=0 ) { m_pair_cont.at(nh).belong_track = true; }
  void QuitTrack( int nh=0 ) { m_pair_cont.at(nh).belong_track = false; }
  bool BelongToTrack( int nh=0 ) const { return m_pair_cont.at(nh).belong_track; }
  bool IsWithinRange( int nh=0 ) const { return m_pair_cont.at(nh).dl_range; }

  void RegisterHits( DCLTrackHit *hit ) const
  { m_register_container.push_back(hit); }

  bool ReCalcDC( bool applyRecursively=false ) { return CalcDCObservables(); }
  bool ReCalcMWPC( bool applyRecursively=false ) { return CalcMWPCObservables(); }

  void TotCut(double min_tot, bool adotp_nan);

  void Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const;

protected:
  void ClearRegisteredHits( void );
};

//_____________________________________________________________________
inline std::ostream&
operator <<( std::ostream& ost, const DCHit& hit )
{
  hit.Print( "", ost );
  return ost;
}

#endif
