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

  ///// for TOF
  double    m_z;

  ///// For E40 Acrylic TOF
  double    m_ofs_dt;

  ///// for CFT
  double m_meanseg;
  double m_maxseg;
  double m_r;
  double m_phi;
  BoolVec m_belong_track;

  mutable std::vector <DCLTrackHit *> m_register_container;

public:
  bool CalcDCObservables( void );
  bool CalcMWPCObservables( void );
  bool CalcFiberObservables( void );
  bool CalcCFTObservables( void );
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

  ///// for TOF
  void SetZ( double z ) { m_z = z; }

  ///// For E40 Acrylic TOF
  void SetOfsdT( double ofs) { m_ofs_dt = ofs;}

  ///// for CFT
  void SetMeanSeg    ( double seg ) { m_meanseg    = seg;   }
  void SetMaxSeg    ( double seg ) { m_maxseg    = seg;   }
  void SetPositionR  ( double r   ) { m_r    = r;   }
  void SetPositionPhi( double phi ) { m_phi  = phi; }
  void SetTdcCFT( int tdc );


  void GateDriftTime(double min, double max, bool select_1st);

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

  ///// for TOF
  double GetZ( void ) const { return m_z; }

  ///// for CFT
  double GetMeanSeg    ( void ) const { return m_meanseg;  } 
  double GetMaxSeg     ( void ) const { return m_maxseg;  } 
  double GetPositionR  ( void ) const { return m_r;        } 
  double GetPositionPhi( void ) const { return m_phi;      } 

  void JoinTrack( int nh=0 ) { m_pair_cont.at(nh).belong_track = true; }
  void QuitTrack( int nh=0 ) { m_pair_cont.at(nh).belong_track = false;}
  bool BelongToTrack( int nh=0 ) const { return m_pair_cont.at(nh).belong_track; }
  bool IsWithinRange( int nh=0 ) const { return m_pair_cont.at(nh).dl_range; }

  void JoinTrackCFT( int nh=0 ) {m_belong_track[nh] = true; }
  void QuitTrackCFT( int nh=0 ) {m_belong_track[nh] = false; }
  bool BelongToTrackCFT( int nh=0 ) const { return m_belong_track[nh]; }

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
