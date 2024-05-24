// -*- C++ -*-

#ifndef DC_HIT_HH
#define DC_HIT_HH

#include <cmath>
#include <vector>
#include <deque>
#include <string>
#include <numeric>

#include <TVector3.h>

#include <std_ostream.hh>

#include "DebugCounter.hh"

//typedef std::vector<bool>   BoolVec;
typedef std::deque<bool>    BoolVec;
typedef std::vector<int>    IntVec;
typedef std::vector<double> DoubleVec;

class DCRawHit;
class DetectorHit;

//_____________________________________________________________________________
class DCHit
{
public:
  DCHit(const DCRawHit* rhit);
  DCHit( int layer );
  DCHit( int layer, int wire );
  DCHit( int cid, int layer, int wire );
  DCHit( DetectorHit *mchit);
  ~DCHit( void );

private:
  DCHit( const DCHit& );
  DCHit& operator =( const DCHit& );

private:
  const DCRawHit* m_raw_hit;
  int       m_cid;
  int       m_layer;
  int       m_wire;
  IntVec    m_tdc;
  IntVec    m_trailing;
  DoubleVec m_dt;
  DoubleVec m_dl;
  DoubleVec m_trailing_time;
  DoubleVec m_tot;

  TVector3  m_wpos;
  TVector3  m_wdir;
  double    m_tilt;
  double    m_rotation;
  int       m_xy;
  BoolVec   m_belong_track;
  BoolVec   m_belong_cluster;
  BoolVec   m_dt_range;
  BoolVec   m_tot_range;
  BoolVec   m_dl_range;

  //  mutable std::vector <DCLTrackHit *> m_register_container;

public:
  bool CalcDCObservables( double retiming=0 );

  void SetDetId( int detid )              { m_cid = detid; }
  void SetLayer( int layer )              { m_layer = layer; }
  void SetWire( int wire )             { m_wire  = wire; }
  void SetTdcVal( int tdc );
  void SetAdcVal( int adc );
  void SetTdcTrailing(int tdc)            { m_trailing.push_back(tdc); }
  void SetDriftTime( double dt )          { m_dt.push_back(dt); }
  void ClearDriftTime( void )             { m_dt.clear(); }
  void SetDriftLength( double dl )        { m_dl.push_back(dl); }
  void SetDriftLength( int ith, double dl ) { m_dl[ith] = dl; }
  void ClearDriftLength( void )           { m_dl.clear();}
  void SetTiltAngle( double angleDegree ) { m_tilt = angleDegree; }
  void SetRotationAngle( double angleDegree ) { m_rotation = angleDegree; }
  void SetTrailingTime( double t )        { m_trailing_time.push_back(t); }
  void SetWirePosition( TVector3 wpos )   { m_wpos     = wpos; }

  int GetDetId( void ) const { return m_cid; }
  int GetLayer( void ) const { return m_layer; }
  int GetWire( void )  const { return m_wire;  }
  int GetXY( void )    const { return m_xy; }

  int GetTdcSize( void ) const { return m_tdc.size(); }
  int GetDriftTimeSize( void ) const { return m_dt.size(); }
  int GetDriftLengthSize( void ) const { return m_dl.size(); }
  int GetTdcTrailingSize( void ) const { return m_trailing.size(); }
  int GetTdcTrailingTimeSize( void ) const { return m_trailing_time.size(); }
  int GetTOTSize( void ) const { return m_tot.size(); }

  int GetTdcVal( int nh=0 ) const { return m_tdc[nh]; }
  int GetTdcTrailing( int nh=0 ) const { return m_trailing[nh]; }
  double GetDriftTime( int nh=0 ) const { return m_dt[nh]; }
  double GetDriftLength( int nh=0 ) const { return m_dl[nh]; }
  double GetTrailingTime( int nh=0 ) const { return m_trailing_time[nh]; }
  double GetTOT( int nh=0 ) const { return m_tot[nh]; }

  double GetResolution( void ) const;
  int GetNHit(bool DT=false, bool TOT=false) const;
  int GetHitID(int i,bool DT=false, bool TOT=false) const;


  double GetTiltAngle( void ) const { return m_tilt; }
  double GetRotationAngle( void ) const { return m_rotation; }
  TVector3 GetWirePosition( void ) const {    return m_wpos;  }
  TVector3 GetWireDirection( void ) const {    return m_wdir;  }

  void JoinTrack( int nh=0 ) { m_belong_track[nh] = true; }
  void QuitTrack( int nh=0 ) { m_belong_track[nh] = false; }
  bool BelongToTrack( int nh=0 ) const { return m_belong_track[nh]; }
  void JoinCluster( int nh=0 ) { m_belong_cluster[nh] = true; }
  void QuitCluster( int nh=0 ) { m_belong_cluster[nh] = false; }
  bool BelongToCluster( int nh=0 ) const { return m_belong_cluster[nh]; }
  bool IsWithinDlRange( int nh=0 ) const { return m_dl_range[nh]; }
  bool IsWithinTotRange( int nh=0 ) const { return m_tot_range[nh]; }
  bool IsWithinDtRange( int nh=0 ) const { return m_dt_range[nh]; }

  bool ReCalcDC( bool applyRecursively=false ) { return CalcDCObservables(); }
  bool CheckRangeHits();

  void Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const;

protected:
  //  void ClearRegisteredHits( void );
};

//_____________________________________________________________________
inline std::ostream&
operator <<( std::ostream& ost, const DCHit& hit )
{
  hit.Print( "", ost );
  return ost;
}

#endif
