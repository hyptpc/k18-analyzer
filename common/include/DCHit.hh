// -*- C++ -*-

#ifndef DC_HIT_HH
#define DC_HIT_HH

#include <cmath>
#include <vector>
#include <deque>
#include <string>
#include <numeric>

#include <TMath.h>
#include <TVector3.h>

#include <std_ostream.hh>

#include "DebugCounter.hh"

class DCRawHit;
class DetectorHit;

//_____________________________________________________________________________
class DCHit
{
public:
  static const TString& ClassName();
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
  int       m_detector_id;
  Int_t     m_plane; // by DIGIT
  int       m_layer; // by DCGEO
  int       m_wire;

  using data_t = std::vector<Double_t>;
  using flag_t = std::deque<Bool_t>;

  data_t m_tdc;
  data_t m_adc;
  data_t m_trailing;

  data_t m_drift_time;
  data_t m_drift_length;
  data_t m_trailing_time;
  data_t m_tot;
  flag_t m_belong_to_track;
  flag_t m_is_good;

  TVector3  m_wpos;
  TVector3  m_wdir;
  double    m_tilt_angle;
  double    m_rotation;
  int       m_xy;
  Double_t m_z;

  flag_t   m_belong_cluster;
  flag_t   m_dt_range;
  flag_t   m_tot_range;
  flag_t   m_dl_range;

  //  mutable std::vector <DCLTrackHit *> m_register_container;

public:
  bool CalcDCObservables(Double_t retiming=0);
  void SetDCData(Double_t dt=0, Double_t dl=0, Double_t tot=TMath::QuietNaN(),
                 Bool_t belong_to_track=false, Bool_t is_good=true);

  void SetDetId( int detid )              { m_detector_id = detid; }
  void SetLayer( int layer )              { m_layer = layer; }
  void SetWire( int wire )             { m_wire  = wire; }
  void SetTdcVal( int tdc );
  void SetAdcVal( int adc );
  void SetTdcTrailing(int tdc)            { m_trailing.push_back(tdc); }
  void SetDriftTime( double dt )          { m_drift_time.push_back(dt); }
  void ClearDriftTime( void )             { m_drift_time.clear(); }
  void SetDriftLength( double dl )        { m_drift_length.push_back(dl); }
  void SetDriftLength( int ith, double dl ) { m_drift_length[ith] = dl; }
  void ClearDriftLength( void )           { m_drift_length.clear();}
  void SetTiltAngle( double angleDegree ) { m_tilt_angle = angleDegree; }
  void SetRotationAngle( double angleDegree ) { m_rotation = angleDegree; }
  void SetTrailingTime( double t )        { m_trailing_time.push_back(t); }
  void SetWirePosition( TVector3 wpos )   { m_wpos     = wpos; }

  int GetDetId( void ) const { return m_detector_id; }
  int GetLayer( void ) const { return m_layer; }
  int GetWire( void )  const { return m_wire;  }
  int GetXY( void )    const { return m_xy; }

  int GetTdcSize( void ) const { return m_tdc.size(); }
  int GetDriftTimeSize( void ) const { return m_drift_time.size(); }
  int GetDriftLengthSize( void ) const { return m_drift_length.size(); }
  int GetTdcTrailingSize( void ) const { return m_trailing.size(); }
  int GetTdcTrailingTimeSize( void ) const { return m_trailing_time.size(); }
  int GetTOTSize( void ) const { return m_tot.size(); }

  int GetTdcVal( int nh=0 ) const { return m_tdc[nh]; }
  int GetTdcTrailing( int nh=0 ) const { return m_trailing[nh]; }
  double GetDriftTime( int nh=0 ) const { return m_drift_time[nh]; }
  double GetDriftLength( int nh=0 ) const { return m_drift_length[nh]; }
  double GetTrailingTime( int nh=0 ) const { return m_trailing_time[nh]; }
  double GetTOT( int nh=0 ) const { return m_tot[nh]; }

  double GetResolution( void ) const;
  int GetNHit(bool DT=false, bool TOT=false) const;
  int GetHitID(int i,bool DT=false, bool TOT=false) const;


  double GetTiltAngle( void ) const { return m_tilt_angle; }
  double GetRotationAngle( void ) const { return m_rotation; }
  TVector3 GetWirePosition( void ) const {    return m_wpos;  }
  TVector3 GetWireDirection( void ) const {    return m_wdir;  }

  void JoinTrack( int nh=0 ) { m_belong_to_track[nh] = true; }
  void QuitTrack( int nh=0 ) { m_belong_to_track[nh] = false; }
  bool BelongToTrack( int nh=0 ) const { return m_belong_to_track[nh]; }
  void JoinCluster( int nh=0 ) { m_belong_cluster[nh] = true; }
  void QuitCluster( int nh=0 ) { m_belong_cluster[nh] = false; }
  bool BelongToCluster( int nh=0 ) const { return m_belong_cluster[nh]; }
  bool IsWithinDlRange( int nh=0 ) const { return m_dl_range[nh]; }
  bool IsWithinTotRange( int nh=0 ) const { return m_tot_range[nh]; }
  bool IsWithinDtRange( int nh=0 ) const { return m_dt_range[nh]; }

  bool ReCalcDC( bool applyRecursively=false ) { return CalcDCObservables(); }
  bool CheckRangeHits();

  void Print(Option_t* arg="") const;

protected:
  //  void ClearRegisteredHits( void );
};

//_____________________________________________________________________
inline const TString&
DCHit::ClassName()
{
  static TString s_name("DCHit");
  return s_name;
}

//_____________________________________________________________________
inline std::ostream&
operator <<(std::ostream& ost, const DCHit& hit)
{
  hit.Print();
  return ost;
}

#endif
