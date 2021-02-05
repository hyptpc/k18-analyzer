// -*- C++ -*-

#ifndef TPC_HIT_HH
#define TPC_HIT_HH

#include <cmath>
#include <vector>
#include <deque>
#include <string>
#include <numeric>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "ThreeVector.hh"
#include <TVector3.h>
#include "DCHit.hh"

// typedef std::deque<bool>    BoolVec;
// typedef std::vector<int>    IntVec;
// typedef std::vector<double> DoubleVec;

//class DCLTrackHit;
class TPCLTrackHit;
class DCHit;

//_____________________________________________________________________________
class TPCHit
{
public:
  static TString ClassName( void );
  TPCHit( Int_t layer, Int_t row );
  TPCHit( int padid, TVector3 pos, double charge);// for cluster hit
  ~TPCHit( void );

private:
  TPCHit( const TPCHit& );
  TPCHit& operator =( const TPCHit& );

protected:
  int       m_pad;
  int       m_layer;
  int       m_row;
  double    m_charge;
  TVector3  m_pos;
  int       m_is_good;
  int       m_is_calculated;

  double  m_wpos;
  double  m_angle;

  ///// for TPC(MWPC)
  int    m_cluster_size;
  double m_mrow;
  bool   m_tpc_flag;

  ///// for TPC
  int    m_hitnum;

  double m_resx;
  double m_resy;
  double m_resz;
  double m_res;

  bool m_belong_track;

  DCHit *m_hit_xz;
  DCHit *m_hit_yz;

  //  mutable std::vector <DCLTrackHit *> m_register_container;
  mutable std::vector <TPCLTrackHit *> m_register_container;

public:
  //bool CalcCFTObservables( void );
  //  bool CalcObservablesSimulation( double dlength);
  bool CalcTPCObservables( void );

  void SetPad( int pad )             { m_pad = pad;          }
  void SetLayer( int layer )         { m_layer = layer;      }
  void SetRow( int row )             { m_row  = row;         }
  void SetCharge( double charge)     { m_charge = charge;   }
  void SetPos ( TVector3 pos)     { m_pos = pos;  }
  void SetWirePosition( double wpos )      { m_wpos  = wpos;}
  void SetAngle( double angle )      { m_angle  = angle;}

  void SetClusterSize( int size )          { m_cluster_size = size; }
  void SetMRow( double mrow )             { m_mrow  = mrow; }
  void SetTPCFlag( bool flag )            { m_tpc_flag = flag; }

  void SetHitNum( int hitnum ) { m_hitnum = hitnum; }

  void SetResX( double resx ) { m_resx = resx; }
  void SetResY( double resy ) { m_resy = resy; }
  void SetResZ( double resz ) { m_resz = resz; }
  void SetRes( double res ) { m_res = res; }

  int GetPad ( void ) const { return m_pad; }
  int GetLayer( void ) const { return m_layer; }
  int GetRow ( void ) const { return m_row; }
  double GetCharge ( void ) const { return m_charge; }
  TVector3 GetPos ( void ) const { return m_pos; }
  double GetWirePosition ( void ) const { return m_wpos; }
  double GetX( void )		const { return m_pos.X(); }
  double GetY( void )		const { return m_pos.Y(); }
  double GetZ( void )		const { return m_pos.Z(); }

  int GetClusterSize ( void ) const { return m_cluster_size; }
  double GetMRow ( void ) const { return m_mrow; }
  bool GetTPCFlag ( void ) const { return m_tpc_flag; }

  int GetHitNum ( void ) const { return m_hitnum; }

  double GetResolutionX( void )  ;
  double GetResolutionY( void )  ;
  double GetResolutionZ( void )  ;
  double GetResolution( void )   ;

  bool IsGoodHit( void ) const { return m_is_good; }


  double GetTiltAngle( void ) const { return m_angle; }

  void JoinTrack( void ) { m_belong_track = true; }
  void QuitTrack( void ) { m_belong_track = false;}
  bool BelongToTrack( void ) const { return m_belong_track; }
  //  bool IsWithinRange( void ) const { return m_pair_cont.at(nh).dl_range; }

  void RegisterHits( TPCLTrackHit *hit ) const
  { m_register_container.push_back(hit); }

  void Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const;

  DCHit *GetHitXZ( void ) const {return m_hit_xz;}
  DCHit *GetHitYZ( void ) const {return m_hit_yz;}


protected:
  void ClearRegisteredHits( void );
};

//_____________________________________________________________________________
inline TString
TPCHit::ClassName( void )
{
  static TString s_name( "TPCHit" );
  return s_name;
}

//_____________________________________________________________________________
inline std::ostream&
operator <<( std::ostream& ost, const TPCHit& hit )
{
  hit.Print( "", ost );
  return ost;
}

#endif
