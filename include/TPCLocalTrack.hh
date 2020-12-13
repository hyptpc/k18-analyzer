// -*- C++ -*-

#ifndef TPC_LOCAL_TRACK_HH
#define TPC_LOCAL_TRACK_HH

#include <vector>
#include <functional>

#include <std_ostream.hh>

#include "TMinuit.h"
#include "TVector3.h"
#include "ThreeVector.hh"
#include "DetectorID.hh"
#include "TPCHit.hh"
#include "TPCCluster.hh"

class TPCHit;
class TPCCluster;
class DCAnalyzer;

//______________________________________________________________________________
class TPCLocalTrack
{
public:
  explicit TPCLocalTrack( void );
  ~TPCLocalTrack( void );

private:
  TPCLocalTrack( const TPCLocalTrack & );
  TPCLocalTrack & operator =( const TPCLocalTrack & );
  TMinuit *minuit;

private:
  bool   m_is_fitted;     // flag of DoFit()
  bool   m_is_calculated; // flag of Calculate()
  std::vector<TPCHit*> m_hit_array;
  std::vector<TPCCluster*> m_cluster_array;

  double m_Ax;
  double m_Ay;
  double m_Au;
  double m_Av;
  // double m_Chix;
  // double m_Chiy;
  // double m_Chiu;
  // double m_Chiv;

  double m_Az;
  double m_Bz;

  double m_x0;
  double m_y0;
  double m_u0;
  double m_v0;
  double m_a;
  double m_b;
  double m_chisqr;
  bool   m_good_for_tracking;
  double m_n_iteration;
  // for SSD
  double m_de;
  

public:
  void         AddTPCHit( TPCHit *hit );
  void         AddTPCCluster( TPCCluster *cluster );
  void         Calculate( void );
  void         DeleteNullHit( void );
  bool         DoLinearFit( void );
  //  bool         DoHelixFit( void );
  bool         DoFit( void );
  int          GetNDF( void ) const;
  int          GetNHit( void ) const { return m_hit_array.size();  }
  TPCHit* GetHit( std::size_t nth ) const;
  bool         IsFitted( void ) const { return m_is_fitted; }
  bool         IsCalculated( void ) const { return m_is_calculated; }
  void         CalcChi2( void);


  void SetAx( double Ax ) { m_Ax = Ax; }
  void SetAy( double Ay ) { m_Ay = Ay; }
  void SetAu( double Au ) { m_Au = Au; }
  void SetAv( double Av ) { m_Av = Av; }
  void SetAz( double Az ){  m_Az = Az; }
  void SetBz( double Bz ){  m_Bz = Bz; }
  // void SetChix( double Chix ) { m_Chix = Chix; }
  // void SetChiy( double Chiy ) { m_Chiy = Chiy; }
  // void SetChiu( double Chiu ) { m_Chiu = Chiu; }
  // void SetChiv( double Chiv ) { m_Chiv = Chiv; }
  void SetDe( double de ) { m_de = de; }

  double GetX0( void ) const { return m_x0; }
  double GetY0( void ) const { return m_y0; }
  double GetU0( void ) const { return m_u0; }
  double GetV0( void ) const { return m_v0; }

  double GetAx( void ) const { return m_Ax; }
  double GetAy( void ) const { return m_Ay; }
  double GetAu( void ) const { return m_Au; }
  double GetAv( void ) const { return m_Av; }

  double GetAz ( void ) const { return m_Az;  }
  double GetBz ( void ) const { return m_Bz;  }


  double GetChiSquare( void ) const { return m_chisqr; }
  // double GetChiX( void ) const { return m_Chix; }
  // double GetChiY( void ) const { return m_Chiy; }
  // double GetChiU( void ) const { return m_Chiu; }
  // double GetChiV( void ) const { return m_Chiv; }
  double GetX( double z ) const { return m_x0+m_u0*z; }
  double GetY( double z ) const { return m_y0+m_v0*z; }
  double GetS( double z, double tilt ) const { return GetX(z)*std::cos(tilt)+GetY(z)*std::sin(tilt); }
  int    GetNIteration( void ) const { return m_n_iteration; }
  double GetTheta( void ) const;
  bool   GoodForTracking( void ) const { return m_good_for_tracking; }
  bool   GoodForTracking( bool status )
  { bool ret = m_good_for_tracking; m_good_for_tracking = status; return ret; }
  double GetDe( void ) const { return m_de; }
  void   Print( const std::string& arg="", std::ostream& ost=hddaq::cout ) const;


};


//______________________________________________________________________________
inline
std::ostream&
operator <<( std::ostream& ost,
	     const TPCLocalTrack& track )
{
  track.Print( "", ost );
  return ost;
}

//______________________________________________________________________________

struct TPCLTrackComp_Nhit
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, bool>
{
  bool operator()( const TPCLocalTrack * const p1,
                   const TPCLocalTrack * const p2 ) const
  {
    int n1=p1->GetNHit(), n2=p2->GetNHit();

    if( n1>=n2 )
      return true;
    else
      return false;
  }
};

//______________________________________________________________________________

struct TPCLTrackComp_Chisqr
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, bool>
{
  bool operator()( const TPCLocalTrack * const p1,
                   const TPCLocalTrack * const p2 ) const
  {
    double chi1=p1->GetChiSquare(),chi2=p2->GetChiSquare();

    if (chi1 <= chi2)
      return true;
    else
      return false;
  }
};

//______________________________________________________________________________
struct TPCLTrackComp // TODO
  : public std::binary_function <TPCLocalTrack *, TPCLocalTrack *, bool>
{
  bool operator()( const TPCLocalTrack * const p1,
		   const TPCLocalTrack * const p2 ) const
  {

    int n1=p1->GetNHit(), n2=p2->GetNHit();
    double chi1=p1->GetChiSquare(), chi2=p2->GetChiSquare();
    if( n1 > n2 ) return true;
    if( n2 > n1 ) return false;
    return ( chi1 <= chi2 );

  }
};

#endif
