/**
 *  file: Hodo2Hit.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO2_HIT_HH
#define HODO2_HIT_HH

#include <cmath>
#include <cstddef>

#include "HodoRawHit.hh"
#include "ThreeVector.hh"

class RawData;

//______________________________________________________________________________
class Hodo2Hit
{
public:
  explicit Hodo2Hit( HodoRawHit *rhit, int index=0 );
  virtual ~Hodo2Hit( void );

private:
  Hodo2Hit( const Hodo2Hit& );
  Hodo2Hit& operator =( const Hodo2Hit& );

protected:
  HodoRawHit          *m_raw;
  bool                 m_is_calculated;
  int                  m_multi_hit;
  std::vector<double>  m_a1;
  std::vector<double>  m_a2;
  std::vector<double>  m_t1;
  std::vector<double>  m_t2;
  std::vector<double>  m_ct1;
  std::vector<double>  m_ct2;
  int                  m_index;

public:
  HodoRawHit* GetRawHit( void )           { return m_raw; }
  int         DetectorId( void )    const { return m_raw->DetectorId(); }
  int         PlaneId( void )       const { return m_raw->PlaneId(); }
  int         SegmentId( void )     const { return m_raw->SegmentId(); }
  bool        Calculate( void );
  bool        IsCalculated( void )  const { return m_is_calculated; }
  int         GetNumOfHit( void )   const { return m_multi_hit; }
  double      GetAUp( int n=0 )     const { return m_a1.at(n); }
  double      GetALeft( int n=0 )   const { return m_a1.at(n); }
  double      GetADown( int n=0 )   const { return m_a2.at(n); }
  double      GetARight( int n=0 )  const { return m_a2.at(n); }
  double      GetTUp( int n=0 )     const { return m_t1.at(n); }
  double      GetTLeft( int n=0 )   const { return m_t1.at(n); }
  double      GetTDown( int n=0 )   const { return m_t2.at(n); }
  double      GetTRight( int n=0 )  const { return m_t2.at(n); }
  double      GetCTUp( int n=0 )    const { return m_ct1.at(n); }
  double      GetCTLeft( int n=0 )  const { return m_ct1.at(n); }
  double      GetCTDown( int n=0 )  const { return m_ct2.at(n); }
  double      GetCTRight( int n=0 ) const { return m_ct2.at(n); }
  double      MeanTime( int n=0 )   const { return 0.5*(m_t1.at(n)+m_t2.at(n)); }
  double      CMeanTime( int n=0 )  const { return 0.5*(m_ct1.at(n)+m_ct2.at(n)); }
  double      DeltaE( int n=0 )     const { return std::sqrt(std::abs(m_a1.at(n)*m_a2.at(n))); }
  double      TimeDiff( int n=0 )   const { return m_ct2.at(n) - m_ct1.at(n); }

  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); }

};


#endif
