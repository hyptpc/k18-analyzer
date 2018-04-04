/**
 *  file: BH2Hit.hh
 *  date: 2017.04.10
 *
 */

#ifndef BH2_HIT_HH
#define BH2_HIT_HH

#include "DebugCounter.hh"
#include "Hodo1Hit.hh"
#include "Hodo2Hit.hh"

//______________________________________________________________________________
class BH2Hit : public Hodo2Hit
{
public:
  BH2Hit( HodoRawHit *rhit, int index=0 );
  ~BH2Hit( void );

private:
  double m_time_offset;

public:
  bool   Calculate( void );
  double UTime0( int n=0 )  const { return m_t1.at(n)  +m_time_offset; }
  double UCTime0( int n=0 ) const { return m_ct1.at(n) +m_time_offset; }
  double DTime0( int n=0 )  const { return m_t2.at(n)  +m_time_offset; }
  double DCTime0( int n=0 ) const { return m_ct2.at(n) +m_time_offset; }
  double Time0( int n=0 )   const { return 0.5*(m_t1.at(n)+m_t2.at(n)) +m_time_offset; }
  double CTime0( int n=0 )  const { return 0.5*(m_ct1.at(n)+m_ct2.at(n)) +m_time_offset; }

  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); };

};

#endif
