/**
 *  file: Hodo1Hit.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO1_HIT_HH
#define HODO1_HIT_HH

#include "Hodo2Hit.hh"

//______________________________________________________________________________
class Hodo1Hit : public Hodo2Hit
{
public:
  explicit Hodo1Hit( HodoRawHit *rhit, int index=0 );
  virtual ~Hodo1Hit();

private:
  Hodo1Hit( const Hodo1Hit& );
  Hodo1Hit& operator =( const Hodo1Hit& );

public:
  virtual double GetA( int n=0 )  const { return m_a1.at(n); }
  virtual double GetCA(int n=0 )  const { return m_ca1.at(n); }
  virtual double GetT( int n=0 )  const { return m_t1.at(n); }
  virtual double GetCT(int n=0 )  const { return m_ct1.at(n); }

  double Time( int n=0 )      const { return GetT(n); }
  double CTime( int n=0 )     const { return GetCT(n); }
  double DeltaE( int n=0 )    const { return GetA(n); }
  double DeltaCE( int n=0 )   const { return GetCA(n); }
  double SumE( int n=0 )      const { return GetA(n); }
  double SumCE( int n=0 )     const { return GetCA(n); }
  double MeanTime( int n=0 )  const { return GetT(n); }
  double CMeanTime( int n=0 ) const { return GetCT(n); }

  bool Calculate ( void );
  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); }

};

#endif
