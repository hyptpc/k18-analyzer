/**
 *  file: Hodo1Hit.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO1_HIT_HH
#define HODO1_HIT_HH

#include "HodoRawHit.hh"

class RawData;

//______________________________________________________________________________
class Hodo1Hit
{
public:
  explicit Hodo1Hit( HodoRawHit *rhit, int index=0 );
  virtual ~Hodo1Hit();

private:
  Hodo1Hit( const Hodo1Hit& );
  Hodo1Hit& operator =( const Hodo1Hit& );

protected:
  HodoRawHit          *m_raw;
  bool                 m_is_calculated;
  int                  m_multi_hit_l;
  int                  m_multi_hit_t;
  std::vector<double>  m_a;
  std::vector<double>  m_t;
  std::vector<double>  m_ct;
  int                  m_index;

public:
  HodoRawHit* GetRawHit( void ) { return m_raw; }
  bool   Calculate( void );
  bool   IsCalculated( void ) const { return m_is_calculated; }
  int    GetNumOfHit(int sel=0) const { return sel==0 ? m_multi_hit_l : m_multi_hit_t; };
  double GetA( int n=0 )   const { return m_a.at(n); }
  double GetT( int n=0 )   const { return m_t.at(n); }
  double GetCT( int n=0 )  const { return m_ct.at(n); }

  double Time( int n=0 )   const { return GetT(n); }
  double CTime( int n=0 )  const { return GetCT(n); }
  double DeltaE( int n=0 ) const { return GetA(n); }

  double GetAUp( int n=0 )     const { return m_a.at(n); }
  double GetALeft( int n=0 )   const { return m_a.at(n); }
  double GetADown( int n=0 )   const { return m_a.at(n); }
  double GetARight( int n=0 )  const { return m_a.at(n); }

  double GetTUp( int n=0 )     const { return m_t.at(n); }
  double GetTLeft( int n=0 )   const { return m_t.at(n); }
  double GetTDown( int n=0 )   const { return m_t.at(n); }
  double GetTRight( int n=0 )  const { return m_t.at(n); }

  double GetCTUp( int n=0 )    const { return m_ct.at(n); }
  double GetCTLeft( int n=0 )  const { return m_ct.at(n); }
  double GetCTDown( int n=0 )  const { return m_ct.at(n); }
  double GetCTRight( int n=0 ) const { return m_ct.at(n); }

  double MeanTime( int n=0 )  const { return GetT(n); }
  double CMeanTime( int n=0 ) const { return GetCT(n); }

  int DetectorId( void ) const { return m_raw->DetectorId(); }
  int PlaneId( void )    const { return m_raw->PlaneId(); }
  int SegmentId( void )  const { return m_raw->SegmentId(); }

  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); }

};

#endif
