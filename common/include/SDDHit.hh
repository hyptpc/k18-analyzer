/**
 *  file: SDDHit.hh
 *  date: 2017.04.10
 *
 */

#ifndef SDD_HIT_HH
#define SDD_HIT_HH

#include "SDDRawHit.hh"

class RawData;

//______________________________________________________________________________
class SDDHit
{
public:
  explicit SDDHit( SDDRawHit *rhit, int index=0 );
  virtual ~SDDHit();

private:
  SDDHit( const SDDHit& );
  SDDHit& operator =( const SDDHit& );

protected:
  SDDRawHit           *m_raw;
  bool                 m_is_calculated;
  int                  m_multi_hit;
  std::vector<double>  m_cadc;
  std::vector<double>  m_a;
  std::vector<double>  m_ca;
  std::vector<double>  m_t;
  std::vector<double>  m_ct;
  int                  m_index;

public:
  SDDRawHit* GetRawHit( void ) { return m_raw; }
  bool   Calculate( void );
  bool   IsCalculated( void ) const { return m_is_calculated; }
  double GetA( int n=0 )   const { return m_a.at(n); }
  double GetCAdc( int n=0 )   const { return m_cadc.at(n); }
  double GetCA( int n=0 )   const { return m_ca.at(n); }
  int  GetASum( )   const { 
    int sum=0; 
    for(int i=0;i<m_raw->SizeAdc();i++) sum+=m_a.at(i); 
    return sum;
  }
  double GetT( int n=0 )   const { return m_t.at(n); }
  double GetCT( int n=0 )  const { return m_ct.at(n); }

  int         GetIndex(void)        const { return m_index; }

  double Time( int n=0 )   const { return GetT(n); }
  double CTime( int n=0 )  const { return GetCT(n); }
  double DeltaE( int n=0 ) const { return GetA(n); }
  double CAdc( int n=0 )   const { return GetCAdc(n); }
  double CDeltaE( int n=0 ) const { return GetCA(n); }

  int DetectorId( void ) const { return m_raw->DetectorId(); }
  int PortId( void )  const { return m_raw->PortId(); }
  int UnitId( void )    const { return m_raw->UnitId(); }
  
  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); }

};

#endif
