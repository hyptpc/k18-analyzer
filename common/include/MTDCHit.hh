// -*- C++ -*-

#ifndef MTDC_HIT_HH
#define MTDC_HIT_HH

#include "MTDCRawHit.hh"

class RawData;

//_____________________________________________________________________________
class MTDCHit
{
public:
  explicit MTDCHit( MTDCRawHit *rhit, int index=0 );
  virtual ~MTDCHit();

private:
  MTDCHit( const MTDCHit& );
  MTDCHit& operator =( const MTDCHit& );

protected:
  MTDCRawHit           *m_raw;
  bool                 m_is_calculated;
  int                  m_multi_hit;
  std::vector<double>  m_t;
  std::vector<double>  m_ct;
  int                  m_index;

public:
  const MTDCRawHit* GetRawHit() const { return m_raw; }
  bool   Calculate();
  bool   IsCalculated() const { return m_is_calculated; }
  double GetT( int n=0 )   const { return m_t.at(n); }
  double GetCT( int n=0 )  const { return m_ct.at(n); }

  double Time( int n=0 )   const { return GetT(n); }
  double CTime( int n=0 )  const { return GetCT(n); }
  int         GetIndex(void)        const { return m_index; }

  int DetectorId() const { return m_raw->DetectorId(); }
  int UnitId()    const { return m_raw->UnitId(); }
  int SegmentId()  const { return m_raw->SegmentId(); }

  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); }

};

#endif
