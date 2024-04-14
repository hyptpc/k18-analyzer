// -*- C++ -*-

#ifndef HODO2_HIT_HH
#define HODO2_HIT_HH

#include <cmath>
#include <cstddef>

#include "HodoRawHit.hh"
#include "ThreeVector.hh"

class DetectorHit;
class RawData;

//_____________________________________________________________________________
class Hodo2Hit
{
public:
  explicit Hodo2Hit( HodoRawHit *rhit, int index=0 );
  explicit Hodo2Hit( DetectorHit *mchit, int index=0 );
  virtual ~Hodo2Hit( void );
  bool operator<(const Hodo2Hit &another) const
  {
    return SegmentId() < another.SegmentId();
  };
private:
  Hodo2Hit( const Hodo2Hit& );
  Hodo2Hit& operator =( const Hodo2Hit& );

protected:
  HodoRawHit          *m_raw;
  bool                 m_is_calculated;
  int                  m_detector_id;
  int                  m_plane_id;
  int                  m_segment_id;
  int                  m_multi_hit;
  std::vector<double>  m_a1;
  std::vector<double>  m_a2;
  std::vector<double>  m_ca1;
  std::vector<double>  m_ca2;
  std::vector<double>  m_t1;
  std::vector<double>  m_t2;
  std::vector<double>  m_ct1;
  std::vector<double>  m_ct2;
  std::vector<double>  m_tm;
  std::vector<double>  m_tsub;
  std::vector<double>  m_ctm;
  std::vector<double>  m_ctsub;
  std::vector<double>  m_de;
  std::vector<double>  m_desum;
  std::vector<double>  m_cde;
  std::vector<double>  m_cdesum;
  int                  m_index;

public:
  HodoRawHit* GetRawHit( void )           { return m_raw; }
  int         DetectorId( void )    const { return m_detector_id; }
  int         PlaneId( void )       const { return m_plane_id;    }
  int         SegmentId( void )     const { return m_segment_id;  }
  bool        IsCalculated( void )  const { return m_is_calculated; }
  int         GetIndex(void)        const { return m_index; }
  int         GetNumOfHit( void )   const { return m_multi_hit; }

  double      GetAUp( int n=0 )     const { return m_a1.size()>n ? m_a1.at(n) : -999; }
  double      GetADown( int n=0 )   const { return m_a2.size()>n ? m_a2.at(n) : -999; }
  double      GetCAUp( int n=0 )    const { return m_ca1.size()>n ? m_ca1.at(n) : -999; }
  double      GetCADown( int n=0 )  const { return m_ca2.size()>n ? m_ca2.at(n) : -999; }
  double      GetTUp( int n=0 )     const { return m_t1.size()>n  ? m_t1.at(n)  : -999; }
  double      GetTDown( int n=0 )   const { return m_t2.size()>n  ? m_t2.at(n)  : -999; }
  double      GetCTUp( int n=0 )    const { return m_ct1.size()>n ? m_ct1.at(n) : -999; }
  double      GetCTDown( int n=0 )  const { return m_ct2.size()>n ? m_ct2.at(n) : -999; }
  virtual double GetA( int ud, int n=0) const { return ud ? GetADown(n) : GetAUp(n); }
  virtual double GetCA(int ud, int n=0) const { return ud ? GetCADown(n) : GetCAUp(n); }
  virtual double GetT( int ud, int n=0) const { return ud ? GetTDown(n) : GetTUp(n); }
  virtual double GetCT(int ud, int n=0) const { return ud ? GetCTDown(n) : GetCTUp(n); }
  double      GetALeft( int n=0 )   const { return GetAUp(n); }
  double      GetARight( int n=0 )  const { return GetADown(n); }
  double      GetCALeft( int n=0 )  const { return GetCAUp(n); }
  double      GetCARight( int n=0 ) const { return GetCADown(n); }
  double      GetTLeft( int n=0 )   const { return GetTUp(n); }
  double      GetTRight( int n=0 )  const { return GetTDown(n); }
  double      GetCTLeft( int n=0 )  const { return GetCTUp(n); }
  double      GetCTRight( int n=0 ) const { return GetCTDown(n); }

  void  SetA1(  double val ) { m_a1.push_back(val); }
  void  SetA2(  double val ) { m_a2.push_back(val); }
  void  SetCA1( double val ) { m_ca1.push_back(val); }
  void  SetCA2( double val ) { m_ca2.push_back(val); }
  void  SetT1(  double val ) { m_t1.push_back(val); }
  void  SetT2(  double val ) { m_t2.push_back(val); }
  void  SetCT1( double val ) { m_ct1.push_back(val); }
  void  SetCT2( double val ) { m_ct2.push_back(val); }
  void  SetA(  int ud, double val ) { if(ud==0) SetA1(val); else SetA2(val); }
  void  SetCA( int ud, double val ) { if(ud==0) SetCA1(val); else SetCA2(val); }
  void  SetT(  int ud, double val ) { if(ud==0) SetT1(val); else SetT2(val); }
  void  SetCT( int ud, double val ) { if(ud==0) SetCT1(val); else SetCT2(val); }
 //   double      MeanTime( int n=0 )   const { return 0.5*(m_t1.at(n)+m_t2.at(n)); }
  //   double      CMeanTime( int n=0 )  const { return 0.5*(m_ct1.at(n)+m_ct2.at(n)); }
  virtual double MeanTime( int n=0 )   const { return m_tm.at(n); }
  virtual double CMeanTime( int n=0 )  const { return m_ctm.at(n); }
  double         TimeDiff( int n=0 )   const { return m_tsub.at(n); }
  double         CTimeDiff( int n=0 )  const { return m_ctsub.at(n); }

  virtual double DeltaE( int n=0 )     const { return m_de.at(n); }
  virtual double DeltaCE( int n=0 )    const { return m_cde.at(n); }
  virtual double SumE( int n=0 )       const { return m_desum.at(n); }
  virtual double SumCE( int n=0 )      const { return m_cdesum.at(n); }
  //  double      TimeDiff( int n=0 )   const { return m_ct2.at(n) - m_ct1.at(n); }

  virtual void Print( int i=0, bool header=true );
  virtual void Print( void );
  virtual bool Calculate( void );
  virtual bool ReCalc( bool applyRecursively=false )
  { return Calculate(); }

};

#endif
