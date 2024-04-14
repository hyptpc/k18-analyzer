/**
 *  file: MTDCRawHit.hh
 *  date: 2017.04.10
 *
 */

#ifndef MTDC_RAW_HIT_HH
#define MTDC_RAW_HIT_HH

#include <cstddef>
#include <vector>
#include <iostream>

//______________________________________________________________________________
class MTDCRawHit
{
public:
  MTDCRawHit( int detector_id, int unit_id, int segment_id );
  ~MTDCRawHit( void );

private:
  int              m_detector_id;
  int              m_unit_id;
  int              m_segment_id;
  std::vector<int> m_leading;
  std::vector<int> m_trailing;

public:
  void SetLeading( int tdc );
  void SetTrailing( int tdc );
  int  GetAdc( int i=0 ) const  {return -1; }
  int  DetectorId( void )      const { return m_detector_id; }
  int  UnitId( void )          const { return m_unit_id;    }
  int  SegmentId( void )       const { return m_segment_id;  }
  // for Multi-hit method
  int  GetLeading( int i=0 )        const { return m_leading.at(i);  }
  int  GetTrailing( int i=0 )        const { return m_trailing.at(i);  }
  int  SizeLeading( void ) const;
  int  SizeTrailing( void ) const;
  int  GetSizeLeading( void )    const { return SizeLeading(); }
  int  GetSizeTrailing( void )  const { return SizeTrailing(); }
  void Clear( void );
  void Print( const std::string& arg="" );
};

#endif
