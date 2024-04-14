/**
 *  file: HodoRawHit.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO_RAW_HIT_HH
#define HODO_RAW_HIT_HH

#include <cstddef>
#include <vector>
#include <iostream>

//______________________________________________________________________________
class HodoRawHit
{
public:
  HodoRawHit( int detector_id, int plane_id, int segment_id );
  ~HodoRawHit( void );

private:
  int              m_detector_id;
  int              m_plane_id;
  int              m_segment_id;
  std::vector<int> m_adc1;
  std::vector<int> m_adc2;
  std::vector<int> m_tdc1;
  std::vector<int> m_tdc2;
  int              m_nhtdc;

public:
  void SetAdc1( int adc );
  void SetAdc2( int adc );
  void SetTdc1( int tdc );
  void SetTdc2( int tdc );
  void SetAdcUp( int adc )    { SetAdc1(adc); }
  void SetAdcLeft( int adc )  { SetAdc1(adc); }
  void SetAdcDown( int adc )  { SetAdc2(adc); }
  void SetAdcRight( int adc ) { SetAdc2(adc); }
  void SetTdcUp( int tdc )    { SetTdc1(tdc); }
  void SetTdcLeft( int tdc )  { SetTdc1(tdc); }
  void SetTdcDown( int tdc )  { SetTdc2(tdc); }
  void SetTdcRight( int tdc ) { SetTdc2(tdc); }
  void SetDetectorId(int val )    { m_detector_id=val; }
  void SetPlaneId(   int val )    { m_plane_id   =val; }
  void SetSegmentId( int val )    { m_segment_id =val; }
  int  DetectorId( void )      const { return m_detector_id; }
  int  PlaneId( void )         const { return m_plane_id;    }
  int  SegmentId( void )       const { return m_segment_id;  }
  // for Multi-hit method
  int  GetNumOfTdcHits( void )   const { return m_nhtdc;       }
  int  GetAdc1( int i=0 )        const { return m_adc1.size()>i ? m_adc1.at(i) : -1;  }
  int  GetAdc2( int i=0 )        const { return m_adc2.size()>i ? m_adc2.at(i) : -1;  }
  int  GetTdc1( int i=0 )        const { return m_tdc1.at(i);  }
  int  GetTdc2( int i=0 )        const { return m_tdc2.at(i);  }
  int  GetAdc( int ud, int i=0)  const { return ud ? GetAdc2(i) : GetAdc1(i); }
  int  GetTdc( int ud, int i=0)  const { return ud ? GetTdc2(i) : GetTdc1(i); }
  int  GetAdcUp( int i=0 )       const { return GetAdc1(i);    }
  int  GetAdcLeft( int i=0 )     const { return GetAdc1(i);    }
  int  GetAdcDown( int i=0 )     const { return GetAdc2(i);    }
  int  GetAdcRight( int i=0 )    const { return GetAdc2(i);    }
  int  GetTdcUp( int i=0 )       const { return GetTdc1(i);    }
  int  GetTdcLeft( int i=0 )     const { return GetTdc1(i);    }
  int  GetTdcDown( int i=0 )     const { return GetTdc2(i);    }
  int  GetTdcRight( int i=0 )    const { return GetTdc2(i);    }
  int  SizeAdc1( void ) const;
  int  SizeAdc2( void ) const;
  int  SizeTdc1( void ) const;
  int  SizeTdc2( void ) const;
  int  GetSizeAdcUp( void )    const { return SizeAdc1(); }
  int  GetSizeAdcLeft( void )  const { return SizeAdc1(); }
  int  GetSizeAdcDown( void )  const { return SizeAdc2(); }
  int  GetSizeAdcRight( void ) const { return SizeAdc2(); }
  int  GetSizeTdcUp( void )    const { return SizeTdc1(); }
  int  GetSizeTdcLeft( void )  const { return SizeTdc1(); }
  int  GetSizeTdcDown( void )  const { return SizeTdc2(); }
  int  GetSizeTdcRight( void ) const { return SizeTdc2(); }
  int  GetSizeAdc( int ud )    const { return ud ? SizeAdc2() : SizeAdc1(); }
  int  GetSizeTdc( int ud )    const { return ud ? SizeTdc2() : SizeTdc1(); }
  void Clear( void );
  void Print( const std::string& arg="" );
};

#endif
