
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
  std::vector<int> m_adc_Hi1;
  std::vector<int> m_adc_Hi2;
  std::vector<int> m_adc_Low1;
  std::vector<int> m_adc_Low2;
  // leading
  std::vector<int> m_tdc1;
  std::vector<int> m_tdc2;
  //trailing
  std::vector<int> m_tdc_t1;
  std::vector<int> m_tdc_t2;
  bool             m_oftdc; // module TDC overflow
  int              m_nhtdc;

public:
  //** Set 1 or 2 **//
  void SetAdc1( int adc );
  void SetAdc2( int adc );
  //** AdcHi  **//
  void SetAdcHi1( int adc );
  void SetAdcHi2( int adc );
  //** AdcLow  **//
  void SetAdcLow1( int adc );
  void SetAdcLow2( int adc );
  //** Leading  **//
  void SetTdc1( int tdc );
  void SetTdc2( int tdc );
  //** Trailing  **//
  void SetTdcT1( int tdc );
  void SetTdcT2( int tdc );

  //** Set Up or Down **//
  void SetAdcUp( int adc )    { SetAdc1(adc); }
  void SetAdcDown( int adc )  { SetAdc2(adc); }
  //** AdcHi  **//
  void SetAdcHiUp( int adc )    { SetAdcHi1(adc); }
  //  void SetAdcLeft( int adc )  { SetAdc1(adc); }
  void SetAdcHiDown( int adc )  { SetAdcHi2(adc); }
  //  void SetAdcRight( int adc ) { SetAdc2(adc); }
  //** AdcLow  **//
  void SetAdcLowUp( int adc )    { SetAdcLow1(adc); }
  void SetAdcLowDown( int adc )  { SetAdcLow2(adc); }
  //** Leading  **//
  void SetTdcUp( int tdc )    { SetTdc1(tdc); }
  //  void SetTdcLeft( int tdc )  { SetTdc1(tdc); }
  void SetTdcDown( int tdc )  { SetTdc2(tdc); }
  //  void SetTdcRight( int tdc ) { SetTdc2(tdc); }
  //** Trailing  **//
  void SetTdcTUp( int tdc )    { SetTdcT1(tdc); }
  //  void SetTdcTLeft( int tdc )  { SetTdcT1(tdc); }
  void SetTdcTDown( int tdc )  { SetTdcT2(tdc); }
  //  void SetTdcTRight( int tdc ) { SetTdcT2(tdc); }

  int  DetectorId( void )      const { return m_detector_id; }
  int  PlaneId( void )         const { return m_plane_id;    }
  int  SegmentId( void )       const { return m_segment_id;  }
  // for Multi-hit method
  void SetTdcOverflow( int fl ) { m_oftdc = static_cast<bool>( fl ); }
  int  GetNumOfTdcHits( void )   const { return m_nhtdc;       }

  //** Get 1 or 2 **//
  int  GetAdc1( int i=0 )        const { return m_adc1.at(i);  }
  int  GetAdc2( int i=0 )        const { return m_adc2.at(i);  }
  //** AdcHi **//
  int  GetAdcHi1( int i=0 )        const { return m_adc_Hi1.at(i);  }
  int  GetAdcHi2( int i=0 )        const { return m_adc_Hi2.at(i);  }
  //** AdcLow **//
  int  GetAdcLow1( int i=0 )        const { return m_adc_Low1.at(i);  }
  int  GetAdcLow2( int i=0 )        const { return m_adc_Low2.at(i);  }
  //** Leading **//
  int  GetTdc1( int i=0 )        const { return m_tdc1.at(i);  }
  int  GetTdc2( int i=0 )        const { return m_tdc2.at(i);  }
  //** Trailing **//
  int  GetTdcT1( int i=0 )       const { return m_tdc_t1.at(i);  }
  int  GetTdcT2( int i=0 )       const { return m_tdc_t2.at(i);  }

  //** Get Up or Down **//
  int  GetAdcUp( int i=0 )       const { return GetAdc1(i);    }
  int  GetAdcDown( int i=0 )     const { return GetAdc2(i);    }
  //** AdcHi **//
  int  GetAdcHiUp( int i=0 )       const { return GetAdcHi1(i);    }
  //  int  GetAdcLeft( int i=0 )     const { return GetAdc1(i);    }
  int  GetAdcHiDown( int i=0 )     const { return GetAdcHi2(i);    }
  //  int  GetAdcRight( int i=0 )    const { return GetAdc2(i);    }
  //** AdcLow **//
  int  GetAdcLowUp( int i=0 )       const { return GetAdcLow1(i);    }
  int  GetAdcLowDown( int i=0 )     const { return GetAdcLow2(i);    }
  //** Leading **//
  int  GetTdcUp( int i=0 )       const { return GetTdc1(i);    }
  //  int  GetTdcLeft( int i=0 )     const { return GetTdc1(i);    }
  int  GetTdcDown( int i=0 )     const { return GetTdc2(i);    }
  //  int  GetTdcRight( int i=0 )    const { return GetTdc2(i);    }
  //** Trailing **//
  int  GetTdcTUp( int i=0 )      const { return GetTdcT1(i);    }
  //  int  GetTdcTLeft( int i=0 )    const { return GetTdcT1(i);    }
  int  GetTdcTDown( int i=0 )    const { return GetTdcT2(i);    }
  //  int  GetTdcTRight( int i=0 )   const { return GetTdcT2(i);    }

  //** Size 1 or 2 **//
  int  SizeAdc1( void ) const;
  int  SizeAdc2( void ) const;
  int  SizeAdcHi1( void ) const;
  int  SizeAdcHi2( void ) const;
  int  SizeAdcLow1( void ) const;
  int  SizeAdcLow2( void ) const;
  int  SizeTdc1( void ) const;
  int  SizeTdc2( void ) const;
  int  SizeTdcT1( void ) const;
  int  SizeTdcT2( void ) const;

  //** Size Up or Down **//
  int  GetSizeAdcUp( void )    const { return SizeAdc1(); }
  int  GetSizeAdcDown( void )  const { return SizeAdc2(); }
  //** AdcHi **//
  int  GetSizeAdcHiUp( void )    const { return SizeAdcHi1(); }
  //  int  GetSizeAdcLeft( void )  const { return SizeAdc1(); }
  int  GetSizeAdcHiDown( void )  const { return SizeAdcHi2(); }
  //  int  GetSizeAdcRight( void ) const { return SizeAdc2(); }
  //** AdcLow **//
  int  GetSizeAdcLowUp( void )    const { return SizeAdcLow1(); }
  int  GetSizeAdcLowDown( void )  const { return SizeAdcLow2(); }
  //** Leading **//
  int  GetSizeTdcUp( void )    const { return SizeTdc1(); }
  //  int  GetSizeTdcLeft( void )  const { return SizeTdc1(); }
  int  GetSizeTdcDown( void )  const { return SizeTdc2(); }
  //  int  GetSizeTdcRight( void ) const { return SizeTdc2(); }
  //** Trailing **//
  int  GetSizeTdcTUp( void )    const { return SizeTdcT1(); }
  //  int  GetSizeTdcTLeft( void )  const { return SizeTdcT1(); }
  int  GetSizeTdcTDown( void )  const { return SizeTdcT2(); }
  //  int  GetSizeTdcTRight( void ) const { return SizeTdcT2(); }

  bool IsTdcOverflow( void )   const { return m_oftdc; }
  void Clear( void );
  void Print( const std::string& arg="" );
};

#endif
