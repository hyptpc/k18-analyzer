/**
 *  file: SDDRawHit.hh
 *  date: 2017.04.10
 *
 */

#ifndef SDD_RAW_HIT_HH
#define SDD_RAW_HIT_HH

#include <cstddef>
#include <vector>
#include <iostream>

//______________________________________________________________________________
class SDDRawHit
{
public:
  SDDRawHit( int detector_id, int port_id, int unit_id );
  ~SDDRawHit( void );

private:
  int              m_detector_id;
  int              m_port_id;
  int              m_unit_id;
  std::vector<int> m_adc;
  std::vector<int> m_leading;
  std::vector<int> m_trailing;
  std::vector<int> m_resetl;
  std::vector<int> m_resett;

public:
  void SetAdc( int i, int adc );
  void SetLeading( int tdc );
  void SetTrailing( int tdc );
  void SetResetLeading( int tdc );
  void SetResetTrailing( int tdc );
  int  DetectorId( void )      const { return m_detector_id; }
  int  UnitId( void )          const { return m_unit_id;    }
  int  PortId( void )       const { return m_port_id;  }
  // for Multi-hit method
  int  GetAdc( int i=0 )        const { return m_adc.at(i);  }
  int  GetAdcSum( )   const { 
    int sum=0; 
    for(int i=0;i<SizeAdc();i++) sum+=m_adc.at(i); 
    return sum;
  }
  int  GetLeading( int i=0 )        const { return m_leading.at(i);  }
  int  GetTrailing( int i=0 )        const { return m_trailing.at(i);  }
  int  GetResetLeading( int i=0 )        const { return m_resetl.at(i);  }
  int  GetResetTrailing( int i=0 )        const { return m_resett.at(i);  }
  int  SizeAdc( void ) const;
  int  SizeLeading( void ) const;
  int  SizeTrailing( void ) const;
  int  GetSizeAdc( void )    const { return SizeAdc(); }
  int  GetSizeLeading( void )    const { return SizeLeading(); }
  int  GetSizeTrailing( void )  const { return SizeTrailing(); }
  int  SizeResetLeading( void ) const;
  int  SizeResetTrailing( void ) const;
  int  GetSizeResetLeading( void )    const { return SizeResetLeading(); }
  int  GetSizeResetTrailing( void )  const { return SizeResetTrailing(); }
  void Clear( void );
  void Print( const std::string& arg="" );
};

#endif
