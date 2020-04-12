/**
 *  file: TPCRawHit.hh
 *  date: 2020.04.11
 *
 */

#ifndef TPC_RAW_HIT_HH
#define TPC_RAW_HIT_HH

#include <cstddef>
#include <string>
#include <vector>

#include <TVector3.h>


//______________________________________________________________________________
class TPCRawHit
{
public:
  TPCRawHit( int padid, double y, double charge );
  ~TPCRawHit( void );

private:
  int      m_pad_id;
  int      m_layer_id;
  int      m_row_id;
  double   m_charge;
  TVector3 m_pos;

public:
  int  PadId( void )		const { return m_pad_id; }
  int  LayerId( void )		const { return m_layer_id; }
  int  RowId( void )		const { return m_row_id; }
  TVector3  Position( void )	const { return m_pos; }
  double X( void )		const { return m_pos.X(); }
  double Y( void )		const { return m_pos.Y(); }
  double Z( void )		const { return m_pos.Z(); }
  double Charge( void )		const { return m_charge; }
  
  void Print( const std::string& arg="" ) const;
};

#endif
