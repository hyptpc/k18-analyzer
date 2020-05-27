/**
 *  file: TPCHit.hh
 *  date: 2020.04.11
 *
 */

#ifndef TPC_HIT_HH
#define TPC_HIT_HH

#include <cstddef>
#include <string>
#include <vector>

#include <TVector3.h>


//______________________________________________________________________________
class TPCHit
{
public:
  TPCHit( int padid, double y, double charge );
  ~TPCHit( void );

private:
  int      m_pad_id;
  int      m_layer_id;
  int      m_row_id;
  double   m_charge;
  TVector3 m_pos;
  int      m_is_good;

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

  bool CalcTPCObservables( void );
  bool IsGoodHit( void ) const { return m_is_good; }
};

#endif
