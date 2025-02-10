#ifndef TKORAW_DATA_HH
#define TKORAW_DATA_HH

#include "RawData.hh"
#include <vector>
//______________________________________________________________________________
class TKORawData : public RawData
{
public:
  TKORawData( void );
  ~TKORawData( void );

public:
  virtual void ClearAll( void );
  virtual bool DecodeHits( void );
  virtual bool DecodeHodoHits( void )  { return false; }
  virtual bool DecodeDCHits( void )    { return false; }
  virtual bool DecodeCDCHits( void )   { return false; }
  virtual bool DecodeMTDCHits( void )  { return false; }
  virtual bool DecodeCalibHits( void ) { return false; }
  virtual void PrintHodo( const int &detid ) const {}
};
#endif
