/**
 *  file: SDDAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef SDD_ANALYZER_HH
#define SDD_ANALYZER_HH

#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"
#include "SDDHit.hh"

class RawData;
class SDDHit;

typedef std::vector <SDDHit*> SDDHitContainer;


//______________________________________________________________________________
class SDDAnalyzer
{
public:
  SDDAnalyzer( void );
  ~SDDAnalyzer( void );

  static SDDAnalyzer& GetInstance( void );

private:
  SDDAnalyzer( const SDDAnalyzer& );
  SDDAnalyzer& operator =( const SDDAnalyzer& );

private:
  SDDHitContainer     m_SDDCont;

public:
  bool DecodeRawHits( RawData* rawData );
  bool DecodeSDDHits( const int &detid, SDDHitContainer &m_Cont, RawData *rawData );

  inline int GetNHits(  int detID )  const;

  inline SDDHit * GetHit( int detid, std::size_t i )  const;  
  bool ReCalcAll( void );

};

inline int
SDDAnalyzer::GetNHits( int detID  )  const
{
  switch(detID){
  case DetIdSDD:
    return m_SDDCont.size();
  default:
    return 0;
  }
}
//______________________________________________________________________________
inline SDDHit*
SDDAnalyzer::GetHit( int detID, std::size_t i ) const
{
  switch(detID){
  case DetIdSDD:
    if( i<m_SDDCont.size() )      return m_SDDCont[i];
  default:
    return 0;
  }  
}
//______________________________________________________________________________
#endif
