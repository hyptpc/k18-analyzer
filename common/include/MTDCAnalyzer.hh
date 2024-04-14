// -*- C++ -*-

#ifndef MTDC_ANALYZER_HH
#define MTDC_ANALYZER_HH

#include <vector>

#include "DetectorID.hh"
#include "RawData.hh"
#include "MTDCHit.hh"

class RawData;
class MTDCHit;

typedef std::vector <MTDCHit*> MTDCHitContainer;


//______________________________________________________________________________
class MTDCAnalyzer
{
public:
  MTDCAnalyzer(const RawData& rawData);
  ~MTDCAnalyzer();

  static MTDCAnalyzer& GetInstance();

private:
  MTDCAnalyzer( const MTDCAnalyzer& );
  MTDCAnalyzer& operator =( const MTDCAnalyzer& );

private:
  const RawData*     m_raw_data;
  MTDCHitContainer     m_TrigCont;
  bool                 m_Flag[32];

public:
  bool DecodeRawHits();
  bool DecodeMTDCHits( const int &detid, MTDCHitContainer &m_Cont);
  inline int GetNHits(  int detID )  const;
  inline const MTDCHit* GetHit( int detid, std::size_t i )  const;
  bool ReCalcAll();
  bool flag(const int &i) { return m_Flag[i]; }
};

inline int
MTDCAnalyzer::GetNHits( int detID  )  const
{
  switch(detID){
  case DetIdTrigFlag:
    return m_TrigCont.size();
  default:
    return 0;
  }
}
//______________________________________________________________________________
inline const MTDCHit*
MTDCAnalyzer::GetHit( int detID, std::size_t i ) const
{
  switch(detID){
  case DetIdTrigFlag:
    if( i<m_TrigCont.size() )      return m_TrigCont[i];
  default:
    return 0;
  }
}
//______________________________________________________________________________
#endif
