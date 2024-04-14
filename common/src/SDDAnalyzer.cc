/**
 *  file: SDDAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "SDDAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "SDDHit.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("SDDAnalyzer");
}

#define Cluster 0

//______________________________________________________________________________
SDDAnalyzer::SDDAnalyzer( void )
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
SDDAnalyzer::~SDDAnalyzer( void )
{
  del::ClearContainer( m_SDDCont );
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
bool
SDDAnalyzer::DecodeRawHits( RawData *rawData )
{
  DecodeSDDHits( DetIdSDD,  m_SDDCont,  rawData );
  return true;  
}
//______________________________________________________________________________
bool
SDDAnalyzer::DecodeSDDHits( const int &detid, SDDHitContainer &m_Cont, RawData *rawData )
{
  del::ClearContainer( m_Cont );
  const SDDRHitContainer &cont = rawData->GetSDDRawHC(detid);
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    SDDRawHit *hit = cont[i];
    if( !hit ) continue;
    //    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    SDDHit *hp = new SDDHit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_Cont.push_back(hp);
    else
      delete hp;
  }//for(i)
  return true;
}

//______________________________________________________________________________
bool
SDDAnalyzer::ReCalcAll( void )
{
  return true;
}

