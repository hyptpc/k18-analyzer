/**
 *  file: MTDCAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "MTDCAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "MTDCHit.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("MTDCAnalyzer");
}

#define Cluster 0

//______________________________________________________________________________
MTDCAnalyzer::MTDCAnalyzer( void )
{
  debug::ObjectCounter::increase(class_name);
  for(int i=0;i<32;i++) m_Flag[i]=false;
}

//______________________________________________________________________________
MTDCAnalyzer::~MTDCAnalyzer( void )
{
  del::ClearContainer( m_TrigCont );
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
bool
MTDCAnalyzer::DecodeRawHits( RawData *rawData )
{
  DecodeMTDCHits( DetIdTrigFlag, m_TrigCont,  rawData );
  return true;  
}
//______________________________________________________________________________
bool
MTDCAnalyzer::DecodeMTDCHits( const int &detid, MTDCHitContainer &m_Cont, RawData *rawData )
{
  del::ClearContainer( m_Cont );
  const MTDCRHitContainer &cont = rawData->GetMTDCRawHC(detid);
  int nh = cont.size();
  //  std::cout<<"DecodeMTDCHits  "<<detid<<"  "<<nh<<std::endl;
  for( int i=0; i<nh; ++i ){
    MTDCRawHit *hit = cont[i];
    if( !hit ) continue;
    //    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    MTDCHit *hp = new MTDCHit( hit );
    if( !hp ) continue;
    if( hp->Calculate() ){
      m_Cont.push_back(hp);
      if(detid==DetIdTrigFlag){
	m_Flag[hit->SegmentId()]=true;
      }
    }
    else
      delete hp;
  }//for(i)
  return true;
}

//______________________________________________________________________________
bool
MTDCAnalyzer::ReCalcAll( void )
{
  return true;
}

