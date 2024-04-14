// -*- C++ -*-

#include "HodoAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

#include "DebugCounter.hh"
#include "DeleteUtility.hh"

namespace
{
  const std::string& class_name("HodoAnalyzer");
  // const double MaxTimeDifBHD =  2.0;
  // const double MaxTimeDifT0 =  2.0;
  // const double MaxTimeDifT0new =  2.0;
}

#define Cluster 0
//______________________________________________________________________________
HodoAnalyzer::HodoAnalyzer( void )
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
HodoAnalyzer::~HodoAnalyzer( void )
{
  del::ClearContainer( m_BHDCont );
  del::ClearContainer( m_T0Cont );
  del::ClearContainer( m_T0newCont );
  del::ClearContainer( m_E0Cont );
  del::ClearContainer( m_DEFCont );
#ifdef E15
  del::ClearContainer( m_BVCCont );
  del::ClearContainer( m_CVCCont );
  del::ClearContainer( m_NCCont );
  del::ClearContainer( m_PCCont );
  del::ClearContainer( m_LBCont );
  del::ClearContainer( m_WVCCont );
  del::ClearContainer( m_BDCont );
  del::ClearContainer( m_BPDCont );
  del::ClearContainer( m_IHCont );
#elif E62
  del::ClearContainer( m_StartCont );
  del::ClearContainer( m_StopCont );
#elif E73
  //  del::ClearContainer( m_FingerCont );
  del::ClearContainer( m_Veto0Cont );
  del::ClearContainer( m_Veto1Cont );
  del::ClearContainer( m_BTCCont );
  del::ClearContainer( m_PbF2Cont );
#elif E73_2024
  del::ClearContainer( m_VetoCont );
  del::ClearContainer( m_BTCCont );
  del::ClearContainer( m_PbF2Cont );
  del::ClearContainer( m_PbGCont );
#elif T98
  del::ClearContainer( m_VetoCont );
  del::ClearContainer( m_BTCCont );
  del::ClearContainer( m_RCCont );
  del::ClearContainer( m_PbF2Cont );
  del::ClearContainer( m_PbGCont );
#endif
#ifdef CDS
  del::ClearContainer( m_CDHCont );
#endif
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
//______________________________________________________________________________
bool
HodoAnalyzer::DecodeRawHits( RawData *rawData )
{
#if T98
  DecodeBHTHits( DetIdBHT, m_BHDCont, rawData );
#elif E73_2024
  DecodeBHTHits( DetIdBHT, m_BHDCont, rawData );
#else
  DecodeHodoHits( DetIdBHD, m_BHDCont, rawData );
#endif
  DecodeHodoHits( DetIdT0,  m_T0Cont,  rawData );
  DecodeHodoHits( DetIdT0new, m_T0newCont, rawData );
  DecodeHodoHits( DetIdDEF, m_DEFCont, rawData );
#ifdef CDS
  DecodeHodoHits( DetIdCDH, m_CDHCont, rawData );
#endif
#ifdef E15
  DecodeHodoHits(  DetIdBVC, m_BVCCont, rawData );
  DecodeHodoHits(  DetIdCVC, m_CVCCont, rawData );
  DecodeHodoHits(  DetIdNC,  m_NCCont,  rawData );
  DecodeHodoHits(  DetIdPC,  m_PCCont,  rawData );
  DecodeHodoHits(  DetIdLB,  m_LBCont,  rawData );
  DecodeHodo1Hits( DetIdWVC, m_WVCCont, rawData );
  DecodeHodoHits(  DetIdBD,  m_BDCont,  rawData );
  DecodeHodoHits(  DetIdBPD, m_BPDCont, rawData );
  DecodeHodo1Hits( DetIdIH,  m_IHCont,  rawData );
#elif E57
  DecodeHodoHits( DetIdE0,  m_E0Cont,  rawData );
#elif E62
  DecodeHodoHits( DetIdE0,  m_E0Cont,  rawData );
  DecodeHodoHits( DetIdStart,  m_StartCont,  rawData );
  DecodeHodoHits( DetIdStop,  m_StopCont,  rawData );
#elif E73
  DecodeHodo1Hits( DetIdPbF2,  m_PbF2Cont, rawData );
  DecodeHodoHits(  DetIdVeto1, m_Veto1Cont,rawData );
  DecodeHodo1Hits( DetIdVeto0, m_Veto0Cont,rawData );
  DecodeHodo1Hits( DetIdBTC,   m_BTCCont,  rawData );
#elif E73_2024
  DecodeHodo1Hits( DetIdPbG,   m_PbGCont,  rawData );
  DecodeHodo1Hits( DetIdPbF2,  m_PbF2Cont, rawData );
  DecodeHodoHits(  DetIdVeto,  m_VetoCont, rawData );
  DecodeHodoHits(  DetIdBTC,   m_BTCCont,  rawData );
#elif T98
  DecodeHodo1Hits( DetIdPbG,   m_PbGCont,  rawData );
  DecodeHodo1Hits( DetIdPbF2,  m_PbF2Cont, rawData );
  DecodeHodoHits(  DetIdVeto,  m_VetoCont, rawData );
  DecodeHodoHits(  DetIdBTC,   m_BTCCont,  rawData );
  DecodeHodoHits(  DetIdRC,    m_RCCont,   rawData );
#endif
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeHodoHits(const int &detid, Hodo2HitContainer &m_Cont, RawData *rawData )
{
  del::ClearContainer( m_Cont );
  const HodoRHitContainer &cont = rawData->GetHodoRawHC(detid);
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    //    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    Hodo2Hit *hp = new Hodo2Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_Cont.push_back(hp);
    else
      delete hp;
  }//for(i)
  //  std::sort(m_Cont.begin(),m_Cont.end());
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeBHTHits(const int &detid, BHTHitContainer &m_Cont, RawData *rawData )
{
  del::ClearContainer( m_Cont );
  const HodoRHitContainer &cont = rawData->GetHodoRawHC(detid);
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    //    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    BHTHit *hp = new BHTHit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_Cont.push_back(hp);
    else
      delete hp;
  }//for(i)
#if 0
  std::cout<<"before:";
  for(auto itr : m_Cont){
    std::cout<<"  "<<itr->SegmentId();
  }
  std::cout<<std::endl;
  std::sort(m_Cont.begin(),m_Cont.end());
  std::cout<<"after:";
  for(auto itr : m_Cont){
    std::cout<<"  "<<itr->SegmentId();
  }
  std::cout<<std::endl;
#endif
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::DecodeHodo1Hits(const int &detid, Hodo2HitContainer &m_Cont, RawData *rawData )
{
  del::ClearContainer( m_Cont );
  const HodoRHitContainer &cont = rawData->GetHodoRawHC(detid);
  int nh = cont.size();
  for( int i=0; i<nh; ++i ){
    HodoRawHit *hit = cont[i];
    if( !hit ) continue;
    //    if( hit->GetTdcUp()<=0 || hit->GetTdcDown()<=0 ) continue;
    Hodo1Hit *hp = new Hodo1Hit( hit );
    if( !hp ) continue;
    if( hp->Calculate() )
      m_Cont.push_back(hp);
    else
      delete hp;
  }//for(i)
  //  std::sort(m_Cont.begin(),m_Cont.end());
  return true;
}

//______________________________________________________________________________
bool
HodoAnalyzer::ReCalcAll( void )
{
  return true;
}
