// -*- C++ -*-

#include "Hodo2Hit.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "ReslMapMan.hh"
#include "RawData.hh"
#include "RootData.hh"

#define PHC 1
#define DE 1
namespace
{
  const std::string& class_name("Hodo2Hit");
  const HodoParamMan& gHodo  = HodoParamMan::GetInstance();
  const HodoPHCMan&   gPHC   = HodoPHCMan::GetInstance();
  const ReslMapMan&   gResol = ReslMapMan::GetInstance();
  const double mc_cm=0.1;
}

//______________________________________________________________________________
Hodo2Hit::Hodo2Hit( HodoRawHit *rhit, int index)
  : m_raw(rhit), m_is_calculated(false),
    m_detector_id(rhit->DetectorId()),
    m_plane_id(   rhit->PlaneId()),
    m_segment_id( rhit->SegmentId()),
    m_index(index)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
Hodo2Hit::Hodo2Hit( DetectorHit *mchit, int index)
  : m_raw(0),
    m_is_calculated(false),
    m_detector_id(mchit->detectorID()),
    m_plane_id(   mchit->layerID()),
    m_segment_id( mchit->channelID()),
    m_index(index)
{
  debug::ObjectCounter::increase(class_name);
#if 0
  TVector3 hitposition;
  if( confMan->GetGeomMapManager()->GetGPos( cid, channel , hitposition ) ){
    hit.SetPos(hitposition);
    double lv=DBL_MIN; TVector3 gpos;
    confMan-> GetGeomMapManager()-> GetGPos(cid, channel, lv, gpos);
    double hitpos;
    if( cid==CID_CDH || cid==CID_IH ) hitpos=pos.Z();
    else if( cid==CID_BVC ) hitpos=pos.X();
    else hitpos=pos.Y();
    hitpos*=mc_cm;
    hit.SetHitPosition(hitpos);
    hit.SetHitPosition(hit.hitpos()+hit.dt()*lv);
    hit.SetTSub(hitpos/lv);
  }else{
    cout<<"!!!!! Hit Position not found !!!!!"<<endl;
    exit(0);
  }
#endif
  m_index++;
  double tresol,eresol;
  gResol.GetResolution(DetectorId(),PlaneId(),SegmentId(),tresol,eresol);
  m_ct1.push_back(mchit->time()+tresol);
  m_ct2.push_back(mchit->time()+tresol);
  m_ctm.push_back(mchit->time()+tresol);
  m_t1.push_back( mchit->time()+tresol);
  m_t2.push_back( mchit->time()+tresol);
  m_tm.push_back( mchit->time()+tresol);
  m_a1.push_back( mchit->de()  +eresol);
  m_a2.push_back( mchit->de()  +eresol);
  m_de.push_back( mchit->de()  +eresol);
  m_cde.push_back(mchit->de()  +eresol);
  m_tsub.push_back(0.);
  m_is_calculated=true;
}
//______________________________________________________________________________
Hodo2Hit::~Hodo2Hit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
bool
Hodo2Hit::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  // std::cout<<func_name<<"  start"<<std::endl;
  if( m_is_calculated ){
    hddaq::cout << func_name << " already calculated" << std::endl;
    return false;
  }
  // if( m_raw->GetNumOfTdcHits()!=2 )
  //   return false;

  int cid  = m_raw->DetectorId();
  int plid = m_raw->PlaneId();
  int seg  = m_raw->SegmentId();
  int adc1 = m_raw->GetAdc1(), adc2=m_raw->GetAdc2();
  //  std::cout<<cid<<"  "<<seg<<std::endl;
  // if( !gPHC.IsReady() ){
  //   hddaq::cout << func_name << " HodoPHCMan must be initialized" << std::endl;
  //   return false;
  // }
  // if( !gHodo.IsReady() ){
  //   hddaq::cout << func_name << " HodoParamMan must be initialized" << std::endl;
  //   return false;
  // }
#if DE
  double dE1 = 0., dE2 = 0.;
  if( adc1>=0 ){
    if( !gHodo.GetDeHighGain( cid, plid, seg, 0, adc1, dE1 ) ){
      return false;
    }
  }
  if( adc2>=0 ){
    if( !gHodo.GetDeHighGain( cid, plid, seg, 1, adc2, dE2 ) ){
      return false;
    }
  }
  m_a1.push_back( dE1 );
  m_a2.push_back( dE2 );
  m_de.push_back(  std::sqrt(std::abs(dE1*dE2)));
  m_cde.push_back( std::sqrt(std::abs(dE1*dE2)));
  m_desum.push_back(  dE1 + dE2 );
  m_cdesum.push_back( dE1 + dE2 );
#endif
  double time,ctime;
  for(int i=0;i<m_raw->GetSizeTdcUp();i++){
    int tdc = m_raw->GetTdc1(i);
    if( tdc<0 ) continue;
    if( !gHodo.GetTime( cid, plid, seg, 0, tdc, time ) ) return false;
    m_t1.push_back( time );
    gPHC.DoCorrection( cid, plid, seg, 0, time, dE1, ctime );
    m_ct1.push_back( ctime );
  }
  for(int i=0;i<m_raw->GetSizeTdcDown();i++){
    int tdc = m_raw->GetTdc2(i);
    if( tdc<0 ) continue;
    if( !gHodo.GetTime( cid, plid, seg, 1, tdc, time ) ) return false;
    m_t2.push_back( time );
    gPHC.DoCorrection(  cid, plid, seg, 1, time, dE2, ctime );
    m_ct2.push_back( ctime );
  }

#if UNIDAQ
  if(m_raw->GetSizeTdcUp()==1&&m_raw->GetSizeTdcDown()==1){
    m_index++;
    m_tm.push_back(0.5*(GetTUp()+GetTDown()));
    m_tsub.push_back(-GetTUp()+GetTDown());
    m_ctm.push_back(0.5*(GetCTUp()+GetCTDown()));
    m_ctsub.push_back(-GetCTUp()+GetCTDown());
  }
#else
  auto cittr1=m_ct1.begin();
  for( auto ittr1=m_t1.begin(); ittr1!=m_t1.end(); ++ittr1 ){
    auto cittr2=m_ct2.begin();
    for( auto ittr2=m_t2.begin(); ittr2!=m_t2.end(); ++ittr2 ){
      if(TMath::Abs(*ittr1-*ittr2)<10){
	m_index++;
	m_tm.push_back(0.5*(*ittr1+*ittr2));
	m_tsub.push_back(*ittr2-*ittr1);
	m_ctm.push_back(0.5*(*cittr1+*cittr2));
	m_ctsub.push_back(*cittr2-*cittr1);
      }
      cittr2++;
    }
    cittr1++;
  }
#endif
  //   std::cout<<func_name<<"  done"<<std::endl;
  m_is_calculated = true;
  return true;
}

//______________________________________________________________________________
void
Hodo2Hit::Print( int i, bool header )
{
  std::cout<<"--------"<<std::endl;
  if(header){
    std::cout<<"[DetId: "<<m_detector_id<<", SegId: "<<m_segment_id<<" ]"<<std::endl;
  }
  std::cout<<i<<" / "<<m_index<<std::endl;
  std::cout<<"MeanTime  = "<<MeanTime(i)<<std::endl;
  std::cout<<"CMeanTime = "<<CMeanTime(i)<<std::endl;
  if(m_detector_id==DetIdBHT){
    std::cout<<"DeltaE    = "<<DeltaE(i)<<std::endl;
    std::cout<<"DeltaCE   = "<<DeltaCE(i)<<std::endl;
  }else{
    std::cout<<"DeltaE    = "<<DeltaE()<<std::endl;
  }
}
//______________________________________________________________________________
void
Hodo2Hit::Print()
{
  std::cout<<"======================"<<std::endl;
  for(int i=0;i<m_index;i++){
    Print(i);
  }
}
