// DCCluster.cc
#include <map>

#include "DCCluster.hh"
#include "DCHit.hh"
#include "DCTimeCorrMan.hh"
#include "DetectorID.hh"
#include "DebugCounter.hh"

namespace{
  const std::string& class_name("DCCluster");
  const DCTimeCorrMan& gDCTC= DCTimeCorrMan::GetInstance();
}
DCCluster::DCCluster()
{
  debug::ObjectCounter::increase(class_name);
}
DCCluster::~DCCluster()
{
  debug::ObjectCounter::decrease(class_name);
}
void DCCluster::Calc( const bool &isMC )
{
  DCHit *hit1=m_hit[0];
  int nth1=m_nthhit[0];
  int cid=hit1->GetDetId();
  if(nhit()==3){
    DCHit *hit2=m_hit[1];
    DCHit *hit3=m_hit[2];
    m_pos=(hit1->GetWirePosition()+hit2->GetWirePosition()+hit3->GetWirePosition())*(1./nhit());
    m_dir=(hit1->GetWireDirection()+hit2->GetWireDirection()+hit3->GetWireDirection())*(1./nhit());
    m_timesub=0;
    m_ctime=-999;
  }else if(nhit()==2){
    DCHit *hit2=m_hit[1];
    int nth2=m_nthhit[1];
    m_pos=(hit1->GetWirePosition()+hit2->GetWirePosition())*(1./nhit());
    m_dir=(hit1->GetWireDirection()+hit2->GetWireDirection())*(1./nhit());
    m_time=(hit1->GetDriftTime(nth1)+hit2->GetDriftTime(nth2))*(1./nhit());
    m_timesub=(hit1->GetDriftTime(nth1)-hit2->GetDriftTime(nth2))/2.;
    //   if(hit1->wx()-hit2->wx()>0){
    //     hit1->SetLeftRight(1);
    //     hit2->SetLeftRight(0);
    //   }else{
    //     hit1->SetLeftRight(0);
    //     hit2->SetLeftRight(1);
    //   }
    int lay1=hit1->GetLayer();
    int lay2=hit2->GetLayer();
    if(cid==DetIdCDC||isMC){
      m_timesub=0;
      m_time=0;
      m_ctime=0;
    }else{
      int lay=std::min(lay1,lay2);
      m_ctime=gDCTC.CalcCValue(cid,lay,0,m_time,m_timesub);
    }
  }else{
    m_pos=hit1->GetWirePosition();
    m_dir=hit1->GetWireDirection();
    m_time=hit1->GetDriftTime(nth1);
    m_timesub=0;
    m_ctime=-999;
  }
  if(cid==DetIdCDC){
    m_pos=m_pos-m_pos.Z()/m_dir.Z()*m_dir;
  }
  //m_pos=mpos
}

