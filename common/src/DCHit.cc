// -*- C++ -*-

#include "DCHit.hh"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "BLDCWireMapMan.hh"
#include "DCTdcCalibMan.hh"
#include "XTMapMan.hh"
#include "UserParamMan.hh"
#include "ReslMapMan.hh"
#include "DetectorID.hh"
#include "DebugCounter.hh"
#include "RootHelper.hh"
#include "RootData.hh"
#ifdef CDS
#include "CDCWireMapMan.hh"
#endif
#include "TRandom.h"

namespace
{
  const std::string& class_name("DCHit");
  const BLDCWireMapMan& gBLDC = BLDCWireMapMan::GetInstance();
#ifdef CDS
  const CDCWireMapMan&  gCDC  = CDCWireMapMan::GetInstance();
#endif
  const UserParamMan&   gUser = UserParamMan::GetInstance();
  const DCTdcCalibMan&  gTdc  = DCTdcCalibMan::GetInstance();
  const XTMapMan&       gXt   = XTMapMan::GetInstance();
  const ReslMapMan&   gResol = ReslMapMan::GetInstance();
  //  const bool SelectTDC1st  = true;
  TVector3 DEFV(-999,-999,-999);
  const double mc_cm=0.1; //mm->cm
}

//______________________________________________________________________________
DCHit::DCHit( void )
  : m_layer(-1), m_wire(-1),
    m_wpos(DEFV), m_wdir(DEFV), m_tilt(0.), m_rotation(0.)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int layer )
  : m_layer( layer ), m_wire(-1),
    m_wpos(DEFV), m_wdir(DEFV), m_tilt(0.), m_rotation(0.)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int layer, int wire )
  : m_layer(layer), m_wire(wire),
    m_wpos(DEFV), m_wdir(DEFV), m_tilt(0.), m_rotation(0.)
{
  debug::ObjectCounter::increase(class_name);
}
//______________________________________________________________________________
DCHit::DCHit( int cid,int layer, int wire )
  : m_cid(cid), m_layer(layer), m_wire(wire),
    m_wpos(DEFV), m_wdir(DEFV), m_tilt(0.), m_rotation(0.)
{
  debug::ObjectCounter::increase(class_name);
}
//______________________________________________________________________________
DCHit::DCHit( DetectorHit *mchit)
  : m_cid(  mchit->detectorID()),
    m_layer(mchit->layerID()   ),
    m_wire( mchit->channelID() ),
    m_wpos(DEFV), m_wdir(DEFV), m_tilt(0.), m_rotation(0.)
{
  debug::ObjectCounter::increase(class_name);
  double tresol,eresol;
  if( m_cid==DetIdCDC ){
#ifdef CDS
    if( !gCDC.IsReady()){
      std::cout<<" CDC wire map is not initialized"<<std::endl;
      return;
    }
    m_wpos  = gCDC.GetWirePos( m_layer,m_wire );
    m_wdir  = gCDC.GetWireDir( m_layer,m_wire );
#endif
  }else if( m_cid==DetIdVFT ){
    m_wpos  = mchit->pos()*mc_cm;
    if(m_wire<32 || (m_wire>63&&m_wire<96) ){
      TVector3 aaa(0,1,1);
      aaa.RotateZ(m_wpos.Phi());
      m_wdir=aaa;
    }else{
      TVector3 aaa(0,-1,1);
      aaa.RotateZ(m_wpos.Phi());
      m_wdir=aaa;
    }
  }else{
    if( !gBLDC.IsReady()) return;
    m_wpos  = gBLDC.CalcWirePosition(  m_cid, m_layer,m_wire );
    m_wdir  = gBLDC.CalcWireDirection( m_cid, m_layer,m_wire );
    m_tilt = gBLDC.GetTiltAngle( m_cid, m_layer );
    m_rotation = gBLDC.GetRotationAngle( m_cid, m_layer );
    m_xy= gBLDC.GetWireMap(m_cid,m_layer)->GetXY();
  }
  gResol.GetResolution(m_cid,m_layer,m_wire,tresol,eresol);
  m_dl.push_back( mchit->dx()*mc_cm + tresol ); //cm
  m_tot.push_back(mchit->de());
  m_belong_track.push_back(  false);
  m_belong_cluster.push_back(false);
  m_dl_range.push_back( true);
  m_tot_range.push_back(true);
  m_dt_range.push_back( true);
  m_dt.push_back(0);
}
//______________________________________________________________________________
DCHit::~DCHit( void )
{
  //  ClearRegisteredHits();
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
DCHit::SetTdcVal( int tdc )
{
  m_tdc.push_back(tdc);
  m_belong_track.push_back(false);
  m_belong_cluster.push_back(false);
  m_dl_range.push_back(false);
  m_tot_range.push_back(false);
  m_dt_range.push_back(false);
}

//______________________________________________________________________________
bool
DCHit::CalcDCObservables( double retiming )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( !gXt.IsReady() ){
    std::cout<<" XT map is not initialized"<<std::endl;
    return false;
  }
  if( !gTdc.IsReady() ){
    std::cout<<" TdcCalib map is not initialized"<<std::endl;
    return false;
  }
  if( m_cid==DetIdCDC ){
#ifdef CDS
    if( !gCDC.IsReady()){
      std::cout<<" CDC wire map is not initialized"<<std::endl;
      return false;
    }
    m_wpos  = gCDC.GetWirePos(  m_layer,m_wire );
    m_wdir  = gCDC.GetWireDir(  m_layer,m_wire );
#endif
  }else{
    if( !gBLDC.IsReady()){
      std::cout<<" BLDC wire map is not initialized"<<std::endl;
      return false;
    }
    m_wpos  = gBLDC.CalcWirePosition( m_cid, m_layer,m_wire );
    m_wdir  = gBLDC.CalcWireDirection( m_cid, m_layer,m_wire );
    m_tilt = gBLDC.GetTiltAngle( m_cid, m_layer );
    m_rotation = gBLDC.GetRotationAngle( m_cid, m_layer );
    m_xy= gBLDC.GetWireMap(m_cid,m_layer)->GetXY();
  }
  if(m_trailing.size()>0&&m_tdc.size()>0&&m_tdc.front()<m_trailing.front()){
    m_trailing.erase(m_trailing.begin());
  }
  bool status = true;
  int nhtrailing = m_trailing.size();
  for ( int i=0; i<nhtrailing; ++i ) {
    double ctime;
    if( !gTdc.GetTime( m_cid, m_layer, m_wire, m_trailing[i], ctime ) )
      return false;
    double dtime=ctime - retiming;
    m_trailing_time.push_back( dtime );
  }
  int  nhtdc = m_tdc.size();
  for ( int i=0; i<nhtdc; ++i ) {
    double ctime;
    if( !gTdc.GetTime( m_cid, m_layer, m_wire, m_tdc[i], ctime ) )
      return false;
    double dtime=ctime - retiming;
    double dlength=gXt.CalcDriftLength( m_cid, m_layer, m_wire, dtime );
    // if( m_cid==DetIdBPC)
    //    std::cout<<i<<" / "<<nhtdc<<"  dt:"<<ctime<<" ,dl:"<<dlength<<std::endl;
    m_dt.push_back( dtime );
    m_dl.push_back( dlength );
    if(i<nhtrailing)   m_tot.push_back( m_trailing_time.at(i) - dtime);
    else m_tot.push_back(0);
  }

  CheckRangeHits();
  return status;
}
bool DCHit::CheckRangeHits(){
  std::string str1=Form("DCDt%d",m_cid);
  std::string str2=Form("DCTot%dL%d",m_cid,m_layer);
  std::string str3=Form("DCTot%d",m_cid);
  //  std::string str3=Form("DCDl%d",m_cid);
  for( int index=m_tdc.size()-1; index!=0-1; --index ){
    if(gUser.IsInRange(str1,m_dt[index])) m_dt_range[index]=true;
    if(gUser.Has(str2)){
      if(gUser.IsInRange(str2,m_tot[index])) m_tot_range[index]=true;
    }else{
      if(gUser.IsInRange(str3,m_tot[index])) m_tot_range[index]=true;
    }
    //    if(gUser.IsInRange(str3,m_dl[index])) m_dl_range[index]=true;
  }
  return true;
}
int DCHit::GetNHit(bool DT,bool TOT) const
{
  int nhit=0;
  for( int index=0; index<m_tdc.size(); ++index ){
    if(DT&&!m_dt_range[index]) continue;
    if(TOT&&!m_tot_range[index]) continue;
    nhit++;
  }
  return nhit;
}
int DCHit::GetHitID(int i,bool DT,bool TOT) const
{
  int nhit=0;
  for( int index=0; index<m_tdc.size(); ++index ){
    if(DT&&!m_dt_range[index]) continue;
    if(TOT&&!m_tot_range[index]) continue;
    if(nhit==i) return index;
    nhit++;
  }
  return -1;
}
//______________________________________________________________________________
double DCHit::GetResolution( void ) const
{
  return 0.1;
  //  return gGeom.GetResolution(m_layer);
}

//______________________________________________________________________________
void
DCHit::Print( const std::string& arg, std::ostream& ost ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  const int w = 16;
  ost << "#D " << func_name << " " << arg << std::endl
      << std::setw(w) << std::left << "layer" << m_layer << std::endl
      << std::setw(w) << std::left << "wire"  << m_wire  << std::endl
      << std::setw(w) << std::left << "wposX"  << m_wpos.X()  << std::endl
      << std::setw(w) << std::left << "wposY"  << m_wpos.Y()  << std::endl
      << std::setw(w) << std::left << "wposZ"  << m_wpos.Z()  << std::endl
      << std::setw(w) << std::left << "tilt" << m_tilt << std::endl
      << std::setw(w) << std::left << "rotation" << m_rotation << std::endl;

  ost << std::setw(w) << std::left << "tdc" << m_tdc.size() << " : ";
  std::copy( m_tdc.begin(), m_tdc.end(),
	     std::ostream_iterator<int>(ost, " ") );
  ost << std::endl;

  ost << std::endl << std::setw(w) << std::left
      << "trailing" << m_trailing.size() << " : ";
  std::copy(m_trailing.begin(), m_trailing.end(),
	    std::ostream_iterator<int>(ost, " "));
  ost << std::endl << std::setw(w) << std::left
      << "drift time" << m_dt.size() << " : ";
  std::copy(m_dt.begin(), m_dt.end(),
	    std::ostream_iterator<double>(ost, " "));
  ost << std::endl << std::setw(w) << std::left
      << "drift length" << m_dl.size() << " : ";
  std::copy(m_dl.begin(), m_dl.end(),
	    std::ostream_iterator<double>(ost, " "));
  ost << std::endl << std::setw(w) << std::left
      << "trailing time" << m_trailing_time.size() << " : ";
  std::copy(m_trailing_time.begin(), m_trailing_time.end(),
	    std::ostream_iterator<double>(ost, " "));
  ost << std::endl;
  ost << std::setw(w) << std::left
      << "belongTrack" << m_belong_track.size() << " : ";
  std::copy(m_belong_track.begin(), m_belong_track.end(),
	    std::ostream_iterator<bool>(ost, " "));
  // ost << std::endl << std::setw(w) << std::left
  //     << "dlRange" << m_dl_range.size() << " : ";
  // std::copy(m_dl_range.begin(), m_dl_range.end(),
  // 	      std::ostream_iterator<bool>(ost, " "));
  //   ost << std::endl;
}
