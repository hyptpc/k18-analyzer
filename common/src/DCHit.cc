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
#include "DCDriftParamMan.hh"
#include "DCRawHit.hh"
#include "DCTdcCalibMan.hh"
#include "FuncName.hh"
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
#include "ThreeVector.hh"

namespace
{
const std::string& class_name("DCHit");
const BLDCWireMapMan& gBLDC = BLDCWireMapMan::GetInstance();
#ifdef CDS
const CDCWireMapMan&  gCDC  = CDCWireMapMan::GetInstance();
#endif
const UserParamMan&   gUser = UserParamMan::GetInstance();
const DCTdcCalibMan&  gTdc  = DCTdcCalibMan::GetInstance();
const auto& gDrift = DCDriftParamMan::GetInstance();
const XTMapMan&       gXt   = XTMapMan::GetInstance();
const ReslMapMan&   gResol = ReslMapMan::GetInstance();
//  const bool SelectTDC1st  = true;
TVector3 DEFV(-999,-999,-999);
const double mc_cm=0.1; //mm->cm
}

//______________________________________________________________________________
DCHit::DCHit(const DCRawHit* rhit)
  : m_raw_hit(rhit),
    m_detector_id(rhit->DetectorId()),
    m_plane(rhit->PlaneId()),
    m_layer(),
    m_wire(rhit->WireId()),
    m_tdc(rhit->GetTdcArray()),
    m_adc(),
    m_trailing(rhit->GetTrailingArray()),
    m_drift_time(),
    m_drift_length(),
    m_tot(),
    m_belong_to_track(),
    m_is_good(),
    m_wpos(DEFV),
    m_wdir(DEFV),
    m_tilt_angle(),
    m_rotation()
{
  std::sort(m_tdc.begin(), m_tdc.end(), std::greater<Double_t>());
  std::sort(m_trailing.begin(), m_trailing.end(), std::greater<Double_t>());
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int layer )
  : m_layer( layer ), m_wire(-1),
    m_wpos(DEFV), m_wdir(DEFV), m_tilt_angle(0.), m_rotation(0.)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::DCHit( int layer, int wire )
  : m_layer(layer), m_wire(wire),
    m_wpos(DEFV), m_wdir(DEFV), m_tilt_angle(0.), m_rotation(0.)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
DCHit::~DCHit( void )
{
  // ClearRegisteredHits();
  debug::ObjectCounter::decrease(class_name);
}

//_____________________________________________________________________________
void
DCHit::SetDCData(Double_t dt, Double_t dl, Double_t tot,
                 Bool_t belong_to_track, Bool_t is_good)
{
  m_drift_time.push_back(dt);
  m_drift_length.push_back(dl);
  m_tot.push_back(tot);
  m_belong_to_track.push_back(belong_to_track);
  m_is_good.push_back(is_good);
}

//______________________________________________________________________________
bool
DCHit::CalcDCObservables( double retiming )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( !gTdc.IsReady() ){
    std::cout<<" DCTdcCalibMan is not initialized"<<std::endl;
    return false;
  }
  if( !gDrift.IsReady() ){
    std::cout<<" DCDriftParamMan is not initialized"<<std::endl;
    return false;
  }

  if( m_detector_id==DetIdCDC ){
#ifdef CDS
    if( !gCDC.IsReady()){
      std::cout<<" CDC wire map is not initialized"<<std::endl;
      return false;
    }
    m_wpos  = gCDC.GetWirePos(  m_layer,m_wire );
    m_wdir  = gCDC.GetWireDir(  m_layer,m_wire );
#endif
  }else{
    // if( !gBLDC.IsReady()){
    //   std::cout<<" BLDC wire map is not initialized"<<std::endl;
    //   return false;
    // }
    // m_wpos  = gBLDC.CalcWirePosition( m_detector_id, m_layer,m_wire );
    // m_wdir  = gBLDC.CalcWireDirection( m_detector_id, m_layer,m_wire );
    // m_tilt_angle = gBLDC.GetTiltAngle( m_detector_id, m_layer );
    // m_rotation = gBLDC.GetRotationAngle( m_detector_id, m_layer );
    // m_xy= gBLDC.GetWireMap(m_detector_id,m_layer)->GetXY();
  }

  std::sort(m_tdc.begin(), m_tdc.end(), std::greater<Double_t>());
  std::sort(m_trailing.begin(), m_trailing.end(), std::greater<Double_t>());

  data_t leading, trailing;
  for(Int_t il=0, nl=m_tdc.size(); il<nl; ++il){
    Double_t l = m_tdc[il];
    Double_t l_next = (il+1) != nl ? m_tdc[il+1] : DBL_MIN;
    Double_t buf = TMath::QuietNaN();
    for(const auto& t: m_trailing){
      if(l_next<t && t<l){
        buf = t;
        break;
      }
    }
    leading.push_back(l);
    trailing.push_back(buf);
    Double_t ctime = TMath::QuietNaN();
    gTdc.GetTime(m_detector_id, m_plane, m_wire, l, ctime);
    Double_t dt = TMath::QuietNaN();
    Double_t dl = TMath::QuietNaN();
    gDrift.CalcDrift(m_detector_id, m_plane, m_wire, ctime, dt, dl);

    Double_t ctime_trailing = TMath::QuietNaN();
    gTdc.GetTime(m_detector_id, m_plane, m_wire, buf, ctime_trailing);
    // Double_t tot = ctime - ctime_trailing;
    Double_t tot = l - buf;

    const Char_t* name = m_raw_hit->DetectorName().Data();
    Bool_t dt_is_good = gUser.IsInRange(Form("%s_DT", name), dt);
    Bool_t dl_is_good = dt_is_good;

    SetDCData(dt, dl, tot, false, dl_is_good);
  }
  m_tdc = leading;
  m_trailing = trailing;
  return true;


  // if(m_trailing.size()>0&&m_tdc.size()>0&&m_tdc.front()<m_trailing.front()){
  //   m_trailing.erase(m_trailing.begin());
  // }
  // bool status = true;
  // int nhtrailing = m_trailing.size();
  // for ( int i=0; i<nhtrailing; ++i ) {
  //   double ctime;
  //   if( !gTdc.GetTime( m_detector_id, m_layer, m_wire, m_trailing[i], ctime ) )
  //     return false;
  //   ctime -= retiming;
  //   m_trailing_time.push_back(ctime);
  // }
  // int  nhtdc = m_tdc.size();
  // for ( int i=0; i<nhtdc; ++i ) {
  //   const Char_t* n = m_raw_hit->DetectorName();
  //   if(!gUser.IsInRange(Form("%s_TDC", n), m_tdc[i]))
  //     return true;

  //   double ctime;
  //   if( !gTdc.GetTime( m_detector_id, m_layer, m_wire, m_tdc[i], ctime ) ){
  //     return false;
  //   }
  //   ctime -= retiming;
  //   double dlength=gXt.CalcDriftLength( m_detector_id, m_layer, m_wire, ctime );
  //   // if( m_detector_id==DetIdBPC)
  //   //    std::cout<<i<<" / "<<nhtdc<<"  dt:"<<ctime<<" ,dl:"<<dlength<<std::endl;
  //   Double_t tot = 0;
  //   if(i<nhtrailing)
  //     tot = m_trailing_time.at(i) - ctime;
  //   m_tot.push_back(tot);

  //   m_drift_time.push_back( ctime );
  //   m_drift_length.push_back( dlength );
  // }

  // CheckRangeHits();
  // return status;
}
bool DCHit::CheckRangeHits()
{
  std::string str1=Form("DCDt%d",m_detector_id);
  std::string str2=Form("DCTot%dL%d",m_detector_id,m_layer);
  std::string str3=Form("DCTot%d",m_detector_id);
  //  std::string str3=Form("DCDl%d",m_detector_id);
  const Char_t* name = m_raw_hit->DetectorName().Data();
  for(Int_t index=m_tdc.size()-1; index!=0-1; --index ){
    if(gUser.IsInRange(str1,m_drift_time[index])) m_dt_range[index]=true;
    m_tot_range[index] = gUser.IsInRange(Form("%s_TOT", name), m_tot[index]);
    // if(gUser.Has(str2)){
    //   if(gUser.IsInRange(str2,m_tot[index])) m_tot_range[index]=true;
    // }else{
    //   if(gUser.IsInRange(str3,m_tot[index])) m_tot_range[index]=true;
    // }
    //    if(gUser.IsInRange(str3,m_drift_length[index])) m_dl_range[index]=true;
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
DCHit::Print(Option_t* arg) const
{
  const Int_t w = 16;
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
              << std::setw(w) << std::left << "detector"
              << m_raw_hit->DetectorName() << std::endl
              << std::setw(w) << std::left << "plane" << m_plane << std::endl
              << std::setw(w) << std::left << "layer" << m_layer << std::endl
              << std::setw(w) << std::left << "wire"  << m_wire  << std::endl
              << std::setw(w) << std::left << "wpos"  << m_wpos  << std::endl
              << std::setw(w) << std::left << "angle" << m_tilt_angle << std::endl
              << std::setw(w) << std::left << "z"     << m_z     << std::endl;

  hddaq::cout << std::setw(w) << std::left << "tdc" << m_tdc.size() << " : ";
  std::copy(m_tdc.begin(), m_tdc.end(),
            std::ostream_iterator<Int_t>(hddaq::cout, " "));
  hddaq::cout << std::endl;
  hddaq::cout << std::setw(w) << std::left << "trailing" << m_trailing.size() << " : ";
  std::copy(m_trailing.begin(), m_trailing.end(),
            std::ostream_iterator<Int_t>(hddaq::cout, " "));
  hddaq::cout << std::endl;
  for(const auto& data_map: std::map<TString, data_t>
        {{"dt", m_drift_time},
         {"dl", m_drift_length},
         {"tot", m_tot},
        }){
    const auto& cont = data_map.second;
    hddaq::cout << std::setw(w) << std::left << data_map.first << cont.size()
                << " : ";
    std::copy(cont.begin(), cont.end(),
              std::ostream_iterator<Double_t>(hddaq::cout, " "));
    hddaq::cout << std::endl;
  }
  for(const auto& data_map: std::map<TString, flag_t>
        {{"belong_to_track", m_belong_to_track},
         {"is_good", m_is_good},
        }){
    const auto& cont = data_map.second;
    hddaq::cout << std::setw(w) << std::left << data_map.first << cont.size()
                << " : ";
    std::copy(cont.begin(), cont.end(),
              std::ostream_iterator<Bool_t>(hddaq::cout, " "));
    hddaq::cout << std::endl;
  }
}
