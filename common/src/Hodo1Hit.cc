/**
 *  file: Hodo1Hit.cc
 *  date: 2017.04.10
 *
 */

#include "Hodo1Hit.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("Hodo2Hit");
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
}

//______________________________________________________________________________
Hodo1Hit::Hodo1Hit( HodoRawHit *rhit, int index )
  : Hodo2Hit(rhit,index)
{
}

//______________________________________________________________________________
Hodo1Hit::~Hodo1Hit( void )
{
}

//______________________________________________________________________________
bool
Hodo1Hit::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_calculated ){
    hddaq::cout << func_name << " already calculated" << std::endl;
    return false;
  }

  // if( m_raw->GetNumOfTdcHits()!=1 )
  //   return false;

  if( !gHodo.IsReady() ){
    hddaq::cout << func_name << " HodoParamMan must be initialized" << std::endl;
    return false;
  }

  if( !gPHC.IsReady() ){
    hddaq::cout << func_name << " HodoPHCMan must be initialized" << std::endl;
    return false;
  }

  int cid  = m_raw->DetectorId();
  int plid = m_raw->PlaneId();
  int seg  = m_raw->SegmentId();

  // hit information
  m_multi_hit = m_raw->SizeTdc1();
  int UorD=0;
  int adc = m_raw->GetAdc1(0);
  if( 0 < m_raw->GetTdc2(0) ){
    UorD=1;
    m_multi_hit = m_raw->SizeTdc2();
    adc = m_raw->GetAdc2(0);
  }

  double dE = 0.;
  if( adc>=0 ){
    if( !gHodo.GetDe( cid, plid, seg, UorD, adc, dE ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetDe("
		  << cid  << ", " << plid << ", " << seg << ", "
		  << UorD << ", " << adc  << ", " << dE  << ")" << std::endl;
      return false;
    }
  }

  m_a1.push_back(dE);

  for( int m=0; m<m_multi_hit; ++m ){
    int    tdc  = -999;
    double time = -999.;
    if(0 == UorD){
      tdc = m_raw->GetTdc1(m);
    }else{
      tdc = m_raw->GetTdc2(m);
    }

    if( tdc<0 ) continue;

    if( !gHodo.GetTime( cid, plid, seg, UorD, tdc, time ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetTime("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << UorD << ", " << tdc  << ", " << time << ")" << std::endl;
      return false;
    }
    m_t1.push_back(time);
    m_index++;
    double ctime = time;
    gPHC.DoCorrection(cid, plid, seg, UorD, time, dE, ctime );
    m_ct1.push_back(ctime);
  }

  m_is_calculated = true;
  return true;
}
