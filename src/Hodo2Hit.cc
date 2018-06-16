/**
 *  file: Hodo2Hit.cc
 *  date: 2017.04.10
 *
 */

#include "Hodo2Hit.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "DebugCounter.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#include <std_ostream.hh>

namespace
{
  const std::string& class_name("Hodo2Hit");
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  const HodoPHCMan&   gPHC = HodoPHCMan::GetInstance();
}

//______________________________________________________________________________
Hodo2Hit::Hodo2Hit( HodoRawHit *rhit, int index)
  : m_raw(rhit), m_is_calculated(false), m_index(index)
{
  debug::ObjectCounter::increase(class_name);
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

  if( m_is_calculated ){
    hddaq::cout << func_name << " already calculated" << std::endl;
    return false;
  }

  if( m_raw->GetNumOfTdcHits()!=2 )
    return false;

  int tdc1 = m_raw->GetTdc1(m_index);
  int tdc2 = m_raw->GetTdc2(m_index);
  if( tdc1<0 || tdc2<0 )
    return false;

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
  int adc1 = m_raw->GetAdc1(), adc2=m_raw->GetAdc2();

  double time1 = 0., time2 = 0.;
  if( !gHodo.GetTime( cid, plid, seg, 0, tdc1, time1 ) ||
      !gHodo.GetTime( cid, plid, seg, 1, tdc2, time2 ) ){
    hddaq::cerr << "#E " << func_name
		<< " something is wrong at GetTime("
		<< cid  << ", " << plid << ", " << seg  << ", "
		<< "U/D, " << tdc1 << "/" << tdc2 << ", " << "time"
		<< ")" << std::endl;
    return false;
  }
  m_t1.push_back( time1 );
  m_t2.push_back( time2 );

  double dE1 = 0., dE2 = 0.;
  if( adc1>=0 ){
    if( !gHodo.GetDe( cid, plid, seg, 0, adc1, dE1 ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetDe("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << "0, " << adc1 << ", " << dE1 << ")" << std::endl;
      return false;
    }
  }
  if( adc2>=0 ){
    if( !gHodo.GetDe( cid, plid, seg, 1, adc2, dE2 ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetDe("
		  << cid  << ", " << plid << ", " << seg  << ", "
		  << "0, " << adc2 << ", " << dE2 << ")" << std::endl;
      return false;
    }
  }

  m_a1.push_back( dE1 );
  m_a2.push_back( dE2 );

  double ctime1 = -999., ctime2 = -999.;
  gPHC.DoCorrection( cid, plid, seg, 0, time1, dE1, ctime1 );
  gPHC.DoCorrection( cid, plid, seg, 1, time2, dE2, ctime2 );

  m_ct1.push_back( ctime1 );
  m_ct2.push_back( ctime2 );

  m_is_calculated = true;
  return true;
}
