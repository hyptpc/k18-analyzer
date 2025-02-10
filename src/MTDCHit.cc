/**
 *  file: MTDCHit.cc
 *  date: 2017.04.10
 *
 */

#include "MTDCHit.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "HodoParamMan.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("MTDCHit");
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
}

//______________________________________________________________________________
MTDCHit::MTDCHit(MTDCRawHit *rhit, int index )
  : m_raw(rhit), m_is_calculated(false), m_index(index)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
MTDCHit::~MTDCHit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
bool
MTDCHit::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_calculated ){
    hddaq::cout << func_name << " already calculated" << std::endl;
    return false;
  }

  if( !gHodo.IsReady() ){
    hddaq::cout << func_name << " HodoParamMan must be initialized" << std::endl;
    return false;
  }

  int cid  = m_raw->DetectorId();
  int uid = m_raw->UnitId();
  int seg  = m_raw->SegmentId();

  // hit information
  m_multi_hit = m_raw->SizeLeading();
  for( int m=0; m<m_multi_hit; ++m ){
    int    tdc  = -999;
    double time = -999.;
    tdc = m_raw->GetLeading(m);
    if( tdc<0 ) continue;
    if( !gHodo.GetTime( cid, uid, seg, 0, tdc, time ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetTime("
		  << cid  << ", " << uid << ", " << seg  << ", "
		  << 0 << ", " << tdc  << ", " << time << ")" << std::endl;
      return false;
    }
    m_t.push_back(time);
    m_index++;
  }
  
  m_is_calculated = true;
  return true;
}
