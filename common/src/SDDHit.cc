/**
 *  file: SDDHit.cc
 *  date: 2017.04.10
 *
 */

#include "SDDHit.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "HodoParamMan.hh"
#include "SDDParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

namespace
{
  const std::string& class_name("SDDHit");
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  const SDDParamMan&  gSDD = SDDParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
}

//______________________________________________________________________________
SDDHit::SDDHit(SDDRawHit *rhit, int index )
  : m_raw(rhit), m_is_calculated(false), m_index(index)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
SDDHit::~SDDHit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
bool
SDDHit::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_is_calculated ){
    hddaq::cout << func_name << " already calculated" << std::endl;
    return false;
  }
  if( !gSDD.IsReady() ){
    hddaq::cout << func_name << " SDDParamMan must be initialized" << std::endl;
    return false;
  }
  // if( !gPHC.IsReady() ){
  //   hddaq::cout << func_name << " HodoPHCMan must be initialized" << std::endl;
  //   return false;
  // }

  int cid  = m_raw->DetectorId();
  int portid  = m_raw->PortId();
  int uid  = m_raw->UnitId();


  // hit informationca
  for( int i = 0; i<m_raw->SizeAdc();i++){
    int adc = m_raw->GetAdc(i);
    int adc1=-1,adc2=-1;
    if(i>0) adc1 = m_raw->GetAdc(i-1);
    if(i<7) adc2 = m_raw->GetAdc(i+1);
    double dE = -999.;
    double CdE = -999.;
    double cadc = 0.;
    if( adc>=0 ){
      if(!gSDD.GetDe( cid, portid, uid, i, adc, dE ) ){
	hddaq::cerr << "#E " << func_name
		    << " something is wrong at GetDe("
		    << cid  << ", " << portid << ", " << uid << ", "
		    << i  << ", " << adc  << ", " << dE  << ")" << std::endl;
	return false;
      }
      if(gSDD.IsCorr() && !gSDD.GetCAdc( cid, portid, uid, i, adc, adc1, adc2, cadc ) ){
	hddaq::cerr << "#E " << func_name
		    << " something is wrong at GetCAdc("
		    << cid  << ", " << portid << ", " << uid << ", "
		    << i  << ", " << adc  << ", " << adc1 <<". "<< adc2 << ")" << std::endl;
	return false;
      }
      if(!gSDD.GetDe( cid, portid, uid, i, cadc, CdE ) ){
	hddaq::cerr << "#E " << func_name
		    << " something is wrong at GetDe("
		    << cid  << ", " << portid << ", " << uid << ", "
		    << i  << ", " << cadc  << ", " << CdE  << ")" << std::endl;
	return false;
      }
    }
    m_a.push_back(dE);
    m_ca.push_back(CdE);
    m_cadc.push_back(cadc);
  }
  
  m_multi_hit = m_raw->SizeLeading();
  for( int m=0; m<m_multi_hit; ++m ){
    int    tdc  = -999;
    double time = -999.;
    tdc = m_raw->GetLeading(m);
    if( tdc<0 ) continue;

    if( !gSDD.GetTime( cid, portid, uid, 0, tdc, time ) ){
      hddaq::cerr << "#E " << func_name
		  << " something is wrong at GetTime("
		  << cid  << ", " << portid << ", " << uid  << ", "
		  << 0 << ", " << tdc  << ", " << time << ")" << std::endl;
      return false;
    }
    m_t.push_back(time);
    //    double ctime = time;
    //    gPHC.DoCorrection(cid, uid, seg, 0, time, dE, ctime );
    //    m_ct.push_back(ctime);
    m_index++;
  }
  
  m_is_calculated = true;
  return true;
}
