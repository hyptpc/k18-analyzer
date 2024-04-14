// -*- C++ -*-
#include "BHTHit.hh"
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

#define PHC 1
#define DE 1

namespace
{
  const std::string& class_name("BHTHit");
  const HodoParamMan& gHodo = HodoParamMan::GetInstance();
  const HodoPHCMan&   gPHC  = HodoPHCMan::GetInstance();
}

//______________________________________________________________________________
BHTHit::BHTHit( HodoRawHit *rhit, int index )
  : Hodo2Hit(rhit,index)
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
BHTHit::~BHTHit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
Bool_t
BHTHit::Calculate( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  //  std::cout<<__FILE__<<"  "<<__LINE__<<std::endl;
  //  if( m_status ){
  if( m_is_calculated ){
    hddaq::cerr << func_name << " already calculated" << std::endl;
    return false;
  }

  // if( m_raw->GetNumOfTdcHits()!=2 )
  //   return false;
  if( !gHodo.IsReady() ){
    hddaq::cerr << func_name << " HodoParamMan must be initialized" << std::endl;
    return false;
  }
#if PHC
  if( !gPHC.IsReady() ){
    hddaq::cerr << func_name << " HodoPHCMan must be initialized" << std::endl;
    return false;
  }
#endif
  // Detector information
  Int_t cid  = m_raw->DetectorId();
  Int_t plid = m_raw->PlaneId();
  Int_t seg  = m_raw->SegmentId();
  
  for(int ud=0;ud<2;ud++){
    int ntdc=m_raw->GetSizeTdc(ud);
    int nadc=m_raw->GetSizeAdc(ud);
    if( ntdc==0 || nadc==0 ) return false; 
    for(int i=0;i<ntdc;i++){ 
      int tdc = m_raw->GetTdc(ud,i);
      if(nadc<i+1) break;
      int adc = m_raw->GetAdc(ud,i);
      if( tdc<0 || adc<0 ){
	if(i==0) return false;
	else break;
      }  
      int j=1; 
      while(tdc<adc){
	if(nadc<i+j+1) break;
	adc=m_raw->GetAdc(ud,i+j);
	j++;
      }
      double tot = tdc - adc;
      double time = 0.;
      if( !gHodo.GetTime( cid, plid, seg, ud, tdc, time )){
	hddaq::cerr << "#E " << func_name
		    << " something is wrong at GetTime("
		    << cid  << ", " << plid << ", " << seg  << ", "
		    << ud << ", " << tdc << ", time" 
		    << ")" << std::endl;
	return false;
      }
      double ctot = 0.;
      if( !gHodo.GetDe( cid, plid, seg, ud, tot, ctot )){
	hddaq::cerr << "#E " << func_name
		    << " something is wrong at GetDe("
		    << cid  << ", " << plid << ", " << seg  << ", "
		    << ud << ", " << tdc << ", time" 
		    << ")" << std::endl;
	return false;
      }
      SetT(ud,time);
      SetCT(ud,time);
      SetA(ud,tot);
      SetCA(ud,ctot);
    }
  }
  auto cittr1=m_ct1.begin();
  for( auto ittr1=m_t1.begin(); ittr1!=m_t1.end(); ++ittr1 ){    
    auto cittr2=m_ct2.begin();
    for( auto ittr2=m_t2.begin(); ittr2!=m_t2.end(); ++ittr2 ){
      if(TMath::Abs(*ittr1-*ittr2)<10){
	m_tm.push_back(0.5*(*ittr1+*ittr2));
	m_tsub.push_back(*ittr2-*ittr1);
	m_ctm.push_back(0.5*(*cittr1+*cittr2));
	m_ctsub.push_back(*cittr2-*cittr1);
	int index1=std::distance(m_t1.begin(), ittr1);
	int index2=std::distance(m_t2.begin(), ittr2);
	m_de.push_back( std::sqrt(std::abs(m_a1.at(index1)*m_a2.at(index2))));
	m_cde.push_back(std::sqrt(std::abs(m_ca1.at(index1)*m_ca2.at(index2))));
	m_desum.push_back( m_a1.at(index1) + m_a2.at(index2) );
	m_cdesum.push_back(m_ca1.at(index1) + m_ca2.at(index2) );
	m_index++;
	break;
      }
      cittr2++;
    }
    cittr1++;
  }
  m_is_calculated = true;
  return true;
}
