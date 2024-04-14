/**
 *  file: HodoParamMan.cc
 *  date: 2017.04.10
 *
 */

#include "HodoParamMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <std_ostream.hh>

#include "TString.h"
#include "DeleteUtility.hh"

namespace
{
  const std::string& class_name("HodoParamMan");
  const int SegMask  = 0x03FF;
  const int CidMask  = 0x00FF;
  const int PlidMask = 0x00FF;
  const int UdMask   = 0x0003;
  const int SegShift  =  0;
  const int CidShift  = 11;
  const int PlidShift = 19;
  const int UdShift   = 27;
}

//______________________________________________________________________________
HodoParamMan::HodoParamMan( void )
  : m_is_ready(false)
{
  m_file_name.clear();
}

//______________________________________________________________________________
HodoParamMan::~HodoParamMan( void )
{
  m_file_name.clear();
  ClearACont(); ClearTCont();
}

//______________________________________________________________________________
void
HodoParamMan::ClearACont( void )
{
  del::ClearMap( m_APContainer );
}

//______________________________________________________________________________
void
HodoParamMan::ClearTCont( void )
{
  del::ClearMap( m_TPContainer );
}

//______________________________________________________________________________
inline int
KEY( int cid, int pl, int seg, int ud )
{
  return ( ( (cid&CidMask) << CidShift  ) |
	   ( (pl&PlidMask) << PlidShift ) |
	   ( (seg&SegMask) << SegShift  ) |
	   ( (ud&UdMask)   << UdShift   ) );
}

//______________________________________________________________________________
bool
HodoParamMan::Initialize( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    hddaq::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }

  ClearACont(); ClearTCont();
  for(int i=0;i<m_file_name.size();i++){
    ReadFile(m_file_name.at(i) );
  }
  m_is_ready = true;
  return true;
}
//______________________________________________________________________________
bool
HodoParamMan::ReadFile( const std::string& file_name )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  std::ifstream ifs( file_name.c_str() );
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << func_name << " file open fail : "
		<< file_name << std::endl;
    return false;
  }

  int invalid=0;
  std::string line;
  while( ifs.good() && std::getline( ifs, line ) ){
    ++invalid;
    if( line[0]=='#' || line.empty() ) continue;
    if( TString(line).Contains("cid",TString::ECaseCompare::kIgnoreCase) ) continue;
    std::istringstream input_line( line );
    int cid=-1, plid=-1, seg=-1, at=-1, ud=-1;
    double p0=-9999., p1=-9999.;
    if( input_line >> cid >> plid >> seg >> at >> ud >> p0 >> p1 ){
      int key = KEY( cid, plid, seg, ud );
      if( at == kAdc ){
	HodoAParam *pre_param = m_APContainer[key];
	HodoAParam *param = new HodoAParam(p0,p1);
	m_APContainer[key] = param;
	if( pre_param ){
	  hddaq::cerr << func_name << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
      }
      else if( at == kTdc ){
	HodoTParam *pre_param = m_TPContainer[key];
	HodoTParam *param = new HodoTParam(p0,p1);
	m_TPContainer[key] = param;
	if( pre_param ){
	  hddaq::cerr << func_name << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
      }else{
	hddaq::cerr << func_name << ": Invalid Input" << std::endl
		    << " ===> (" << invalid << "a)" << line << " " << std::endl;
      } /* if(at) */
    }
    else {
      hddaq::cerr << func_name << ": Invalid Input" << std::endl
		  << " ===> (" << invalid << "b)" << line << " " << std::endl;
    } /* if( input_line >> ) */
  } /* while( std::getline ) */

  return true;
}

//______________________________________________________________________________
bool
HodoParamMan::Initialize( const std::string& file_name )
{
  m_file_name.clear();
  m_file_name.push_back(file_name);
  return Initialize();
}
//______________________________________________________________________________
bool
HodoParamMan::Initialize( const std::string& file_name1, const std::string& file_name2 )
{
  m_file_name.clear();
  m_file_name.push_back(file_name1);
  m_file_name.push_back(file_name2);
  return Initialize();
}

//______________________________________________________________________________
bool
HodoParamMan::GetTime( int cid, int plid, int seg, int ud, int tdc, double &time ) const
{
  HodoTParam* map = GetTmap( cid, plid, seg, ud );
  if(!map){
    hddaq::cerr << "#E "
		<< " something is wrong at HodoParamMan::GetTime( Cid="
		<< cid  << ", Plane=" << plid << ", Seg=" << seg  << ", UD="
		<< ud << ", TDC=" << tdc << ", " << "time"
		<< ")" << std::endl;    
    return false;
  }
  time = map->Time( tdc );
  return true;
}

//______________________________________________________________________________
bool
HodoParamMan::GetDe( int cid, int plid, int seg, int ud, int adc, double &de ) const
{
  HodoAParam* map = GetAmap( cid, plid, seg, ud );
  if(!map){
    hddaq::cerr << "#E "
		<< " something is wrong at HodoParamMan::GetDe( Cid="
		<< cid  << ", Plane=" << plid << ", Seg=" << seg  << ", UD="
		<< ud << ", ADC=" << adc << ", " << "de"
		<< ")" << std::endl;    

    return false;
  }
  de = map->DeltaE( adc );
  return true;
}

//______________________________________________________________________________
HodoTParam*
HodoParamMan::GetTmap( int cid, int plid, int seg, int ud ) const
{
  int key = KEY(cid,plid,seg,ud);
  HodoTParam* map     = 0;
  TIterator   itr     = m_TPContainer.find(key);
  TIterator   itr_end = m_TPContainer.end();
  if( itr!=itr_end ) map = itr->second;
  return map;
}

//______________________________________________________________________________
HodoAParam*
HodoParamMan::GetAmap( int cid, int plid, int seg, int ud ) const
{
  int key = KEY(cid,plid,seg,ud);
  HodoAParam* map     = 0;
  AIterator   itr     = m_APContainer.find(key);
  AIterator   itr_end = m_APContainer.end();
  if( itr!=itr_end ) map = itr->second;
  return map;
}
