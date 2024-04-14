/**
 *  file: SDDParamMan.cc
 *  date: 2017.04.10
 *
 */

#include "SDDParamMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <std_ostream.hh>

#include "DeleteUtility.hh"

namespace
{
  const std::string& class_name("SDDParamMan");
  const int SegMask  = 0x03FF;
  const int CidMask  = 0x00FF;
  const int PlidMask = 0x00FF;
  const int UdMask   = 0x000F;
  const int SegShift  =  0;
  const int CidShift  = 11;
  const int PlidShift = 19;
  const int UdShift   = 27;
}

//______________________________________________________________________________
SDDParamMan::SDDParamMan( void )
  : m_is_ready(false),
    m_file_name1(""),
    m_file_name2("")
{}

//______________________________________________________________________________
SDDParamMan::~SDDParamMan( void )
{
  ClearACont(); ClearTCont(); ClearCCont();
}

//______________________________________________________________________________
void
SDDParamMan::ClearACont( void )
{
  del::ClearMap( m_APContainer );
}

//______________________________________________________________________________
void
SDDParamMan::ClearTCont( void )
{
  del::ClearMap( m_TPContainer );
}

//______________________________________________________________________________
void
SDDParamMan::ClearCCont( void )
{
  del::ClearMap( m_CPContainer );
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
SDDParamMan::Initialize( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    hddaq::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }

  ClearACont(); ClearTCont(); ClearCCont();

  if( ReadFile(m_file_name1) )   m_is_ready = true;
  if( ReadFile(m_file_name2) )   m_is_corr = true;

  return m_is_ready;
}
bool
SDDParamMan::ReadFile( const std::string& file_name )
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
    if( line[0]=='#' || line.empty() || TString(line).Contains("cid",TString::ECaseCompare::kIgnoreCase) ) continue;
    std::istringstream input_line( line );
    int cid=-1, plid=-1, seg=-1, at=-1, ud=-1;
    double p0=-9999., p1=-9999.;
    if( input_line >> cid >> plid >> seg >> at >> ud >> p0 >> p1 ){
      int key = KEY( cid, plid, seg, ud );
      if( at == kAdc ){
	SDDAParam *pre_param = m_APContainer[key];
	SDDAParam *param = new SDDAParam(p0,p1);
	m_APContainer[key] = param;
	if( pre_param ){
	  hddaq::cerr << func_name << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key  = " << key 
		      << " cid  = "<< cid
		      << " plid = "<< plid
		      << " seg  = "<< seg
		      << " ud   = "<< ud
		      <<std::endl;
	  delete pre_param;
	}
      }
      else if( at == kTdc ){
	SDDTParam *pre_param = m_TPContainer[key];
	SDDTParam *param = new SDDTParam(p0,p1);
	m_TPContainer[key] = param;
	if( pre_param ){
	  hddaq::cerr << func_name << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
      }
      else if( at == kPre || at == kPost ){
	SDDCParam *param = m_CPContainer[key];
	if( !param )
	  param = new SDDCParam();
	if( at==kPre )	param->SetParameter1(p0,p1);
	if( at==kPost )	param->SetParameter2(p0,p1);
	m_CPContainer[key] = param;
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
SDDParamMan::Initialize( const std::string& file_name )
{
  m_file_name1 = file_name;
  return Initialize();
}

//______________________________________________________________________________
bool
SDDParamMan::Initialize( const std::string& file_name1, const std::string& file_name2 )
{
  m_file_name1 = file_name1;
  m_file_name2 = file_name2;
  return Initialize();
}

//______________________________________________________________________________
bool
SDDParamMan::GetTime( int cid, int plid, int seg, int ud, int tdc, double &time ) const
{
  SDDTParam* map = GetTmap( cid, plid, seg, ud );
  if(!map) return false;
  time = map->Time( tdc );
  return true;
}

//______________________________________________________________________________
bool
SDDParamMan::GetDe( int cid, int plid, int seg, int ud, int adc, double &de ) const
{
  SDDAParam* map = GetAmap( cid, plid, seg, ud );
  if(!map) return false;
  de = map->DeltaE( adc );
  return true;
}

//______________________________________________________________________________
bool
SDDParamMan::GetCAdc( int cid, int plid, int seg, int ud, int adc0, int adc1, int adc2, double &corr ) const
{
  SDDCParam* map = GetCmap( cid, plid, seg, ud );
  if(!map) return false;
  corr = map->CAdc( adc0, adc1, adc2 );
  return true;
}

//______________________________________________________________________________
SDDTParam*
SDDParamMan::GetTmap( int cid, int plid, int seg, int ud ) const
{
  int key = KEY(cid,plid,seg,ud);
  SDDTParam* map     = 0;
  TIterator   itr     = m_TPContainer.find(key);
  TIterator   itr_end = m_TPContainer.end();
  if( itr!=itr_end ) map = itr->second;
  return map;
}

//______________________________________________________________________________
SDDAParam*
SDDParamMan::GetAmap( int cid, int plid, int seg, int ud ) const
{
  int key = KEY(cid,plid,seg,ud);
  SDDAParam* map     = 0;
  AIterator   itr     = m_APContainer.find(key);
  AIterator   itr_end = m_APContainer.end();
  if( itr!=itr_end ) map = itr->second;
  return map;
}

//______________________________________________________________________________
SDDCParam*
SDDParamMan::GetCmap( int cid, int plid, int seg, int ud ) const
{
  int key = KEY(cid,plid,seg,ud);
  SDDCParam* map     = 0;
  CIterator   itr     = m_CPContainer.find(key);
  CIterator   itr_end = m_CPContainer.end();
  if( itr!=itr_end ) map = itr->second;
  return map;
}
