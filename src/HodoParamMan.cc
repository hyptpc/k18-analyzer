// -*- C++ -*-

#include "HodoParamMan.hh"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "FuncName.hh"

namespace
{
const Int_t SegMask  = 0x03FF;
const Int_t CidMask  = 0x00FF;
const Int_t PlidMask = 0x00FF;
const Int_t UdMask   = 0x0003;
const Int_t SegShift  =  0;
const Int_t CidShift  = 11;
const Int_t PlidShift = 19;
const Int_t UdShift   = 27;
}

//______________________________________________________________________________
HodoParamMan::HodoParamMan( void )
  : m_is_ready( false ),
    m_file_name()
{}

//______________________________________________________________________________
HodoParamMan::~HodoParamMan( void )
{
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
inline Int_t
KEY( Int_t cid, Int_t pl, Int_t seg, Int_t ud )
{
  return ( ( (cid&CidMask) << CidShift  ) |
	   ( (pl&PlidMask) << PlidShift ) |
	   ( (seg&SegMask) << SegShift  ) |
	   ( (ud&UdMask)   << UdShift   ) );
}

//______________________________________________________________________________
Bool_t
HodoParamMan::Initialize( void )
{
  if( m_is_ready ){
    hddaq::cerr << "#W " << FUNC_NAME
		<< " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    hddaq::cerr << FUNC_NAME << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  ClearACont(); ClearTCont();

  Int_t invalid=0;
  TString line;
  while( ifs.good() && line.ReadLine( ifs ) ){
    ++invalid;
    if( line.IsNull() || line[0]=='#' ) continue;
    std::istringstream input_line( line.Data() );
    Int_t cid=-1, plid=-1, seg=-1, at=-1, ud=-1;
    Double_t p0=-9999., p1=-9999.;
    Double_t p2=-9999., p3=-9999., p4=-9999., p5=-9999.;
    if( input_line >> cid >> plid >> seg >> at >> ud >> p0 >> p1 ){
      Int_t key = KEY( cid, plid, seg, ud );
      if( at == kAdc ){
	HodoAParam* pre_param = m_APContainer[key];
	HodoAParam* param = new HodoAParam( p0, p1 );
	m_APContainer[key] = param;
	if( pre_param ){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
      }
      else if( at == kTdc ){
	HodoTParam *pre_param = m_TPContainer[key];
	HodoTParam *param = new HodoTParam( p0, p1 );
	m_TPContainer[key] = param;
	if( pre_param ){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " key = " << key << std::endl;
	  delete pre_param;
	}
      }else if(at == 3){// for fiber position correction
	if(input_line  >> p2 >> p3>> p4 >> p5 ){
	  HodoFParam *pre_param = m_FPContainer[key];
	  HodoFParam *param = new HodoFParam( p0, p1, p2, p3, p4, p5 );
	  m_FPContainer[key] = param;
	  if( pre_param ){
	    hddaq::cerr << FUNC_NAME << ": duplicated key "
			<< " following record is deleted." << std::endl
			<< " key = " << key << std::endl;
	    delete pre_param;
	  }
	}
      }else{
	hddaq::cerr << FUNC_NAME << ": Invalid Input" << std::endl
		    << " ===> (" << invalid << "a)" << line << " " << std::endl;
      } /* if(at) */
    }
    else {
      hddaq::cerr << FUNC_NAME << ": Invalid Input" << std::endl
		  << " ===> (" << invalid << "b)" << line << " " << std::endl;
    } /* if( input_line >> ) */
  } /* while( std::getline ) */

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
Bool_t
HodoParamMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//______________________________________________________________________________
Bool_t
HodoParamMan::GetTime( Int_t cid, Int_t plid, Int_t seg, Int_t ud, Int_t tdc,
                       Double_t &time ) const
{
  HodoTParam* map = GetTmap( cid, plid, seg, ud );
  if( !map ) return false;
  time = map->Time( tdc );
  return true;
}

//______________________________________________________________________________
Bool_t
HodoParamMan::GetDe( Int_t cid, Int_t plid, Int_t seg, Int_t ud, Int_t adc,
                     Double_t &de ) const
{
  HodoAParam* map = GetAmap( cid, plid, seg, ud );
  if( !map ) return false;
  de = map->DeltaE( adc );
  return true;
}

//______________________________________________________________________________
Double_t
HodoParamMan::GetP0( Int_t cid, Int_t plid, Int_t seg, Int_t ud ) const
{
  HodoAParam* map = GetAmap( cid, plid, seg, ud );
  if( !map ) return -1;
  Double_t p0 = map->Pedestal();
  return p0;
}

//______________________________________________________________________________
Double_t
HodoParamMan::GetP1( Int_t cid, Int_t plid, Int_t seg, Int_t ud ) const
{
  HodoAParam* map = GetAmap( cid, plid, seg, ud );
  if( !map ) return -1;
  Double_t p1 = map->Gain();
  return p1;
}

//______________________________________________________________________________
Double_t
HodoParamMan::GetPar( Int_t cid, Int_t plid, Int_t seg, Int_t ud, Int_t i ) const
{
  HodoFParam *map = GetFmap( cid, plid, seg, ud );
  if( !map ) return -1;
  Double_t par=0;
  if( i==0 ) par=map->par0();
  else if( i==1 ) par=map->par1();
  else if( i==2 ) par=map->par2();
  else if( i==3 ) par=map->par3();
  else if( i==4 ) par=map->par4();
  else if( i==5 ) par=map->par5();
  return par;
}

//______________________________________________________________________________
Double_t
HodoParamMan::GetOffset( Int_t cid, Int_t plid, Int_t seg, Int_t ud ) const
{
  HodoTParam* map = GetTmap( cid, plid, seg, ud );
  if( !map ) return -9999.;
  return map->Offset();
}

//______________________________________________________________________________
Double_t
HodoParamMan::GetGain( Int_t cid, Int_t plid, Int_t seg, Int_t ud ) const
{
  HodoTParam* map = GetTmap( cid, plid, seg, ud );
  if( !map ) return -9999.;
  return map->Gain();
}

//______________________________________________________________________________
HodoTParam*
HodoParamMan::GetTmap( Int_t cid, Int_t plid, Int_t seg, Int_t ud ) const
{
  Int_t key = KEY( cid, plid, seg, ud );
  TIterator itr = m_TPContainer.find( key );
  if( itr != m_TPContainer.end() )
    return itr->second;
  else
    return nullptr;
}

//______________________________________________________________________________
HodoAParam*
HodoParamMan::GetAmap( Int_t cid, Int_t plid, Int_t seg, Int_t ud ) const
{
  Int_t key = KEY( cid, plid, seg, ud );
  AIterator itr = m_APContainer.find( key );
  if( itr != m_APContainer.end() )
    return itr->second;
  else
    return nullptr;
}

//______________________________________________________________________________
HodoFParam*
HodoParamMan::GetFmap( Int_t cid, Int_t plid, Int_t seg, Int_t ud ) const
{
  Int_t key = KEY( cid, plid, seg, ud );
  FIterator itr = m_FPContainer.find( key );
  if( itr != m_FPContainer.end() )
    return itr->second;
  else
    return nullptr;
}
