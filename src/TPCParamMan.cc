// -*- C++ -*-

#include "TPCParamMan.hh"

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
const Int_t RowMask = 0xff;
const Int_t LayerMask = 0xff;
const Int_t RowShift = 0;
const Int_t LayerShift = 8;

inline Int_t
MakeKey( Int_t layer, Int_t row )
{
  return ( ( ( layer & LayerMask) << LayerShift ) |
	   ( ( row & RowMask) << RowShift ) );
}

}

//_____________________________________________________________________________
TPCParamMan::TPCParamMan( void )
  : m_is_ready( false ),
    m_file_name()
{
}

//_____________________________________________________________________________
TPCParamMan::~TPCParamMan( void )
{
  ClearACont();
  ClearTCont();
  ClearYCont();
  ClearCoboCont();
}

//_____________________________________________________________________________
void
TPCParamMan::ClearACont( void )
{
  del::ClearMap( m_APContainer );
}

//_____________________________________________________________________________
void
TPCParamMan::ClearTCont( void )
{
  del::ClearMap( m_TPContainer );
}

//_____________________________________________________________________________
void
TPCParamMan::ClearYCont( void )
{
  del::ClearMap( m_YPContainer );
}

//_____________________________________________________________________________
void
TPCParamMan::ClearCoboCont( void )
{
  del::ClearMap( m_CoboContainer );
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::Initialize( void )
{
  if( m_is_ready ){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    hddaq::cerr << FUNC_NAME << " file open fail : "
		<< m_file_name << std::endl;
    return false;
  }

  ClearACont();
  ClearTCont();
  ClearYCont();
  ClearCoboCont();

  Int_t line_number = 0;
  TString line;
  while( ifs.good() && line.ReadLine( ifs ) ){
    ++line_number;
    if( line.IsNull() || line[0]=='#' ) continue;
    std::istringstream input_line( line.Data() );
    Int_t layer, row, aty;
    Double_t p0, p1;
    if( input_line >> layer >> row >> aty >> p0 >> p1 ){
      Int_t key = MakeKey( layer, row );
      switch( aty ){
      case kAdc: {
	TPCAParam* pre_param = m_APContainer[key];
	TPCAParam* param = new TPCAParam( p0, p1 );
	m_APContainer[key] = param;
	if( pre_param ){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " layer = " << layer << ","
                      << " row = " << row
                      << " aty = " << aty
                      << std::endl;
	  delete pre_param;
	}
        break;
      }
      case kTdc: {
	TPCTParam *pre_param = m_TPContainer[key];
	TPCTParam *param = new TPCTParam( p0, p1 );
	m_TPContainer[key] = param;
	if( pre_param ){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " layer = " << layer << ","
                      << " row = " << row
                      << " aty = " << aty
                      << std::endl;
	  delete pre_param;
	}
        break;
      }
      case kY: {
        TPCYParam *pre_param = m_YPContainer[key];
        TPCYParam *param = new TPCYParam( p0, p1 );
        m_YPContainer[key] = param;
        if( pre_param ){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " layer = " << layer << ","
                      << " row = " << row
                      << " aty = " << aty
                      << std::endl;
          delete pre_param;
        }
        break;
      }
      case kCobo: {
	//layer is Cobo number (not real layer)
	Int_t Cobo = layer;
	
        TPCCoboParam *pre_param = m_CoboContainer[Cobo];
        TPCCoboParam *param = new TPCCoboParam( p0 );
        m_CoboContainer[Cobo] = param;
        if( pre_param ){
	  hddaq::cerr << FUNC_NAME << ": duplicated key "
		      << " following record is deleted." << std::endl
		      << " layer = " << layer << ","
                      << " row = " << row
                      << " aty = " << aty
                      << std::endl;
          delete pre_param;
        }
        break;
      }
      default:
        hddaq::cerr << FUNC_NAME << ": Invalid Input" << std::endl
                    << " ===> L" << line_number << " " << line << std::endl;
        break;
      } // switch
    } else {
      hddaq::cerr << FUNC_NAME << ": Invalid Input" << std::endl
		  << " ===> L" << line_number << " " << line << std::endl;
    }
  } // while

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::GetCDe( Int_t layer, Int_t row, Double_t adc, Double_t& de ) const
{
  TPCAParam* map = GetAmap( layer, row );
  if( map ){
    de = map->CDeltaE( adc );
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record for"
		<< " LayerId=" << std::setw(3) << std::dec << layer
		<< " RowId="   << std::setw(3) << std::dec << row
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::GetCTime( Int_t layer, Int_t row, Double_t tdc,
                       Double_t& time ) const
{
  TPCTParam* map = GetTmap( layer, row );
  if( map ){
    time = map->Time( tdc );
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record for"
		<< " LayerId=" << std::setw(3) << std::dec << layer
		<< " RowId="   << std::setw(3) << std::dec << row
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
Bool_t
TPCParamMan::GetC_Clock( Int_t layer, Int_t row, Double_t clock, Double_t& clock_cor ) const
{
  Int_t Cobo;
  if(layer<4||layer>5)
    Cobo = layer/4; 
  else if(layer==4){
    if(60<=row&&row<=99)
      Cobo = 1;
    else
      Cobo = 0;
  }
  else if(layer==5){
    if(72<=row&&row<=119)
      Cobo = 1;
    else
      Cobo = 0;
  }
  // const Double_t syn_t[8]={-36., 14., 5., 5.,
  // 			   4.,-16., 34., 13.};

  Double_t sync_t;
  TPCCoboParam* map = GetCobomap( Cobo );
  if( map ){
    sync_t = map->Clk_param( );
    
    clock_cor = clock;
    if(clock > sync_t)
      clock_cor -= 80.;

    //    std::cout<<"Cobo"<<Cobo<<", sync_t="<<sync_t<<std::endl;
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record for"
		<< " LayerId=" << std::setw(3) << std::dec << layer
		<< " RowId="   << std::setw(3) << std::dec << row
		<< std::endl;
    return false;
  }

}

//_____________________________________________________________________________
Bool_t
TPCParamMan::GetY( Int_t layer, Int_t row, Double_t time, Double_t& y ) const
{
  TPCYParam* map = GetYmap( layer, row );
  if( map ){
    y = map->Y( time );
    return true;
  } else {
    hddaq::cerr << FUNC_NAME << ": No record for"
		<< " LayerId=" << std::setw(3) << std::dec << layer
		<< " RowId="   << std::setw(3) << std::dec << row
		<< std::endl;
    return false;
  }
}

//_____________________________________________________________________________
TPCAParam*
TPCParamMan::GetAmap( Int_t layer, Int_t row ) const
{
  Int_t key = MakeKey( layer, row );
  AIterator itr = m_APContainer.find( key );
  if( itr != m_APContainer.end() )
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
TPCTParam*
TPCParamMan::GetTmap( Int_t layer, Int_t row ) const
{
  Int_t key = MakeKey( layer, row );
  TIterator itr = m_TPContainer.find( key );
  if( itr != m_TPContainer.end() )
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
TPCYParam*
TPCParamMan::GetYmap( Int_t layer, Int_t row ) const
{
  Int_t key = MakeKey( layer, row );
  YIterator itr = m_YPContainer.find( key );
  if( itr != m_YPContainer.end() )
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
TPCCoboParam*
TPCParamMan::GetCobomap( Int_t Cobo ) const
{
  CoboIterator itr = m_CoboContainer.find( Cobo );
  if( itr != m_CoboContainer.end() )
    return itr->second;
  else
    return nullptr;
}
