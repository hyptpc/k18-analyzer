// -*- C++ -*-

#include "DCDriftParamMan.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "Exception.hh"
#include "FuncName.hh"

namespace
{
inline void
DecodeKey( Int_t key, Int_t& plane, Int_t& wire )
{
  wire  = key & 0x3ff;
  plane = key >> 10;
}

inline Int_t
MakeKey( Int_t plane, Double_t wire )
{
  return (plane<<10) | Int_t(wire);
}
}

//_____________________________________________________________________________
DCDriftParamMan::DCDriftParamMan( void )
  : m_is_ready(false),
    m_file_name(),
    m_container()
{
}

//_____________________________________________________________________________
DCDriftParamMan::~DCDriftParamMan( void )
{
  ClearElements();
}

//_____________________________________________________________________________
void
DCDriftParamMan::ClearElements( void )
{
  del::ClearMap( m_container );
}

//_____________________________________________________________________________
Bool_t
DCDriftParamMan::Initialize( void )
{
  if( m_is_ready ){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    hddaq::cerr << FUNC_NAME
                << " file open fail : " << m_file_name << std::endl;
    return false;
  }

  ClearElements();

  TString line;
  while( ifs.good() && line.ReadLine( ifs ) ){
    if( line.IsNull() || line[0]=='#' ) continue;
    std::istringstream iss( line.Data() );
    // wid is not used at present
    // reserved for future updates
    Int_t pid, wid, type, np;
    Double_t param;
    std::vector<Double_t> q;
    if( iss >> pid >> wid >> type >> np ){
      while( iss >> param )
	q.push_back( param );
    }
    if( q.size() < 2 ){
      hddaq::cerr << FUNC_NAME << " format is wrong : " << line << std::endl;
      continue;
    }
    Int_t key = MakeKey( pid, 0 );
    DCDriftParamRecord *record = new DCDriftParamRecord( type, np, q );
    if( m_container[key] ){
      hddaq::cerr << FUNC_NAME << " "
		  << "duplicated key is deleted : " << key << std::endl;
      delete m_container[key];
    }
    m_container[key] = record;
  }

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
DCDriftParamMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
DCDriftParamRecord*
DCDriftParamMan::GetParameter( Int_t PlaneId, Double_t WireId ) const
{
  WireId = 0;
  Int_t key = MakeKey( PlaneId, WireId );
  DCDriftIterator itr = m_container.find( key );
  if( itr != m_container.end() )
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
Double_t
DCDriftParamMan::DriftLength1( Double_t dt, Double_t vel )
{
  return dt*vel;
}

//_____________________________________________________________________________
Double_t
DCDriftParamMan::DriftLength2( Double_t dt,
                               Double_t p1, Double_t p2, Double_t p3,
			       Double_t st, Double_t p5, Double_t p6 )
{
  Double_t dtmax=10.+p2+1./p6;
  Double_t dl;
  Double_t alph=-0.5*p5*st+0.5*st*p3*p5*(p1-st)
    +0.5*p3*p5*(p1-st)*(p1-st);
  if( dt<-10. || dt>dtmax+10. )
    dl = -500.;
  else if( dt<st )
    dl = 0.5*(p5+p3*p5*(p1-st))/st*dt*dt;
  else if( dt<p1 )
    dl = alph+p5*dt-0.5*p3*p5*(p1-dt)*(p1-dt);
  else if( dt<p2 )
    dl = alph+p5*dt;
  else
    dl = alph+p5*dt-0.5*p6*p5*(dt-p2)*(dt-p2);
  return dl;
}

//_____________________________________________________________________________
//DL = a0 + a1*Dt + a2*Dt^2
Double_t
DCDriftParamMan::DriftLength3( Double_t dt, Double_t p1, Double_t p2,
                               Int_t PlaneId )
{
  if ( PlaneId>=1 && PlaneId <=4 ) {
    if( dt > 130. )
      return 999.9;
  } else if (PlaneId == 5) {
    if (dt > 500.)
      ;
  } else if (PlaneId>=6 && PlaneId <=11) {
    if (dt > 80.)
      return 999.9;
  } else if (PlaneId>=101 && PlaneId <=124) {
    if (dt > 80.)
      return 999.9;
  }

  return dt*p1+dt*dt*p2;
}

//_____________________________________________________________________________
Double_t
DCDriftParamMan::DriftLength4( Double_t dt,
                               Double_t p1, Double_t p2, Double_t p3 )
{
  if (dt < p2)
    return dt*p1;
  else
    return p3*(dt-p2) + p1*p2;
}

//_____________________________________________________________________________
Double_t
DCDriftParamMan::DriftLength5( Double_t dt, Double_t p1, Double_t p2,
			       Double_t p3, Double_t p4, Double_t p5 )
{
  if (dt < p2)
    return dt*p1;
  else if (dt >= p2 && dt < p4)
    return p3*(dt-p2) + p1*p2;
  else
    return p5*(dt-p4) + p3*(p4-p2) + p1*p2;
}

//_____________________________________________________________________________
////////////////Now using
// DL = a0 + a1*Dt + a2*Dt^2 + a3*Dt^3 + a4*Dt^4 + a5*Dt^5
Double_t
DCDriftParamMan::DriftLength6( Int_t PlaneId, Double_t dt,
			       Double_t p1, Double_t p2, Double_t p3,
                               Double_t p4, Double_t p5 )
{
  Double_t dl = dt*p1+dt*dt*p2+p3*pow(dt, 3.0)+p4*pow(dt, 4.0)+p5*pow(dt, 5.0);

  switch( PlaneId ){
    // BC3&4
  case 113: case 114: case 115: case 116: case 117: case 118:
  case 119: case 120: case 121: case 122: case 123: case 124:
    //if( dt<-10. || dt>50. ) // Tight drift time selection
    if (dt<-10 || dt>80) // Loose drift time selection
      return 999.9;
    if( PlaneId==123 || PlaneId==124 ){
      if( dt>35 ) dt=35.;
      dl = dt*p1+dt*dt*p2+p3*pow(dt, 3.0)+p4*pow(dt, 4.0)+p5*pow(dt, 5.0);
    }else if( dt>32. ){
      dt = 32.;
    }
    if( dl>1.5 )
      return 1.5;
    if( dl<0. )
      return 0.;
    else
      return dl;
    break;
    // SDC1
  case 1: case 2: case 3: case 4: case 5: case 6:
    if( dt<-10 || 150<dt ) // Loose drift time selection
      return 999.9;
    // if( dt>120. ){
    //   dt = 120.;
    //   dl = dt*p1+dt*dt*p2+p3*pow(dt, 3.0)+p4*pow(dt, 4.0)+p5*pow(dt, 5.0);
    // }
    if( dl>3.0 || dt>120. )
      return 3.0;
    if( dl<0. )
      return 0.;
    else
      return dl;
    break;
  case 7: case 8: case 9:
    return 0.;
    // SDC2
  case 31: case 32: case 33: case 34:
    if( dt<-20. || dt>150. )
      return 999.9;
    if( dl>4.5 || dt>120. )
      return 4.5;
    if( dl<0. )
      return 0.;
    else
      return dl;
    break;
    //For SDC3
  case 35: case 36: case 37: case 38:
    if( dt<-20. || dt>300. )
      return 999.9;
    if( dl>10. || dt>250. )
      return 10.0;
    if( dl<0. )
      return 0.;
    else
      return dl;
    break;
  default:
    throw Exception( FUNC_NAME+" invalid plane id : "
                     +TString::Itoa(PlaneId, 10) );
  }
}

//_____________________________________________________________________________
Bool_t
DCDriftParamMan::CalcDrift( Int_t PlaneId, Double_t WireId, Double_t ctime,
			    Double_t & dt, Double_t & dl ) const
{
  DCDriftParamRecord *record = GetParameter( PlaneId, WireId );
  if( !record ){
    hddaq::cerr << FUNC_NAME << " No record. "
		<< " PlaneId=" << std::setw(3) << PlaneId
		<< " WireId="  << std::setw(3) << WireId << std::endl;
    return false;
  }

  Int_t type = record->type;
  // Int_t np   = record->np;
  std::vector<Double_t> p = record->param;

  dt = p[0]-ctime;

  switch( type ){
  case 1:
    dl=DriftLength1( dt, p[1] );
    return true;
  case 2:
    dl=DriftLength2( dt, p[1], p[2], p[3], p[4],
		     p[5], p[6] );
    return true;
  case 3:
    dl=DriftLength3( dt, p[1], p[2], PlaneId );
    return true;
  case 4:
    dl=DriftLength4( dt, p[1], p[2], p[3] );
    return true;
  case 5:
    dl=DriftLength5( dt, p[1], p[2], p[3], p[4], p[5] );
    return true;
  case 6:
    dl=DriftLength6( PlaneId, dt, p[1], p[2], p[3], p[4], p[5] );
    return true;
  default:
    hddaq::cerr << FUNC_NAME << " invalid type : " << type << std::endl;
    return false;
  }
}
