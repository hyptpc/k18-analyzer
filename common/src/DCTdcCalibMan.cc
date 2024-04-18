// -*- C++ -*-

#include "DCTdcCalibMan.hh"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <TRandom.h>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "FuncName.hh"

//______________________________________________________________________________
struct DCTdcCalMap
{
  double p0, p1;
  DCTdcCalMap( double q0, double q1 )
    : p0(q0), p1(q1)
  {}
};

//______________________________________________________________________________
DCTdcCalibMan::DCTdcCalibMan()
  : m_is_ready(false)
{
}

//______________________________________________________________________________
DCTdcCalibMan::~DCTdcCalibMan()
{
  ClearElements();
}

//______________________________________________________________________________
void
DCTdcCalibMan::ClearElements()
{
  del::ClearMap( m_container );
}

//______________________________________________________________________________
inline unsigned int
MakeKey( int det_id, int plane_id, int wire_id )
{
  return (det_id<<20) | (plane_id<<10) | int(wire_id);
}

//______________________________________________________________________________
bool
DCTdcCalibMan::Initialize()
{
  if( m_is_ready ){
    hddaq::cerr << "#W " << FUNC_NAME
		<< " already initialied" << std::endl;
    return false;
  }
  ClearElements();
  for(int i=0;i<m_file_names.size();i++){
    ReadFile(m_file_names.at(i) );
  }
  m_is_ready = true;
  return m_is_ready;
}
//______________________________________________________________________________
bool
DCTdcCalibMan::ReadFile( const TString& file_name )
{
  if(file_name=="") return true;
  std::ifstream ifs(file_name);
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << FUNC_NAME << " "
		<< "file open fail : " << file_name << std::endl;
    return false;
  }

  std::string line;
  while( ifs.good() && std::getline(ifs, line) ){
    if( line.empty() || line[0]=='#' ) continue;
    if(TString(line).Contains("cid",TString::ECaseCompare::kIgnoreCase) ) continue;
    std::istringstream iss( line );
    int    det_id=-1, plane_id=-1, wire_id=-1;
    double p0=-9999., p1=-9999.;
    if( iss >> det_id >> plane_id >> wire_id >> p0 >> p1 ){
      unsigned int key = MakeKey( det_id, plane_id, wire_id );
      DCTdcCalMap *tdc_calib = new DCTdcCalMap( p0, p1 );
      if( m_container[key] ) delete m_container[key];
      m_container[key] = tdc_calib;
    }
    else{
      hddaq::cerr << FUNC_NAME << ": Bad format => "
		  << line << std::endl;
    }
  }
  return true;
}

//______________________________________________________________________________
bool
DCTdcCalibMan::Initialize( const TString& file_name )
{
  m_file_names.push_back(file_name);
  return Initialize();
}

//______________________________________________________________________________
bool
DCTdcCalibMan::Initialize( const TString& file_name,const TString& file_name2 )
{
  m_file_names.push_back(file_name);
  m_file_names.push_back(file_name2);
  return Initialize();
}

//______________________________________________________________________________
DCTdcCalMap*
DCTdcCalibMan::GetMap( int det_id, int plane_id, int wire_id ) const
{
  unsigned int key = MakeKey( det_id, plane_id, wire_id );
  DCTdcCalMap *map = 0;
  DCTdcIterator itr = m_container.find(key);
  if( itr!=m_container.end() ) map = itr->second;
  return map;
}

//______________________________________________________________________________
bool
DCTdcCalibMan::GetTime( int det_id, int plane_id, int wire_id,
			int tdc, double& time ) const
{
  DCTdcCalMap *tdc_calib = GetMap( det_id, plane_id, wire_id );
  if( tdc_calib ){
    time = ( tdc + gRandom->Uniform(-0.5,0.5) ) * (tdc_calib->p1)  + (tdc_calib->p0);
    return true;
  }
  tdc_calib = GetMap( det_id, plane_id, 0 );
  if( tdc_calib ){
    time = ( tdc + gRandom->Uniform(-0.5,0.5) ) * (tdc_calib->p1)  + (tdc_calib->p0);
    return true;
  }
  else{
    hddaq::cerr << FUNC_NAME << ": No record. "
		<< " DetectorId=" << std::setw(3) << std::dec << det_id
		<< " PlaneId=" << std::setw(3) << std::dec << plane_id
		<< " WireId="  << std::setw(3) << std::dec << wire_id
		<< std::endl;
    return false;
  }
}

//______________________________________________________________________________
bool
DCTdcCalibMan::GetTdc( int det_id, int plane_id, int wire_id,
		       double time, int &tdc ) const
{
  DCTdcCalMap *tdc_calib = GetMap( det_id, plane_id, wire_id );
  if( tdc_calib ){
    tdc = int((time-(tdc_calib->p0))/(tdc_calib->p1));
    return true;
  }
  tdc_calib = GetMap( det_id, plane_id, 0 );
  if( tdc_calib ){
    tdc = int((time-(tdc_calib->p0))/(tdc_calib->p1));
    return true;
  }
  else{
    hddaq::cerr << FUNC_NAME << ": No record. "
		<< " DetectorId=" << std::setw(3) << std::dec << det_id
		<< " PlaneId=" << std::setw(3) << std::dec << plane_id
		<< " WireId="  << std::setw(3) << std::dec << wire_id
		<< std::endl;
    return false;
  }
}
