/**
 *  file: SsdParamMan.cc
 *  date: 2017.04.10
 *
 */

#include "SsdParamMan.hh"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include <std_ostream.hh>

#include "DeleteUtility.hh"

namespace
{
  const std::string& class_name("SsdParamMan");
  const int PlidMask  = 0x00ff;
  const int SegMask   = 0x0fff;
  const int PlidShift = 19;
  const int SegShift  =  0;
}

//______________________________________________________________________________
SsdParamMan::SsdParamMan( void )
  : m_is_ready(false),
    m_file_name("")
{
}

//______________________________________________________________________________
SsdParamMan::~SsdParamMan( void )
{
  ClearACont();
}

//______________________________________________________________________________
void
SsdParamMan::ClearACont( void )
{
  del::ClearMap( m_APContainer );
}

//______________________________________________________________________________
inline int
MakeKey( int pl, int seg )
{
  return ( ((pl&PlidMask)<<PlidShift) | ((seg&SegMask)<<SegShift) );
}

//______________________________________________________________________________
bool
SsdParamMan::Initialize( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  ClearACont();

  std::ifstream f( m_file_name.c_str() );
  if( !f.is_open() ){
    hddaq::cerr << func_name << ": file open fail" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string line;
  while( f.good() && std::getline( f, line ) ){
    if( line.empty() || line[0]=='#' ) continue;
    std::stringstream input_line( line );
    int plid; int seg;
    double p0[NumOfSampleSSD];
    double p1[NumOfSampleSSD];
    if( input_line >> plid >> seg
	>> p0[0] >> p1[0] >> p0[1] >> p1[1] >> p0[2] >> p1[2] >> p0[3] >> p1[3]
	>> p0[4] >> p1[4] >> p0[5] >> p1[5] >> p0[6] >> p1[6] >> p0[7] >> p1[7] ){
      int key = MakeKey( plid, seg );
      std::vector<double> pedestal;
      std::vector<double> rms;
      for(int i=0; i<NumOfSampleSSD; ++i){
	pedestal.push_back( p0[i] );
	rms.push_back( p1[i] );
      }
      SsdAParam *param = new SsdAParam( pedestal, rms );
      SsdAParam *pre_param = m_APContainer[key];
      m_APContainer[key] = param;
      if( pre_param ){
	hddaq::cerr << func_name << ": duplicated id number. "
		    << " following record is deleted." << std::endl;
	hddaq::cerr << " key = " << key << std::endl;
	delete pre_param;
      }
    }else{
      hddaq::cerr << "#E " << func_name << " invalid format "
		  << std::endl << line << std::endl;
    }
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
bool
SsdParamMan::Initialize( const std::string& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//______________________________________________________________________________
bool
SsdParamMan::GetDe( int plid, int seg, int isample, int adc, double &de ) const
{
  SsdAParam *map = GetAmap( plid, seg );
  if( !map ) return false;
  de = map->de( isample, adc );
  return true;
}

//______________________________________________________________________________
bool
SsdParamMan::GetRms( int plid, int seg, int isample, double &rms ) const
{
  SsdAParam *map = GetAmap( plid, seg );
  if( !map ) return false;
  rms = map->rms( isample );
  return true;
}

//______________________________________________________________________________
SsdAParam*
SsdParamMan::GetAmap( int plid, int seg ) const
{
  int key = MakeKey( plid, seg );
  SsdAParam *map = 0;
  AIterator i = m_APContainer.find( key );
  if( i!=m_APContainer.end() ) map = i->second;
  return map;
}
