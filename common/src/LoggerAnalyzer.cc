#include "LoggerAnalyzer.hh"

#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

#include <std_ostream.hh>

#include "DetectorID.hh"
#include "UnpackerManager.hh"
#include "TString.h"

namespace
{
  const std::string& class_name("LoggerAnalyzer");
  using namespace hddaq::unpacker;
  const UnpackerManager& gUnpacker =GUnpacker::get_instance();
  static const int ndata      = 3;
  static const int xstep      = 25750;
  static const int ystep      = 25900;
}

//______________________________________________________________________________
LoggerAnalyzer::LoggerAnalyzer( void )
  : m_serial(0)
{
}

//______________________________________________________________________________
LoggerAnalyzer::~LoggerAnalyzer( void )
{
}

//______________________________________________________________________________
void
LoggerAnalyzer::ClearAll( void )
{
  m_serial     = 0;
  m_raw_x      = 9999999;
  m_raw_y      = 9999999;
  m_segment    = -1;
  m_x          = -9999999.;
  m_y          = -9999999.;
  for(int i=0;i<10;i++){
    m_raw_data[i]=9999999;       
    m_data[i]=-9999999.;    
  }
}

//______________________________________________________________________________
bool
LoggerAnalyzer::Decode( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  ClearAll();

  static const int device_id = gUnpacker.get_device_id("Logger");
  static const int serial_id = gUnpacker.get_data_id("Logger", "serial");
  static const int xpos_id   = gUnpacker.get_data_id("Logger", "xpos");
  static const int ypos_id   = gUnpacker.get_data_id("Logger", "ypos");

  // serial
  if( gUnpacker.get_entries( device_id, 0, 0, 0, serial_id ) ){
    m_serial = gUnpacker.get( device_id, 0, 0, 0, serial_id );
  }
  // xpos
  if( gUnpacker.get_entries( device_id, 0, 0, 0, xpos_id ) ){
    m_raw_x = gUnpacker.get( device_id, 0, 0, 0, xpos_id );
    m_x     = m_raw_x/1000.;
  }
  // ypos
  if( gUnpacker.get_entries( device_id, 0, 0, 0, ypos_id ) ){
    int tmp = gUnpacker.get( device_id, 0, 0, 0, ypos_id );
    m_raw_y = - ( tmp & 0x800000 ) | ( tmp & 0x7ffffff);
    //    if(m_raw_y>2e8) m_raw_y=0;
    m_y     = m_raw_y/1000.;
  }

  for(int i=0;i<ndata;i++){
    const int tmpid   = gUnpacker.get_data_id("Logger", Form("data%d",i+1) );
    if( gUnpacker.get_entries( device_id, 0, 0, 0, tmpid ) ){
      m_raw_data[i] = gUnpacker.get( device_id, 0, 0, 0, tmpid );
      m_data[i]=m_raw_data[i]/100.;
    }
  }
  int tmpx=-1,tmpy=-1;
  for(int i=0;i<8;i++){
    if(fabs(i*xstep-m_raw_x)<xstep/2){
      tmpx=i;
      break;
    }
  }
  for(int i=0;i<5;i++){
    if(fabs(i*ystep-m_raw_y)<ystep/2){
      tmpy=i;
      break;
    }
  }
  if(tmpx>-1&&tmpy>-1){
    m_segment=tmpx+tmpy*8;
  }
  return true;
}
//______________________________________________________________________________
void
LoggerAnalyzer::Print( const std::string& arg ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  hddaq::cout << "#D " << func_name << " " << arg << std::endl
	      << "Serial       " << m_serial << std::endl
	      << std::setprecision(0) << std::fixed
	      << "(RawX,RawY)  ( "
	      << std::setw(8) << m_raw_x << ", "
	      << std::setw(8) << m_raw_y << " )" << std::endl
	      << std::setprecision(3)
	      << "(X,Y)        ( "
	      << std::setw(8) << m_x << ", "
	      << std::setw(8) << m_y << " )" << std::endl
	      << "Data1(T in)   " << m_data[0] << std::endl
	      << "Data2(T out)  " << m_data[1] << std::endl
	      << "Data3(Mag A)  " << m_data[2] << std::endl;
}
