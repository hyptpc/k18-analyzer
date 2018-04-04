/**
 *  file: EMCAnalyzer.cc
 *  date: 2017.04.10
 *
 */

#include "EMCAnalyzer.hh"

#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>

#include <std_ostream.hh>

#include "DetectorID.hh"
#include "UnpackerManager.hh"

namespace
{
  const std::string& class_name("EMCAnalyzer");
  using namespace hddaq::unpacker;
  const UnpackerManager& gUnpacker =GUnpacker::get_instance();
}

//______________________________________________________________________________
EMCAnalyzer::EMCAnalyzer( void )
  : m_serial(0),
    m_raw_x(-9999999.),
    m_raw_y(-9999999.),
    m_x(-9999999.),
    m_y(-9999999.),
    m_state(-1),
    m_utime(0),
    m_ltime(0),
    m_time(0)
{
}

//______________________________________________________________________________
EMCAnalyzer::~EMCAnalyzer( void )
{
}

//______________________________________________________________________________
void
EMCAnalyzer::ClearAll( void )
{
  m_serial     = 0;
  m_raw_x      = -9999999.;
  m_raw_y      = -9999999.;
  m_x          = -9999999.;
  m_y          = -9999999.;
  m_state      = -1;
  m_utime      = 0;
  m_ltime      = 0;
  m_time       = 0;
}

//______________________________________________________________________________
bool
EMCAnalyzer::Decode( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  ClearAll();

  static const int device_id = gUnpacker.get_device_id("EMC");
  static const int serial_id = gUnpacker.get_data_id("EMC", "serial");
  static const int xpos_id   = gUnpacker.get_data_id("EMC", "xpos");
  static const int ypos_id   = gUnpacker.get_data_id("EMC", "ypos");
  static const int state_id  = gUnpacker.get_data_id("EMC", "state");
  static const int utime_id  = gUnpacker.get_data_id("EMC", "utime");
  static const int ltime_id  = gUnpacker.get_data_id("EMC", "ltime");

  // serial
  if( gUnpacker.get_entries( device_id, 0, 0, 0, serial_id ) ){
    m_serial = gUnpacker.get( device_id, 0, 0, 0, serial_id );
  }
  // xpos
  if( gUnpacker.get_entries( device_id, 0, 0, 0, xpos_id ) ){
    m_raw_x = gUnpacker.get( device_id, 0, 0, 0, xpos_id );
    m_x     = 500. - ( m_raw_x/1000. );
  }
  // ypos
  if( gUnpacker.get_entries( device_id, 0, 0, 0, ypos_id ) ){
    m_raw_y = gUnpacker.get( device_id, 0, 0, 0, ypos_id );
    m_y     = 500. - ( m_raw_y/1000. );
  }
  // state
  if( gUnpacker.get_entries( device_id, 0, 0, 0, state_id ) ){
    m_state = gUnpacker.get( device_id, 0, 0, 0, state_id );
  }
  // utime
  if( gUnpacker.get_entries( device_id, 0, 0, 0, utime_id ) ){
    m_utime = gUnpacker.get( device_id, 0, 0, 0, utime_id );
  }
  // ltime
  if( gUnpacker.get_entries( device_id, 0, 0, 0, ltime_id ) ){
    m_ltime = gUnpacker.get( device_id, 0, 0, 0, ltime_id );
  }

  unsigned long long u = (unsigned long long)(m_utime)<<32;
  unsigned long long l = (unsigned long long)(m_ltime)<<0;
  m_time = ( u + l )/1000000;

  return true;
}

//______________________________________________________________________________
void
EMCAnalyzer::Print( const std::string& arg ) const
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
	      << "Time         " << m_time << std::endl;
}
