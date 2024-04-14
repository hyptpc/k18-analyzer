/**
 *  file: SDDRawHit.cc
 *  date: 2017.04.10
 *
 */

#include "SDDRawHit.hh"

#include <iterator>

#include "DebugCounter.hh"

namespace
{
  const std::string class_name("SDDRawHit");
}

//______________________________________________________________________________
SDDRawHit::SDDRawHit( int detector_id, int port_id, int unit_id )
  : m_detector_id(detector_id),
    m_port_id(port_id),
    m_unit_id(unit_id),
    m_adc(8, 0),
    m_leading(1, -1), m_trailing(1, -1),
    m_resetl(1, -1), m_resett(1, -1)    
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
SDDRawHit::~SDDRawHit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
SDDRawHit::SetAdc( int i, int adc )
{
  m_adc.at(i) = adc;
}

//______________________________________________________________________________
void
SDDRawHit::SetLeading( int tdc )
{
  if(-1 == m_leading.at(0)){
    m_leading.at(0) = tdc;
  }else{
    m_leading.push_back(tdc);
  }
}

//______________________________________________________________________________
void
SDDRawHit::SetResetLeading( int tdc )
{
  if(-1 == m_resetl.at(0)){
    m_resetl.at(0) = tdc;
  }else{
    m_resetl.push_back(tdc);
  }
}

//______________________________________________________________________________
void
SDDRawHit::SetTrailing( int tdc )
{
  if(-1 == m_trailing.at(0)){
    m_trailing.at(0) = tdc;
  }else{
    m_trailing.push_back(tdc);
  }
}
//______________________________________________________________________________
void
SDDRawHit::SetResetTrailing( int tdc )
{
  if(-1 == m_resett.at(0)){
    m_resett.at(0) = tdc;
  }else{
    m_resett.push_back(tdc);
  }
}

//______________________________________________________________________________
int
SDDRawHit::SizeAdc( void ) const
{
  if(-1 == m_adc.at(0))
    return 0;
  else
    return m_adc.size();
}

//______________________________________________________________________________
int SDDRawHit::SizeLeading( void ) const
{
  if(-1 == m_leading.at(0))
    return 0;
  else
    return m_leading.size();
}

//______________________________________________________________________________
int
SDDRawHit::SizeTrailing( void ) const
{
  if(-1 == m_trailing.at(0))
    return 0;
  else
    return m_trailing.size();
}

//______________________________________________________________________________
int SDDRawHit::SizeResetLeading( void ) const
{
  if(-1 == m_resetl.at(0))
    return 0;
  else
    return m_resetl.size();
}

//______________________________________________________________________________
int
SDDRawHit::SizeResetTrailing( void ) const
{
  if(-1 == m_resett.at(0))
    return 0;
  else
    return m_resett.size();
}

//______________________________________________________________________________
void
SDDRawHit::Clear( void )
{
  m_adc.clear();
  m_leading.clear();
  m_trailing.clear();
  m_resetl.clear();
  m_resett.clear();
  m_adc.push_back(-1);
  m_leading.push_back(-1);
  m_trailing.push_back(-1);
  m_resetl.push_back(-1);
  m_resett.push_back(-1);
}

//______________________________________________________________________________
void
SDDRawHit::Print( const std::string& arg )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  hddaq::cerr << func_name << " " << arg << std::endl
	      << "detector_id = " << m_detector_id << std::endl
	      << "port_id    = " << m_port_id    << std::endl
	      << "unit_id    = " << m_unit_id    << std::endl;

  
  std::vector<int>::const_iterator itr, end;
  hddaq::cout << "adc         = " << m_adc.size() << " ";
  std::copy( m_adc.begin(), m_adc.end(),
	     std::ostream_iterator<int>(hddaq::cout," ") );
  hddaq::cout << std::endl
	      << "leading        = " << m_leading.size() << " ";
  std::copy( m_leading.begin(), m_leading.end(),
	     std::ostream_iterator<int>(hddaq::cout," ") );
  hddaq::cout << std::endl
	      << "trailing        = " << m_trailing.size() << " ";
  std::copy( m_trailing.begin(), m_trailing.end(),
	     std::ostream_iterator<int>(hddaq::cout," ") );
  hddaq::cout << std::endl;
}
