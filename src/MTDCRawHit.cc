/**
 *  file: SDDRawHit.cc
 *  date: 2017.04.10
 *
 */

#include "MTDCRawHit.hh"

#include <iterator>

#include "DebugCounter.hh"

namespace
{
  const std::string class_name("MTDCRawHit");
}

//______________________________________________________________________________
MTDCRawHit::MTDCRawHit( int detector_id, int unit_id, int segment_id )
  : m_detector_id(detector_id),
    m_unit_id(unit_id),
    m_segment_id(segment_id),
    m_leading(1, -1), m_trailing(1, -1)
    
{
  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
MTDCRawHit::~MTDCRawHit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
MTDCRawHit::SetLeading( int tdc )
{
  if(-1 == m_leading.at(0)){
    m_leading.at(0) = tdc;
  }else{
    m_leading.push_back(tdc);
  }
}

//______________________________________________________________________________
void
MTDCRawHit::SetTrailing( int tdc )
{
  if(-1 == m_trailing.at(0)){
    m_trailing.at(0) = tdc;
  }else{
    m_trailing.push_back(tdc);
  }
}

//______________________________________________________________________________
int MTDCRawHit::SizeLeading( void ) const
{
  if(-1 == m_leading.at(0))
    return 0;
  else
    return m_leading.size();
}

//______________________________________________________________________________
int
MTDCRawHit::SizeTrailing( void ) const
{
  if(-1 == m_trailing.at(0))
    return 0;
  else
    return m_trailing.size();
}

//______________________________________________________________________________
void
MTDCRawHit::Clear( void )
{
  m_leading.clear();
  m_trailing.clear();
  m_leading.push_back(-1);
  m_trailing.push_back(-1);
}

//______________________________________________________________________________
void
MTDCRawHit::Print( const std::string& arg )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  hddaq::cerr << func_name << " " << arg << std::endl
	      << "detector_id = " << m_detector_id << std::endl
	      << "unit_id    = " << m_unit_id    << std::endl
	      << "segment_id  = " << m_segment_id   << std::endl;
  
  std::vector<int>::const_iterator itr, end;
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
