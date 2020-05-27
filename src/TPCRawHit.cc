/**
 *  file: TPCRawHit.cc
 *  date: 2020.04.11
 *
 */

#include "TPCRawHit.hh"

#include <iostream>
#include <iterator>

#include "std_ostream.hh"

#include "DebugCounter.hh"
#include "TPCPadHelper.hh"

namespace
{
  const std::string& class_name("TPCRawHit");
}

//______________________________________________________________________________
TPCRawHit::TPCRawHit( int padid, double y, double charge )
  : m_pad_id(padid),
    m_charge(charge)
{
  m_pos = tpc::getPosition(padid);
  m_pos.SetY(y);
  m_layer_id = tpc::getLayerID(padid);
  m_row_id   = tpc::getRowID(padid);

  debug::ObjectCounter::increase(class_name);
}

//______________________________________________________________________________
TPCRawHit::~TPCRawHit( void )
{
  debug::ObjectCounter::decrease(class_name);
}

//______________________________________________________________________________
void
TPCRawHit::Print( const std::string& arg ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  std::cerr << func_name << " " << arg << std::endl
	      << "padID  = "  << m_pad_id  << std::endl
	      << "(x,y,z)= (" << m_pos.X() <<
	      		   ","<< m_pos.Y() <<
			   ","<< m_pos.Z()<<")" << std::endl
	      << "charge = "  << m_charge  << std::endl;
  std::cerr << std::endl;
}
