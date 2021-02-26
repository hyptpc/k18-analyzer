// -*- C++ -*-

#include "TPCRawHit.hh"

#include <iostream>
#include <iterator>

#include <std_ostream.hh>

#include "DebugCounter.hh"
#include "FuncName.hh"
#include "TPCPadHelper.hh"

//_____________________________________________________________________________
TPCRawHit::TPCRawHit( Int_t layer, Int_t row )
  : m_layer_id( layer ),
    m_row_id( row )
{
  debug::ObjectCounter::increase( ClassName() );
}

//_____________________________________________________________________________
TPCRawHit::~TPCRawHit( void )
{
  debug::ObjectCounter::decrease( ClassName() );
}

//_____________________________________________________________________________
void
TPCRawHit::AddFadc( Int_t adc )
{
  m_fadc.push_back( adc );
}

//_____________________________________________________________________________
void
TPCRawHit::Print( Option_t* ) const
{
  hddaq::cout << FUNC_NAME << " " << std::endl
              << "   layer = " << m_layer_id  << std::endl
              << "   row   = " << m_row_id << std::endl
              << std::endl;
}
