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
Double_t
TPCRawHit::MaxAdc( Int_t min_t, Int_t max_t )
{
  int max_adc = 0;
  for(int i=0; i<m_fadc.size(); ++i){
    if(i>=min_t&&i<=max_t){
      if(max_adc<m_fadc[i])
	max_adc = m_fadc[i];
    } 
  }
  return max_adc;
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
