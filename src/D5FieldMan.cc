// -*- C++ -*-

#include "D5FieldMan.hh"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

#include "FuncName.hh"

#include <std_ostream.hh>

//_____________________________________________________________________________
D5FieldMan::D5FieldMan( void )
  : m_table(),
    m_is_ready( false ),
    m_file_name()
{
}

//_____________________________________________________________________________
D5FieldMan::~D5FieldMan( void )
{
}

//_____________________________________________________________________________
Bool_t
D5FieldMan::Initialize( void )
{
  if( m_is_ready ){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    hddaq::cerr << FUNC_NAME << " file open fail : "
                << m_file_name << std::endl;
    return false;
  }

  // Each line: "run D5[T] rms[T]"; '#' starts a comment.
  std::string line;
  while( std::getline( ifs, line ) ){
    const auto comment = line.find( '#' );
    if( comment != std::string::npos ) line.erase( comment );
    std::istringstream iss( line );
    UInt_t run;
    Entry entry;
    if( !( iss >> run ) ) continue;
    if( !( iss >> entry.field >> entry.rms ) ){
      hddaq::cerr << FUNC_NAME << " invalid line : " << line << std::endl;
      return false;
    }
    m_table[run] = entry;
  }

  hddaq::cout << FUNC_NAME << " " << m_table.size() << " runs from "
              << m_file_name << std::endl;
  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
D5FieldMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
D5FieldMan::GetField( UInt_t run, Double_t& field, Double_t& rms ) const
{
  const auto itr = m_table.find( run );
  if( itr == m_table.end() ) return false;
  field = itr->second.field;
  rms   = itr->second.rms;
  return true;
}

//_____________________________________________________________________________
Double_t
D5FieldMan::CentralMomentum( Double_t field )
{
  return 0.299792458 * std::abs( field ) * kRadius * kFieldFactor;
}
