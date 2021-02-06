// -*- C++ -*-

#include "UserParamMan.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

#include <std_ostream.hh>

#include "FuncName.hh"

namespace
{
const Double_t default_value = -9999.;
}

// if no parameter,
//   0: throw exception
//   1: return default value
#define ReturnDefaultValue 0

//______________________________________________________________________________
UserParamMan::UserParamMan( void )
  : m_is_ready(false), m_file_name("")
{
}

//______________________________________________________________________________
UserParamMan::~UserParamMan( void )
{
}

//______________________________________________________________________________
bool
UserParamMan::Initialize( void )
{
  std::ifstream ifs( m_file_name.c_str() );
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << FUNC_NAME << " "
		<< "No such parameter file : " << m_file_name << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line.empty() || line[0]=='#' ) continue;
    std::istringstream input_line(line);

    std::string first_param;
    input_line >> first_param;

    std::string& key = first_param;
    ParamArray   param_array;
    double       param;
    while( input_line >> param ){
      param_array.push_back( param );
    }

    m_param_map[key] = param_array;
  }

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
bool
UserParamMan::Initialize( const std::string& filename )
{
  m_file_name = filename;
  return Initialize();
};

//______________________________________________________________________________
int
UserParamMan::GetSize( const std::string& key ) const
{
  PIterator itr = m_param_map.find(key);
  if( itr==m_param_map.end() ){
    Print(m_file_name);
    hddaq::cerr << "#E " << FUNC_NAME << " "
		<< "No such key : " << key << std::endl;
    return 0;
  }

  return itr->second.size();
}

//______________________________________________________________________________
double
UserParamMan::GetParameter( const std::string& key, std::size_t i ) const
{
  std::stringstream param;
  param << key << "(" << i << ")";

  PIterator itr = m_param_map.find(key);

  if( itr==m_param_map.end() ){
    Print( m_file_name );
#if ReturnDefaultValue
    hddaq::cerr << "#E " << FUNC_NAME
		<< "set default value : "       << param.str() << " -> "
		<< default_value << std::endl;
    return default_value;
#else
    throw std::out_of_range( FUNC_NAME+" No such key : "+key );
#endif
  }

  if( i+1 > itr->second.size() ){
    Print( m_file_name );
#if ReturnDefaultValue
    hddaq::cerr << "#E " << FUNC_NAME
		<< "set default value : "        << param.str() << " -> "
		<< default_value << std::endl;
    return default_value;
#else
    throw std::out_of_range( FUNC_NAME+" No such key : "+key );
#endif
  }

  return itr->second.at(i);
}

//______________________________________________________________________________
void
UserParamMan::Print( const std::string& arg ) const
{
  hddaq::cout << "#D " << FUNC_NAME << " " << arg << std::endl;

  const int w = 20;
  PIterator itr, end=m_param_map.end();
  for( itr=m_param_map.begin(); itr!=end; ++itr){
    hddaq::cout << " key = " << std::setw(w) << std::left
		<< itr->first << itr->second.size() << " : ";
    const std::size_t size = itr->second.size();
    for( std::size_t i=0; i<size; ++i ){
      hddaq::cout << std::setw(5) << std::right
		  << itr->second.at(i) << " ";
    }
    hddaq::cout << std::endl;
  }
}
