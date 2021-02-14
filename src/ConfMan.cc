// -*- C++ -*-

#include "ConfMan.hh"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <vector>

#include <lexical_cast.hh>
#include <filesystem_util.hh>
#include <replace_string.hh>
#include <std_ostream.hh>

// #include "BH1Filter.hh"
#include "BH1Match.hh"
#include "BH2Filter.hh"
#include "DCGeomMan.hh"
#include "DCTdcCalibMan.hh"
#include "DCDriftParamMan.hh"
#include "EventDisplay.hh"
#include "FieldMan.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "K18TransMatrix.hh"
#include "MatrixParamMan.hh"
#include "MsTParamMan.hh"
#include "UnpackerManager.hh"
#include "UserParamMan.hh"

namespace
{
using hddaq::unpacker::GUnpacker;
const TString kConfFile( "CONF" );
TString       sConfDir;
auto&       gUnpacker = GUnpacker::get_instance();
const auto& gMatrix   = MatrixParamMan::GetInstance();
const auto& gUser     = UserParamMan::GetInstance();
}

//_____________________________________________________________________________
ConfMan::ConfMan( void )
  : m_is_ready(false)
{
}

//_____________________________________________________________________________
ConfMan::~ConfMan( void )
{
}

//_____________________________________________________________________________
Bool_t
ConfMan::Initialize( void )
{
  if( m_is_ready ){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file[kConfFile] );
  if( !ifs.is_open() ){
    hddaq::cerr << FUNC_NAME << " cannot open file : "
                << m_file[kConfFile] << std::endl;
    return false;
  }

  hddaq::cout << FUNC_NAME << " " << m_file[kConfFile] << std::endl;
  sConfDir = hddaq::dirname( m_file[kConfFile].Data() );

  TString line;
  while( ifs.good() && line.ReadLine( ifs ) ){
    if( line.IsNull() || line[0]=='#' ) continue;

    line.ReplaceAll( ",",  "" ); // remove ,
    line.ReplaceAll( ":",  "" ); // remove :
    line.ReplaceAll( "\"",  "" ); // remove "

    std::istringstream iss( line.Data() );
    std::istream_iterator<std::string> begin( iss );
    std::istream_iterator<std::string> end;
    std::vector<TString> v( begin, end );
    if( v.size()<2 ) continue;

    TString key = v[0];
    TString val = v[1];
    hddaq::cout << " key = "   << std::setw(10) << std::left << key
		<< " value = " << std::setw(30) << std::left << val
		<< std::endl;

    m_file[key]   = FilePath( val );
    m_string[key] = val;
    m_double[key] = val.Atof();
    m_int[key]    = val.Atoi();
    m_bool[key]   = ( val.Atoi() == 1 );
  }

  if ( !InitializeParameterFiles() || !InitializeHistograms() )
    return false;

  if( gMatrix.IsReady() )
    gMatrix.Print2D();
  if( gUser.IsReady() )
    gUser.Print();

  m_is_ready = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::Initialize( const TString& file_name )
{
  m_file[kConfFile] = file_name;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
ConfMan::InitializeUnpacker( void )
{
  gUnpacker.set_config_file( m_file["UNPACK"].Data(),
			     m_file["DIGIT"].Data(),
			     m_file["CMAP"].Data() );
  return true;
}

//_____________________________________________________________________________
Bool_t
ConfMan::Finalize( void )
{
  return FinalizeProcess();
}

//_____________________________________________________________________________
TString
ConfMan::FilePath( const TString& src ) const
{
  std::ifstream tmp( src );
  if( tmp.good() )
    return src;
  else
    return sConfDir + "/" + src;
}
