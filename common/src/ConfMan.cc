/**
 *  file: ConfMan.cc
 *  date: 2017.04.10
 *
 */

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

#include <TMacro.h>
#include <TFile.h>
#include <TSystem.h>

#include "UnpackerManager.hh"
#include "UserParamMan.hh"
#include "ReslMapMan.hh"

namespace
{
  using namespace hddaq::unpacker;
  const std::string& class_name("ConfMan");
  std::string sConfDir;
  UnpackerManager&      gUnpacker = GUnpacker::get_instance();
  const UserParamMan&   gUser     = UserParamMan::GetInstance();
  enum EArg
    {
      kArgProcess,
      kArgConfFile,
      kArgInFile,
      kArgOutFile,
      kArgRunNum,
      kArgc
    };
}

//______________________________________________________________________________
ConfMan::ConfMan( void )
  : m_is_ready(false), m_runnum(-1), m_maxloop(-1), m_skip(-1), m_resolution_cdc(-1)
{
}

//______________________________________________________________________________
ConfMan::~ConfMan( void )
{
}

//______________________________________________________________________________
bool
ConfMan::Initialize( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_is_ready ){
    hddaq::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file[kConfFile].c_str() );
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << func_name
		<< " cannot open file : " << m_file[kConfFile] << std::endl;
    return false;
  }
  m_keys_used.push_back(kConfFile);
  hddaq::cout << "#D " << func_name << std::endl;

  sConfDir = hddaq::dirname(m_file[kConfFile]);

  std::string line;
  while( ifs.good() && std::getline( ifs, line ) ){
    if( line.empty() || line[0]=='#' ) continue;

    hddaq::replace_all(line, ",",  ""); // remove ,
    hddaq::replace_all(line, ":",  ""); // remove :
    hddaq::replace_all(line, "\"", ""); // remove "

    std::istringstream iss( line );
    std::istream_iterator<std::string> begin( iss );
    std::istream_iterator<std::string> end;
    std::vector<std::string> v( begin, end );
    if( v.size()<2 ) continue;

    std::string key = v[0];
    hddaq::cout << " key = "   << std::setw(10) << std::left << key
		<< " value = " << std::setw(30) << std::left << v[1]
		<< std::endl;

    std::string tmpstr=v[1];
    if(v[1].find("_list")!=std::string::npos){
      tmpstr="none";
      std::ifstream ifs2(FilePath(v[1]));
      int run1, run2;
      std::string line2;
      while( ifs2.good() && std::getline( ifs2, line2 ) ){		
	if( line2[0]=='#' || line2.empty() ) continue;
	std::istringstream iss2( line2 );
	if(iss2>>run1>>run2>>tmpstr){
	  if(m_runnum>=run1&&m_runnum<=run2){
	    hddaq::cout << "       " << std::setw(10);
	    hddaq::cout <<v[1]<<" : "<<tmpstr<<std::endl;
	    break;	  
	  }
	  tmpstr="none";
	}
      }
      ifs2.close();
    }
    if(tmpstr.find("none")!=std::string::npos) continue;
    //    std::cout<<key<<"  "<<tmpstr<<std::endl;
    m_file[key]   = FilePath(tmpstr);
    m_string[key] = tmpstr;
    m_double[key] = hddaq::lexical_cast<double>(tmpstr);
    m_int[key]    = hddaq::lexical_cast<int>(tmpstr);
    m_bool[key]   = hddaq::lexical_cast<bool>(tmpstr);
  }

  if ( !InitializeParameterFiles() || !InitializeHistograms() )
    return false;

  if( gUser.IsReady() )
    gUser.Print();

  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
bool
ConfMan::Initialize( const std::string& file_name, const int &runnum )
{
  m_file[kConfFile] = file_name;
  m_runnum          = runnum;
  return Initialize();
}

//______________________________________________________________________________
bool
ConfMan::Initialize( int argc, char **argv )
{
  return Initialize();
}

//______________________________________________________________________________
bool
ConfMan::ParseCommand( int argc, char **argv )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  std::vector<std::string> arg( argv, argv + argc );
  const std::string& process = arg[kArgProcess];
  if( argc<kArgc ){
    hddaq::cout << "#D Usage: " << hddaq::basename(process)
  		<< " --conf=[analyzer config file]"
  		<< " --infile=[data input stream]"
  		<< " --outfile=[output root file]"
  		<< " --run=[run number]"
  		<< "(--=[skip])"
  		<< "(--skip=[skip])"
  		<< "(--maxloop=[max loop])"
		<< "(--cdctree=[tracking file name])"
		<< "(--cdcresol=[cdc resolution])"
  		<< std::endl;
    return EXIT_SUCCESS;
  }

  std::istringstream iss;
  std::string tmp;
  m_file[kCDCFile]="None";
  m_is_mc=false;
  m_is_cosmic=false;
  for (int i = 1 ; i < argc ; i++) {
    std::string arg = argv[i];
    std::cout<<arg<<std::endl;
    iss.str("");
    iss.clear();
    if (arg.substr(0, 7) == "--conf=") {
      iss.str(arg.substr(7));
      iss >> tmp;
      m_file[kConfFile] = tmp;
    }
    if (arg.substr(0, 9) == "--infile=") {
      iss.str(arg.substr(9));
      iss >> tmp;
      m_file[kInFile] = tmp;
    }
    if (arg.substr(0, 10) == "--outfile=") {
      iss.str(arg.substr(10));
      iss >> tmp;
      m_file[kOutFile] = tmp;
    }
    if (arg.substr(0, 6) == "--run=") {
      iss.str(arg.substr(6));
      iss >> m_runnum;
    }
    if (arg.substr(0, 10) == "--maxloop=") {
      iss.str(arg.substr(10));
      int tmpint;
      iss >> tmpint;
      if(m_maxloop<0 || tmpint<m_maxloop) m_maxloop=tmpint;
    }
    if (arg.substr(0, 7) == "--skip=") {
      iss.str(arg.substr(7));
      iss >> m_skip;
    }
    if (arg.substr(0, 10) == "--cdctree=") {
      iss.str(arg.substr(10));
      iss >> tmp;
      m_file[kCDCFile] = tmp;
    }
    if (arg.substr(0, 11) == "--cdcresol=") {
      iss.str(arg.substr(11));
      iss >> m_resolution_cdc;
    }
    if (arg.substr(0, 8) == "--cosmic") {
      m_is_cosmic=true;
    }
    if (arg.substr(0, 4) == "--mc") {
      m_is_mc=true;
    }
  }
  hddaq::cout << "#D " << func_name
	      << " runnumber           : " << m_runnum<<std::endl;
  hddaq::cout << "#D " << func_name
	      << " configuration file  : " << m_file[kConfFile] << std::endl;
  hddaq::cout << "#D " << func_name
	      << " input raw data file : " << m_file[kInFile] << std::endl;
  hddaq::cout << "#D " << func_name
	      << " output root file    : " << m_file[kOutFile] << std::endl;
  hddaq::cout << "#D " << func_name
	      << " cdc track root file : " << m_file[kCDCFile] << std::endl;
  hddaq::cout << "#D " << func_name
	      << " max loop            : " << m_maxloop << std::endl;
  hddaq::cout << "#D " << func_name
	      << " skip                : " << m_skip << std::endl;
  hddaq::cout << "#D " << func_name
	      << " cosmic flag         : " << m_is_cosmic << std::endl;
  hddaq::cout << "#D " << func_name
	      << " MC flag             : " << m_is_mc << std::endl;
  return true;
}
//______________________/________________________________________________________
bool
ConfMan::InitializeUnpacker( void )
{
  gUnpacker.set_config_file( m_file["UNPACK"],
			     m_file["DIGIT"],
			     m_file["CMAP"] );
  return true;
}

//______________________________________________________________________________
bool
ConfMan::Finalize( void )
{
  WriteParams();
  return FinalizeProcess();
}
//______________________________________________________________________________
void
ConfMan::WriteParams( void )
{
  gFile->mkdir("param");
  gFile->cd("param");
  for(auto itr = m_keys_used.begin(); itr != m_keys_used.end(); ++itr) {
    TString tmp=m_file[*itr];
    tmp=tmp(tmp.Last('/')+1,tmp.Sizeof());
    tmp=tmp.ReplaceAll('.','_');
    TMacro paramfile;
    paramfile.SetName(tmp.Data());
    paramfile.SetTitle(m_file[*itr].c_str());
    paramfile.ReadFile(m_file[*itr].c_str());
    paramfile.Write();
    std::cout<<std::left<<std::setw(10)<<*itr
	     <<std::left<<std::setw(30)<<tmp
	     <<m_file[*itr].c_str()<<std::endl;
  }
  gFile->cd();
}
//______________________________________________________________________________
bool
ConfMan::BeginRun()
{
  if(IsMC()){
    ReslMapMan&  gResol = ReslMapMan::GetInstance();
    double dlres=get_resolution_cdc();
    if(gResol.IsReady()&&dlres>0){
      gResol.SetParam(100,0,0,dlres,0);//cdc
    }
    gResol.PrintMap();
  }
  return BeginRunProcess();
}

//______________________________________________________________________________
std::string
ConfMan::FilePath( const std::string& src ) const
{
  return sConfDir + "/" + src;
  std::ifstream tmp( src.c_str() );
  if ( tmp.good() )
    return src;
  else
    return sConfDir + "/" + src;
}
