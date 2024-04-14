// -*- C++ -*-

#include "HodoPHCMan.hh"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include <std_ostream.hh>

#include "TString.h"
#include "DeleteUtility.hh"
#include "MathTools.hh"

namespace
{
  const std::string& class_name("HodoPHCMan");
  const int SegMask  = 0x03FF;
  const int CidMask  = 0x00FF;
  const int PlidMask = 0x00FF;
  const int UdMask   = 0x0003;
  const int SegShift  =  0;
  const int CidShift  = 11;
  const int PlidShift = 19;
  const int UdShift   = 27;
}

//______________________________________________________________________________
HodoPHCParam::HodoPHCParam( int type, int n_param, std::vector<double> parlist )
  : m_type(type),
    m_n_param(n_param),
    m_param_list(parlist)
{
}

//______________________________________________________________________________
HodoPHCParam::~HodoPHCParam( void )
{
}

//______________________________________________________________________________
inline int
MakeKey( int cid, int pl, int seg, int ud )
{
  return ( ( (cid&CidMask) << CidShift  ) |
	   ( (pl&PlidMask) << PlidShift ) |
	   ( (seg&SegMask) << SegShift  ) |
	   ( (ud&UdMask)   << UdShift   ) );
}

//______________________________________________________________________________
double
HodoPHCParam::DoPHC( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  double ctime = time;

  switch(m_type){
  case 0:
    ctime = time; break;
  case 1:
    ctime = type1Correction(time,de); break;
  case 2:
    ctime = type2Correction(time,de); break;
  case 3:
    ctime = type3Correction(time,de); break; // SDD
  case 4:
    ctime = type4Correction(time,de); break; // CVC/NC/PC
  case 5:
    ctime = type5Correction(time,de); break; // CDC in E15 (Run65)
  case 6:
    ctime = type6Correction(time,de); break; // PbF2 non-linearity
  default:
    hddaq::cerr << func_name << ": No Correction Method. type="
		<< m_type << std::endl;
  }
  return ctime;
}

//______________________________________________________________________________
double
HodoPHCParam::DoRPHC( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  double ctime = time;

  switch(m_type){
  case 0:
    ctime = time; break;
  case 1:
    ctime = type1RCorrection(time,de); break;
  default:
    hddaq::cerr << func_name << ": No Correction Method. type="
		<< m_type << std::endl;
  }
  return ctime;
}

//______________________________________________________________________________
double
HodoPHCParam::type1Correction( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_param_list.size()<2 )
    throw std::out_of_range(func_name+" invalid parameter");

  if( std::abs(de)<math::Epsilon() )
    de = math::Epsilon();

  return time-(m_param_list[0]/sqrt(std::abs(de)) + m_param_list[1]);
}
//______________________________________________________________________________
double
HodoPHCParam::type2Correction( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_param_list.size()<3 )
    throw std::out_of_range(func_name+" invalid parameter");

  // Correction function for fiber is quadratic function
  return time-(m_param_list[0]/sqrt(std::abs(de)) + m_param_list[1] + m_param_list[2]*de);
}

double
HodoPHCParam::type3Correction( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_param_list.size()<6 )
    throw std::out_of_range(func_name+" invalid parameter");

  return time-(m_param_list[0]/sqrt(std::abs(de)) + m_param_list[1] + m_param_list[2]*de + m_param_list[3]*de*de+ m_param_list[4]*de*de*de+ m_param_list[5]*de*de*de*de);
  // Correction function for fiber is quadratic function
  //return time-(m_param_list[0]*w*w*w*w + m_param_list[1]*w*w*w + m_param_list[2]*w*w + m_param_list[3]*w + m_param_list[4]);
}

double
HodoPHCParam::type4Correction( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_param_list.size()<3 )
    throw std::out_of_range(func_name+" invalid parameter");

  // Correction function for fiber is quadratic function
  return time-(m_param_list[0]/sqrt(std::abs(de)) + m_param_list[1] + m_param_list[2]/de);
}

double
HodoPHCParam::type5Correction( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_param_list.size()<4 )
    throw std::out_of_range(func_name+" invalid parameter");

  // Correction function for fiber is quadratic function
  double corr=m_param_list[0]/sqrt(std::abs(de)) + m_param_list[1] + m_param_list[2]*de + m_param_list[3]/de;
  return time-corr;
}

double
HodoPHCParam::type6Correction( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_param_list.size()<3 )
    throw std::out_of_range(func_name+" invalid parameter");

  // Correction function for fiber is quadratic function
  if(de<1) return de;
  double corr=de+m_param_list[0]*pow(de-1,2);
  return corr;
}

//______________________________________________________________________________
double
HodoPHCParam::type1RCorrection( double time, double de ) const
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");

  if( m_param_list.size()<3 )
    throw std::out_of_range(func_name+" invalid parameter");
  if( std::abs(de-m_param_list[1])<math::Epsilon() )
    de = m_param_list[1] + math::Epsilon();

  return time+m_param_list[0]/sqrt( std::abs(de-m_param_list[1]) )-m_param_list[2];
}

//______________________________________________________________________________
HodoPHCMan::HodoPHCMan( void )
  : m_is_ready(false)
{
}

//______________________________________________________________________________
HodoPHCMan::~HodoPHCMan()
{
  ClearElements();
}

//______________________________________________________________________________
void HodoPHCMan::ClearElements( void )
{
  del::ClearMap( m_container );
}

//______________________________________________________________________________
bool
HodoPHCMan::Initialize( void )
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_is_ready ){
    hddaq::cerr << "#W " << func_name
		<< " already initialied" << std::endl;
    return false;
  }
  ClearElements();

  for(int i=0;i<m_file_names.size();i++){
    if(!ReadFile(m_file_names[i])) return false;
  }
  m_is_ready = true;
  return true;
}

//______________________________________________________________________________
bool
HodoPHCMan::Initialize( const std::string& file_name1, const std::string& file_name2, const std::string& file_name3 )
{
  m_file_names.push_back(file_name1);
  m_file_names.push_back(file_name2);
  m_file_names.push_back(file_name3);
  return Initialize();
}
//______________________________________________________________________________
bool
HodoPHCMan::ReadFile( const std::string& file_name )
{
  if(file_name=="") return true;
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  std::ifstream ifs( file_name.c_str() );
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << func_name << " file open fail : "
		<< file_name << std::endl;
    return false;
  }
  std::string line;
  while( ifs.good() && std::getline(ifs,line) ){
    if( line.empty() || line[0]=='#' ) continue;
    if( TString(line).Contains("cid",TString::ECaseCompare::kIgnoreCase) ) continue;
    std::istringstream input_line( line );
    int cid=-1, plid=-1, seg=-1, ud=-1, type=-1, np=-1;
    std::vector<double> par;
    if( input_line >> cid >> plid >> seg >> ud >> type >> np ){
      double p = 0.;
      while( input_line >> p ) par.push_back(p);
      int key = MakeKey(cid,plid,seg,ud);
      HodoPHCParam *param = new HodoPHCParam(type,np,par);
      HodoPHCParam *pre_param = m_container[key];
      m_container[key] = param;
      if( pre_param ){
	hddaq::cerr << func_name << ": duplicated key "
		    << " following record is deleted." << std::endl
		    << " key = " << key << std::endl;
	delete pre_param;
      }
    }
    else{
      hddaq::cerr << func_name << ": Invalid format" << std::endl
		  << " ===> " << line << std::endl;
    }
  }
  return true;
}
//______________________________________________________________________________
bool
HodoPHCMan::DoCorrection( int cid, int plid, int seg, int ud,
			  double time, double de, double & ctime ) const
{
  ctime = time;
  HodoPHCParam* map = GetMap(cid,plid,seg,ud);
  if(!map) return false;
  ctime = map->DoPHC(time,de);
  return true;
}

//______________________________________________________________________________
bool
HodoPHCMan::DoRCorrection( int cid, int plid, int seg, int ud,
			   double time, double de, double & ctime ) const
{
  ctime = time;
  HodoPHCParam* map = GetMap(cid,plid,seg,ud);
  if(!map) return false;
  ctime = map->DoRPHC(time,de);
  return true;
}

//______________________________________________________________________________
HodoPHCParam*
HodoPHCMan::GetMap( int cid, int plid, int seg, int ud ) const
{
  int key = MakeKey(cid,plid,seg,ud);
  //  HodoPHCParam* map = 0;
  PhcPIterator  itr = m_container.find(key);
  if( itr != m_container.end() ) return itr->second;
  key = MakeKey(cid,plid,-1,ud);
  itr = m_container.find(key);
  if( itr != m_container.end() ) return itr->second;
  return 0;
}
