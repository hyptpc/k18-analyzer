// ReslMapMan.cpp
#include "ReslMapMan.hh"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include <TRandom.h>
namespace{
  const std::string& class_name("ReslMapMan");
  const unsigned int KEYMASK  = 0x000F;
  // |0111|1111|0001|1111|0000|1111|1111|0011|
  const unsigned int WMASK    = 0x00FF;      /* Wire Mask 6 Bits (0-255) */
  const unsigned int LMASK    = 0x001F;      /* Layer Mask 5 Bits (0-31) */
  const unsigned int CMASK    = 0x007F;      /* CID Mask 7 Bits (0-31) */
  const int          WSHIFT   =  4;
  const int          LSHIFT   = 16;
  const int          CSHIFT   = 24;
  const unsigned int KEYFLAG  = 0x0003;
}
inline int 
KEY(const int &cid,const int &layer,const int &wire)
{
  return ( ( (cid&CMASK)   << CSHIFT ) |
	   ( (layer&LMASK) << LSHIFT ) |
	   ( (wire&WMASK)  << WSHIFT ) |
	   KEYFLAG );
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
ReslMapMan::ReslMapMan()
  :m_is_ready(false)
{
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ReslMapMan::Initialize( const char *file_name )
{
  FileName=file_name;
  return Initialize();
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ReslMapMan::Initialize( const std::string& file_name )
{
  FileName=file_name;
  return Initialize();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ReslMapMan::Initialize()
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  if( m_is_ready ){
    std::cerr << "#W " << func_name
	      << " already initialied" << std::endl;
    return false;
  }
  std::ifstream ifs(FileName.Data());
  if(!ifs.good()){
    std::cerr << " File open fail. [" << FileName << "]" << std::endl;
    return false;
  }
  reslContainer.clear();

  int invalid=0;
  std::string line;
  while( ifs.good() && std::getline( ifs, line ) ){
    ++invalid;
    if( line[0]=='#' || line.empty() ) continue;
    if( TString(line).Contains("cid",TString::ECaseCompare::kIgnoreCase) ) continue;
    std::istringstream input_line( line );
    int cid=-1, plid=-1, seg=-1;
    double tresl=-1, eresl=-1;
    if( input_line >> cid >> plid >> seg >> tresl >> eresl ){
      unsigned int key = KEY(cid,plid,seg);
      ReslMap *rmap = new ReslMap();
      rmap->SetTResl(tresl);
      rmap->SetEResl(eresl);
      reslContainer[key] = *rmap;
      delete rmap;
    }
    else {
      std::cerr << func_name << ": Invalid Input" << std::endl
		<< " ===> (" << invalid << "b)" << line << " " << std::endl;
    } /* if( input_line >> ) */
  }
  //  PrintMap();
  m_is_ready = true;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ReslMapMan::SetParam( const int &cid, const int &layer, const int &wire, const double &tresl, const double &eresl )
{
  static const TString funcname = "ReslMapMan::SetParam";
  unsigned int key;

  key = KEY(cid,layer,wire);
  ReslContainer::iterator ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){ // replace parameter
    (ir->second).SetTResl(tresl);
    (ir->second).SetEResl(eresl);
    //(ir->second) = resl;
    return true;
  }
  else{ // add parameter
      ReslMap *rmap = new ReslMap();
      rmap->SetTResl(tresl); rmap->SetEResl(eresl);
      reslContainer[key] = *rmap;
      delete rmap;
    //reslContainer[key] = resl;
  }
  return 0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ReslMapMan::GetResolution( const int &cid, const int &layer, const int &wire, double &tresl, double &eresl ) const
{
  static const TString funcname = "ReslMapMan::CalcResolution";
  unsigned int key;
  ReslContainer::const_iterator ir;

  // wire# != 0
  key = KEY(cid,layer,wire);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = gRandom->Gaus(0,(ir->second).tresl() );
    eresl = gRandom->Gaus(0,(ir->second).eresl() );
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return;
  }

  // wire# = 0
  key = KEY(cid,layer,0);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = gRandom->Gaus(0,(ir->second).tresl() );
    eresl = gRandom->Gaus(0,(ir->second).eresl() );
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return;
  }

  // layer# = 0
  key = KEY(cid,0,0);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = gRandom->Gaus(0,(ir->second).tresl() );
    eresl = gRandom->Gaus(0,(ir->second).eresl() );
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return;
  }

  // no param
  tresl = eresl = 0;
  return;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool ReslMapMan::GetParam( const int &cid, const int &layer, const int &wire, double &tresl, double &eresl )
{
  static const TString funcname = "ReslMapMan::GetParam";
  unsigned int key;
  ReslContainer::iterator ir;
  //  std::cout<<cid<<"  "<<layer<<"  "<<wire<<std::endl;
  // wire# != 0
  key = KEY(cid,layer,wire);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = (ir->second).tresl();
    eresl = (ir->second).eresl();
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return true;
  }

  // wire# = 0
  key = KEY(cid,layer,0);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = (ir->second).tresl();
    eresl = (ir->second).eresl();
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return true;
  }

  // layer# = 0
  key = KEY(cid,0,0);
  ir = reslContainer.find(key);
  if( ir != reslContainer.end() ){
    tresl = (ir->second).tresl();
    eresl = (ir->second).eresl();
    //    std::cout<<"tresl eresl :"<<tresl<<" "<<eresl<<std::endl;
    return true;
  }

  // no param
  tresl = eresl = 0;
  return false;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void ReslMapMan::PrintMap( std::ostream &p_out ) const 
{
  static const TString funcname = "ReslMapMan::PrintMap";
  unsigned int key;
  ReslMap rmap;
  std::cout << " ---- " << funcname << " ---- " << std::endl;
  for( ReslContainer::const_iterator ir=reslContainer.begin();
       ir!=reslContainer.end(); ir++ ){
    key = ir->first;
    rmap = ir->second;
    p_out << std::setw(5) << ((key>>CSHIFT)&CMASK)
	  << std::setw(8) << ((key>>LSHIFT)&LMASK)
	  << std::setw(8) << ((key>>WSHIFT)&WMASK)
	  << std::setw(10) << rmap.tresl()
	  << std::setw(10) << rmap.eresl()
	  << std::endl;
  }
}
