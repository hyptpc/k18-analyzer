// XTMapMan.cpp

#include <string>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <new>

#include <std_ostream.hh>

#include "DetectorID.hh"
#include "XTMapMan.hh"
#include "TFile.h"

namespace
{
  const std::string& class_name("XTMapMan");
  const unsigned int KEYMASK  = 0x000F;
  // |0111|1111|0001|1111|0000|1111|1111|0011|
  const unsigned int WMASK    = 0x00FF;      /* Wire Mask 6 Bits (0-255) */
  const unsigned int LMASK    = 0x001F;      /* Layer Mask 5 Bits (0-31) */
  const unsigned int CMASK    = 0x007F;      /* CID Mask 7 Bits (0-31) */
  const int          WSHIFT   =  4;
  const int          LSHIFT   = 16;
  const int          CSHIFT   = 24;
  const unsigned int KEYFLAG  = 0x0003;
  const int MAXCHAR = 512;
  const int MAXTOKEN = 30;
  const double unit=0.1; //mm->cm
}

inline int 
KEY( int cid,int layer, int wire)
{
  return ( ( (cid  &CMASK) << CSHIFT) | 
	   ( (layer&LMASK) << LSHIFT) | 
	   ( (wire &WMASK) << WSHIFT) |
	   KEYFLAG );
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
XTMap::XTMap()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double XTMap::GetParam( int i ) const 
{
  if( 0<=i && i<(int)Par.size() )
    return Par[i];
  else
    return -999.;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMap::Clear()
{
  Par.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
XTMapMan::XTMapMan()
  :m_isready(false)
{
  FileNameBLDC = "";
  FileNameCDC = "";
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
XTMapMan::XTMapMan( const XTMapMan &right )
{
  FileNameBLDC  = right.FileNameBLDC;
  FileNameCDC   = right.FileNameCDC;
  for( XTMapContainer::const_iterator i=right.xtContainer.begin();
       i!=right.xtContainer.end(); i++ ){
    xtContainer[i->first] = i->second;
  }
  for( TGraphContainer::const_iterator i=right.tgContainer.begin();
       i!=right.tgContainer.end(); i++ ){
    tgContainer[i->first] = i->second;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
XTMapMan::~XTMapMan()
{
  xtContainer.clear();
  tgContainer.clear();
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::SetFileNames( const TString & filename1, const TString & filename2 )
{
  FileNameBLDC = filename1;
  FileNameCDC = filename2;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool XTMapMan::Initialize( const char *filename1, const char *filename2 )
{
  FileNameBLDC=filename1;
  FileNameCDC=filename2;
  return Initialize();
}
bool XTMapMan::Initialize( const std::string& filename1, const std::string& filename2 )
{
  FileNameBLDC=filename1;
  FileNameCDC=filename2;
  return Initialize();
}
bool XTMapMan::Initialize( const char *filename1 )
{
  FileNameBLDC=filename1;
  return Initialize();
}
bool XTMapMan::Initialize( const std::string& filename1 )
{
  FileNameBLDC=filename1;
  return Initialize();
}

bool XTMapMan::Initialize()
{
  InitializeBLDC();
  InitializeCDC();
  m_isready=true;
  return true;
}
bool XTMapMan::InitializeBLDC()
{
  static const TString funcname = "XTMapMan::InitializeBLDC";
  //  std::cout << "[" << funcname << "] Initialization start ... "<< std::endl;
  int cid;
  int layer,wire,npar;
  double cellsize;
  int nd;
  unsigned int key;
  char str[MAXCHAR];
  FILE *fp;

  tgContainer.clear();
  std::vector<double> threshold;
  const char* DELIMITER = " ";
  if( FileNameBLDC!="" ){
    if( (fp=fopen(FileNameBLDC.Data(), "r"))==0 ){
      std::cerr << " File open fail. [" << FileNameBLDC << "]" << std::endl;
      exit(-1);
    }    
    while( fgets(str,MAXCHAR,fp)!=0 ){
      if( str[0]=='#' ) continue;
      if( (nd=sscanf(str,"XTParam: %d %lf %d ",&cid,&cellsize,&npar))==3 ){
	cellsize*=unit;
      }else{
	const char* token[MAXTOKEN] = {};
	int n=0;
	token[0] = strtok(str, DELIMITER);
	if (token[0] == NULL) continue;
	if( !strcmp(token[0], "#") ) continue;
	for (n = 1; n < MAXTOKEN; n++){
	  token[n] = strtok( NULL, DELIMITER);
	  if (token[n] == NULL ) break;
	}
	if( !strcmp(token[0], "Threshold:") ){
	  threshold.clear();
	  for(int i=0;i<npar;i++){
	    double tmpth=atof(token[i+1])*cellsize/100.;
	    threshold.push_back(tmpth);
	  }
	}else if(n==npar+3){
	  if(atoi(token[0])!=cid){
	    std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
	    std::cerr << std::string(str) << std::endl;
	    exit(-1);
	  }
	  layer=atoi(token[1]);
	  wire=atoi(token[2]);
	  if(wire==0) wire=-1; // temporal 20211212. parameter file should be changed.
	  key = KEY(cid,layer,wire);
	  TGraph* gr=new TGraph(npar);
	  for (int i = 0; i < npar; i++){
	    gr->SetPoint(i,atof(token[3+i]),threshold[i]);
	  }
	  tgContainer[key] = *gr;
	  gr->SetName(Form("xt_%d_%d_%d",cid,layer,wire));
	  gROOT->Append(gr);
	}
      }
    }
    fclose(fp);
  }
  //  std::cout <</* "[" << funcname << "] Initialization*/" finish." << std::endl;
  return true;
}

bool XTMapMan::InitializeCDC()
{
  static const std::string func_name("["+class_name+"::"+__func__+"()]");
  //  std::cout << func_name << " Initialization start ..." << std::endl;
  std::ifstream ifs( FileNameCDC.Data() );
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << func_name << " file open fail : "
		<< FileNameCDC << std::endl;
    return false;
  }
  
  xtContainer.clear();
  std::string line;
  while( ifs.good() && std::getline( ifs, line ) ){
    if( line[0]=='#' || line.empty() ) continue;
    if( TString(line).Contains("cid",TString::ECaseCompare::kIgnoreCase) ) continue;
    std::istringstream input_line( line );
    int cid=-1,layer=-1,wire=-1,npar=-1;
    double par[10];
    if( input_line >> cid >> layer >> wire >> npar ){
      int key = KEY(cid,layer,wire);
      XTMap *xtmap = new XTMap();
      //std::cout<<cid<<"  "<<layer<<"  "<<wire<<" "<<npar<<", params= ";
      for( int i=0; i<npar; i++ ){
	input_line >> par[i];
	//std::cout<<par[i]<<"  ";
	xtmap->SetParam(par[i]);
      }
      //      std::cout<<std::endl;
      xtContainer[key] = *xtmap;
      delete xtmap;
    }
  }
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
// bool XTMapMan::SetParam( const int &cid, const int &layer, const int &wire, const int &npar, double *par )
// {
//   static const TString funcname = "XTMapMan::SetParam";
//   unsigned int key;
//   key = KEY(cid,layer,wire);
//   XTMapContainer::iterator ix = xtContainer.find(key);
//   return false;
// }

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double XTMapMan::CalcDriftLength( const int &cid, const int &layer, const int &wire, const double &dt ) const
{
  static const TString funcname = "XTMapMan::CalcDriftLength";
  unsigned int key;
  double dl = 0.;
  if(cid==DetIdCDC){
    key = KEY(cid,layer,wire);
    XTMapContainer::const_iterator ix=xtContainer.find(key);
    if( ix == xtContainer.end() ){
      key = KEY(cid,layer,-1);
      ix = xtContainer.find(key);
    }
    if( ix != xtContainer.end() ){
      for( int i=0; i<(ix->second).nParam(); i++ ){
	double t = 1.;
	for( int j=0; j<i; j++ ){ t *= dt; }
	dl += ((ix->second).GetParam(i))*t;
      }
      return dl;
    }
  }else{
    key = KEY(cid,layer,wire);
    TGraphContainer::const_iterator ix=tgContainer.find(key);
    if( ix == tgContainer.end() ){
      key = KEY(cid,layer,-1);
      ix = tgContainer.find(key);
    }
    if( ix != tgContainer.end() ){
      dl = (ix->second).Eval(dt);
      return dl;
    }
  }
  hddaq::cerr << "#E "
	      << " no parameter found at XTMapMan::CalcDriftLength( Cid="
	      << cid  << ", Layer=" << layer << ", Wire=" << wire
	      << " )" << std::endl;     
  return dl;
}

// // + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double XTMapMan::CalcDxDt( const int &cid, const int &layer, const int &wire, const double &dt ) const
{
  // differential dt - dx function at dt
  static const TString funcname = "XTMapMan::CalcDxDt";
  unsigned int key;
  double dxdt = 0.;
  if(cid==DetIdCDC){
    key = KEY(cid,layer,wire);
    XTMapContainer::const_iterator ix=xtContainer.find(key);
    if( ix == xtContainer.end() ){
      key = KEY(cid,layer,-1);
      ix = xtContainer.find(key);
    }
    if( ix != xtContainer.end() ){
      for( int i=0; i<(ix->second).nParam()-1; i++ ){
	double t = 1.;
	for( int j=0; j<i; j++ ){ t *= dt; }
	dxdt += (i+1)*((ix->second).GetParam(i+1))*t;
      }
      return dxdt;
    }
  }else{
    key = KEY(cid,layer,wire);
    TGraphContainer::const_iterator ix=tgContainer.find(key);
    if( ix == tgContainer.end() ){
      key = KEY(cid,layer,-1);
      ix = tgContainer.find(key);
    }
    if( ix != tgContainer.end() ){
      dxdt = (ix->second).Eval(dt);
      return dxdt;
    }
  }
  return dxdt;
}

// double XTMapMan::CalcDriftTime( const int &cid, const int &layer, const int &wire, const double &dl ) const
// {
//   return dl;
// }


// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TGraph* XTMapMan::GetGraph( const int &cid, const int &layer, const int &wire )
{
  unsigned int key;
  TGraphContainer::iterator itg;
  key = KEY(cid,layer,wire);
  itg = tgContainer.find(key);
  if( itg != tgContainer.end() ){
    return &(itg->second);
  }
  return 0;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
// int XTMapMan::nparam( const int &cid, const int &layer, const int &wire ) const 
// {
//   static const TString funcname = "XTMapMan::nparam";
//   unsigned int key;

//   // // wire# != 0
//   // key = KEY(cid,layer,wire);
//   // XTMapContainer::iterator ix = xtContainer.find(key);
//   // if( ix != xtContainer.end() ){
//   //   return (ix->second).nParam();
//   // }

//   // // wire# = 0
//   // key = KEY(cid,layer,0);
//   // ix = xtContainer.find(key);
//   // if( ix != xtContainer.end() ){
//   //   return (ix->second).nParam();
//   // }
//   return 0;
// }
// // + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
// double XTMapMan::param( const int &cid, const int &layer, const int &wire, const int &i )
// {
//   static const TString funcname = "XTMapMan::param";
//   unsigned int key;

//   // // wire# != 0
//   // key = KEY(cid,layer,wire);
//   // XTMapContainer::iterator ix = xtContainer.find(key);
//   // if( ix != xtContainer.end() ){
//   //   return (ix->second).GetParam(i);
//   // }

//   // // wire# = 0
//   // key = KEY(cid,layer,0);
//   // ix = xtContainer.find(key);
//   // if( ix != xtContainer.end() ){
//   //   return (ix->second).GetParam(i);
//   // }
//   return 0;
// }

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::PrintMapHeader( std::ostream &p_out )
{
  static const TString funcname = "XTMapMan::PrintMapHeader";

  std::cout << " ---- " << funcname << " ---- " << std::endl;
  p_out << "#" << std::endl
	<< "# XT curve for Drift chamers." << std::endl
	<< "#" << std::endl
	<< "# units : [nsec] --> [mm]" << std::endl
	<< "# At first, XT relations are represented by n-th polynomials." << std::endl
	<< "# dx = p0 + p1*dt + p2*dt**2 + p3*dt**3 + p4*dt**4 + p5*dt**5" << std::endl
	<< "#" << std::endl
	<< "# The wire#=0 means that the value is applied to all wires." << std::endl
	<< "# If wire#!=0 exists after wire#=0, the value is uploaded. " << std::endl
	<< "#" << std::endl;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void XTMapMan::PrintMap( const int &Cid, std::ostream &p_out )
{
  // static const TString funcname = "XTMapMan::PrintMap";
  // unsigned int key;
  // TGraph gr;
  // int cid, cid_old=-1;

  // std::cout << " ---- " << funcname << " ---- " << std::endl;
  // //  p_out<<"XTParam: "<<Cid<<" 
  // for( XTMapContainer2::const_iterator i=xtContainer2.begin();
  //      i!=xtContainer2.end(); i++ ){
  //   key = i->first;
  //   gr = i->second;
  //   cid = ((key>>CSHIFT)&CMASK);
  //   if( !( cid==Cid || Cid==-1 ) ) continue;
  //   if( cid!=cid_old ){
  //     p_out	<< "# CID  Layer    Wire    npar            p0             p1             p2             p3             p4             p5" << std::endl;
  //     cid_old = cid;
  //   }
  //   p_out << std::setw(5) << cid
  // 	  << std::setw(8) << ((key>>LSHIFT)&LMASK)
  // 	  << std::setw(8) << ((key>>WSHIFT)&WMASK)
  //     	  << std::setw(7) << xtmap.nParam();
  //   p_out.setf(std::ios::showpoint);
  //   p_out.setf(std::ios::scientific, std::ios::floatfield);
  //   for( int j=0; j<xtmap.nParam(); j++ ){
  //     p_out  << std::setprecision(5) << std::setw(15) << xtmap.GetParam(j);
  //   }
  //   p_out << std::endl;
  // }
}

