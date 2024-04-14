// CounterMapMan.cc

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>

#include "CounterMapMan.hh"
#include "DetectorID.hh"
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
namespace{
  const unsigned int KEYMASK  = 0x0007;
  // |0000|0000|1111|0001|1111|0011|1111|0011|
  const unsigned int AMASK    = 0x003F;      /* A Mask 6 Bits (0-63) */
  const unsigned int NMASK    = 0x001F;      /* N Mask 5 Bits (0-31) */
  const unsigned int CMASK    = 0x000F;      /* C Mask 4 Bits (0-15) */
  const int          ASHIFT   =  4;
  const int          NSHIFT   = 12;
  const int          CSHIFT   = 20;
  const unsigned int KEYFLAG  = 0x0003;
  
  // |1111|1110|1 111|0111|1111|0111|1111|1 100|
  const unsigned int SEGMASK  = 0x00FF;    /* Segment  (0-255) */
  const unsigned int CIDMASK  = 0x007F;    /* CId      (0-127) */
  const unsigned int UDMASK   = 0x0007;    /* UD       (0-7)   */
  const unsigned int ATMASK   = 0x0003;    /* AT       (0-1)   */
  const unsigned int LAYMASK  = 0x007F;    /* Layer    (0-127)  */
  const int          SEGSHIFT =  3;
  const int          CIDSHIFT = 12;
  const int          UDSHIFT  = 20;
  const int          ATSHIFT  = 23;
  const int          LAYSHIFT = 25;
  const unsigned int RKEYFLAG = 0x0004;
  
  int KEY(const int &c, const int &n, const int &a)
  {  return ((((c)&CMASK)<<CSHIFT) | 
	     (((n)&NMASK)<<NSHIFT) | 
	     (((a)&AMASK)<<ASHIFT)| KEYFLAG ); 
  }
  int RKEY(const int &at,const int &ud,const int &cid,const int &seg,const int &lay)
  { return ((((at)&ATMASK)<<ATSHIFT) | 
	    (((ud)&UDMASK)<<UDSHIFT) | 
	    (((cid)&CIDMASK)<<CIDSHIFT) |
	    (((seg)&SEGMASK)<<SEGSHIFT) | 
	    (((lay)&LAYMASK)<<LAYSHIFT) | RKEYFLAG );
  } 
  const int MAXCHAR = 144;
}

CounterMapMan::CounterMapMan()
  :m_isready(false)
{
  FileName.clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CounterMapMan::Clear()
{
  fContainer.clear();
  bContainer.clear();
  fCrateDef.clear();
  bCrateDef.clear();
  nameCNAContainer.clear();
  nameCounterContainer.clear();
  NSMP = NCH_SCA = 0;
  for(int i=0;i<2;i++)
    for(int j=0;j<23;j++)
      CRATE_TYPE[i][j]=DRT;
  for(int j=0;j<23;j++)
    CRATE_TYPE[2][j]=NORMAL;
  for(int i=3;i<6;i++)
    for(int j=0;j<23;j++)
      CRATE_TYPE[i][j]=DRT;
  for(int i=6;i<10;i++)
    for(int j=0;j<23;j++)
      CRATE_TYPE[i][j]=NORMAL;
  for(int i=10;i<20;i++)
    for(int j=0;j<23;j++)
      CRATE_TYPE[i][j]=DRT;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CounterMapMan::Initialize( const char *file_name )
{
  FileName.clear();
  FileName.push_back(file_name);
  return Initialize();
}
bool CounterMapMan::Initialize( const std::string& file_name )
{
  FileName.clear();
  FileName.push_back(file_name);
  return Initialize();
}
bool CounterMapMan::Initialize( const std::string& file_name1, const std::string& file_name2 )
{
  FileName.clear();
  FileName.push_back(file_name1);
  FileName.push_back(file_name2);
  return Initialize();
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CounterMapMan::Initialize()
{
  static const TString funcname = "CounterMapMan::Initialize";
  //  std::cout << "[" << funcname << "] Initialization start ...";
  Clear();
  // at first, fill CID_TEMP
  unsigned int key, rkey;
  for( int cr=0; cr<15; cr++ ){
    for( int sl=1; sl<24; sl++ ){
      for( int ch=0; ch<32; ch++ ){
	key = KEY(cr,sl,ch);
	rkey = RKEY(0,0,127,0,0);
	fContainer[key] = rkey;
      }
    }
  }
  for(int i=0;i<nFiles();i++){
    ReadFile(FileName.at(i));
  }
#if 0
  for(int i=0;i<NSMP;i++){
    std::cout<<"CrateType"<<i<<"\t";
    for(int j=0;j<23;j++) std::cout<<CRATE_TYPE[i][j]<<" ";
    std::cout<<std::endl;
  }
#endif 
  m_isready=true;
  return true;
}

bool CounterMapMan::ReadFile( const TString filename )
{
  static const TString funcname = "CounterMapMan::ReadFile";
  int c,n,a,at,ud,cid,seg;
  int lay,wire;
  int address;
  int nd;
  int drtmin,drtmax,normmin,normmax;
  unsigned int key, rkey;
  char name[MAXCHAR];
  char str[MAXCHAR];
  FILE *fp;
  
  if( (fp=fopen(filename.Data(), "r"))==0 ){
    std::cerr << " File open fail. [" << filename << "]" << std::endl;
    exit(-1);
  }
  
  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;
    c=n=a=at=ud=cid=seg=lay=wire=0;
    
    if( (nd=sscanf(str,"CrateDef: %d %x", &c, &address)) == 2 ) {
#if 0
      std::cout << std::hex << " address: 0x" << address
		<< " crate: " << c << std::dec << std::endl;
#endif
      fCrateDef[c] = address;
      bCrateDef[address] = c;
    }
    else if( (nd=sscanf(str,"CrateType: %d MIX DRT %d %d NORMAL %d %d", &n, &drtmin, &drtmax, &normmin, &normmax)) == 5 ) {
      for(int i=drtmin;i<=drtmax;i++)
	CRATE_TYPE[n][i-1]=DRT;
      for(int i=normmin;i<=normmax;i++)
	CRATE_TYPE[n][i-1]=NORMAL;
    }
    else if( (nd=sscanf(str,"CrateType: %d %s", &n, name)) == 2 ) {
      if(strcmp(name,"DRT")==0)
	for(int i=0;i<23;i++)
	  CRATE_TYPE[n][i]=DRT;
      else if(strcmp(name,"NORMAL")==0)
	for(int i=0;i<23;i++)
	  CRATE_TYPE[n][i]=NORMAL;
    }
    else if( (nd=sscanf(str,"NumSMP: %d", &n)) == 1 ) {
      NSMP = n;
    }
    else if( (nd=sscanf(str,"NumScaler: %d", &n)) == 1 ) {
      NCH_SCA = n;
    }
    else if( (nd=sscanf(str,"%d %d %d %d %d %d %d %d %s", &c, &n, &a, &cid, &lay, &wire, &at, &ud, name)) == 9 ) {
      // drift chamber
#if 0
      if(c==10)
      std::cout << c << "  " << n << "  " << a << "  " << cid << "  "
		<< lay << "  " << wire << "  " << at << "  " << ud << "  " << name << std::endl;
#endif
      if( cid != DetIdCDC ){
	key = KEY(c,n,a);
	rkey = RKEY(at,ud,cid,wire,lay);
	fContainer[key] = rkey;
	bContainer[rkey] = key;
	nameCNAContainer[key] = name;
	nameCounterContainer[rkey] = name;
      }else{
	for(int ich=0;ich<16;ich++){
	  key = KEY(c,n,a+ich);
	  rkey = RKEY(at,ud,cid,ich,lay);
	  fContainer[key] = rkey;
	  bContainer[rkey] = key;
	  nameCNAContainer[key] = name;
	  nameCounterContainer[rkey] = name;	  
	}
      }
    }
    else if( (nd=sscanf(str,"%d %d %d %d %d %d %d %s", &c, &n, &a, &cid, &seg, &at, &ud, name)) == 8 ) {
      // counter
#if 0
      if(c==2&&n==1)
      std::cout << c << "  " << n << "  " << a << "  " << cid << "  "
		<< seg << "  " << at << "  " << ud << "  " << name << std::endl;
#endif
      key = KEY(c,n,a);
      rkey = RKEY(at,ud,cid,seg,0);
      fContainer[key] = rkey;
      bContainer[rkey] = key;
      nameCNAContainer[key] = name;
      nameCounterContainer[rkey] = name;
    }
    else if( (nd=sscanf(str,"%d %d %d %d %d %d %d", &c, &n, &a, &cid, &seg, &at, &ud)) == 7 ) {
#if 0
      std::cout << c << "  " << n << "  " << a << "  " << cid << "  "
		<< seg << "  " << at << "  " << ud << std::endl;
#endif
      // counter without name
      key = KEY(c,n,a);
      rkey = RKEY(at,ud,cid,seg,0);
      fContainer[key] = rkey;
      bContainer[rkey] = key;
      nameCNAContainer[key] = "";
      nameCounterContainer[rkey] = "";
    }
    else{
      std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
      std::cerr << std::string(str) << std::endl;
    }
  }


  fclose(fp);

  //PrintMap();
  //std::cout <</* "[" << funcname << "] Initialization*/" finish." << std::endl;
#if 0
  for(int i=0;i<NSMP;i++){
    std::cout<<"CrateType"<<i<<"\t";
    for(int j=0;j<23;j++) std::cout<<CRATE_TYPE[i][j]<<" ";
    std::cout<<std::endl;
  }
#endif 
  m_isready=true;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TString CounterMapMan::GetName( const int &c, const int &n, const int &a )
{
  unsigned int key = KEY(c,n,a);
  nameCNAMapContainer::const_iterator ic = nameCNAContainer.find(key);
  if( ic != nameCNAContainer.end() ) return ic->second;
  else                               return "";
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TString CounterMapMan::GetName( const int &cid, const int &seg, const int &at, const int &ud )
{
  unsigned int rkey = RKEY(at,ud,cid,seg,0);
  nameCounterMapContainer::const_iterator ic = nameCounterContainer.find(rkey);
  if( ic != nameCounterContainer.end() ) return ic->second;
  else                                   return "";
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TString CounterMapMan::GetName( const int &cid, const int &lay, const int &wire, const int &at, const int &ud )
{
  unsigned int rkey = RKEY(at,ud,cid,wire,lay);
  nameCounterMapContainer::const_iterator ic = nameCounterContainer.find(rkey);
  if( ic != nameCounterContainer.end() ) return ic->second;
  else                                   return "";
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int CounterMapMan::GetCrateNum( int address )
{
  return bCrateDef[address];
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
int CounterMapMan::GetSMPAddress( int c )
{
  return fCrateDef[c];
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CounterMapMan::GetCNA( int cid, int seg, int at, int ud, int &c, int &n, int &a )
{
  static const TString funcname = "CMapMan::GetCNA";
  unsigned int key, rkey;

  rkey = RKEY(at,ud,cid,seg,0);
  key  = bContainer[rkey];
  if ( (key&KEYMASK) == KEYFLAG ) {

    c = (key>>CSHIFT)&CMASK;
    n = (key>>NSHIFT)&NMASK;
    a = (key>>ASHIFT)&AMASK;
    return true;

  } else {
#if 0
    std::cerr << "[" << funcname << "]: No address found ..."
	      " CID=" << cid << " SEG=" << seg
              << " AT=" << at << " UD=" << ud << std::endl;
#endif
    return false;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CounterMapMan::GetInfo( int c, int n, int a, int &cid, int &lay, int &seg, int &at, int &ud )
{
  static const TString funcname = "CounterMapMan::GetInfo";
  unsigned int key, rkey;

  key = KEY(c,n,a);
  rkey = fContainer[key];
  if( (rkey&KEYMASK)==RKEYFLAG ){
    cid = (rkey>>CIDSHIFT)&CIDMASK;
    seg = (rkey>>SEGSHIFT)&SEGMASK;
    lay = (rkey>>LAYSHIFT)&LAYMASK;
    at  = (rkey>>ATSHIFT )&ATMASK;
    ud  = (rkey>>UDSHIFT )&UDMASK;
    return true;
  }
  else{
    std::cerr << "[" << funcname << "]: No address found ... C=" 
	      << c << " N=" << n << " A=" << a << std::endl;
    return false;
  }
}

int CounterMapMan::GetCID( int c, int n, int a )
{
  static const TString funcname = "CounterMapMan::GetCID";
  unsigned int key, rkey;

  key = KEY(c,n,a);
  rkey = fContainer[key];
  if( (rkey&KEYMASK)==RKEYFLAG ){
    return ((rkey>>CIDSHIFT)&CIDMASK);
  }
  else{
    std::cerr << "[" << funcname << "]: No address found ... C=" 
	      << c << " N=" << n << " A=" << a << std::endl;
    return -1;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool CounterMapMan::GetCNA( int cid, int lay, int wire, int at, int ud, int &c, int &n, int &a )
{
  static const TString funcname = "CMapMan::GetCNA";
  unsigned int key, rkey;

  rkey = RKEY(at,ud,cid,wire,lay);
  key  = bContainer[rkey];
  if ( (key&KEYMASK) == KEYFLAG ) {

    c = (key>>CSHIFT)&CMASK;
    n = (key>>NSHIFT)&NMASK;
    a = (key>>ASHIFT)&AMASK;
    return true;

  } else {
    std::cerr << "[" << funcname << "]: No address found ..."
	      << " CID=" << cid << " LAYER=" << lay
              << " WIRE=" << wire << " AT:" << at << " UD:" << ud << std::endl;
    c = n = a = -1;
    return false;

  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CounterMapMan::PrintSimpleMap( std::ostream &p_out )
{
  static const TString funcname = "CMapMan::PrintMap";
  unsigned int key, rkey;

  std::cout << " ---- " << funcname << " ---- " << std::endl;
  for( fCounterMapContainer::const_iterator i=fContainer.begin();
       i!=fContainer.end(); i++ ){
    key = i->first; rkey = i->second;
    int cr = ((key>>CSHIFT)&CMASK);
    if( cr==0 || cr==1 || cr==2 || cr==6 ){
      p_out.setf(std::ios::showpoint);
      p_out  << std::setw(5) << cr
	     << std::setw(5) << ((key>>NSHIFT)&NMASK)
	     << std::setw(5) << ((key>>ASHIFT)&AMASK)
	     << std::setw(5) << ((rkey>>CIDSHIFT)&CIDMASK)
	     << std::setw(5) << ((rkey>>LAYSHIFT)&LAYMASK)
	     << std::setw(5) << ((rkey>>SEGSHIFT)&SEGMASK)
	     << std::setw(5) << ((rkey>>ATSHIFT )&ATMASK)
	     << std::setw(5) << ((rkey>>UDSHIFT )&UDMASK)
	     << std::endl;
    }
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void CounterMapMan::PrintMap()
{
  static const TString funcname = "CMapMan::PrintMap";
  unsigned int key, rkey;

  std::cout << " ---- " << funcname << " ---- " << std::endl;
  for( fCounterMapContainer::const_iterator i=fContainer.begin();
       i!=fContainer.end(); i++ ){
    key = i->first; rkey = i->second;
    std::cout << key  << "  "
	      << ((key>>CSHIFT)&CMASK) << " "
	      << ((key>>NSHIFT)&NMASK) << " "
	      << ((key>>ASHIFT)&AMASK) << " "
	      << rkey << "  "
	      << ((rkey>>CIDSHIFT)&CIDMASK) << " "
	      << ((rkey>>LAYSHIFT)&LAYMASK) << " "
	      << ((rkey>>SEGSHIFT)&SEGMASK) << " "
	      << ((rkey>>ATSHIFT )&ATMASK) << " "
	      << ((rkey>>UDSHIFT )&UDMASK)
	      << std::endl;
  }
}

