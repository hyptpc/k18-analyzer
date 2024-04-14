// CounterMapMan.h
#ifndef CounterMapMan_h
#define CounterMapMan_h 1

#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "TString.h"

class CounterMapMan
{
private:
  CounterMapMan( void );
public:
  static CounterMapMan& GetInstance(void);
  static const std::string& ClassName( void );
  ~CounterMapMan(void) {}
  bool Initialize();
  bool Initialize( const char *file_name );
  bool Initialize( const std::string& file_name );
  bool Initialize( const std::string& file_name,const std::string& file_name2 );
  bool ReadFile(const TString filename);
  
 private:
  static const int MMAXSMP = 20;
  static const int MAXSLOT = 23;
  static const int NORMAL  = 1;
  static const int DRT     = 0;

  std::vector<TString> FileName;
  int NSMP;
  int NCH_SCA;
  int CRATE_TYPE[MMAXSMP][MAXSLOT];

  typedef std::map <unsigned int, unsigned int> fCounterMapContainer;
  typedef std::map <unsigned int, unsigned int> bCounterMapContainer;
  fCounterMapContainer fContainer;
  bCounterMapContainer bContainer;

  typedef std::map <unsigned int, unsigned int> fCrateDefContainer;
  typedef std::map <unsigned int, unsigned int> bCrateDefContainer;
  fCrateDefContainer fCrateDef;
  bCrateDefContainer bCrateDef;

  typedef std::map <unsigned int, std::string> nameCNAMapContainer;
  typedef std::map <unsigned int, std::string> nameCounterMapContainer;
  nameCNAMapContainer nameCNAContainer;
  nameCounterMapContainer nameCounterContainer;

  bool m_isready;

 public:
  bool IsReady( void ) const { return m_isready; }
  int nFiles() { return FileName.size(); }
  TString GetFileName( const int &i ) { return FileName.at(i); }
  bool GetInfo( int c, int n, int a, int &cid, int &lay, int &seg, int &at, int &ud );
  int GetCID( int c, int n, int a );
  TString GetName( const int &c, const int &n, const int &a );

  // for Hodoscope or Cherenkov
  bool GetCNA( int cid, int seg, int at, int ud, int &c, int &n, int &a );
  TString GetName( const int &cid, const int &seg, const int &at, const int &ud);
  // for DC
  bool GetCNA( int cid, int layer, int wire, int at, int ud, int &c, int &n, int &a );
  TString GetName( const int &cid, const int &lay, const int &wire, const int &at, const int &ud );

  int GetCrateNum( int address );

  int GetSMPAddress( int c );
  int GetNumSMP() { return NSMP; }
  int GetNumScaler() { return NCH_SCA; }
  int GetCrateType( const int &cr, const int &sl ) { return CRATE_TYPE[cr][sl-1]; } //sl 1 origin
  void PrintSimpleMap( std::ostream &p_out = std::cout );
  void PrintMap();

  void Clear(); 
};
//______________________________________________________________________________
inline CounterMapMan&
CounterMapMan::GetInstance( void )
{
  static CounterMapMan g_instance;
  return g_instance;
}
inline const std::string&
CounterMapMan::ClassName( void )
{
  static std::string g_name("CounterMapMan");
  return g_name;
}
#endif
