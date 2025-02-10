// XTMapMan.h

#ifndef XTMapMan_h
#define XTMapMan_h 1

#include <map>
#include <string>
#include <vector>
#include <iostream>
#include <stdlib.h>

#include "TH1.h"
#include "TGraph.h"
#include "TROOT.h"

class XTMap
{
 public:
  XTMap();
  ~XTMap() {};
 private:
  std::vector <double> Par;
 public:
  void SetParam( const double &par ) { Par.push_back(par); }
  int nParam() const { return Par.size(); }
  std::vector <double> *GetParam() { return &Par; }
  double GetParam( int i ) const; 
  void Clear();
};

class XTMapMan
{
public:
  static XTMapMan& GetInstance(void);
  static const std::string& ClassName( void );
  ~XTMapMan();

public:
  void SetFileNames( const TString & filename1, const TString & filename2 );
  bool Initialize();
  bool Initialize( const char *filename1, const char *filename2 );
  bool Initialize( const std::string& filename1, const std::string& filename2 );
  bool Initialize( const char *filename1 );
  bool Initialize( const std::string& filename1 );
  bool InitializeBLDC( void );
  bool InitializeCDC( void );
  
private:
  XTMapMan();
  XTMapMan( const XTMapMan &right );
  
private:
  TString FileNameBLDC;
  TString FileNameCDC;
  typedef std::map <unsigned int, TGraph > TGraphContainer;
  typedef std::map <unsigned int, XTMap>   XTMapContainer;
  TGraphContainer  tgContainer;
  XTMapContainer   xtContainer;
  
  bool m_isready;

 public:
  bool IsReady() const {return m_isready; }
  void SetGainMapMan( const XTMapContainer container )  { xtContainer = container; }  
  TString GetFileNameBLDC() { return FileNameBLDC; }
  TString GetFileNameCDC() { return FileNameCDC; }
  double CalcDriftLength( const int &cid, const int &layer, const int &wire, const double &dt ) const;

  //  bool SetParam( const int &cid, const int &layer, const int &wire, const int &npar, double *par );
  TGraph* GetGraph( const int &cid, const int &layer, const int &wire);

  void PrintMapHeader( std::ostream &p_out = std::cout );
  void PrintMap( const int &Cid = -1, std::ostream &p_out = std::cout );
 
  // int nparam( const int &cid, const int &layer, const int &wire ) const;
  // double param( const int &cid, const int &layer, const int &wire, const int &i ) const;

  // simulation
  double CalcDxDt( const int &cid, const int &layer, const int &wire, const double &dt ) const;
  // double CalcDriftTime( const int &cid, const int &layer, const int &wire, const double &dl );
};
inline XTMapMan&
XTMapMan::GetInstance( void )
{
  static XTMapMan g_instance;
  return g_instance;
}
inline const std::string&
XTMapMan::ClassName( void )
{
  static std::string g_name("XTMapMan");
  return g_name;
}

#endif
