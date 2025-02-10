// DCTimeCorrMan.h

#ifndef DCTimeCorrMan_h
#define DCTimeCorrMan_h 1

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>

#include "TGraph.h"
#include "TFile.h"
#include "TROOT.h"

class DCTimeCorrMan
{
 public:
  static DCTimeCorrMan& GetInstance(void);
  static const std::string& ClassName(void);
  ~DCTimeCorrMan();

  void SetFileName( const TString & filename1) { FileNameBLDC=filename1; }
  void SetFileNames( const TString & filename1, const TString & filename2 );
  bool Initialize();
  bool Initialize( const char *filename1, const char *filename2 );
  bool Initialize( const std::string& filename1, const std::string& filename2 );
  bool Initialize( const char *filename1 );
  bool Initialize( const std::string& filename1 );

  void SetFileNameCDC( const TString & filename );
  void SetFileNameBLDC( const TString & filename );
  
 private:
  DCTimeCorrMan();
  DCTimeCorrMan( const DCTimeCorrMan &right );
  TString FileNameCDC;
  TString FileNameBLDC;
  typedef std::map < int, TGraph> DCTimeCorrHistContainer;
  DCTimeCorrHistContainer dctimecorrContainer;
  bool m_isready;
  
 public:
  bool IsReady() const {return m_isready; }
  void SetDCTimeCorrMan( const DCTimeCorrHistContainer container )  { dctimecorrContainer = container; }

  TString GetFileNameCDC() { return FileNameCDC; }
  TString GetFileNameBLDC() { return FileNameBLDC; }

  double CalcCValue( const int &cid, const int &layer, const int &wire,
		     const double &timemean, const double &timesub ) const;
  double CalcDATValue( const int &cid, const int &layer, const int &wire,
		       const double &timemean, const double &timesub ) const;
};
inline DCTimeCorrMan&
DCTimeCorrMan::GetInstance( void )
{
  static DCTimeCorrMan g_instance;
  return g_instance;
}
inline const std::string&
DCTimeCorrMan::ClassName( void )
{
  static std::string g_name("DCTimeCorrMan");
  return g_name;
}

#endif
