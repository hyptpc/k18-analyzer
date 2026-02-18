// DCTimeCorrMan.h

#ifndef DCTimeCorrMan_h
#define DCTimeCorrMan_h 1

#include <map>
#include <string>
#include <vector>

#include <TGraph.h>
#include <TString.h>

class DCTimeCorrMan
{
 public:
  static DCTimeCorrMan& GetInstance(void);
  static const TString& ClassName(void);
  ~DCTimeCorrMan();

  void SetFileName( const TString & filename1) { FileNameBLDC=filename1; }
  void SetFileNames( const TString & filename1, const TString & filename2 );
  Bool_t Initialize();
  Bool_t Initialize( const char *filename1, const char *filename2 );
  Bool_t Initialize( const TString& filename1, const TString& filename2 );
  Bool_t Initialize( const char *filename1 );
  Bool_t Initialize( const TString& filename1 );

  void SetFileNameCDC( const TString & filename );
  void SetFileNameBLDC( const TString & filename );
  
 private:
  DCTimeCorrMan();
  DCTimeCorrMan( const DCTimeCorrMan &right );
  TString FileNameCDC;
  TString FileNameBLDC;
  typedef std::map < Int_t, TGraph> DCTimeCorrHistContainer;
  DCTimeCorrHistContainer dctimecorrContainer;
  Bool_t m_isready;
  
 public:
  Bool_t IsReady() const {return m_isready; }
  void SetDCTimeCorrMan( const DCTimeCorrHistContainer container )  { dctimecorrContainer = container; }

  TString GetFileNameCDC() { return FileNameCDC; }
  TString GetFileNameBLDC() { return FileNameBLDC; }

  Double_t CalcCValue( const Int_t &cid, const Int_t &layer, const Int_t &wire,
		     const Double_t &timemean, const Double_t &timesub ) const;
  Double_t CalcDATValue( const Int_t &cid, const Int_t &layer, const Int_t &wire,
		       const Double_t &timemean, const Double_t &timesub ) const;
};
inline DCTimeCorrMan&
DCTimeCorrMan::GetInstance( void )
{
  static DCTimeCorrMan g_instance;
  return g_instance;
}
inline const TString&
DCTimeCorrMan::ClassName( void )
{
  static TString g_name("DCTimeCorrMan");
  return g_name;
}

#endif
