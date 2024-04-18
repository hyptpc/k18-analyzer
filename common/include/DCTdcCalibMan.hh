// -*- C++ -*-

#ifndef DC_TDC_CALIB_MAN_HH
#define DC_TDC_CALIB_MAN_HH

#include <map>
#include <vector>

#include <TString.h>

struct DCTdcCalMap;

//______________________________________________________________________________
class DCTdcCalibMan
{
public:
  static DCTdcCalibMan& GetInstance();
  static const TString& ClassName();
  ~DCTdcCalibMan();

private:
  DCTdcCalibMan();
  DCTdcCalibMan( const DCTdcCalibMan & );
  DCTdcCalibMan& operator =( const DCTdcCalibMan & );

private:
  typedef std::map <unsigned int, DCTdcCalMap*> DCTdcContainer;
  typedef DCTdcContainer::const_iterator DCTdcIterator;
  bool           m_is_ready;
  std::vector<TString>    m_file_names;
  DCTdcContainer m_container;

public:
  bool Initialize();
  bool Initialize( const TString& file_name );
  bool Initialize( const TString& file_name, const TString& file_name2 );
  bool ReadFile( const TString& file_name );
  bool IsReady() const { return m_is_ready; }
  bool GetTime( int det_id, int plane_id, int wire_id, int tdc, double& time ) const;
  bool GetTdc( int det_id, int plane_id, int wire_id, double time, int& tdc ) const;
  //  void SetFileName( const TString& file_name ) { m_file_name = file_name; }

private:
  DCTdcCalMap* GetMap( int det_id, int plane_id, int wire_id ) const;
  void         ClearElements();
};

//______________________________________________________________________________
inline DCTdcCalibMan&
DCTdcCalibMan::GetInstance()
{
  static DCTdcCalibMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const TString&
DCTdcCalibMan::ClassName()
{
  static TString g_name("DCTdcCalibMan");
  return g_name;
}

#endif
