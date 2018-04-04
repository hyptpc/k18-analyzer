/**
 *  file: SsdParamMan.hh
 *  date: 2017.04.10
 *
 */

#ifndef SSD_PARAM_MAN_HH
#define SSD_PARAM_MAN_HH

#include <string>
#include <map>
#include <vector>
#include "DetectorID.hh"

//______________________________________________________________________________
class SsdAParam
{
public:
  inline SsdAParam( std::vector<double> pedestal, std::vector<double> rms )
    : m_pedestal(pedestal), m_rms(rms)
  {}
  ~SsdAParam( void )
  {}

private:
  SsdAParam( void );
  SsdAParam( const SsdAParam& );
  SsdAParam& operator =( const SsdAParam& );

private:
  std::vector<double> m_pedestal;
  std::vector<double> m_rms;

public:
  double pedestal( int i ) const { return m_pedestal[i]; }
  double rms( int i )      const { return m_rms[i]; }
  double de( int i, int adc ) const
  { return ((double)adc-m_pedestal[i]); }
  int    adc( int i, double de ) const
  { return (int)(de+m_pedestal[i]); }
};

//______________________________________________________________________________
class SsdParamMan
{
public:
  static SsdParamMan&       GetInstance( void );
  static const std::string& ClassName( void );
  ~SsdParamMan( void );

private:
  SsdParamMan( void );
  SsdParamMan( const SsdParamMan& );
  SsdParamMan& operator =( const SsdParamMan& );

  typedef std::map <int, SsdAParam *> AContainer;
  typedef std::map <int, SsdAParam *>::const_iterator AIterator;

private:
  bool        m_is_ready;
  std::string m_file_name;
  AContainer  m_APContainer;

public:
  void SetFileName( const std::string & file_name ) { m_file_name = file_name; }
  bool Initialize( void );
  bool Initialize( const std::string& file_name );
  bool IsReady( void ) const { return m_is_ready; }
  bool GetDe( int plid, int seg, int isample, int adc, double &de ) const;
  bool GetRms( int plid, int seg, int isample, double &rms ) const;
  bool GetAdc(int plid, int seg, int isample, double de, int &adc ) const;

private:
  void       ClearACont( void );
  SsdAParam* GetAmap( int plid, int seg ) const;

};

//______________________________________________________________________________
inline SsdParamMan&
SsdParamMan::GetInstance( void )
{
  static SsdParamMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
SsdParamMan::ClassName( void )
{
  static std::string g_name("SsdParamMan");
  return g_name;
}

#endif
