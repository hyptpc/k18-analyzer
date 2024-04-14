/**
 *  file: SDDParamMan.hh
 *  date: 2017.04.10
 *
 */

#ifndef SDD_PARAM_MAN_HH
#define SDD_PARAM_MAN_HH

#include <string>
#include <map>
#include "TRandom3.h"

//______________________________________________________________________________
//SDD TDC to Time
class SDDTParam
{
public:
  inline SDDTParam( double offset, double gain )
    : m_offset(offset), m_gain(gain)
  {}
  ~SDDTParam( void )
  {}

private:
  SDDTParam( void );
  SDDTParam( const SDDTParam& );
  SDDTParam& operator =( const SDDTParam& );

private:
  double m_offset;
  double m_gain;

public:
  double Offset( void )  const { return m_offset; }
  double Gain( void )    const { return m_gain; }
  double Time( int tdc ) const
  { return ( (double)tdc - m_offset + gRandom->Uniform(-0.5,0.5) ) * m_gain; }
  int    Tdc( double time ) const
  { return (int)( time/m_gain + m_offset ); }
};

//______________________________________________________________________________
//SDD ADC to DeltaE
class SDDAParam
{
public:
  inline SDDAParam( double pedestal, double gain )
    : m_pedestal(pedestal), m_gain(gain)
  {}
  ~SDDAParam( void )
  {}

private:
  SDDAParam( void );
  SDDAParam( const SDDAParam& );
  SDDAParam& operator =( const SDDAParam& );

private:
  double m_pedestal;
  double m_gain;

public:
  double Pedestal( void )  const { return m_pedestal; }
  double Gain( void )      const { return m_gain; }
  double DeltaE( int adc ) const
  { return ( (double)adc + gRandom->Uniform(-0.5,0.5) - m_pedestal ) * (m_gain); }
  int    Adc( double de ) const
  { return (int)( m_gain * de + m_pedestal * ( 1. - de ) ); }
};

//______________________________________________________________________________
//SDD ADC with a correction by neighbouring pixels
class SDDCParam
{
public:
  inline SDDCParam( )
    : m_pedestal1(-1), m_pedestal2(-1), m_gain1(-1), m_gain2(-1)
  {}
  ~SDDCParam( void )
  {}

private:
  SDDCParam( const SDDCParam& );
  SDDCParam& operator =( const SDDCParam& );

private:
  double m_pedestal1;
  double m_pedestal2;
  double m_gain1;
  double m_gain2;

public:
  void SetParameter1(double ped1, double gain1){ m_pedestal1=ped1; m_gain1=gain1; }
  void SetParameter2(double ped2, double gain2){ m_pedestal2=ped2; m_gain2=gain2; }
  double CAdc( int adc, int adc1, int adc2 ) const
  { double  corr =0;
    if(m_gain1>0 ) corr += ( adc1 - m_pedestal1 ) * (m_gain1);
    if(m_gain2>0 ) corr += ( adc2 - m_pedestal2 ) * (m_gain2);
    return (double)adc-corr;
  }
};

//______________________________________________________________________________
//SDDParam Main Class
class SDDParamMan
{
public:
  static SDDParamMan&      GetInstance( void );
  static const std::string& ClassName( void );
  ~SDDParamMan( void );

private:
  SDDParamMan( void );
  SDDParamMan( const SDDParamMan& );
  SDDParamMan& operator =( const SDDParamMan& );

private:
  enum eAorT { kAdc, kTdc, kPre, kPost, kType };
  typedef std::map<int, SDDTParam*> TContainer;
  typedef std::map<int, SDDAParam*> AContainer;
  typedef std::map<int, SDDCParam*> CContainer;
  typedef TContainer::const_iterator TIterator;
  typedef AContainer::const_iterator AIterator;
  typedef CContainer::const_iterator CIterator;
  bool        m_is_ready;
  bool        m_is_corr;
  std::string m_file_name1;
  std::string m_file_name2;
  TContainer  m_TPContainer;
  AContainer  m_APContainer;
  CContainer  m_CPContainer;

public:
  bool Initialize( void );
  bool Initialize( const std::string& file_name );
  bool Initialize( const std::string& file_name1, const std::string &file_name2 );
  bool ReadFile( const std::string& file_name );
  bool IsReady( void ) const { return m_is_ready; }
  bool IsCorr( void ) const { return m_is_corr; }
  bool GetTime( int cid, int plid, int seg, int ud, int tdc, double &time ) const;
  bool GetDe( int cid, int plid, int seg, int ud, int adc, double &de ) const;
  bool GetTdc( int cid, int plid, int seg, int ud, double time, int &tdc ) const;
  bool GetAdc( int cid, int plid, int seg, int ud, double de, int &adc ) const;
  bool GetCAdc( int cid, int plid, int seg, int ud, int adc0, int adc1, int adc2, double &corr ) const;
  void SetFileName1( const std::string& file_name ) { m_file_name1 = file_name; }
  void SetFileName2( const std::string& file_name ) { m_file_name2 = file_name; }

private:
  SDDTParam *GetTmap( int cid, int plid, int seg, int ud ) const;
  SDDAParam *GetAmap( int cid, int plid, int seg, int ud ) const;
  SDDCParam *GetCmap( int cid, int plid, int seg, int ud ) const;
  void ClearACont( void );
  void ClearTCont( void );
  void ClearCCont( void );
};

//______________________________________________________________________________
inline SDDParamMan&
SDDParamMan::GetInstance( void )
{
  static SDDParamMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
SDDParamMan::ClassName( void )
{
  static std::string g_name("SDDParamMan");
  return g_name;
}

#endif
