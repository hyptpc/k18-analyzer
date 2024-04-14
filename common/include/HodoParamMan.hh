/**
 *  file: HodoParamMan.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO_PARAM_MAN_HH
#define HODO_PARAM_MAN_HH

#include <string>
#include <map>
#include "TRandom3.h"

//______________________________________________________________________________
//Hodo TDC to Time
class HodoTParam
{
public:
  inline HodoTParam( double offset, double gain )
    : m_offset(offset), m_gain(gain)
  {}
  ~HodoTParam( void )
  {}

private:
  HodoTParam( void );
  HodoTParam( const HodoTParam& );
  HodoTParam& operator =( const HodoTParam& );

private:
  double m_offset;
  double m_gain;

public:
  double Offset( void )  const { return m_offset; }
  double Gain( void )    const { return m_gain; }
  double Time( int tdc, bool RAND=false ) const
  { double tmp = RAND ? gRandom->Uniform(-0.5,0.5) : 0;
    return ( (double)tdc + tmp ) * m_gain + m_offset; }
  int    Tdc( double time ) const
  { return (int)( (time-m_offset) /m_gain  ); }
};

//______________________________________________________________________________
//Hodo ADC to DeltaE
class HodoAParam
{
public:
  inline HodoAParam( double pedestal, double gain )
    : m_pedestal(pedestal), m_gain(gain)
  {}
  ~HodoAParam( void )
  {}

private:
  HodoAParam( void );
  HodoAParam( const HodoAParam& );
  HodoAParam& operator =( const HodoAParam& );

private:
  double m_pedestal;
  double m_gain;

public:
  double Pedestal( void )  const { return m_pedestal; }
  double Gain( void )      const { return m_gain; }
  double DeltaE( int adc, bool RAND=false ) const
  { double tmp = RAND ? gRandom->Uniform(-0.5,0.5) : 0;
    return ( (double)adc - m_pedestal +tmp ) * (m_gain); }
  int    Adc( double de ) const
  { return (int)( m_gain * de + m_pedestal * ( 1. - de ) ); }
};

//______________________________________________________________________________
//HodoParam Main Class
class HodoParamMan
{
public:
  static HodoParamMan&      GetInstance( void );
  static const std::string& ClassName( void );
  ~HodoParamMan( void );

private:
  HodoParamMan( void );
  HodoParamMan( const HodoParamMan& );
  HodoParamMan& operator =( const HodoParamMan& );

private:
  enum eAorT { kAdc, kTdc, kAorT };
  typedef std::map<int, HodoTParam*> TContainer;
  typedef std::map<int, HodoAParam*> AContainer;
  typedef TContainer::const_iterator TIterator;
  typedef AContainer::const_iterator AIterator;
  bool        m_is_ready;
  std::vector<std::string> m_file_name;
  TContainer  m_TPContainer;
  AContainer  m_APContainer;

public:
  bool Initialize( void );
  bool Initialize( const std::string& file_name );
  bool Initialize( const std::string& file_name, const std::string& file_name2 );
  bool ReadFile( const std::string& file_name );
  bool IsReady( void ) const { return m_is_ready; }
  bool GetTime( int cid, int plid, int seg, int ud, int tdc, double &time ) const;
  bool GetDe( int cid, int plid, int seg, int ud, int adc, double &de ) const;
  bool GetTdc( int cid, int plid, int seg, int ud, double time, int &tdc ) const;
  bool GetAdc( int cid, int plid, int seg, int ud, double de, int &adc ) const;

private:
  HodoTParam *GetTmap( int cid, int plid, int seg, int ud ) const;
  HodoAParam *GetAmap( int cid, int plid, int seg, int ud ) const;
  void ClearACont( void );
  void ClearTCont( void );
};

//______________________________________________________________________________
inline HodoParamMan&
HodoParamMan::GetInstance( void )
{
  static HodoParamMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
HodoParamMan::ClassName( void )
{
  static std::string g_name("HodoParamMan");
  return g_name;
}

#endif
