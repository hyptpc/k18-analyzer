/**
 *  file: HodoParamMan.hh
 *  date: 2017.04.10
 *
 */

#ifndef HODO_PARAM_MAN_HH
#define HODO_PARAM_MAN_HH

#include <string>
#include <map>

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
  double Time( int tdc ) const
  { return ( (double)tdc - m_offset ) * m_gain; }
  int    Tdc( double time ) const
  { return (int)( time/m_gain + m_offset ); }
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
  double DeltaE( int adc ) const
  { return ( (double)adc - m_pedestal ) / (m_gain - m_pedestal ); }
  int    Adc( double de ) const
  { return (int)( m_gain * de + m_pedestal * ( 1. - de ) ); }
};

//______________________________________________________________________________
//correction parameter for CFT fiber position
class HodoFParam
{
public:
  inline HodoFParam( double par0, double par1, double par2, double par3, double par4, double par5 )
    : Par0(par0), Par1(par1), Par2(par2), Par3(par3), Par4(par4), Par5(par5)
  {}
  ~HodoFParam() {}
private:
  HodoFParam();
  HodoFParam( const HodoFParam & );
  HodoFParam & operator = ( const HodoFParam & );
private:
  double Par0, Par1, Par2, Par3, Par4, Par5;
public:
  double par0( void ) const { return Par0; }
  double par1( void ) const { return Par1; }
  double par2( void ) const { return Par2; }
  double par3( void ) const { return Par3; }
  double par4( void ) const { return Par4; }
  double par5( void ) const { return Par5; }
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
  typedef std::map<int, HodoFParam*> FContainer;
  typedef TContainer::const_iterator TIterator;
  typedef AContainer::const_iterator AIterator;
  typedef FContainer::const_iterator FIterator;
  bool        m_is_ready;
  std::string m_file_name;
  TContainer  m_TPContainer;
  AContainer  m_APContainer;
  FContainer  m_FPContainer;

public:
  bool Initialize( void );
  bool Initialize( const std::string& file_name );
  bool IsReady( void ) const { return m_is_ready; }
  bool GetTime( int cid, int plid, int seg, int ud, int tdc, double &time ) const;
  bool GetDe( int cid, int plid, int seg, int ud, int adc, double &de ) const;
  bool GetTdc( int cid, int plid, int seg, int ud, double time, int &tdc ) const;
  bool GetAdc( int cid, int plid, int seg, int ud, double de, int &adc ) const;
  void SetFileName( const std::string& file_name ) { m_file_name = file_name; }

  double GetP0( int cid, int plid, int seg, int ud) const;
  double GetP1( int cid, int plid, int seg, int ud) const;
  double GetPar( int cid, int plid, int seg, int ud ,int i) const; // i=0~5

  double  GetOffset(int cid, int plid, int seg, int ud) const;
  double  GetGain(int cid, int plid, int seg, int ud) const;

private:
  HodoTParam *GetTmap( int cid, int plid, int seg, int ud ) const;
  HodoAParam *GetAmap( int cid, int plid, int seg, int ud ) const;
  HodoFParam *GetFmap( int cid, int plid, int seg, int ud ) const;
  void ClearACont( void );
  void ClearTCont( void );
  void ClearFCont( void );
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
