/**
 *  file: ScalerAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef SCALER_ANALYZER_HH
#define SCALER_ANALYZER_HH

#include "DetectorID.hh"

#include <iomanip>
#include <iostream>
#include <vector>

typedef unsigned long long Scaler;

//______________________________________________________________________________
class ScalerAnalyzer
{
public:
  static ScalerAnalyzer& GetInstance( void );
  ~ScalerAnalyzer( void );

private:
  ScalerAnalyzer( void );
  ScalerAnalyzer( const ScalerAnalyzer& );
  ScalerAnalyzer& operator=( const ScalerAnalyzer& );

public:
  struct ScalerInfo
  {
    std::string name;
    int         module_id;
    int         channel;
    bool        flag_disp;
    Scaler      data; // integral value in a spill
    Scaler      curr; // current value
    Scaler      prev; // previous value

    ScalerInfo( void ) {}
    ScalerInfo( std::string n, int m, int c, bool f=true )
      : name(n), module_id(m), channel(c), flag_disp(f),
	data(0), curr(0), prev(0)
    {
    }
    inline void Print( void ) const
    {
      std::cout << std::setw(8)  << name
		<< std::setw(2)  << module_id
		<< std::setw(2)  << channel
		<< std::setw(2)  << flag_disp
		<< std::setw(10) << data << std::endl;
    }
  };

private:
  static const std::size_t MaxRow = 32;
  enum eDisp   { kLeft, kRight, MaxColumn };
  // for Phase1
  enum eScaler { kScaler1, kScaler2, kScaler3, nScaler };

  typedef std::vector< std::vector<ScalerInfo> > ScalerList;

  std::ostream& m_ost;
  ScalerList    m_info;
  bool          m_spill_increment;
  bool          m_separate_comma;

public:
  void        Clear( void );
  bool        Decode( void );
  bool        FlagDisp( int i, int j ) const { return m_info[i][j].flag_disp; }
  Scaler      Get( int i, int j ) const { return m_info[i][j].data; }
  Scaler      Get( const std::string& name ) const;
  std::string GetName( int i, int j ) const { return m_info[i][j].name; }
  ScalerInfo  GetScalerInfo( int i, int j ) const { return m_info[i][j]; }
  bool        Initialize( void );
  int         ModuleId( int i, int j ) const { return m_info[i][j].module_id; }
  void        Print( const std::string& arg="" ) const;
  bool        SpillIncrement( void ) const { return m_spill_increment; }
  std::string SeparateComma( Scaler number ) const;
  void        SetComma( bool flag=true ) { m_separate_comma = flag; }

};

//______________________________________________________________________________
inline ScalerAnalyzer&
ScalerAnalyzer::GetInstance( void )
{
  static ScalerAnalyzer g_instance;
  return g_instance;
}

#endif
