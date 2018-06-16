/**
 *  file: EMCAnalyzer.hh
 *  date: 2017.04.10
 *
 */

#ifndef EMC_DATA_HH
#define EMC_DATA_HH

#include <ctime>
#include <string>
#include <vector>

//______________________________________________________________________________
class EMCAnalyzer
{
public:
  static EMCAnalyzer& GetInstance( void );
  ~EMCAnalyzer( void );

private:
  EMCAnalyzer( void );
  EMCAnalyzer( const EMCAnalyzer& );
  EMCAnalyzer& operator=( const EMCAnalyzer& );

private:
  unsigned int m_serial;
  double       m_raw_x;
  double       m_raw_y;
  double       m_x;
  double       m_y;
  int          m_state;
  unsigned int m_utime;
  unsigned int m_ltime;
  std::time_t  m_time;

public:
  void         ClearAll( void );
  bool         Decode( void );
  unsigned int Serial( void ) const { return m_serial; }
  double       RawX( void ) const { return m_raw_x; }
  double       RawY( void ) const { return m_raw_y; }
  double       X( void ) const { return m_x; }
  double       Y( void ) const { return m_y; }
  int          State( void ) const { return m_state; }
  unsigned int UTime( void ) const { return m_utime; }
  unsigned int LTime( void ) const { return m_ltime; }
  std::time_t  Time( void )  const { return m_time;  }
  void         Print( const std::string& arg="" ) const;

};

//______________________________________________________________________________
inline EMCAnalyzer&
EMCAnalyzer::GetInstance( void )
{
  static EMCAnalyzer g_instance;
  return g_instance;
}

#endif
