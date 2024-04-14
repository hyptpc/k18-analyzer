#ifndef LOGGER_DATA_HH
#define LOGGER_DATA_HH

#include <ctime>
#include <string>
#include <vector>

//______________________________________________________________________________
class LoggerAnalyzer
{
public:
  static LoggerAnalyzer& GetInstance( void );
  ~LoggerAnalyzer( void );

private:
  LoggerAnalyzer( void );
  LoggerAnalyzer( const LoggerAnalyzer& );
  LoggerAnalyzer& operator=( const LoggerAnalyzer& );

private:
  int m_serial;
  int m_raw_x;
  int m_raw_y;
  int m_raw_data[10];  
  int          m_segment;
  double       m_x;
  double       m_y;
  double       m_data[10];

public:
  void         ClearAll( void );
  bool         Decode( void );
  int Serial( void ) const { return m_serial; }
  int RawX( void )   const { return m_raw_x; }
  int RawY( void )   const { return m_raw_y; }
  int RawData( const int &i )   const { return m_raw_data[i]; }
  double       X( void ) const { return m_x; }
  double       Y( void ) const { return m_y; }
  double       Data( const int &i )   const { return m_data[i]; }
  double       GetDetectorTemperature()   const { return Data(0); }
  double       GetRoomTemperature()       const { return Data(1); }
  double       GetMagnetCurrent()         const { return Data(2); }
  int          GetSegment()               const { return m_segment; }
  void         Print( const std::string& arg="" ) const;
};

//______________________________________________________________________________
inline LoggerAnalyzer&
LoggerAnalyzer::GetInstance( void )
{
  static LoggerAnalyzer g_instance;
  return g_instance;
}

#endif
