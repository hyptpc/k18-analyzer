// -*- C++ -*-

#ifndef RESLMAN_HH
#define RESLMAN_HH

#include <map>
#include <iostream>
#include <TString.h>

class ReslMap
{
 public:
  ReslMap(){};
  ~ReslMap(){};

 private:
  double TResl, EResl;

 public:
  double tresl() const { return TResl; }
  double eresl() const { return EResl; }
  void SetTResl( const double &t ) { TResl = t; }
  void SetEResl( const double &e ) { EResl = e; }
};

class ReslMapMan
{
 public:
  static ReslMapMan& GetInstance(void);
  static const std::string& ClassName( void );
  ~ReslMapMan(void) {}
  void SetFileName( const TString & filename ) { FileName = filename; }
  bool Initialize();
  bool Initialize( const char *file_name );
  bool Initialize( const std::string& file_name );

private:
  ReslMapMan(void);

private:
  TString FileName;
  typedef std::map <unsigned int, ReslMap> ReslContainer;
  ReslContainer reslContainer;
  bool m_is_ready;

 public:
  bool IsReady( void ) const { return m_is_ready; }
  TString GetFileName() { return FileName; }
  void GetResolution( const int &cid, const int &layer, const int &wire, double &tresl, double &eresl ) const;
  bool SetParam( const int &cid, const int &layer, const int &wire, const double &tresl, const double &eresl );
  bool GetParam( const int &cid, const int &layer, const int &wire, double &tresl, double &eresl );

  // void PrintMap( const int &Cid = -1, std::ostream &p_out = std::cout ) const;
  void PrintMap( std::ostream &p_out = std::cout ) const;
};

inline ReslMapMan&
ReslMapMan::GetInstance( void )
{
  static ReslMapMan g_instance;
  return g_instance;
}
inline const std::string&
ReslMapMan::ClassName( void )
{
  static std::string g_name("ReslMan");
  return g_name;
}
#endif
