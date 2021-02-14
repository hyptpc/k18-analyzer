// -*- C++ -*-

#ifndef CONF_MAN_HH
#define CONF_MAN_HH

#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include <TString.h>

#include <std_ostream.hh>

class VEvent;

//_____________________________________________________________________________
class ConfMan
{
public:
  static const TString& ClassName( void );
  static ConfMan&       GetInstance( void );
  ~ConfMan( void );

private:
  ConfMan( void );
  ConfMan( const ConfMan& );
  ConfMan& operator=( const ConfMan& );

private:
  typedef std::map<TString, TString>  StrList;
  typedef std::map<TString, Double_t> DoubleList;
  typedef std::map<TString, Int_t>    IntList;
  typedef std::map<TString, Bool_t>   BoolList;
  typedef StrList::const_iterator     StrIterator;
  typedef DoubleList::const_iterator  DoubleIterator;
  typedef IntList::const_iterator     IntIterator;
  typedef BoolList::const_iterator    BoolIterator;
  Bool_t     m_is_ready;
  StrList    m_file;
  StrList    m_string;
  DoubleList m_double;
  IntList    m_int;
  BoolList   m_bool;

public:
  VEvent* EventAllocator( void );
  Bool_t  Initialize( void );
  Bool_t  Initialize( const TString& file_name );
  Bool_t  InitializeHistograms( void );
  Bool_t  InitializeParameterFiles( void );
  Bool_t  InitializeUnpacker( void );
  Bool_t  IsReady( void ) const { return m_is_ready; }
  Bool_t  Finalize( void );
  Bool_t  FinalizeProcess( void );
  // templates
  template <typename T>
  static const T& Get( const TString& key ) { return T(); }
  template <typename T>
  Bool_t    InitializeParameter( void );
  template <typename T>
  Bool_t    InitializeParameter( const TString& key );
  template <typename T>
  Bool_t    InitializeParameter( const TString& key1,
                                 const TString& key2 );

private:
  TString FilePath( const TString& src ) const;
  Bool_t  ShowResult( Bool_t s, const TString& name ) const;
};

//_____________________________________________________________________________
inline const TString&
ConfMan::ClassName( void )
{
  static TString s_name( "ConfMan" );
  return s_name;
}

//_____________________________________________________________________________
inline ConfMan&
ConfMan::GetInstance( void )
{
  static ConfMan s_instance;
  return s_instance;
}

//_____________________________________________________________________________
template <>
inline const TString&
ConfMan::Get<TString>( const TString& key )
{
  return GetInstance().m_string[key];
}

//_____________________________________________________________________________
template <>
inline const Double_t&
ConfMan::Get<Double_t>( const TString& key )
{
  return GetInstance().m_double[key];
}

//_____________________________________________________________________________
template <>
inline const Int_t&
ConfMan::Get<Int_t>( const TString& key )
{
  return GetInstance().m_int[key];
}

//_____________________________________________________________________________
template <>
inline const Bool_t&
ConfMan::Get<Bool_t>( const TString& key )
{
  return GetInstance().m_bool[key];
}

//_____________________________________________________________________________
inline Bool_t
ConfMan::ShowResult( Bool_t s, const TString& name ) const
{
  if( s )
    hddaq::cout << std::setw(20) << std::left
		<< " ["+name+"]"
		<< "-> Initialized" << std::endl;
  else
    hddaq::cout << std::setw(20) << std::left
		<< " ["+name+"]"
		<< "-> Failed" << std::endl;
  return s;
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
ConfMan::InitializeParameter( void )
{
  return
    ShowResult( T::GetInstance().Initialize(),
		T::ClassName() );
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
ConfMan::InitializeParameter( const TString& key )
{
  return
    ShowResult( T::GetInstance().Initialize(m_file[key]),
		T::ClassName() );
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
ConfMan::InitializeParameter( const TString& key1,
			      const TString& key2 )
{
  return
    ShowResult( T::GetInstance().Initialize(m_file[key1],
					    m_file[key2]),
		T::ClassName() );
}

#endif
