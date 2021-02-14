// -*- C++ -*-

#ifndef FIELD_MAN_HH
#define FIELD_MAN_HH

#include "ThreeVector.hh"
#include <string>
#include <vector>

class KuramaFieldMap;
class FieldElements;

typedef std::vector<FieldElements*> FEContainer;
typedef FEContainer::const_iterator FEIterator;

//_____________________________________________________________________________
class FieldMan
{
public:
  static const TString& ClassName( void );
  static FieldMan&      GetInstance( void );
  ~FieldMan( void );

private:
  FieldMan( void );
  FieldMan( const FieldMan& );
  FieldMan & operator =( const FieldMan& );

private:
  Bool_t          m_is_ready;
  TString         m_file_name;
  KuramaFieldMap* m_kurama_map;
  FEContainer     m_element_list;

public:
  Bool_t      Initialize( void );
  Bool_t      Initialize( const TString& file_name );
  Bool_t      IsReady( void ) const { return m_is_ready; }
  ThreeVector GetField( const ThreeVector& position ) const;
  ThreeVector GetdBdX( const ThreeVector& position ) const;
  ThreeVector GetdBdY( const ThreeVector& position ) const;
  ThreeVector GetdBdZ( const ThreeVector& position ) const;
  void        ClearElementsList( void );
  void        AddElement( FieldElements* element );
  void        SetFileName( const TString& file_name )
    { m_file_name = file_name; }
  Double_t    StepSize( const ThreeVector& position,
                        Double_t default_step_size,
                        Double_t min_step_size ) const;
};

//_____________________________________________________________________________
inline const TString&
FieldMan::ClassName( void )
{
  static TString s_name( "FieldMan" );
  return s_name;
}

//_____________________________________________________________________________
inline FieldMan&
FieldMan::GetInstance( void )
{
  static FieldMan s_instance;
  return s_instance;
}

#endif
