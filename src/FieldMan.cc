// -*- C++ -*-

#include "FieldMan.hh"

#include <cmath>

#include <std_ostream.hh>

#include "FieldElements.hh"
#include "FuncName.hh"
#include "KuramaFieldMap.hh"


namespace
{
const Double_t Delta = 0.1;
}

//_____________________________________________________________________________
FieldMan::FieldMan( void )
  : m_is_ready( false ),
    m_kurama_map( nullptr )
{
}

//_____________________________________________________________________________
FieldMan::~FieldMan( void )
{
  if( m_kurama_map )
    delete m_kurama_map;
}

//_____________________________________________________________________________
Bool_t
FieldMan::Initialize( void )
{
  if( m_is_ready ){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  if( m_kurama_map )
    delete m_kurama_map;

  m_kurama_map = new KuramaFieldMap( m_file_name );
  if( m_kurama_map )
    m_is_ready = m_kurama_map->Initialize();
  else
    m_is_ready = false;

  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
FieldMan::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
ThreeVector
FieldMan::GetField( const ThreeVector& position ) const
{
  ThreeVector field( 0., 0., 0. );
  if( m_kurama_map ){
    Double_t p[3], b[3];
    p[0] = position.x()*0.1;
    p[1] = position.y()*0.1;
    p[2] = position.z()*0.1;
    if( m_kurama_map->GetFieldValue( p, b ) ){
      field.SetX(b[0]);
      field.SetY(b[1]);
      field.SetZ(b[2]);
    }
  }

#if 1
  FEIterator itr, itr_end = m_element_list.end();
  for( itr=m_element_list.begin(); itr!=itr_end; ++itr ){
    if( (*itr)->ExistField( position ) )
      field += (*itr)->GetField( position );
  }
#endif

  return field;
}

//_____________________________________________________________________________
ThreeVector
FieldMan::GetdBdX( const ThreeVector& position ) const
{
  ThreeVector p1 = position + ThreeVector( Delta, 0., 0. );
  ThreeVector p2 = position - ThreeVector( Delta, 0., 0. );
  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );
  return 0.5/Delta*(B1-B2);
}

//_____________________________________________________________________________
ThreeVector
FieldMan::GetdBdY( const ThreeVector& position ) const
{
  ThreeVector p1 = position + ThreeVector( 0., Delta, 0. );
  ThreeVector p2 = position - ThreeVector( 0., Delta, 0. );
  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );
  return 0.5/Delta*(B1-B2);
}

//_____________________________________________________________________________
ThreeVector
FieldMan::GetdBdZ( const ThreeVector& position ) const
{
  ThreeVector p1 = position + ThreeVector( 0., 0., Delta );
  ThreeVector p2 = position - ThreeVector( 0., 0., Delta );
  ThreeVector B1 = GetField( p1 );
  ThreeVector B2 = GetField( p2 );
  return 0.5/Delta*(B1-B2);
}

//_____________________________________________________________________________
void
FieldMan::ClearElementsList( void )
{
  m_element_list.clear();
}

//_____________________________________________________________________________
void
FieldMan::AddElement( FieldElements *element )
{
  m_element_list.push_back( element );
}

//_____________________________________________________________________________
Double_t
FieldMan::StepSize( const ThreeVector &position,
		    Double_t def_step_size, Double_t min_step_size ) const
{
  Double_t d = std::abs( def_step_size );
  Double_t s = def_step_size/d;
  min_step_size = std::abs( min_step_size );

  Bool_t flag = true;
  FEIterator itr, itr_end = m_element_list.end();
  while ( flag && d>min_step_size ){
    for( itr=m_element_list.begin(); itr!=itr_end; ++itr ){
      if( (*itr)->CheckRegion( position, d ) != FieldElements::FERSurface() ){
        flag=false;
        break;
      }
    }
    d *= 0.5;
  }
  if( flag )
    return s*min_step_size;
  else
    return s*d;
}
