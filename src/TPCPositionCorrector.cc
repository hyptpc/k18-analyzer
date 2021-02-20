// -*- C++ -*-

#include "TPCPositionCorrector.hh"

#include <fstream>

#include <std_ostream.hh>

// #include "FieldElements.hh"
#include "FuncName.hh"
// #include "KuramaFieldMap.hh"


namespace
{
//const Double_t Delta = 0.1;
}

//_____________________________________________________________________________
TPCPositionCorrector::TPCPositionCorrector( void )
  : m_is_ready( false ),
    m_file_name(),
    m_correction_map(),
    m_n_x(),
    m_n_y(),
    m_n_z(),
    m_x0(),
    m_y0(),
    m_z0(),
    m_dx(),
    m_dy(),
    m_dz()
{
}

//_____________________________________________________________________________
TPCPositionCorrector::~TPCPositionCorrector( void )
{
}

//_____________________________________________________________________________
TVector3
TPCPositionCorrector::Correct( const TVector3& pos ) const
{
  return pos + GetCorrectionVector(pos);
}

//_____________________________________________________________________________
TVector3
TPCPositionCorrector::GetCorrectionVector( const TVector3& pos ) const
{
  Int_t ix1, ix2, iy1, iy2, iz1, iz2;
  ix1 = Int_t( ( pos.X() - m_x0 )/m_dx );
  iy1 = Int_t( ( pos.Y() - m_y0 )/m_dy );
  iz1 = Int_t( ( pos.Z() - m_z0 )/m_dz );

  Double_t wx1, wx2, wy1, wy2, wz1, wz2;
  if( ix1<0 ) { ix1=ix2=0; wx1=1.; wx2=0.; }
  else if( ix1>=m_n_x-1 ) { ix1=ix2=m_n_x-1; wx1=1.; wx2=0.; }
  else { ix2=ix1+1; wx1=(m_x0+m_dx*ix2-pos.X())/m_dx; wx2=1.-wx1; }

  if( iy1<0 ) { iy1=iy2=0; wy1=1.; wy2=0.; }
  else if( iy1>=m_n_y-1 ) { iy1=iy2=m_n_y-1; wy1=1.; wy2=0.; }
  else { iy2=iy1+1; wy1=(m_y0+m_dy*iy2-pos.Y())/m_dy; wy2=1.-wy1; }

  if( iz1<0 ) { iz1=iz2=0; wz1=1.; wz2=0.; }
  else if( iz1>=m_n_z-1 ) { iz1=iz2=m_n_z-1; wz1=1.; wz2=0.; }
  else { iz2=iz1+1; wz1=(m_z0+m_dz*iz2-pos.Z())/m_dz; wz2=1.-wz1; }

  Double_t cvx1 = wx1*wy1*m_correction_map[ix1][iy1][iz1].X()
    + wx1*wy2*m_correction_map[ix1][iy2][iz1].X()
    + wx2*wy1*m_correction_map[ix2][iy1][iz1].X()
    + wx2*wy2*m_correction_map[ix2][iy2][iz1].X();
  Double_t cvx2 = wx1*wy1*m_correction_map[ix1][iy1][iz2].X()
    + wx1*wy2*m_correction_map[ix1][iy2][iz2].X()
    + wx2*wy1*m_correction_map[ix2][iy1][iz2].X()
    + wx2*wy2*m_correction_map[ix2][iy2][iz2].X();
  Double_t cvx  = wz1*cvx1 + wz2*cvx2;

  Double_t cvy1 = wx1*wy1*m_correction_map[ix1][iy1][iz1].Y()
    + wx1*wy2*m_correction_map[ix1][iy2][iz1].Y()
    + wx2*wy1*m_correction_map[ix2][iy1][iz1].Y()
    + wx2*wy2*m_correction_map[ix2][iy2][iz1].Y();
  Double_t cvy2 = wx1*wy1*m_correction_map[ix1][iy1][iz2].Y()
    + wx1*wy2*m_correction_map[ix1][iy2][iz2].Y()
    + wx2*wy1*m_correction_map[ix2][iy1][iz2].Y()
    + wx2*wy2*m_correction_map[ix2][iy2][iz2].Y();
  Double_t cvy  = wz1*cvy1 + wz2*cvy2;

  Double_t cvz1 = wx1*wy1*m_correction_map[ix1][iy1][iz1].Z()
    + wx1*wy2*m_correction_map[ix1][iy2][iz1].Z()
    + wx2*wy1*m_correction_map[ix2][iy1][iz1].Z()
    + wx2*wy2*m_correction_map[ix2][iy2][iz1].Z();
  Double_t cvz2 = wx1*wy1*m_correction_map[ix1][iy1][iz2].Z()
    + wx1*wy2*m_correction_map[ix1][iy2][iz2].Z()
    + wx2*wy1*m_correction_map[ix2][iy1][iz2].Z()
    + wx2*wy2*m_correction_map[ix2][iy2][iz2].Z();
  Double_t cvz  = wz1*cvz1 + wz2*cvz2;

  return TVector3( cvx, cvy, cvz );
}

//_____________________________________________________________________________
Bool_t
TPCPositionCorrector::Initialize( void )
{
  if( m_is_ready ){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  std::ifstream ifs( m_file_name );
  if( !ifs.is_open() ){
    hddaq::cerr << FUNC_NAME << " file open fail : "
                << m_file_name << std::endl;
    return false;
  }

  if( !( ifs >> m_n_x >> m_n_y >> m_n_z >>
         m_x0 >> m_y0 >> m_z0 >> m_dx >> m_dy >> m_dz ) ){
    hddaq::cerr << FUNC_NAME << " invalid format" << std::endl;
    return false;
  }

  if( m_n_x <=0 || m_n_y <= 0 || m_n_z <= 0 ){
    hddaq::cerr << FUNC_NAME << " n_x, n_y, n_z must be positive." << std::endl;
    return false;
  }

  m_correction_map.resize( m_n_x );
  for( auto& x : m_correction_map ){
    x.resize( m_n_y );
    for( auto& y : x ){
      y.resize( m_n_z );
    }
  }

  hddaq::cout << " reading correction map " << std::flush;

  Int_t line = 0;
  Double_t x, y, z, vx, vy, vz;
  while( ifs.good() ){
    if( line++%1000==0 ) hddaq::cout << "." << std::flush;
    ifs >> x >> y >> z >> vx >> vy >> vz;
    Int_t ix = Int_t( ( x - m_x0 + m_dx )/m_dx );
    Int_t iy = Int_t( ( y - m_y0 + m_dy )/m_dy );
    Int_t iz = Int_t( ( z - m_z0 + m_dz )/m_dz );
    if( ix >= 0 && ix < m_n_x &&
        iy >= 0 && iy < m_n_y &&
        iz >= 0 && iz < m_n_z ){
      m_correction_map[ix][iy][iz].SetX( vx );
      m_correction_map[ix][iy][iz].SetY( vy );
      m_correction_map[ix][iy][iz].SetZ( vz );
    }
  }

  hddaq::cout << " done" << std::endl;

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
TPCPositionCorrector::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
// ThreeVector
// TPCPositionCorrector::GetField( const ThreeVector& position ) const
// {
//   ThreeVector field( 0., 0., 0. );
//   if( m_kurama_map ){
//     Double_t p[3], b[3];
//     p[0] = position.x()*0.1;
//     p[1] = position.y()*0.1;
//     p[2] = position.z()*0.1;
//     if( m_kurama_map->GetFieldValue( p, b ) ){
//       field.SetX(b[0]);
//       field.SetY(b[1]);
//       field.SetZ(b[2]);
//     }
//   }

// #if 1
//   FEIterator itr, itr_end = m_element_list.end();
//   for( itr=m_element_list.begin(); itr!=itr_end; ++itr ){
//     if( (*itr)->ExistField( position ) )
//       field += (*itr)->GetField( position );
//   }
// #endif

//   return field;
// }

//_____________________________________________________________________________
// ThreeVector
// TPCPositionCorrector::GetdBdX( const ThreeVector& position ) const
// {
//   ThreeVector p1 = position + ThreeVector( Delta, 0., 0. );
//   ThreeVector p2 = position - ThreeVector( Delta, 0., 0. );
//   ThreeVector B1 = GetField( p1 );
//   ThreeVector B2 = GetField( p2 );
//   return 0.5/Delta*(B1-B2);
// }

//_____________________________________________________________________________
// ThreeVector
// TPCPositionCorrector::GetdBdY( const ThreeVector& position ) const
// {
//   ThreeVector p1 = position + ThreeVector( 0., Delta, 0. );
//   ThreeVector p2 = position - ThreeVector( 0., Delta, 0. );
//   ThreeVector B1 = GetField( p1 );
//   ThreeVector B2 = GetField( p2 );
//   return 0.5/Delta*(B1-B2);
// }

//_____________________________________________________________________________
// ThreeVector
// TPCPositionCorrector::GetdBdZ( const ThreeVector& position ) const
// {
//   ThreeVector p1 = position + ThreeVector( 0., 0., Delta );
//   ThreeVector p2 = position - ThreeVector( 0., 0., Delta );
//   ThreeVector B1 = GetField( p1 );
//   ThreeVector B2 = GetField( p2 );
//   return 0.5/Delta*(B1-B2);
// }

//_____________________________________________________________________________
// void
// TPCPositionCorrector::ClearElementsList( void )
// {
//   m_element_list.clear();
// }

//_____________________________________________________________________________
// void
// TPCPositionCorrector::AddElement( FieldElements *element )
// {
//   m_element_list.push_back( element );
// }

//_____________________________________________________________________________
// Double_t
// TPCPositionCorrector::StepSize( const ThreeVector &position,
// 		    Double_t def_step_size, Double_t min_step_size ) const
// {
//   Double_t d = std::abs( def_step_size );
//   Double_t s = def_step_size/d;
//   min_step_size = std::abs( min_step_size );

//   Bool_t flag = true;
//   FEIterator itr, itr_end = m_element_list.end();
//   while ( flag && d>min_step_size ){
//     for( itr=m_element_list.begin(); itr!=itr_end; ++itr ){
//       if( (*itr)->CheckRegion( position, d ) != FieldElements::FERSurface() ){
//         flag=false;
//         break;
//       }
//     }
//     d *= 0.5;
//   }
//   if( flag )
//     return s*min_step_size;
//   else
//     return s*d;
// }
