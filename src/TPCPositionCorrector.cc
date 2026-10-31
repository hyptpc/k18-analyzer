// -*- C++ -*-

#include "TPCPositionCorrector.hh"

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "FuncName.hh"
#include "TPCPadHelper.hh"

#include <std_ostream.hh>

//_____________________________________________________________________________
TPCPositionCorrector::TPCPositionCorrector( void )
  : m_map(NumOfLayersTPC),
    m_is_ready( false ),
    m_file_name(),
    m_n_y(),
    m_y0(),
    m_dy()
{
}

//_____________________________________________________________________________
TPCPositionCorrector::~TPCPositionCorrector( void )
{
}

//_____________________________________________________________________________
TVector3
TPCPositionCorrector::Correct( const TVector3& pos, Int_t layer, Int_t row  ) const
{
  return pos + GetCorrectionVector(pos, layer, row);
}

//_____________________________________________________________________________
TVector3
TPCPositionCorrector::GetCorrectionVector( const TVector3& pos, Int_t layer, Int_t row ) const
{
  // Clamp to endpoint nodes outside [y0, y_last] before any index math.
  // Comparing real y avoids Int_t truncation that extrapolated in (y0-dy, y0).
  const Double_t y = pos.Y();
  const Double_t y_last = m_y0 + m_dy * static_cast<Double_t>(m_n_y - 1);

  Int_t iy1 = 0;
  Int_t iy2 = 0;
  Double_t wy1 = 1.;
  Double_t wy2 = 0.;

  if( y <= m_y0 || m_n_y == 1 ){
    iy1 = iy2 = 0;
    wy1 = 1.;
    wy2 = 0.;
  } else if( y >= y_last ){
    iy1 = iy2 = m_n_y - 1;
    wy1 = 1.;
    wy2 = 0.;
  } else {
    iy1 = static_cast<Int_t>((y - m_y0) / m_dy);
    if( iy1 < 0 ){
      iy1 = 0;
    }
    if( iy1 >= m_n_y - 1 ){
      iy1 = m_n_y - 2;
    }
    iy2 = iy1 + 1;
    wy1 = (m_y0 + m_dy * static_cast<Double_t>(iy2) - y) / m_dy;
    wy2 = 1. - wy1;
  }

  return wy1 * m_map[layer][row][iy1] + wy2 * m_map[layer][row][iy2];
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

  // Line-oriented read. Skip blank lines and '#' comment lines
  // (same convention as ConfMan / TPCParamMan: line[0] == '#').
  auto skip_blank_or_comment = [](const std::string& line) -> Bool_t {
    if( line.empty() || line[0] == '#' ){
      return true;
    }
    std::istringstream blank_iss(line);
    std::string tok;
    return !(blank_iss >> tok);
  };

  hddaq::cout << " reading correction map " << std::flush;

  constexpr Double_t k_y_tol = 1.0e-6;
  Int_t line_no = 0;
  std::string line;
  Bool_t header_ok = false;
  while( std::getline(ifs, line) ){
    ++line_no;
    if( skip_blank_or_comment(line) ){
      continue;
    }
    std::istringstream iss(line);
    if( !(iss >> m_n_y >> m_y0 >> m_dy) ){
      hddaq::cerr << FUNC_NAME << " parse error in " << m_file_name
                  << " at line " << line_no
                  << " (expect n_y y0 dy): " << line << std::endl;
      return false;
    }
    std::string extra;
    if( iss >> extra ){
      hddaq::cerr << FUNC_NAME << " parse error in " << m_file_name
                  << " at line " << line_no
                  << " (extra token '" << extra << "' in header): "
                  << line << std::endl;
      return false;
    }
    header_ok = true;
    break;
  }
  if( !header_ok ){
    hddaq::cerr << FUNC_NAME << " invalid format in " << m_file_name
                << " (missing header n_y y0 dy)." << std::endl;
    return false;
  }

  if( m_n_y <= 0 ){
    hddaq::cerr << FUNC_NAME << " n_y must be positive." << std::endl;
    return false;
  }
  if( !(std::isfinite(m_y0) && std::isfinite(m_dy)) || m_dy <= 0. ){
    hddaq::cerr << FUNC_NAME << " y0/dy must be finite and dy > 0." << std::endl;
    return false;
  }

  for( Int_t l = 0; l < NumOfLayersTPC; l++ ){
    Int_t nrow = tpc::padParameter[l][1];
    m_map[l].resize( nrow );
    for( auto& x : m_map[l] ){
      x.resize( m_n_y );
    }
  }

  // Track which (layer,row,iy) cells were filled for a gap/dup warning.
  std::vector<std::vector<std::vector<char>>> filled(NumOfLayersTPC);
  for( Int_t l = 0; l < NumOfLayersTPC; l++ ){
    Int_t nrow = tpc::padParameter[l][1];
    filled[l].assign(nrow, std::vector<char>(m_n_y, 0));
  }

  Int_t n_data = 0;
  Int_t n_dup = 0;
  Int_t n_bad = 0;
  Int_t n_filled = 0;
  while( std::getline(ifs, line) ){
    ++line_no;
    if( skip_blank_or_comment(line) ){
      continue;
    }

    std::istringstream iss(line);
    Int_t layer = 0;
    Int_t row = 0;
    Double_t y = 0.;
    Double_t vx = 0.;
    Double_t vz = 0.;
    if( !(iss >> layer >> row >> y >> vx >> vz) ){
      hddaq::cerr << FUNC_NAME << " parse error in " << m_file_name
                  << " at line " << line_no
                  << " (expect layer row y vx vz): " << line << std::endl;
      return false;
    }
    std::string extra;
    if( iss >> extra ){
      hddaq::cerr << FUNC_NAME << " parse error in " << m_file_name
                  << " at line " << line_no
                  << " (extra token '" << extra << "'): " << line << std::endl;
      return false;
    }

    if( n_data++ % 1000 == 0 ) hddaq::cout << "." << std::flush;

    if( layer < 0 || layer >= NumOfLayersTPC ){
      ++n_bad;
      continue;
    }
    const Int_t nrow = tpc::padParameter[layer][1];
    if( row < 0 || row >= nrow ){
      ++n_bad;
      continue;
    }
    if( !(std::isfinite(y) && std::isfinite(vx) && std::isfinite(vz)) ){
      ++n_bad;
      continue;
    }

    const Double_t iy_f = (y - m_y0) / m_dy;
    const Int_t iy = static_cast<Int_t>(std::lround(iy_f));
    if( iy < 0 || iy >= m_n_y ){
      ++n_bad;
      continue;
    }
    const Double_t y_grid = m_y0 + m_dy * static_cast<Double_t>(iy);
    if( std::fabs(y - y_grid) > k_y_tol * std::max(1., std::fabs(m_dy)) ){
      ++n_bad;
      continue;
    }

    if( filled[layer][row][iy] ){
      ++n_dup;
    }
    filled[layer][row][iy] = 1;
    m_map[layer][row][iy].SetX(vx);
    m_map[layer][row][iy].SetY(0.);
    m_map[layer][row][iy].SetZ(vz);
    ++n_filled;
  }

  if( !ifs.eof() ){
    hddaq::cerr << FUNC_NAME << " read error in " << m_file_name
                << " near line " << line_no << std::endl;
    return false;
  }

  if( n_filled == 0 ){
    hddaq::cerr << FUNC_NAME << " no valid map cells in " << m_file_name
                << " (empty body is not a zero map)." << std::endl;
    return false;
  }

  Int_t n_miss = 0;
  for( Int_t l = 0; l < NumOfLayersTPC; l++ ){
    const Int_t nrow = tpc::padParameter[l][1];
    for( Int_t r = 0; r < nrow; r++ ){
      for( Int_t iy = 0; iy < m_n_y; iy++ ){
        if( !filled[l][r][iy] ){
          ++n_miss;
        }
      }
    }
  }

  hddaq::cout << " done" << std::endl;
  if( n_bad > 0 ){
    hddaq::cerr << FUNC_NAME << " skipped " << n_bad
                << " invalid map rows (layer/row/y/NaN)." << std::endl;
  }
  if( n_dup > 0 ){
    hddaq::cerr << FUNC_NAME << " warning: " << n_dup
                << " duplicate (layer,row,y) cells overwritten." << std::endl;
  }
  if( n_miss > 0 ){
    hddaq::cerr << FUNC_NAME << " warning: " << n_miss
                << " missing (layer,row,y) cells left as zero." << std::endl;
  }

  m_is_ready = true;
  return m_is_ready;
}

/*
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
*/
//_____________________________________________________________________________
Bool_t
TPCPositionCorrector::Initialize( const TString& file_name )
{
  m_file_name = file_name;
  return Initialize();
}
