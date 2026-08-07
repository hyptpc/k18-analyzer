// -*- C++ -*-

#include "D5TransMatrix.hh"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include "FuncName.hh"

//_____________________________________________________________________________
D5TransMatrix&
D5TransMatrix::GetInstance()
{
  static D5TransMatrix s_instance;
  return s_instance;
}

//_____________________________________________________________________________
D5TransMatrix::D5TransMatrix()
  : m_is_ready(false),
    m_z_ref_ready(false),
    m_file_name(""),
    m_central_momentum(0.0),
    m_d5_z_in(RefZIn),
    m_d5_z_out(RefZOut),
    m_matrix_1st(6, 6)  // R: 6x6
{
  // T[i]: 6x6 quadratic block; i=0..4 = x,u,y,v,delta (D5Matrix: 21..25)
  for (Int_t i=0; i<5; ++i)
    m_matrix_2nd[i].ResizeTo(6, 6);
}

//_____________________________________________________________________________
D5TransMatrix::~D5TransMatrix()
{
}

//_____________________________________________________________________________
void
D5TransMatrix::Clear()
{
  m_is_ready = false;
  m_z_ref_ready = false;
  m_central_momentum = 0.0;
  ResetRefPlanesToDefault();
  m_matrix_1st.Zero();
  for(Int_t i=0; i<5; ++i)
    m_matrix_2nd[i].Zero();
}

//_____________________________________________________________________________
Bool_t
D5TransMatrix::Initialize()
{
  if (m_is_ready) return true;
  if (m_file_name.IsNull()) return false;
  
  std::ifstream ifs(m_file_name.Data());
  if (!ifs.is_open()) {
    std::cerr << FUNC_NAME << " file open fail: " << m_file_name << std::endl;
    return false;
  }

  Clear();
  ResetRefPlanesToDefault();

  std::string line;
  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream iss(line);
    std::string key;
    iss >> key;

    if (key == "CentralMomentum:") {
      iss >> m_central_momentum;
    } else if (key == "D5ZIn:") {
      iss >> m_d5_z_in;
    } else if (key == "D5ZOut:") {
      iss >> m_d5_z_out;
    } else if (key == "D5Matrix:") {
      Int_t order;
      iss >> order;
      if (order == 1) {
        // R-matrix: 6 rows x 6 cols
        for (Int_t i=0; i<6; ++i) {
          for (Int_t j=0; j<6; ++j) {
            ifs >> m_matrix_1st[i][j];
          }
        }
      } else if (order >= 21 && order <= 25) {
        // T-block for out[order-21]; 6x6 in par-index (j,k)
        const Int_t idx = order - 21;
        for (Int_t i=0; i<6; ++i) {
          for (Int_t j=0; j<6; ++j) {
            ifs >> m_matrix_2nd[idx][i][j];
          }
        }
      }
    }
  }

  m_z_ref_ready = true;
  m_is_ready = true;
  std::cout << FUNC_NAME << " initialized with " << m_file_name
            << " d5_z_in=" << m_d5_z_in
            << " d5_z_out=" << m_d5_z_out << std::endl;
  return true;
}

//_____________________________________________________________________________
Bool_t
D5TransMatrix::Initialize(const TString& filename)
{
  if (m_is_ready && m_file_name == filename) return true;
  SetFileName(filename);
  m_is_ready = false;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
D5TransMatrix::Transport(const Double_t* in, Double_t* out) const
{
  if (!m_is_ready) return false;

  // parin[6] = (x, u, y, v, 0, delta) in matrix units (x,y: cm; u,v: mrad; delta: %)
  const Double_t unit = 10.0; // mm to cm
  Double_t parin[6] = { in[0]/unit, in[1], in[2]/unit, in[3], 0.0, in[4] };
  TMatrixD m_in(6, 1);
  for(Int_t i=0; i<6; ++i) m_in[i][0] = parin[i];

  // p = parin (6x1).  out_i = (R*p)_i + p^T T_i p  (i=0..3 used)
  TMatrixD m_out1st(6, 1);
  m_out1st.Mult(m_matrix_1st, m_in);  // R * p

  Double_t out2nd[6] = {0.,0.,0.,0.,0.,0.};
  for (Int_t i=0; i<5; ++i) {  // out2nd[i] = p^T T_i p
    for (Int_t j=0; j<6; ++j) {
      for (Int_t k=0; k<6; ++k) {
        out2nd[i] += m_matrix_2nd[i][j][k] * parin[j] * parin[k];
      }
    }
  }

  out[0] = (m_out1st[0][0] + out2nd[0]) * unit; // cm to mm
  out[1] = (m_out1st[1][0] + out2nd[1]);
  out[2] = (m_out1st[2][0] + out2nd[2]) * unit; // cm to mm
  out[3] = (m_out1st[3][0] + out2nd[3]);

  return true;
}
