// -*- C++ -*-

#ifndef D5_TRANS_MATRIX_HH
#define D5_TRANS_MATRIX_HH

#include <TMatrixD.h>
#include <TString.h>

class D5TransMatrix
{
public:
  // Default reference planes [mm] in FF (overridden by D5TransferMatrix.param)
  static constexpr Double_t RefZIn  = 0.0;
  static constexpr Double_t RefZOut = -1300.885;

  static const TString& ClassName();
  static D5TransMatrix& GetInstance();
  ~D5TransMatrix();

private:
  D5TransMatrix();
  D5TransMatrix(const D5TransMatrix&);
  D5TransMatrix& operator=(const D5TransMatrix&);

private:
  Bool_t   m_is_ready;
  Bool_t   m_z_ref_ready;
  TString  m_file_name;
  Double_t m_central_momentum;
  Double_t m_d5_z_in;
  Double_t m_d5_z_out;
  // 1st order R: 6x6, state (x, u, y, v, _, delta) -> out (same 6-vector)
  TMatrixD m_matrix_1st;
  // 2nd order T_i (6x6): component i of aberration is quadratic form p^T T_i p
  TMatrixD m_matrix_2nd[5];

  void ResetRefPlanesToDefault()
  {
    m_d5_z_in = RefZIn;
    m_d5_z_out = RefZOut;
  }

public:
  void     Clear();
  Double_t GetCentralMomentum() const { return m_central_momentum; }
  Bool_t   Initialize();
  Bool_t   Initialize(const TString& filename);
  Bool_t   IsReady() const { return m_is_ready; }
  Bool_t   IsRefPlaneReady() const { return m_z_ref_ready; }
  Double_t GetD5ZIn() const { return m_d5_z_in; }
  Double_t GetD5ZOut() const { return m_d5_z_out; }
  void     SetFileName(const TString& filename) { m_file_name = filename; }

  // Transport upstream track parameters to downstream
  // in: (x, u, y, v, delta) [mm, mrad, mm, mrad, %]
  // out: (x, u, y, v) [mm, mrad, mm, mrad]
  Bool_t   Transport(const Double_t* in, Double_t* out) const;
};

//_____________________________________________________________________________
inline const TString&
D5TransMatrix::ClassName()
{
  static TString s_name("D5TransMatrix");
  return s_name;
}

#endif
