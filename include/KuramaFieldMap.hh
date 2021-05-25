// -*- C++ -*-

#ifndef KURAMA_FIELD_MAP_HH
#define KURAMA_FIELD_MAP_HH

#include <vector>
#include <TString.h>

//_____________________________________________________________________________
class KuramaFieldMap
{
public:
  static const TString& ClassName();
  KuramaFieldMap(const TString& file_name);
  ~KuramaFieldMap();

private:
  KuramaFieldMap(const KuramaFieldMap&);
  KuramaFieldMap& operator =(const KuramaFieldMap&);

private:
  struct XYZ { Double_t x, y, z; };
  typedef std::vector<std::vector<std::vector<XYZ>>> Field;
  Bool_t   m_is_ready;
  TString  m_file_name;
  Field    B;
  Int_t    Nx;
  Int_t    Ny;
  Int_t    Nz;
  Double_t X0;
  Double_t Y0;
  Double_t Z0;
  Double_t dX;
  Double_t dY;
  Double_t dZ;

public:
  Bool_t Initialize();
  Bool_t IsReady() const { return m_is_ready; }
  Bool_t GetFieldValue(const Double_t pointCM[3], Double_t* BfieldTesla) const;

private:
  void ClearField();
};

//_____________________________________________________________________________
inline const TString&
KuramaFieldMap::ClassName()
{
  static TString s_name("KuramaFieldMap");
  return s_name;
}

#endif
