// -*- C++ -*-

#ifndef FIELD_MAP_HH
#define FIELD_MAP_HH

#include <vector>
#include <TString.h>

//_____________________________________________________________________________
class FieldMap
{
public:
  static const TString& ClassName();
  FieldMap(const TString& file_name);
  FieldMap(const TString& file_name, const Double_t measure, const Double_t calc);
  ~FieldMap();

private:
  FieldMap(const FieldMap&);
  FieldMap& operator =(const FieldMap&);

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

  Double_t valueMeasure;
  Double_t valueCalc;
  void SetValueMeasure(const Double_t val) { valueMeasure = val; }
  void SetValueCalc(const Double_t val) { valueCalc = val; }

public:
  Bool_t Initialize();
  Bool_t IsReady() const { return m_is_ready; }
  Bool_t GetFieldValue(const Double_t pointCM[3], Double_t* BfieldTesla) const;

private:
  void ClearField();
};

//_____________________________________________________________________________
inline const TString&
FieldMap::ClassName()
{
  static TString s_name("FieldMap");
  return s_name;
}

#endif
