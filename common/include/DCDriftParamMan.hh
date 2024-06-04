// -*- C++ -*-

#ifndef DC_DRIFT_PARAM_MAN_HH
#define DC_DRIFT_PARAM_MAN_HH

#include <map>
#include <vector>

#include <TGraph.h>
#include <TString.h>

//_____________________________________________________________________________
class DCDriftParamMan
{
public:
  static const TString&   ClassName();
  static DCDriftParamMan& GetInstance();
  ~DCDriftParamMan();

private:
  DCDriftParamMan();
  DCDriftParamMan(const DCDriftParamMan&);
  DCDriftParamMan& operator =(const DCDriftParamMan&);

private:
  using DCDriftFunctionContainer = std::map<TString, TGraph*>;
  using DCDFC = DCDriftFunctionContainer;
  // typedef DCDriftContainer::const_iterator DCDriftIterator;
  Bool_t  m_is_ready;
  TString m_file_name;
  DCDFC   m_container;

public:
  Bool_t CalcDrift(const TString& detector_name, Int_t plane_id,
                   Double_t wire_id, Double_t ctime,
                   Double_t& dt, Double_t& dl) const;
  Bool_t Initialize();
  Bool_t Initialize(const TString& file_name);
  Bool_t IsReady() const { return m_is_ready; }
  void   SetFileName(const TString& file_name) { m_file_name = file_name; }

private:
  void                ClearElements();
  const TGraph* GetParameter(const TString& detector_name,
                             Int_t plane_id, Int_t wire_id) const;
};

//_____________________________________________________________________________
inline const TString&
DCDriftParamMan::ClassName()
{
  static TString s_name("DCDriftParamMan");
  return s_name;
}

//_____________________________________________________________________________
inline DCDriftParamMan&
DCDriftParamMan::GetInstance()
{
  static DCDriftParamMan s_instance;
  return s_instance;
}

#endif
