// -*- C++ -*-

#include "DCDriftParamMan.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <TFile.h>
#include <TKey.h>
#include <TMath.h>
#include <TObjString.h>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"

namespace
{
const auto qnan = TMath::QuietNaN();
}

//_____________________________________________________________________________
DCDriftParamMan::DCDriftParamMan()
  : m_is_ready(false),
    m_file_name(),
    m_container()
{
}

//_____________________________________________________________________________
DCDriftParamMan::~DCDriftParamMan()
{
  ClearElements();
}

//_____________________________________________________________________________
void
DCDriftParamMan::ClearElements()
{
  del::ClearMap(m_container);
}

//_____________________________________________________________________________
Bool_t
DCDriftParamMan::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  auto prev_file = gFile;

  auto f = TFile::Open(m_file_name);
  if(!f || !f->IsOpen()){
    hddaq::cerr << FUNC_NAME
                << " file open fail : " << m_file_name << std::endl;
    return false;
  }

  ClearElements();

  TIter itr(f->GetListOfKeys());
  for(TKey* key=(TKey*)itr(); itr!=TIter::End(); key=(TKey*)itr()){
    auto array = TString(key->GetName()).Tokenize("_");
    TString k;
    k += dynamic_cast<TObjString*>(array->First())->GetString();
    k += "_";
    k += dynamic_cast<TObjString*>(array->Last())->GetString();
    m_container[k] = dynamic_cast<TGraph*>(key->ReadObj());
  }

  f->Close();
  prev_file->cd();

  m_is_ready = true;
  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
DCDriftParamMan::Initialize(const TString& file_name)
{
  m_file_name = file_name;
  return Initialize();
}

//_____________________________________________________________________________
const TGraph*
DCDriftParamMan::GetParameter(const TString& detector_name,
                              Int_t plane_id, Int_t /* wire_id */) const
{
  return m_container.at(Form("%s_plane%d", detector_name.Data(), plane_id));
}

//_____________________________________________________________________________
Bool_t
DCDriftParamMan::CalcDrift(const TString& detector_name, Int_t plane_id,
                           Double_t wire_id, Double_t ctime,
                           Double_t& dt, Double_t& dl) const
{
  auto g1 = GetParameter(detector_name, plane_id, wire_id);
  dt = ctime;
  dl = g1->Eval(dt, nullptr, "S");
  return true;
}
