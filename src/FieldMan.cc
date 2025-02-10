// -*- C++ -*-

#include "FieldMan.hh"

#include <cmath>

#include <std_ostream.hh>

#include "FieldElements.hh"
#include "FuncName.hh"
#include "S2sFieldMap.hh"
#include "ConfMan.hh"

namespace
{
const auto& gConf = ConfMan::GetInstance();
const auto& valueNMR  = ConfMan::Get<Double_t>("FLDNMR");
const auto& valueCalc = ConfMan::Get<Double_t>("FLDCALC");
const auto& valueHSHall = ConfMan::Get<Double_t>("HSFLDHALL");
const auto& valueHSCalc = ConfMan::Get<Double_t>("HSFLDCALC");
}

namespace
{
const Double_t Delta = 0.1;
}

//_____________________________________________________________________________
FieldMan::FieldMan()
  : m_is_ready(false),
    m_s2s_map(nullptr),
    m_shs_map(nullptr)
{
}

//_____________________________________________________________________________
FieldMan::~FieldMan()
{
  if(m_s2s_map) delete m_s2s_map;
  if(m_shs_map)    delete m_shs_map;
}

//_____________________________________________________________________________
Bool_t
FieldMan::Initialize()
{
  if(m_is_ready){
    hddaq::cerr << FUNC_NAME << " already initialied" << std::endl;
    return false;
  }

  if(m_s2s_map)
    delete m_s2s_map;

  if(!m_file_name_s2s.IsNull()){
    m_s2s_map = new S2sFieldMap(m_file_name_s2s, valueNMR, valueCalc);
    if(!m_s2s_map->Initialize())
      return false;
  }

  if(m_shs_map)
    delete m_shs_map;

  std::cout << m_file_name_shs << " " << m_file_name_shs.IsNull() << std::endl;
  if(!m_file_name_shs.IsNull()){
    std::cout << m_file_name_shs << " " << m_file_name_shs.IsNull() << std::endl;
    m_shs_map = new S2sFieldMap(m_file_name_shs, valueHSHall, valueHSCalc);
    if(!m_shs_map->Initialize())
      return false;
  }

  m_is_ready = (m_s2s_map || m_shs_map);
  return m_is_ready;
}

//_____________________________________________________________________________
Bool_t
FieldMan::Initialize(const TString& file_name_s2s)
{
  m_file_name_s2s = file_name_s2s;
  return Initialize();
}

//_____________________________________________________________________________
Bool_t
FieldMan::Initialize(const TString& file_name_s2s,
                     const TString& file_name_shs)
{
  m_file_name_s2s = file_name_s2s;
  m_file_name_shs    = file_name_shs;
  return Initialize();
}

//_____________________________________________________________________________
TVector3
FieldMan::GetField(const TVector3& position) const
{
  TVector3 field(0., 0., 0.);
  if(m_s2s_map // && m_shs_map
    ){
    Double_t p[3], b_s2s[3]; //, b_shs[3];
    p[0] = position.x()*0.1;
    p[1] = position.y()*0.1;
    p[2] = position.z()*0.1;
    if(m_s2s_map->GetFieldValue(p, b_s2s)//  &&
       // m_shs_map->GetFieldValue(p, b_shs)
      ){
      field.SetX(b_s2s[0]);
      field.SetY(b_s2s[1]);
      field.SetZ(b_s2s[2]);
    }
  }

#if 1
  FEIterator itr, itr_end = m_element_list.end();
  for(itr=m_element_list.begin(); itr!=itr_end; ++itr){
    if((*itr)->ExistField(position))
      field += (*itr)->GetField(position);
  }
#endif

  return field;
}

//_____________________________________________________________________________
TVector3
FieldMan::GetdBdX(const TVector3& position) const
{
  TVector3 p1 = position + TVector3(Delta, 0., 0.);
  TVector3 p2 = position - TVector3(Delta, 0., 0.);
  TVector3 B1 = GetField(p1);
  TVector3 B2 = GetField(p2);
  return 0.5/Delta*(B1-B2);
}

//_____________________________________________________________________________
TVector3
FieldMan::GetdBdY(const TVector3& position) const
{
  TVector3 p1 = position + TVector3(0., Delta, 0.);
  TVector3 p2 = position - TVector3(0., Delta, 0.);
  TVector3 B1 = GetField(p1);
  TVector3 B2 = GetField(p2);
  return 0.5/Delta*(B1-B2);
}

//_____________________________________________________________________________
TVector3
FieldMan::GetdBdZ(const TVector3& position) const
{
  TVector3 p1 = position + TVector3(0., 0., Delta);
  TVector3 p2 = position - TVector3(0., 0., Delta);
  TVector3 B1 = GetField(p1);
  TVector3 B2 = GetField(p2);
  return 0.5/Delta*(B1-B2);
}

//_____________________________________________________________________________
void
FieldMan::ClearElementsList()
{
  m_element_list.clear();
}

//_____________________________________________________________________________
void
FieldMan::AddElement(FieldElements *element)
{
  m_element_list.push_back(element);
}

//_____________________________________________________________________________
Double_t
FieldMan::StepSize(const TVector3 &position,
                   Double_t def_step_size, Double_t min_step_size) const
{
  Double_t d = TMath::Abs(def_step_size);
  Double_t s = def_step_size/d;
  min_step_size = TMath::Abs(min_step_size);

  Bool_t flag = true;
  FEIterator itr, itr_end = m_element_list.end();
  while (flag && d>min_step_size){
    for(itr=m_element_list.begin(); itr!=itr_end; ++itr){
      if((*itr)->CheckRegion(position, d) != FieldElements::FERSurface()){
        flag=false;
        break;
      }
    }
    d *= 0.5;
  }
  if(flag)
    return s*min_step_size;
  else
    return s*d;
}
