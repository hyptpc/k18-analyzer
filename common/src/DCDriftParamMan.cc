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

#include <TMath.h>

#include <std_ostream.hh>

#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"

namespace
{
const auto qnan = TMath::QuietNaN();

inline void
DecodeKey(Int_t key, Int_t& detector_id, Int_t& plane_id, Int_t& wire_id)
{
  wire_id  = key & 0xfff;
  plane_id = (key >> 12) && 0xfff;
  detector_id = (key >> 24);
}

inline Int_t
MakeKey(Int_t detector_id, Int_t plane_id, Double_t wire_id)
{
  return (detector_id << 24) | (plane_id << 12) | Int_t(wire_id);
}
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

  std::ifstream ifs(m_file_name);
  if(!ifs.is_open()){
    hddaq::cerr << FUNC_NAME
                << " file open fail : " << m_file_name << std::endl;
    return false;
  }

  ClearElements();

  TString line;
  while(ifs.good() && line.ReadLine(ifs)){
    if(line.IsNull() || line[0]=='#') continue;
    std::istringstream iss(line.Data());
    // wid is not used at present
    // reserved for future updates
    Int_t cid, pid, wid, type, np;
    Double_t param;
    std::vector<Double_t> q;
    if(iss >> cid >> pid >> wid >> type >> np){
      while(iss >> param)
	q.push_back(param);
    }
    if(q.size() < 2){
      hddaq::cerr << FUNC_NAME << " format is wrong : " << line << std::endl;
      continue;
    }
    Int_t key = MakeKey(cid, pid, 0);
    auto record = new DCDriftParamRecord(type, np, q);
    if(m_container[key]){
      hddaq::cerr << FUNC_NAME << " "
		  << "duplicated key is deleted : " << key << std::endl;
      delete m_container[key];
    }
    m_container[key] = record;
  }

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
DCDriftParamRecord*
DCDriftParamMan::GetParameter(Int_t detector_id,
                              Int_t plane_id, Double_t wire_id) const
{
  wire_id = 0;
  Int_t key = MakeKey(detector_id, plane_id, wire_id);
  auto itr = m_container.find(key);
  if(itr != m_container.end())
    return itr->second;
  else
    return nullptr;
}

//_____________________________________________________________________________
Double_t
DCDriftParamMan::DriftLength1(Double_t dt, Double_t vel)
{
  return dt*vel;
}

//_____________________________________________________________________________
Double_t
DCDriftParamMan::DriftLength2(Double_t dt,
                              Double_t p1, Double_t p2, Double_t p3,
                              Double_t st, Double_t p5, Double_t p6)
{
  Double_t dtmax=10.+p2+1./p6;
  Double_t dl;
  Double_t alph=-0.5*p5*st+0.5*st*p3*p5*(p1-st)
    +0.5*p3*p5*(p1-st)*(p1-st);
  if(dt<-10. || dt>dtmax+10.)
    dl = -500.;
  else if(dt<st)
    dl = 0.5*(p5+p3*p5*(p1-st))/st*dt*dt;
  else if(dt<p1)
    dl = alph+p5*dt-0.5*p3*p5*(p1-dt)*(p1-dt);
  else if(dt<p2)
    dl = alph+p5*dt;
  else
    dl = alph+p5*dt-0.5*p6*p5*(dt-p2)*(dt-p2);
  return dl;
}

//_____________________________________________________________________________
//DL = a0 + a1*Dt + a2*Dt^2
Double_t
DCDriftParamMan::DriftLength3(Double_t dt, Double_t p1, Double_t p2,
                              Int_t plane_id)
{
  if (plane_id>=1 && plane_id <=4) {
    if(dt > 130.)
      return qnan;
  } else if (plane_id == 5) {
    if (dt > 500.)
      ;
  } else if (plane_id>=6 && plane_id <=11) {
    if (dt > 80.)
      return qnan;
  } else if (plane_id>=101 && plane_id <=124) {
    if (dt > 80.)
      return qnan;
  }

  return dt*p1+dt*dt*p2;
}

//_____________________________________________________________________________
Double_t
DCDriftParamMan::DriftLength4(Double_t dt,
                              Double_t p1, Double_t p2, Double_t p3)
{
  if (dt < p2)
    return dt*p1;
  else
    return p3*(dt-p2) + p1*p2;
}

//_____________________________________________________________________________
Double_t
DCDriftParamMan::DriftLength5(Double_t dt, Double_t p1, Double_t p2,
                              Double_t p3, Double_t p4, Double_t p5)
{
  if (dt < p2)
    return dt*p1;
  else if (dt >= p2 && dt < p4)
    return p3*(dt-p2) + p1*p2;
  else
    return p5*(dt-p4) + p3*(p4-p2) + p1*p2;
}

//_____________________________________________________________________________
////////////////Now using
// DL = a0 + a1*Dt + a2*Dt^2 + a3*Dt^3 + a4*Dt^4 + a5*Dt^5
Double_t
DCDriftParamMan::DriftLength6(Int_t plane_id, Double_t dt,
                              Double_t p1, Double_t p2, Double_t p3,
                              Double_t p4, Double_t p5)
{
  Double_t dl = (p1*dt
                 + p2*dt*dt
                 + p3*dt*dt*dt
                 + p4*dt*dt*dt*dt
                 + p5*dt*dt*dt*dt*dt);
  return dl;
}

//_____________________________________________________________________________
Bool_t
DCDriftParamMan::CalcDrift(Int_t detector_id, Int_t plane_id,
                           Double_t wire_id, Double_t ctime,
                           Double_t& dt, Double_t& dl) const
{
  auto record = GetParameter(detector_id, plane_id, wire_id);
  if(!record){
    hddaq::cerr << FUNC_NAME << " No record. "
		<< " detector_id=" << std::setw(3) << detector_id
		<< " plane_id=" << std::setw(3) << plane_id
		<< " wire_id="  << std::setw(3) << wire_id << std::endl;
    return false;
  }

  Int_t type = record->type;
  // Int_t np   = record->np;
  std::vector<Double_t> p = record->param;

  dt = ctime + p[0];

  switch(type){
  case 1:
    dl=DriftLength1(dt, p[1]);
    return true;
  case 2:
    dl=DriftLength2(dt, p[1], p[2], p[3], p[4],
                    p[5], p[6]);
    return true;
  case 3:
    dl=DriftLength3(dt, p[1], p[2], plane_id);
    return true;
  case 4:
    dl=DriftLength4(dt, p[1], p[2], p[3]);
    return true;
  case 5:
    dl=DriftLength5(dt, p[1], p[2], p[3], p[4], p[5]);
    return true;
  case 6:
    dl=DriftLength6(plane_id, dt, p[1], p[2], p[3], p[4], p[5]);
    return true;
  default:
    hddaq::cerr << FUNC_NAME << " invalid type : " << type << std::endl;
    return false;
  }
}
