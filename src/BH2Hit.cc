// -*- C++ -*-

#include "BH2Hit.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <std_ostream.hh>

#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "RawData.hh"

namespace
{
const auto& gHodo = HodoParamMan::GetInstance();
// const auto& gPHC = HodoPHCMan::GetInstance();
}

//_____________________________________________________________________________
BH2Hit::BH2Hit(HodoRawHit *rhit, double max_time_diff)
  : HodoHit(rhit, max_time_diff)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
BH2Hit::~BH2Hit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
bool
BH2Hit::Calculate()
{
  if(!HodoHit::Calculate())
    return false;

  Int_t cid  = m_raw->DetectorId();
  Int_t plid = m_raw->PlaneId();
  Int_t seg  = m_raw->SegmentId();
  gHodo.GetTime(cid, plid, seg, HodoRawHit::kExtra, 0., m_time_offset);

  return true;
}

//_____________________________________________________________________________
Double_t
BH2Hit::UTime0(Int_t i) const
{
  return m_time_leading.at(HodoRawHit::kUp).at(i) + m_time_offset;
}

//_____________________________________________________________________________
Double_t
BH2Hit::UCTime0(Int_t i) const
{
  return m_ctime_leading.at(HodoRawHit::kUp).at(i) + m_time_offset;
}

//_____________________________________________________________________________
Double_t
BH2Hit::DTime0(Int_t i) const
{
  return m_time_leading.at(HodoRawHit::kDown).at(i) + m_time_offset;
}

//_____________________________________________________________________________
Double_t
BH2Hit::DCTime0(Int_t i) const
{
  return m_ctime_leading.at(HodoRawHit::kDown).at(i) + m_time_offset;
}
