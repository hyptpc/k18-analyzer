// -*- C++ -*-

#include "FiberHit.hh"

#include <limits>
#include <stdexcept>

#include <TMath.h>
#include <TSystem.h>

#include "ConfMan.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "PrintHelper.hh"
#include "RawData.hh"

#include <std_ostream.hh>

namespace
{
const auto qnan = TMath::QuietNaN();
const auto& gHodo = HodoParamMan::GetInstance();
const auto& gPHC  = HodoPHCMan::GetInstance();
}

//_____________________________________________________________________________
FiberHit::FiberHit(HodoRawHit *rhit)
  : HodoHit(rhit),
    m_position(qnan)
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
FiberHit::~FiberHit()
{
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
bool
FiberHit::Calculate()
{
  if(!HodoHit::Calculate())
    return false;

  m_is_clustered.clear();

  Int_t id    = m_raw->DetectorId();
  Int_t plane = m_raw->PlaneId();
  Int_t seg   = m_raw->SegmentId();

  data_t trailing(m_n_ch);
  for(Int_t ch=0; ch<m_n_ch; ++ch){
    const auto& l_cont = m_time_leading[ch];
    m_ctime_leading[ch].clear();
    m_ctime_trailing[ch].clear();
    for(Int_t il=0, nl=l_cont.size(); il<nl; ++il){
      Double_t l = l_cont[il];
      Double_t l_next = (il+1) != nl ? l_cont[il+1] : DBL_MAX;
      Double_t buf = qnan;
      for(const auto& t: m_time_trailing[ch]){
        if(l<t && t<l_next){
          buf = t;
          break;
        }
      }
      trailing[ch].push_back(buf);
      Double_t ctime = qnan;
      Double_t tot = buf - l;
      // for BHT
      Double_t de = TMath::QuietNaN();
      gHodo.GetDeHighGain(id, plane, seg, ch, tot, de);
      m_de_high.at(ch).push_back(de);
      gPHC.DoCorrection(id, plane, seg, ch, l, de, ctime);
      m_ctime_leading[ch].push_back(ctime);
      m_ctime_trailing[ch].push_back(ctime + tot); // no use
      m_is_clustered.push_back(false);
    }
  }
  m_time_trailing = trailing;
  m_is_calculated = true;

  return true;
}

//_____________________________________________________________________________
Double_t
FiberHit::MeanTimeOverThreshold(Int_t j) const
{
  try {
    if(m_n_ch == 1){
      return TimeOverThreshold(HodoRawHit::kUp, j);
    }else{
      return TMath::Sqrt(
        TMath::Abs(TimeOverThreshold(HodoRawHit::kUp, j) *
                   TimeOverThreshold(HodoRawHit::kDown, j)));
    }
  }catch(const std::out_of_range&){
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
void
FiberHit::Print(Option_t* arg) const
{
  PrintHelper helper(3, std::ios::fixed);
  hddaq::cout << FUNC_NAME << " " << arg << std::endl
              << "detector_name = " << m_raw->DetectorName() << std::endl
              << "detector_id   = " << m_raw->DetectorId() << std::endl
              << "plane_name    = " << m_raw->PlaneName()  << std::endl
              << "plane_id      = " << m_raw->PlaneId()    << std::endl
              << "segment_id    = " << m_raw->SegmentId()  << std::endl
              << "n_ch          = " << m_n_ch              << std::endl
              << "de            = " << DeltaE() << std::endl
              << "mt/cmt        = " << MeanTime()
              << " / " << CMeanTime() << std::endl
              << "mtot          = " << MeanTOT() << std::endl
              << "tdiff/ctdiff  = " << TimeDiff()
              << " / " << CTimeDiff() << std::endl;
  for(const auto& data_map: std::map<TString, data_t>
        {{"de-hi  ", m_de_high},      {"de-lo  ", m_de_low},
         {"time-l ", m_time_leading}, {"time-t ", m_time_trailing},
         {"ctime-l", m_ctime_leading}, {"ctime-t", m_ctime_trailing}
        }){
    for(const auto& cont: data_map.second){
      hddaq::cout << " " << data_map.first << ":" << cont.size()
                  << " ";
      for(const auto& val: cont){
        hddaq::cout << val << " ";
      }
    }
    hddaq::cout << std::endl;
  }
  hddaq::cout << " tot    :";
  for(Int_t ch=0; ch<m_n_ch; ++ch){
    for(Int_t j=0, n=GetEntries(ch); j<n; ++j){
      hddaq::cout << " " << TOT(ch, j);
    }
  }
  hddaq::cout << std::endl;
}
