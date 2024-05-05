// -*- C++ -*-

#include "FiberHit.hh"

#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <limits>

#include <std_ostream.hh>

#include "ConfMan.hh"
#include "DCGeomMan.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FuncName.hh"
#include "HodoParamMan.hh"
#include "HodoPHCMan.hh"
#include "PrintHelper.hh"
#include "RawData.hh"

namespace
{
const auto qnan = TMath::QuietNaN();
const auto& gGeom = DCGeomMan::GetInstance();
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
      gPHC.DoCorrection(id, plane, seg, ch, l, tot, ctime);
      m_ctime_leading[ch].push_back(ctime);
      m_ctime_trailing[ch].push_back(ctime + tot); // no use
      m_is_clustered.push_back(false);
    }
  }
  m_time_trailing = trailing;

  // Int_t layer = gGeom.GetDetectorId(DetectorName()+"-"+PlaneName());
  // m_position = gGeom.CalcWirePosition(layer, seg);
  // m_dxdw     = gGeom.dXdW(layer);

  return true;

#if 0
  if(id==113){// CFT ADC
    Double_t nhit_adc = m_raw->GetSizeAdc1();
    if(nhit_adc>0){
      Double_t hi  =  m_raw->GetAdc1();
      Double_t low =  m_raw->GetAdc2();
      Double_t pedeHi  = gHodo.GetP0(id, plane, seg, 0);
      Double_t pedeLow = gHodo.GetP0(id, plane, seg, 1);
      Double_t gainHi  = gHodo.GetP1(id, plane, seg, 0);// pedestal+mip(or peak)
      Double_t gainLow = gHodo.GetP1(id, plane, seg, 1);// pedestal+mip(or peak)
      Double_t Alow = gHodo.GetP0(id,plane,0/*seg*/,3);// same value for the same layer
      Double_t Blow = gHodo.GetP1(id,plane,0/*seg*/,3);// same value for the same layer

      if (hi>0)
	m_adc_hg  = hi  - pedeHi;
      if (low>0)
	m_adc_lg = low - pedeLow;
      //m_mip_hg  = (hi  - pedeHi)/(gainHi  - pedeHi);
      //m_mip_lg = (low - pedeLow)/(gainLow - pedeLow);
      if (m_pedcor_hg>-2000 && hi >0)
	m_adc_hg  = hi  + m_pedcor_hg;
      if (m_pedcor_lg>-2000 && low >0)
	m_adc_lg  = low  + m_pedcor_lg;

      if (m_adc_hg>0/* && gainHi > 0*/)
	m_mip_hg  = m_adc_hg/gainHi ;
      if (m_adc_lg>0/* && gainLow >0*/)
	m_mip_lg = m_adc_lg/gainLow;
      if(m_mip_lg>0){
	m_dE_lg = -(Alow/Blow) * log(1. - m_mip_lg/Alow);// [MeV]
	if(1-m_mip_lg/Alow<0){ // when pe is too big
	  m_dE_lg = -(Alow/Blow) * log(1. - (Alow-0.001)/Alow);// [MeV] almost max
	}
      }else{
	m_dE_lg = 0;// [MeV]
      }
#if 0
      if(m_dE_lg > 10){
        std::cout << "layer = " << plane << ", seg = " << seg << ", adcLow = "
                  << m_adc_lg << ", dE = " << m_dE_lg << ", mip_lg = "
                  << m_mip_lg << ", gainLow = " << gainLow << ", gainHi = "
                  << gainHi << std::endl;
      }
#endif
      for(auto& pair: m_pair_cont){
	Double_t time= pair.time_l;
	Double_t ctime = -100;
	if (m_adc_hg>20) {
	  gPHC.DoCorrection(id, plane, seg, m_ud, time, m_adc_hg, ctime);
	  pair.ctime_l = ctime;
	} else
	  pair.ctime_l = time;
      }
    }
  }
#endif
  m_is_calculated = true;

  HodoHit::Print();
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
  for(const auto data_map: std::map<TString, data_t>
        {{"de-hi  ", m_de_high},      {"de-lo  ", m_de_low},
         {"time-l ", m_time_leading}, {"time-t ", m_time_trailing},
         {"ctime-l", m_ctime_leading}, {"ctime-t", m_ctime_trailing}
        }){
    for(const auto& cont: data_map.second){
      hddaq::cout << " " << data_map.first << ":" << cont.size()
                  << " ";
      std::copy(cont.begin(), cont.end(),
                std::ostream_iterator<Double_t>(hddaq::cout," "));
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
