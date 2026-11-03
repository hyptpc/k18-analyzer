// -*- C++ -*-

#include "HodoAnalyzer.hh"

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

#include "BH2Hit.hh"
#include "CherenkovHelper.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "FiberCluster.hh"
#include "FiberHit.hh"
#include "FuncName.hh"
#include "HodoCluster.hh"
#include "HodoHit.hh"
#include "HodoParamMan.hh"
#include "HodoRawHit.hh"
#include "RawData.hh"
#include "UserParamMan.hh"

#include <std_ostream.hh>

#define REQDE   0

namespace
{
const auto& gUser = UserParamMan::GetInstance();
const auto& gHodo = HodoParamMan::GetInstance();
const Int_t kMinAdcCmChannels = 3;
}

//_____________________________________________________________________________
HodoAnalyzer::HodoAnalyzer(const RawData& raw_data)
  : m_raw_data(&raw_data),
    m_hodo_hit_collection(),
    m_hodo_cluster_collection(),
    m_adc_cm_det(),
    m_adc_cm()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
void
HodoAnalyzer::SetAdcCommonModeDetector(const TString& name)
{
  if (std::find(m_adc_cm_det.begin(), m_adc_cm_det.end(), name)
      == m_adc_cm_det.end())
    m_adc_cm_det.push_back(name);
}

//_____________________________________________________________________________
const HodoAnalyzer::AdcCommonMode&
HodoAnalyzer::GetAdcCommonMode(const TString& name) const
{
  static const AdcCommonMode zero;
  const auto it = m_adc_cm.find(name);
  return (it != m_adc_cm.end()) ? it->second : zero;
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::AdcCommonModeU(const TString& name) const
{
  return GetAdcCommonMode(name).u;
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::AdcCommonModeD(const TString& name) const
{
  return GetAdcCommonMode(name).d;
}

//_____________________________________________________________________________
Int_t
HodoAnalyzer::AdcCommonModeNChU(const TString& name) const
{
  return GetAdcCommonMode(name).nch_u;
}

//_____________________________________________________________________________
Int_t
HodoAnalyzer::AdcCommonModeNChD(const TString& name) const
{
  return GetAdcCommonMode(name).nch_d;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::UpdateAdcCommonMode(const TString& name)
{
  const TString tdc_key = name + "_TDC";
  if (!gUser.Has(tdc_key)) {
    m_adc_cm.erase(name);
    static std::set<TString> s_warned;
    if (s_warned.insert(name).second) {
      hddaq::cerr << FUNC_NAME
                  << " ADC CM enabled but missing UserParam \""
                  << tdc_key << "\"" << std::endl;
    }
    return false;
  }

  AdcCommonMode cm;
  std::vector<Double_t> adc_u;
  std::vector<Double_t> adc_d;
  for (const auto& rhit : m_raw_data->GetHodoRawHitContainer(name)) {
    if (!rhit)
      continue;
    const Int_t det_id = rhit->DetectorId();
    const Int_t plane = rhit->PlaneId();
    const Int_t seg = rhit->SegmentId();
    for (const Int_t ch : {HodoRawHit::kUp, HodoRawHit::kDown}) {
      Bool_t has_signal = false;
      for (const auto& tdc : rhit->GetArrayTdcLeading(ch)) {
        if (gUser.IsInRange(tdc_key, tdc)) {
          has_signal = true;
          break;
        }
      }
      if (has_signal)
        continue;
      const Double_t ped = gHodo.GetAdcPedestal(det_id, plane, seg, ch);
      if (!std::isfinite(ped))
        continue;
      const auto& adcs = rhit->GetArrayAdcHigh(ch);
      if (adcs.empty())
        continue;
      const Double_t ped_sub = adcs.at(0) - ped;
      if (ch == HodoRawHit::kUp)
        adc_u.push_back(ped_sub);
      else
        adc_d.push_back(ped_sub);
    }
  }

  // Sample counts before the min-channel cut. A short side may still get
  // the other side's CM via fallback; nch stays the raw count.
  cm.nch_u = static_cast<Int_t>(adc_u.size());
  cm.nch_d = static_cast<Int_t>(adc_d.size());

  auto median = [](std::vector<Double_t>& v) -> Double_t {
    if (v.empty())
      return 0.;
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    if (n % 2)
      return v[n / 2];
    return 0.5 * (v[n / 2 - 1] + v[n / 2]);
  };

  const Bool_t ok_u = cm.nch_u >= kMinAdcCmChannels;
  const Bool_t ok_d = cm.nch_d >= kMinAdcCmChannels;
  const Double_t cm_u = ok_u ? median(adc_u) : 0.;
  const Double_t cm_d = ok_d ? median(adc_d) : 0.;
  cm.u = ok_u ? cm_u : (ok_d ? cm_d : 0.);
  cm.d = ok_d ? cm_d : (ok_u ? cm_u : 0.);
  m_adc_cm[name] = cm;
  return true;
}

//_____________________________________________________________________________
HodoAnalyzer::~HodoAnalyzer()
{
  for(auto& elem: m_hodo_hit_collection)
    del::ClearContainer(elem.second);
  for(auto& elem: m_hodo_cluster_collection)
    del::ClearContainer(elem.second);
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
std::vector<Double_t>
HodoAnalyzer::GetOfflineNpe(const TString& name) const
{
  auto it = m_offline_npe.find(name);
  if(it == m_offline_npe.end())
    it = m_offline_npe.emplace(name, CherenkovHelper::Compute(*m_raw_data, name)).first;
  return it->second;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcHit(const TString& name, Bool_t applyRecursively)
{
  for(auto& cl: GetHitContainer(name)){
    cl->ReCalc(applyRecursively);
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcCluster(const TString& name, Bool_t applyRecursively)
{
  if (applyRecursively)
    ReCalcHit(name);
  for (auto& cl : GetClusterContainer(name))
    cl->ReCalc(false);
  return true;
}

//_____________________________________________________________________________
Bool_t
HodoAnalyzer::ReCalcAll()
{
  return true;
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCut(const TString& name)
{
  const Double_t min = gUser.Get(name + "_Time", 0);
  const Double_t max = gUser.Get(name + "_Time", 1);
  TimeCut(m_hodo_hit_collection[name], min, max);
  TimeCut(m_hodo_cluster_collection[name], min, max);
}

//_____________________________________________________________________________
void
HodoAnalyzer::TimeCut(const TString& name, Double_t min, Double_t max)
{
  TimeCut(m_hodo_hit_collection[name], min, max);
  TimeCut(m_hodo_cluster_collection[name], min, max);
}

//_____________________________________________________________________________
// Implementation of Time cut for the cluster container
template <typename T>
void
HodoAnalyzer::TimeCut(std::vector<T>& cont,
                      Double_t min, Double_t max)
{
  std::vector<T> DeleteCand;
  std::vector<T> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t ctime = cont.at(i)->CMeanTime();
    if(min < ctime && ctime < max){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
void
HodoAnalyzer::TotCut(const TString& name, Double_t min, Double_t max,
                     Bool_t adopt_nan)
{
  TotCut(m_hodo_cluster_collection[name], min, max, adopt_nan);
}

//_____________________________________________________________________________
// Implementation of tot cut for the cluster container
template <typename T>
void
HodoAnalyzer::TotCut(std::vector<T>& cont,
                     Double_t min, Double_t max, Bool_t adopt_nan)
{
  std::vector<T> DeleteCand;
  std::vector<T> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t tot = cont.at(i)->TOT();
    if(TMath::IsNaN(tot) && adopt_nan){
      ValidCand.push_back(cont.at(i));
    }else if(min < tot && tot < max){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
// Implementation of tot cut for the cluster container
// template <typename T>
// void
// HodoAnalyzer::WidthCutR(std::vector<T>& cont,
//                         Double_t min, Double_t max, Bool_t adopt_nan)
// {
//   std::vector<T> DeleteCand;
//   std::vector<T> ValidCand;
//   for(Int_t i=0, n=cont.size(); i<n; ++i){
//     Double_t width = cont.at(i)->Width();
//     //Double_t width = -cont.at(i)->minWidth();//reverse

//     if(TMath::IsNaN(width) && adopt_nan){
//       ValidCand.push_back(cont.at(i));
//     }else if(min < width && width < max){
//       ValidCand.push_back(cont.at(i));
//     }else{
//       DeleteCand.push_back(cont.at(i));
//     }
//   }

//   del::ClearContainer(DeleteCand);

//   cont.clear();
//   cont.resize(ValidCand.size());
//   std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
//   ValidCand.clear();
// }

//_____________________________________________________________________________
// Implementation of ADC cut for the cluster container
template <typename T>
void
HodoAnalyzer::AdcCut(std::vector<T>& cont,
                     Double_t amin, Double_t amax)
{
  std::vector<T> DeleteCand;
  std::vector<T> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t adc = cont.at(i)->MaxAdcLow();
    if(amin < adc && adc < amax){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}
//_____________________________________________________________________________
void
HodoAnalyzer::DeCut(const TString& name, Double_t min, Double_t max)
{
  DeCut(m_hodo_cluster_collection[name], min, max);
}

//_____________________________________________________________________________
// Implementation of Time cut for the cluster container
template <typename T>
void
HodoAnalyzer::DeCut(std::vector<T>& cont,
                      Double_t min, Double_t max)
{
  std::vector<T> DeleteCand;
  std::vector<T> ValidCand;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    Double_t de = cont.at(i)->DeltaE();
    if(min < de && de < max){
      ValidCand.push_back(cont.at(i));
    }else{
      DeleteCand.push_back(cont.at(i));
    }
  }

  del::ClearContainer(DeleteCand);

  cont.clear();
  cont.resize(ValidCand.size());
  std::copy(ValidCand.begin(), ValidCand.end(), cont.begin());
  ValidCand.clear();
}

//_____________________________________________________________________________
const HodoCluster*
HodoAnalyzer::GetTime0Cluster() const
{
  const HodoCluster* cl_time0 = nullptr;
  Double_t min_mt = DBL_MAX;
  for(const auto& cluster : GetClusterContainer("BH2")){
    for(Int_t j=0, m=cluster->ClusterSize(); j<m; ++j){
      Double_t mt = cluster->MeanTimeAt(j);
      if(true
         && TMath::Abs(mt) < TMath::Abs(min_mt)
         && gUser.IsInRange("T0_Time", mt)
#if REQDE
         && gUser.IsInRange("T0_DeltaE", cluster->DeltaE())
#endif
         ){
        min_mt = mt;
        cl_time0 = cluster;
      }
    }
  }
  return cl_time0;
}

//_____________________________________________________________________________
const HodoCluster*
HodoAnalyzer::GetBtof0Cluster() const
{
  const HodoCluster* cl_btof0 = nullptr;
  Double_t time0 = Time0();
  if(TMath::IsNaN(time0)) return nullptr;
  Double_t min_btof = DBL_MAX;
  for(Int_t i=0, n=GetNClusters("BHT"); i<n; ++i){
    const auto& cluster = GetCluster("BHT", i);
    for(Int_t j=0, m=cluster->ClusterSize(); j<m; ++j){
      Double_t cmt  = cluster->CMeanTimeAt(j);
      Double_t btof = cmt - time0;
      if(true
         && TMath::Abs(btof) < TMath::Abs(min_btof)
         && gUser.IsInRange("BTOF", btof)
#if REQDE
         && gUser.IsInRange("BHT_DeltaE", cluster->DeltaE())
#endif
         ){
        min_btof = btof;
        cl_btof0 = cluster;
      }
    }
  }
  return cl_btof0;
}

//_____________________________________________________________________________
const HodoCluster*
HodoAnalyzer::GetFtof0Cluster() const
{
  const HodoCluster* cl_ftof0 = nullptr;
  Double_t time0 = Time0();
  if(TMath::IsNaN(time0)) return nullptr;
  Double_t min_ftof = DBL_MAX;
  for(Int_t i=0, n=GetNClusters("CVC"); i<n; ++i){
    const auto& cluster = GetCluster("CVC", i);
    for(Int_t j=0, m=cluster->ClusterSize(); j<m; ++j){
      Double_t cmt  = cluster->CMeanTimeAt(j);
      Double_t ftof = cmt - time0;
      if(true
         && TMath::Abs(ftof) < TMath::Abs(min_ftof)
         && gUser.IsInRange("FTOF", ftof)
         && gUser.IsInRange("CVC_TimeDiff", cluster->TimeDiffAt(j))
#if REQDE
         && gUser.IsInRange("CVC_DeltaE", cluster->DeltaE())
#endif
         ){
        min_ftof = ftof;
        cl_ftof0 = cluster;
      }
    }
  }
  return cl_ftof0;
}

//_____________________________________________________________________________
const HodoHC&
HodoAnalyzer::GetHitContainer(const TString& name) const
{
  auto itr = m_hodo_hit_collection.find(name);
  if(itr == m_hodo_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoHC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
const HodoCC&
HodoAnalyzer::GetClusterContainer(const TString& name) const
{
  auto itr = m_hodo_cluster_collection.find(name);
  if(itr == m_hodo_cluster_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoCC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Time0() const
{
  const auto& cl = GetTime0Cluster();
  if(cl) return cl->CTime0();
  else   return TMath::QuietNaN();
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Btof0() const
{
  const auto& cl = GetBtof0Cluster();
  if(cl){
    Double_t btof0 = DBL_MAX;
    for(Int_t j=0, m=cl->ClusterSize(); j<m; ++j){
      Double_t cmt  = cl->CMeanTimeAt(j);
      Double_t btof = cmt - Time0();
      if(true
         && TMath::Abs(btof) < TMath::Abs(btof0)
         ){
        btof0 = btof;
      }
    }
    return btof0;
  }else{
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Ftof0() const
{
  const auto& cl = GetFtof0Cluster();
  if(cl){
    Double_t ftof0 = DBL_MAX;
    for(Int_t j=0, m=cl->ClusterSize(); j<m; ++j){
      Double_t cmt  = cl->CMeanTimeAt(j);
      Double_t ftof = cmt - Time0();
      if(true
         && TMath::Abs(ftof) < TMath::Abs(ftof0)
         ){
        ftof0 = ftof;
      }
    }
    return ftof0;
  }else{
    return TMath::QuietNaN();
  }
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Time0Seg() const
{
  const auto& cl = GetTime0Cluster();
  if(cl) return cl->MeanSeg();
  else   return TMath::QuietNaN();
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Btof0Seg() const
{
  const auto& cl = GetBtof0Cluster();
  if(cl) return cl->MeanSeg();
  else   return TMath::QuietNaN();
}

//_____________________________________________________________________________
Double_t
HodoAnalyzer::Ftof0Seg() const
{
  const auto& cl = GetFtof0Cluster();
  if(cl) return cl->MeanSeg();
  else   return TMath::QuietNaN();
}
