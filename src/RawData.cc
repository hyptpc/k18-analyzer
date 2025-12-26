// -*- C++ -*-

#include "RawData.hh"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <std_ostream.hh>
#include <UnpackerConfig.hh>
#include <UnpackerManager.hh>
#include <UnpackerXMLReadDigit.hh>

#include "ConfMan.hh"
#include "DCRawHit.hh"
#include "DebugCounter.hh"
#include "DeleteUtility.hh"
#include "DetectorID.hh"
#include "Exception.hh"
#include "FuncName.hh"
#include "HodoRawHit.hh"
#include "TPCRawHit.hh"
#include "MathTools.hh"
#include "TPCPadHelper.hh"
#include "UserParamMan.hh"

namespace
{
using namespace hddaq::unpacker;
const auto& gUnpacker     = GUnpacker::get_instance();
const auto& gUser         = UserParamMan::GetInstance();

  // for TPC baseline corrected HC
  const TString nameCorTPC = "CorTPC";

  ///// for CorrectBaselineTPC()
  TH1D* h_baseline = nullptr;
  Double_t f_baseline(Double_t* x, Double_t* par)
  {
    // par[0]: adc offset, par[1]: scale, par[2]: time offset
    if(!h_baseline){
      throw Exception("something is wrong in [RawData::CorrectBaselineTPC()]");
      // return TMath::QuietNaN();
    }
    Int_t floor = TMath::FloorNint(par[2]);
    Double_t frac = par[2] - floor;
    Int_t bin_left = h_baseline->GetXaxis()->FindBin(x[0] + floor);
    Double_t val_left = h_baseline->GetBinContent(bin_left);
    Double_t val_right = h_baseline->GetBinContent(bin_left+1);
    return
      par[0] + par[1]*((1-frac)*val_left + frac*val_right);
  }
}

#define DebugEvDisp    0


//_____________________________________________________________________________
RawData::RawData()
  : m_is_decoded(),
    m_hodo_raw_hit_collection(),
    m_dc_raw_hit_collection(),
    m_tpc_raw_hit_collection(),
    m_baseline()
{
  debug::ObjectCounter::increase(ClassName());
}

//_____________________________________________________________________________
RawData::~RawData()
{
  Clear();
  debug::ObjectCounter::decrease(ClassName());
}

//_____________________________________________________________________________
void
RawData::Clear(const TString& name)
{
  if(name.IsNull()){
    for(auto& elem: m_hodo_raw_hit_collection)
      del::ClearContainer(elem.second);
    for(auto& elem: m_dc_raw_hit_collection)
      del::ClearContainer(elem.second);
    for(auto& elem: m_tpc_raw_hit_collection)
      del::ClearContainer(elem.second);
    m_hodo_raw_hit_collection.clear();
    m_dc_raw_hit_collection.clear();
    m_tpc_raw_hit_collection.clear();
  }else{
    del::ClearContainer(m_hodo_raw_hit_collection[name]);
    del::ClearContainer(m_dc_raw_hit_collection[name]);
    del::ClearContainer(m_tpc_raw_hit_collection[name]);
  }
}

//_____________________________________________________________________________
Bool_t
RawData::DecodeHits(const TString& name)
{
  static const auto& digit_info = GConfig::get_instance().get_digit_info();

  if(m_is_decoded[name]){
    hddaq::cerr << FUNC_NAME << " " << name << " is already decoded"
                << std::endl;
    return false;
  }

  if(name.IsNull()){
    Bool_t ret = true;
    for(const auto& n: digit_info.get_name_list()){
      if(!n.empty())
        ret &= DecodeHits(n);
    }
    return ret;
  }

  auto id = digit_info.get_device_id(name.Data());
  const TString type = digit_info.get_device_type(id);

#if 0
  hddaq::cout << FUNC_NAME << std::endl
              << id << " " << name << " " << type << std::endl;
#endif

  if(type.IsNull())
    return false;

  Clear(name);

  if(name == "VFT")
    return true; // ignore

  Bool_t is_hodo  = type.Contains("Hodo", TString::kIgnoreCase);
  Bool_t is_fiber = type.Contains("Fiber", TString::kIgnoreCase);
  Bool_t is_dc    = type.Contains("DC", TString::kIgnoreCase);
  Bool_t is_tpc   = type.BeginsWith("TPC", TString::kIgnoreCase);
  Bool_t is_dummy = type.Contains("dummy", TString::kIgnoreCase);

  if(is_dummy) return false;
  if(!is_hodo && !is_fiber && !is_dc && !is_tpc) return false;

  for(Int_t plane=0, n_plane=gUnpacker.get_n_plane(id);
      plane<n_plane; ++plane){
    for(Int_t seg=0, n_seg=gUnpacker.get_n_segment(id, plane);
        seg<n_seg; ++seg){
      for(Int_t ch=0, n_ch=gUnpacker.get_n_ch(id, plane, seg);
          ch<n_ch; ++ch){
        for(Int_t data=0, n_data=gUnpacker.get_n_data(id, plane, seg, ch);
            data<n_data; ++data){
          for(Int_t i=0, n=gUnpacker.get_entries(id, plane, seg, ch, data);
              i<n; ++i){
            UInt_t val = gUnpacker.get(id, plane, seg, ch, data, i);
            if(is_hodo)  AddHodoRawHit(name, plane, seg, ch, data, val);
            if(is_fiber) AddFiberRawHit(name, plane, seg, ch, data, val);
            if(is_dc)    AddDCRawHit(name, plane, seg, ch, data, val);
	    if(is_tpc){
	      AddTPCRawHit(name, plane, seg, ch, data, val, nullptr);
	      //Corrected TPC HC before baseline correction
	      AddTPCRawHit(nameCorTPC, plane, seg, ch, data, val, nullptr);
	    }
          }
        }
      }
    }
  }

#if 0
  hddaq::cout << FUNC_NAME << std::endl
              << id << " " << name << " " << type <<  " decoded" << std::endl;
#endif

  // For AC-SUM
  if(name == "AC"){
    Double_t suma = 0;
    for(const auto& hit: GetHodoRawHC(name)){
      if(hit->SegmentId() != 0) suma += hit->GetAdc();
    }
    AddHodoRawHit(name, 0, 0, 0, 0, suma);
  }

  m_is_decoded[name] = true;
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::DecodeTPCHits()
{

  static const Bool_t BaselineCorrectionTPC
    = (gUser.GetParameter("BaselineCorrectionTPC") == 1);
  if(m_is_decoded["TPC"]){
    hddaq::cout << FUNC_NAME << " " << "already decoded!" << std::endl;
    return false;
  }

  if(!DecodeHits("TPC"))
    return false;

  bool is_baselinecorrected = true;
  if(BaselineCorrectionTPC) is_baselinecorrected = CorrectBaselineTPC();

  /*
   * if correction is skipped or null baseline is found,
   * The TPCCorHC is deeply copied from the TPCRawHC.
   * (give up the baseline correction)
   * So, in any case, the TPCCorHC will be used in TPCAnalyzer not the TPCRawHC.
   */

  if(!m_baseline&&is_baselinecorrected){
    //del::ClearContainerAll(m_TPCCorHC);
    Clear(nameCorTPC);
    auto& cont = m_tpc_raw_hit_collection["TPC"];
    for(Int_t i=0, n=cont.size(); i<n; ++i){
      TPCRawHit *hit = cont[i];
      double raw_rms = hit->RawRMS();
      for(const auto& adc: hit->Fadc()){
	double datatype = 0; //not used for TPC
	AddTPCRawHit(nameCorTPC,
		     hit->LayerId(),
		     0,
		     hit->RowId(),
		     datatype, adc,
		     nullptr, raw_rms);
      }
    }
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::CorrectBaselineTPC()
{
  static const Double_t MinRms = gUser.GetParameter("MinBaseRmsTPC");
  static const Int_t MinTimeBucket = gUser.GetParameter("TimeBucketTPC", 0);
  static const Int_t MaxTimeBucket = gUser.GetParameter("TimeBucketTPC", 1);
  static const Int_t NumOfTimeBucket = gUser.GetParameter("NumOfTimeBucket");

  if(!m_is_decoded["TPC"]){
    hddaq::cerr << FUNC_NAME << " DecodeTPCHits() must be done!" << std::endl;
    return false;
  }

  //del::ClearContainerAll(m_TPCCorHC);
  Clear(nameCorTPC);

  TH1D h1("baseline", "Baseline", NumOfTimeBucket, 0, NumOfTimeBucket);
  h_baseline = &h1;

#if DebugEvDisp
  gStyle->SetOptStat(0);
  // gStyle->SetOptStat(1110);
  // gStyle->SetOptFit(1);
  static TCanvas c1("c"+FUNC_NAME, FUNC_NAME, 1200, 900);
  c1.cd();
  TH1D h2(FUNC_NAME+"-h2", "Corrected FADC",
	  NumOfTimeBucket, 0, NumOfTimeBucket);
#endif

  m_baseline = nullptr;
  Double_t min_ref = 1e5;
  auto& cont = m_tpc_raw_hit_collection["TPC"];
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    TPCRawHit *hit = cont[i];
    if(hit->FadcSize() != NumOfTimeBucket){
      hddaq::cout << FUNC_NAME << " " << "FadcSize!=NumOfTimeBucket!" << hit->FadcSize()<<" "<<NumOfTimeBucket<<std::endl;
      //hit->Print();
      continue;
    }
    ///// Minimum RMS method (unused)
    // auto ref = hit->RMS(MinTimeBucket, MaxTimeBucket);
    // if(hit->RMS(MinTimeBucket, MaxTimeBucket) < 20.) continue;

    ///// Minimum Amplitude method
    auto ref = hit->MaxAdc(MinTimeBucket, MaxTimeBucket)
      - hit->Mean(MinTimeBucket, MaxTimeBucket);
    double rms = hit -> RMS(MinTimeBucket, MaxTimeBucket);
    if(ref < min_ref && rms > MinRms){
      min_ref = ref;
      m_baseline = hit;
    }
  }
  if(!m_baseline || !h_baseline){
    hddaq::cerr << FUNC_NAME << " no reference was found." << std::endl;
    return false;
  }

  const auto& base_fadc = m_baseline->Fadc();
  const Double_t base_ped = m_baseline->Mean(0, MinTimeBucket);
  for(Int_t i=0; i<NumOfTimeBucket; ++i){
    h_baseline->SetBinContent(i+1, base_fadc.at(i) - base_ped);
  }
  h_baseline->SetTitle(Form("Baseline Layer#%d Row#%d",
			    m_baseline->LayerId(),
			    m_baseline->RowId()));

#if DebugEvDisp
  {
    m_baseline->Print();
  }
#endif

  for(Int_t i=0, n=cont.size(); i<n; ++i){
    TPCRawHit *hit = cont[i];
    double raw_rms = hit->RMS(0, NumOfTimeBucket);
    TH1D h_fadc("h_fadc", Form("FADC Layer#%d Row#%d;sample# ;ADC ch",
			       hit->LayerId(), hit->RowId()),
		NumOfTimeBucket, 0, NumOfTimeBucket);
    h_fadc.SetLineWidth(2);
    const auto& fadc = hit->Fadc();
    for(Int_t i=0, n=fadc.size(); i<n; ++i){
      h_fadc.SetBinContent(i+1, fadc.at(i));
    }
    TF1 f1("f1", f_baseline, 0, NumOfTimeBucket, 3);
    //f1.SetParameter(0, hit->Mean(0, MinTimeBucket));
    f1.SetParameter(0, hit->Mean(MaxTimeBucket, NumOfTimeBucket));
    f1.SetParameter(1, 1.);
    f1.SetParameter(2, 0.);
    f1.SetParLimits(0, 0, 4000.);
    f1.SetParLimits(1, -5, 5.);
    f1.SetParLimits(2, -10, 10);
    //h_fadc.Fit("f1", "Q", "", 0, MinTimeBucket);
    h_fadc.Fit("f1", "Q", "", MaxTimeBucket, NumOfTimeBucket);
    Double_t max_cadc = -1e10;
    for(Int_t i=0, n=fadc.size(); i<n; ++i){
      Double_t cadc = fadc.at(i) - f1.Eval(i);
      double datatype = 0; //not used for TPC
	AddTPCRawHit(nameCorTPC,
		     hit->LayerId(),
		     0,
		     hit->RowId(),
		     datatype, cadc,
		     f1.GetParameters(),
		     raw_rms);
      if(cadc > max_cadc && i > 25){
	max_cadc = cadc;
      }
    }
#if DebugEvDisp
    if(max_cadc > 100 && max_cadc<1000)
      {
	h2.Reset();
	h2.SetLineWidth(2);
	Double_t max_adc = -1e10;
	for(Int_t i=0, n=fadc.size(); i<n; ++i){
	  h2.SetBinContent(i+1, fadc.at(i) - f1.Eval(i));
	  if(fadc.at(i) > max_adc) max_adc = fadc.at(i);
	}
	h_fadc.Draw();
	// h_fadc.SetMinimum(min_adc - 100);
	h_fadc.SetMinimum(-100);
	h_fadc.SetMaximum(max_adc + 100);
	// h_fadc.SetMaximum(1000);
	// h_fadc.SetMaximum(2000);
	//TF1 f2("f2", f_baseline, MinTimeBucket, NumOfTimeBucket, 3);
	TF1 f2("f2", f_baseline, 0, NumOfTimeBucket, 3);
	f2.SetParameters(f1.GetParameters());
	f2.Draw("same");
	h2.SetLineColor(kGreen+1);
	h2.SetLineWidth(2);
	h2.Draw("same");
	// h2.SetMinimum(-500);
	// h2.SetMaximum( 500);
	gPad->Modified();
	gPad->Update();
	c1.Print("c1.pdf");
	getchar();
      }
#endif
  }

  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddHodoRawHit(const TString& name, Int_t plane, Int_t seg,
		       Int_t ch, Int_t data, Double_t val)
{
  auto& cont = m_hodo_raw_hit_collection[name];
  HodoRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    HodoRawHit* q = cont[i];
    if(true
       && q->DetectorName() == name
       && q->PlaneId() == plane
       && q->SegmentId() == seg){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit(name, plane, seg);
    cont.push_back(p);
  }

  if(data == gUnpacker.get_data_id(name, "adc")){
    p->SetAdc(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "leading")){
    p->SetTdcLeading(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "trailing")){
    p->SetTdcTrailing(ch, val);
  // }else if(data == gUnpacker.get_data_id(name, "cstop")){
  //   ;
  }else if(data == gUnpacker.get_data_id(name, "overflow")){
    p->SetTdcOverflow(ch, val);
  }
  else{
    hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
                << " Detector   = " << name  << std::endl
                << " Plane      = " << plane << std::endl
                << " Segment    = " << seg   << std::endl
                << " Channel    = " << ch    << std::endl
                << " Data       = " << data  << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddFiberRawHit(const TString& name, Int_t plane, Int_t seg,
		       Int_t ch, Int_t data, Double_t val)
{
  auto& cont = m_hodo_raw_hit_collection[name];
  HodoRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    HodoRawHit* q = cont[i];
    if(true
       && q->DetectorName() == name
       && q->PlaneId() == plane
       && q->SegmentId() == seg){
      p=q; break;
    }
  }
  if(!p){
    p = new HodoRawHit(name, plane, seg);
    cont.push_back(p);
  }

  if(data == gUnpacker.get_data_id(name, "leading")){
    p->SetTdcLeading(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "trailing")){
    p->SetTdcTrailing(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "highgain")){
    p->SetAdcHigh(ch, val);
  }else if(data == gUnpacker.get_data_id(name, "lowgain")){
    p->SetAdcLow(ch, val);
  }else{
    hddaq::cerr << FUNC_NAME << " wrong data type " << std::endl
                << " Detector   = " << name  << std::endl
                << " Plane      = " << plane << std::endl
                << " Segment    = " << seg   << std::endl
                << " Channel    = " << ch    << std::endl
                << " Data       = " << data  << std::endl;
    return false;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddDCRawHit(const TString& name, Int_t plane, Int_t seg,
                     Int_t ch, Int_t data, Double_t val)
{
  auto& cont = m_dc_raw_hit_collection[name];
  Int_t wire = ch;
  DCRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    DCRawHit* q = cont[i];
    if(true
       && q->DetectorName() == name
       && q->PlaneId() == plane
       && q->WireId() == wire){
      p=q; break;
    }
  }
  if(!p){
    p = new DCRawHit(name, plane, wire);
    cont.push_back(p);
  }

  if(data == gUnpacker.get_data_id(name, "leading")){
      p->SetTdc(val);
  }else if(data == gUnpacker.get_data_id(name, "trailing")){
      p->SetTrailing(val);
  }else if(data == gUnpacker.get_data_id(name, "overflow")){
    p->SetTdcOverflow(val);
  }else{
    hddaq::cerr << FUNC_NAME << " unknown data type " << std::endl
                << " Detector = " << name  << std::endl
		<< " PlaneId  = " << plane << std::endl
		<< " WireId   = " << wire  << std::endl
		<< " DataType = " << data  << std::endl
		<< " Value    = " << val   << std::endl;
  }
  return true;
}

//_____________________________________________________________________________
Bool_t
RawData::AddTPCRawHit(const TString& name,
		      Int_t layer/* plane */,
		      Int_t seg/*not used for TPC*/,
		      Int_t row/*ch*/,
		      Int_t data/* not used for TPC */,
		      Double_t val/*adc*/,
		      Double_t* pars,
		      Double_t raw_rms)
{
  if(data!=0) return false;
  auto& cont = m_tpc_raw_hit_collection[name];
  TPCRawHit* p = nullptr;
  for(Int_t i=0, n=cont.size(); i<n; ++i){
    TPCRawHit* q = cont[i];
    if(q->LayerId() == layer && q->RowId() == row){
      p=q; break;
    }
  }
  if(!p){
    p = new TPCRawHit(layer, row, pars);
    p->SetRawRMS(raw_rms);
    cont.push_back(p);
  }
  p->AddFadc(val);
  p->SetDetectorName(name);
  return true;
}

//_____________________________________________________________________________
const HodoRHC&
RawData::GetHodoRawHitContainer(const TString& name) const
{
  auto itr = m_hodo_raw_hit_collection.find(name);
  if(itr == m_hodo_raw_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static HodoRHC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
const HodoRHC&
RawData::GetHodoRawHitContainer(Int_t det_id) const
{
  for(const auto& [key, value]: m_hodo_raw_hit_collection){
    if(gUnpacker.get_device_id(key) == det_id){
      return value;
    }
  }
  static HodoRHC null_container;
  return null_container;
}

//_____________________________________________________________________________
const DCRHC&
RawData::GetDCRawHitContainer(const TString& name) const
{
  auto itr = m_dc_raw_hit_collection.find(name);
  if(itr == m_dc_raw_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static DCRHC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}

//_____________________________________________________________________________
const DCRHC&
RawData::GetDCRawHitContainer(Int_t det_id) const
{
  for(const auto& [key, value]: m_dc_raw_hit_collection){
    if(gUnpacker.get_device_id(key) == det_id){
      return value;
    }
  }
  static DCRHC null_container;
  return null_container;
}

//_____________________________________________________________________________
const DCRHC&
RawData::GetDCRawHitContainer(Int_t det_id, Int_t plane) const
{
  const auto& cont_all = GetDCRawHitContainer(det_id);
  static DCRHC cont_plane;
  cont_plane.clear();
  for(const auto& hit: cont_all){
    if(hit->PlaneId() == plane){
      cont_plane.push_back(hit);
    }
  }
  return cont_plane;
}

//_____________________________________________________________________________
const TPCRHC&
RawData::GetTPCRawHitContainer(const TString& name) const
{
  auto itr = m_tpc_raw_hit_collection.find(name);
  if(itr == m_tpc_raw_hit_collection.end()){
    // throw Exception(FUNC_NAME + " No such detector: " + name);
    static TPCRHC null_container;
    return null_container;
  }else{
    return itr->second;
  }
}
/*
//_____________________________________________________________________________
const TPCRHC&
RawData::GetTPCRawHitContainer(Int_t det_id) const
{
  for(const auto& [key, value]: m_tpc_raw_hit_collection){
    if(gUnpacker.get_device_id(key) == det_id){
      return value;
    }
  }
  static TPCRHC null_container;
  return null_container;
}
*/
//_____________________________________________________________________________
const TPCRHC&
RawData::GetTPCRawHitContainer(Int_t layer) const
{
  TString name = "TPC";
  const auto& cont_all = GetTPCRawHitContainer(name);
  static TPCRHC cont_layer;
  cont_layer.clear();
  for(const auto& hit: cont_all){
    if(hit->LayerId() == layer){
      cont_layer.push_back(hit);
    }
  }
  return cont_layer;
}

//_____________________________________________________________________________
const TPCRHC&
RawData::GetTPCCorHitContainer(Int_t layer) const
{
  TString name = "CorTPC";
  const auto& cont_all = GetTPCRawHitContainer(name);
  static TPCRHC cont_layer;
  cont_layer.clear();
  for(const auto& hit: cont_all){
    if(hit->LayerId() == layer){
      cont_layer.push_back(hit);
    }
  }
  return cont_layer;
}

//_____________________________________________________________________________
void
RawData::Print(Option_t* arg) const
{
  for(const auto& p: m_hodo_raw_hit_collection){
    for(const auto& hit: p.second){
      if(!arg || hit->DetectorName().EqualTo(arg))
        hit->Print();
    }
  }
  for(const auto& p: m_dc_raw_hit_collection){
    for(const auto& hit: p.second){
      if(!arg || hit->DetectorName().EqualTo(arg))
        hit->Print();
    }
  }
  for(const auto& p: m_tpc_raw_hit_collection){
    for(const auto& hit: p.second){
      if(!arg || hit->DetectorName().EqualTo(arg))
        hit->Print();
    }
  }

}
