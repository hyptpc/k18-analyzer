// -*- C++ -*-

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

//_____________________________________________________________________________
void
ExampleHodoscope(const TString& root_path="tmp.root")
{
  using seg_t = std::vector<Int_t>;
  using adc_t = std::vector<Double_t>;
  using tdc_t = std::vector<std::vector<Double_t>>;

  TFile::Open(root_path);
  TTreeReader reader("hodo", gFile);
  TTreeReaderValue<UInt_t> evnum(reader, "event_number");
  TTreeReaderValue<seg_t> t0_raw_seg(reader, "t0_raw_seg");
  TTreeReaderValue<adc_t> t0_adc_u(reader, "t0_adc_u");
  TTreeReaderValue<tdc_t> t0_tdc_u(reader, "t0_tdc_u");

  while(reader.Next()){
    Int_t ievent = reader.GetCurrentEntry();
    std::cout << TString('=', 80) << std::endl;
    for(Int_t i=0, n=(*t0_raw_seg).size(); i<n; ++i){
      auto seg = (*t0_raw_seg).at(i);
      auto adc_u = (*t0_adc_u).at(i);
      std::cout << seg << " " << adc_u << " ";
      for(const auto& tdc_u: (*t0_tdc_u).at(i))
        std::cout << tdc_u << " ";
      std::cout << std::endl;
    }
  }
}
