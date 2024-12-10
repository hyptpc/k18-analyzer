// -*- C++ -*-

void htof_tree_example(TString file_path="tmp.root") {
  TFile::Open(file_path);
  if (!gFile || !gFile->IsOpen()) {
    std::cerr << "Failed to open " << file_path << std::endl;
    return;
  }
  auto tree = gFile->Get<TTree>("tree");
  if (!tree) {
    std::cerr << "Failed to find tree in " << file_path << std::endl;
    return;
  }

  const Int_t n_seg = 4;
  TH1* h_htof_adc[n_seg];
  TH1* h_htof_tdc[n_seg];
  for (Int_t seg=0; seg<n_seg; ++seg) {
    h_htof_adc[seg] = new TH1D(Form("HTOF_ADC_%d", seg+1),
                               Form("HTOF_ADC_%d", seg+1),
                               0x1000, 0, 0x1000);
    h_htof_tdc[seg] = new TH1D(Form("HTOF_TDC_%d", seg+1),
                               Form("HTOF_TDC_%d", seg+1),
                               10000, 0, 200000);
  }

  TTreeReader reader(tree);
  TTreeReaderValue<Int_t> evnum(reader, "evnum");
  TTreeReaderValue<Int_t> spill(reader, "spill");
  TTreeReaderValue<std::vector<Double_t>> adc(reader, "adc");
  TTreeReaderValue<std::vector<std::vector<Double_t>>> tdc(reader, "tdc");

  while (reader.Next()) {
    for (Int_t seg=0; seg<n_seg; ++seg) {
      auto a = adc->at(seg);
      h_htof_adc[seg]->Fill(a);
      for (const auto& t : tdc->at(seg)) {
        h_htof_tdc[seg]->Fill(t);
      }
    }
  }

  auto c_adc = new TCanvas("c_adc", "c_adc");
  auto c_tdc = new TCanvas("c_tdc", "c_tdc");
  c_adc->Divide(2, 2);
  c_tdc->Divide(2, 2);
  for (Int_t seg=0; seg<n_seg; ++seg) {
    c_adc->cd(seg+1);
    h_htof_adc[seg]->Draw();
    c_tdc->cd(seg+1);
    h_htof_tdc[seg]->Draw();
    h_htof_tdc[seg]->GetXaxis()->SetRangeUser(70000, 110000);
  }
}
