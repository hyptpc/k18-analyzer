// -*- C++ -*-

#ifndef DST_HELPER_HH
#define DST_HELPER_HH

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>
#include <type_traits>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <spdlog/spdlog.h>

#include <filesystem_util.hh>

#include "DCAnalyzer.hh"

// if event number mismatch is found, exit process.
#define CheckEventNumberMismatch 1

//______________________________________________________________________________
namespace dst
{
// implemented in each dst
extern std::vector<TString> ArgName;
extern std::vector<TString> TreeName;
extern std::vector<TFile*>  TFileCont;
extern std::vector<TTree*>  TTreeCont;
extern std::vector<TTreeReader*> TTreeReaderCont;
Bool_t InitializeEvent();
Bool_t DstOpen(std::vector<std::string> arg);
Bool_t DstRead();
Bool_t DstRead(Int_t ievent);
Bool_t DstRead(Int_t ievent, DCAnalyzer *DCAna);
Bool_t DstClose();

//______________________________________________________________________________
inline Bool_t
CheckArg(const std::vector<std::string>& arg)
{
  const Int_t n = arg.size();
  Bool_t status = (n == ArgName.size() && n == TreeName.size());

  if(!status){
    std::string usage = "Usage : " + std::string(hddaq::basename(arg[0]));
    for(Int_t i=1; i<ArgName.size(); ++i){
      usage += " " + std::string(ArgName[i]);
    }
    spdlog::info(usage);
    spdlog::error("Argument count mismatch: provided = {}, expected = {} (ArgName.size() = {}, TreeName.size() = {})",
                  n, ArgName.size(), ArgName.size(), TreeName.size());
    return false;
  }

  for(Int_t i=0; i<n; ++i){
    spdlog::info("key = {:<18} arg[{}] = {}", ArgName[i].Data(), i, arg[i]);
  }

  TFileCont.resize(n); TTreeCont.resize(n); TTreeReaderCont.resize(n);
  return (ArgName.size()==TreeName.size());
}

//______________________________________________________________________________
inline Bool_t
OpenFile(TFile*& file, const TString& name)
{
  file = new TFile(name);
  if(!file || file->IsZombie()){
    spdlog::error("failed to open TFile : {}", name.Data());
    return false;
  }
  return true;
}

//______________________________________________________________________________
inline Bool_t
OpenTree(TFile* file, TTree*& tree, const TString& name)
{
  if(!file || !file->IsOpen()) return false;
  tree = (TTree*)file->Get(name);
  if(!tree){
    spdlog::error("failed to open TTree : {}", name.Data());
    return false;
  }
  return true;
}

//______________________________________________________________________________
inline Bool_t
CheckEntries(const std::vector<TTree*>& TTreeCont)
{
  const Int_t n = TTreeCont.size();
  Bool_t status = true;
  std::vector<Int_t> entries(n, -1);
  for(Int_t i=0; i<n; ++i){
    if(!TTreeCont[i]) continue;
    entries[i] = TTreeCont[i]->GetEntries();
    if(i>0 && entries[i]!=entries[i-1] && entries[i-1]!=-1){
      status = false;
    }
  }
  if(!status){
#if CheckEventNumberMismatch
    spdlog::error("Entries Mismatch");
#else
    spdlog::warn("Entries Mismatch");
#endif
    for(Int_t i=0; i<n; ++i){
      if(!TTreeCont[i]) continue;
#if CheckEventNumberMismatch
      spdlog::error("   {:8} {}", TTreeCont[i]->GetName(), entries[i]);
#else
      spdlog::warn("   {:8} {}", TTreeCont[i]->GetName(), entries[i]);
#endif
    }
  }
#if CheckEventNumberMismatch
  return status;
#else
  return true;
#endif
}

//______________________________________________________________________________
inline Int_t
GetEntries(const std::vector<TTree*>& TTreeCont)
{
  std::vector<Int_t> nevent;
  for(Int_t i=0, n=TTreeCont.size(); i<n; ++i){
    if(TTreeCont[i]){
#if CheckEventNumberMismatch
      return TTreeCont[i]->GetEntries();
#else
      nevent.push_back(TTreeCont[i]->GetEntries());
#endif
    }
  }
#if CheckEventNumberMismatch
  return 0;
#else
  return *std::min_element(nevent.begin(), nevent.end());
#endif
}

//______________________________________________________________________________
inline Bool_t
GetEntry(Int_t ievent)
{
  for(Int_t i=0, n=TTreeCont.size(); i<n; ++i){
    if(TTreeCont[i]){
      TTreeCont[i]->GetEntry(ievent);
      if(TTreeReaderCont[i]){
        TTreeReaderCont[i]->SetEntry(ievent);
      }
    }
  }
  return true;
}

//______________________________________________________________________________
template <class... Vecs>
inline void clear_all(Vecs&... vecs) {
  (vecs.clear(), ...);
}

//______________________________________________________________________________
template <class SizeT, class... Vecs>
inline void resize_all(SizeT n, Vecs&... vecs) {
  (vecs.resize(n), ...);
}

//______________________________________________________________________________
inline Bool_t
SetupReader(Int_t index, const std::string& label)
{
  if (!TFileCont[index] || TFileCont[index]->IsZombie()) {
    spdlog::error("Failed to open TFileCont[{}] (Null or Zombie).", label);
    return false;
  }
  TTreeReaderCont[index] = new TTreeReader(TreeName[index], TFileCont[index]);
  return true;
}

//______________________________________________________________________________
template <class T>
inline void SetBranch(TTreeReader* reader, const Char_t* name, T*& ptr)
{
  using ValueType = typename std::decay<decltype(**ptr)>::type;
  ptr = new TTreeReaderValue<ValueType>(*reader, name);
}

} // namespace dst

#endif
