// -*- C++ -*-

#ifndef COUNTER_MAP_MAN_HH
#define COUNTER_MAP_MAN_HH

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <Rtypes.h>
#include <TString.h>


class CounterMapMan
{
public:
  static CounterMapMan& GetInstance();
  static const std::string& ClassName();
  ~CounterMapMan();

private:
  CounterMapMan();
  CounterMapMan(const CounterMapMan&);
  CounterMapMan& operator=(const CounterMapMan&);

public:
  Bool_t Initialize();
  Bool_t Initialize(const char* file_name);
  Bool_t Initialize(const std::string& file_name);
  Bool_t Initialize(const std::string& file_name1, const std::string& file_name2);
  Bool_t ReadFile(const TString& filename);

private:
  static const Int_t MMAXSMP = 20;
  static const Int_t MAXSLOT = 23;
  static const Int_t NORMAL  = 1;
  static const Int_t DRT     = 0;

  std::vector<TString> FileName;
  Int_t NSMP;
  Int_t NCH_SCA;
  Int_t CRATE_TYPE[MMAXSMP][MAXSLOT];

  typedef std::map <UInt_t, UInt_t> fCounterMapContainer;
  typedef std::map <UInt_t, UInt_t> bCounterMapContainer;
  fCounterMapContainer fContainer;
  bCounterMapContainer bContainer;

  typedef std::map <UInt_t, UInt_t> fCrateDefContainer;
  typedef std::map <UInt_t, UInt_t> bCrateDefContainer;
  fCrateDefContainer fCrateDef;
  bCrateDefContainer bCrateDef;

  typedef std::map <UInt_t, TString> nameCNAMapContainer;
  typedef std::map <UInt_t, TString> nameCounterMapContainer;
  nameCNAMapContainer nameCNAContainer;
  nameCounterMapContainer nameCounterContainer;

  Bool_t m_isready;

public:
  Bool_t  IsReady() const { return m_isready; }
  Int_t   nFiles() const { return FileName.size(); }
  TString GetFileName(Int_t i) const { return FileName.at(i); }
  Bool_t  GetInfo(Int_t c, Int_t n, Int_t a, Int_t& cid, Int_t& lay, Int_t& seg, Int_t& at, Int_t& ud);
  Int_t   GetCID(Int_t c, Int_t n, Int_t a);
  TString GetName(Int_t c, Int_t n, Int_t a);

  // for Hodoscope or Cherenkov
  Bool_t  GetCNA(Int_t cid, Int_t seg, Int_t at, Int_t ud, Int_t& c, Int_t& n, Int_t& a);
  TString GetName(Int_t cid, Int_t seg, Int_t at, Int_t ud);
  // for DC
  Bool_t  GetCNA(Int_t cid, Int_t layer, Int_t wire, Int_t at, Int_t ud, Int_t& c, Int_t& n, Int_t& a);
  TString GetName(Int_t cid, Int_t lay, Int_t wire, Int_t at, Int_t ud);

  Int_t   GetCrateNum(Int_t address);
  Int_t   GetSMPAddress(Int_t c);
  Int_t   GetNumSMP() const { return NSMP; }
  Int_t   GetNumScaler() const { return NCH_SCA; }
  Int_t   GetCrateType(Int_t cr, Int_t sl) const { return CRATE_TYPE[cr][sl-1]; } //sl 1 origin
  void    PrintSimpleMap(std::ostream& p_out = std::cout);
  void    PrintMap();
  void    Clear();
};

//______________________________________________________________________________
inline CounterMapMan&
CounterMapMan::GetInstance()
{
  static CounterMapMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
inline const std::string&
CounterMapMan::ClassName()
{
  static std::string g_name("CounterMapMan");
  return g_name;
}

#endif
