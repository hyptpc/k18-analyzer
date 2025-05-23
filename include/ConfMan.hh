// -*- C++ -*-

#ifndef CONF_MAN_HH
#define CONF_MAN_HH

#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include <TString.h>

#include <spdlog/spdlog.h>

class TNamed;

//_____________________________________________________________________________
class ConfMan
{
public:
  static const TString& ClassName();
  static ConfMan&       GetInstance();
  ~ConfMan();

private:
  ConfMan();
  ConfMan(const ConfMan&);
  ConfMan& operator=(const ConfMan&);

private:
  using StrList        = std::map<TString, TString>;
  using DoubleList     = std::map<TString, Double_t>;
  using IntList        = std::map<TString, Int_t>;
  using BoolList       = std::map<TString, Bool_t>;
  using StrIterator    = StrList::const_iterator;
  using DoubleIterator = DoubleList::const_iterator;
  using IntIterator    = IntList::const_iterator;
  using BoolIterator   = BoolList::const_iterator;
  Bool_t     m_is_ready;
  StrList    m_file;
  StrList    m_string;
  DoubleList m_double;
  IntList    m_int;
  BoolList   m_bool;
  TString    m_buf;
  TNamed*    m_object;

public:
  Bool_t  Finalize();
  Bool_t  FinalizeProcess();
  Bool_t  Initialize();
  Bool_t  Initialize(const TString& file_name);
  Bool_t  InitializeHistograms();
  Bool_t  InitializeParameterFiles();
  Bool_t  InitializeUnpacker();
  Bool_t  IsReady() const { return m_is_ready; }
  void    WriteParameters();
  // templates
  template <typename T>
  static const T& Get(const TString& key) { return T(); }
  template <typename T>
  Bool_t InitializeParameter();
  template <typename T>
  Bool_t InitializeParameter(const TString& key);
  template <typename T>
  Bool_t InitializeParameter(const TString& key1,
                             const TString& key2);

private:
  TString FilePath(const TString& src) const;
  Bool_t  ShowResult(Bool_t s, const TString& name) const;
};

//_____________________________________________________________________________
inline const TString&
ConfMan::ClassName()
{
  static TString s_name("ConfMan");
  return s_name;
}

//_____________________________________________________________________________
inline ConfMan&
ConfMan::GetInstance()
{
  static ConfMan s_instance;
  return s_instance;
}

//_____________________________________________________________________________
template <>
inline const TString&
ConfMan::Get<TString>(const TString& key)
{
  return GetInstance().m_string[key];
}

//_____________________________________________________________________________
template <>
inline const Double_t&
ConfMan::Get<Double_t>(const TString& key)
{
  return GetInstance().m_double[key];
}

//_____________________________________________________________________________
template <>
inline const Int_t&
ConfMan::Get<Int_t>(const TString& key)
{
  return GetInstance().m_int[key];
}

//_____________________________________________________________________________
template <>
inline const Bool_t&
ConfMan::Get<Bool_t>(const TString& key)
{
  return GetInstance().m_bool[key];
}

//_____________________________________________________________________________
inline Bool_t
ConfMan::ShowResult(Bool_t s, const TString& name) const
{
  if(s)
    spdlog::info(" [{}] -> initialized", name.Data());
  else
    spdlog::error(" [{}] -> failed", name.Data());
  return s;
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
ConfMan::InitializeParameter()
{
  return
    ShowResult(T::GetInstance().Initialize(),
               T::ClassName());
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
ConfMan::InitializeParameter(const TString& key)
{
  return
    ShowResult(T::GetInstance().Initialize(m_file[key]),
               T::ClassName());
}

//_____________________________________________________________________________
template <typename T>
inline Bool_t
ConfMan::InitializeParameter(const TString& key1,
                             const TString& key2)
{
  return
    ShowResult(T::GetInstance().Initialize(m_file[key1], m_file[key2]),
               T::ClassName());
}

#endif
