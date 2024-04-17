// -*- C++ -*-

#ifndef CONF_MAN_HH
#define CONF_MAN_HH

#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include <std_ostream.hh>

class VEvent;

namespace
{
const std::string& kConfFile("CONF");
const std::string& kInFile("IN");
const std::string& kOutFile("OUT");
const std::string& kCDCFile("CDC");
}
//______________________________________________________________________________
class ConfMan
{
public:
  static ConfMan& GetInstance();
  ~ConfMan();

private:
  ConfMan();
  ConfMan(const ConfMan&);
  ConfMan& operator=(const ConfMan&);

private:
  typedef std::map<std::string, std::string> StrList;
  typedef std::map<std::string, double>      DoubleList;
  typedef std::map<std::string, int>         IntList;
  typedef std::map<std::string, bool>        BoolList;
  typedef StrList::const_iterator    StrIterator;
  typedef DoubleList::const_iterator DoubleIterator;
  typedef IntList::const_iterator    IntIterator;
  typedef BoolList::const_iterator   BoolIterator;

  bool        m_is_ready;
  StrList     m_file;
  StrList     m_string;
  DoubleList  m_double;
  IntList     m_int;
  BoolList    m_bool;

  int         m_runnum;
  int         m_maxloop;
  int         m_skip;
  double      m_resolution_cdc;
  bool        m_is_mc;
  bool        m_is_cosmic;

  std::vector<std::string> m_keys_used;

public:
  VEvent* EventAllocator();
  template <typename T>
  static const T& Get(const std::string& key) { return T(); }
  bool    Initialize();
  bool    Initialize(int argc, char **argv);
  bool    Initialize(const std::string& file_name, const int &runnum=-1);
  bool    InitializeHistograms();
  bool    InitializeParameterFiles();
  bool    InitializeUnpacker();
  bool    ParseCommand(int argc, char **argv);
  bool    IsReady() const { return m_is_ready; }
  bool    IsCosmic() const { return m_is_cosmic; }
  bool    IsMC() const { return m_is_mc; }
  bool    Finalize();
  bool    FinalizeProcess();
  bool    BeginRun();
  bool    BeginRunProcess();
  // Initialize Parameter
  template <typename T>
  bool    InitializeParameter();
  template <typename T>
  bool    InitializeParameter(const std::string& key);
  template <typename T>
  bool    InitializeParameter(const std::string& key1,
                              const std::string& key2);
  template <typename T>
  bool    InitializeParameter(const std::string& key1,
                              const std::string& key2,
                              const std::string& key3);

  int     get_run_number()     const { return m_runnum; }
  int     get_maxloop()        const { return m_maxloop; }
  int     get_skip()           const { return m_skip; }
  double  get_resolution_cdc() const { return m_resolution_cdc; }

  std::string get_conffilename()  { return m_file[kConfFile]; }
  std::string get_infilename()    { return m_file[kInFile]; }
  std::string get_outfilename()   { return m_file[kOutFile]; }
  std::string get_cdcfilename()   { return m_file[kCDCFile]; }

  void WriteParams();

private:
  std::string FilePath(const std::string& src) const;
  bool        ShowResult(bool s, const std::string& name) const;
};

//______________________________________________________________________________
inline ConfMan&
ConfMan::GetInstance()
{
  static ConfMan g_instance;
  return g_instance;
}

//______________________________________________________________________________
template <>
inline const std::string&
ConfMan::Get<std::string>(const std::string& key)
{
  return GetInstance().m_string[key];
}

//______________________________________________________________________________
template <>
inline const double&
ConfMan::Get<double>(const std::string& key)
{
  return GetInstance().m_double[key];
}

//______________________________________________________________________________
template <>
inline const int&
ConfMan::Get<int>(const std::string& key)
{
  return GetInstance().m_int[key];
}

//______________________________________________________________________________
template <>
inline const bool&
ConfMan::Get<bool>(const std::string& key)
{
  return GetInstance().m_bool[key];
}

//______________________________________________________________________________
inline bool
ConfMan::ShowResult(bool s, const std::string& name) const
{
  if(s)
    hddaq::cout << std::setw(20) << std::left
		<< " ["+name+"]"
		<< "-> Initialized" << std::endl;
  else
    hddaq::cout << std::setw(20) << std::left
		<< " ["+name+"]"
		<< "-> Failed" << std::endl;
  return s;
}

//______________________________________________________________________________
template <typename T>
inline bool
ConfMan::InitializeParameter()
{
  return
    ShowResult(T::GetInstance().Initialize(),
               T::ClassName());
}

//______________________________________________________________________________
template <typename T>
inline bool
ConfMan::InitializeParameter(const std::string& key)
{
  m_keys_used.push_back(key);
  return
    ShowResult(T::GetInstance().Initialize(m_file[key]),
               T::ClassName());
}

//______________________________________________________________________________
template <typename T>
inline bool
ConfMan::InitializeParameter(const std::string& key1,
                             const std::string& key2)
{
  m_keys_used.push_back(key1);
  m_keys_used.push_back(key2);
  return
    ShowResult(T::GetInstance().Initialize(m_file[key1],
                                           m_file[key2]),
               T::ClassName());
}

//______________________________________________________________________________
template <typename T>
inline bool
ConfMan::InitializeParameter(const std::string& key1,
                             const std::string& key2,
                             const std::string& key3)
{
  m_keys_used.push_back(key1);
  m_keys_used.push_back(key2);
  m_keys_used.push_back(key3);
  return
    ShowResult(T::GetInstance().Initialize(m_file[key1],
                                           m_file[key2],
                                           m_file[key3]),
               T::ClassName());
}

#endif
