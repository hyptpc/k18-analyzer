// -*- C++ -*-

#ifndef DEBUG_COUNTER_HH
#define DEBUG_COUNTER_HH

#include <map>
#include <string>
#include <vector>

#include <TString.h>

#include <std_ostream.hh>

//_____________________________________________________________________________
namespace debug
{
class ObjectCounter
{
public:
  static ObjectCounter& GetInstance();
  static TString& ClassName();
  ~ObjectCounter();

private:
  ObjectCounter();
  ObjectCounter(const ObjectCounter&);
  ObjectCounter& operator =(const ObjectCounter&);

private:
  typedef std::map<TString,int> ObjectMap;
  typedef ObjectMap::const_iterator ObjectIter;
  ObjectMap m_map;

public:
  void check(const TString& arg="") const;
  void print(const TString& arg="") const;

public:
  static void decrease(const TString& key);
  static void increase(const TString& key);
};

//_____________________________________________________________________
inline ObjectCounter&
ObjectCounter::GetInstance()
{
  static ObjectCounter g_instance;
  return g_instance;
}

//_____________________________________________________________________
inline TString&
ObjectCounter::ClassName()
{
  static TString s_name("ObjectCounter");
  return s_name;
}

//_____________________________________________________________________
inline void
ObjectCounter::decrease(const TString& key)
{
#ifdef MemoryLeak
  --(GetInstance().m_map[key]);
#endif
}

//_____________________________________________________________________
inline void
ObjectCounter::increase(const TString& key)
{
#ifdef MemoryLeak
  ++(GetInstance().m_map[key]);
#endif
}

}

#endif
