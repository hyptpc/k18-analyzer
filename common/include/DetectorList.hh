// DetectorList.h
#ifndef DetectorList_h
#define DetectorList_h 1

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <TString.h>

class DetectorList 
{
public:
  static DetectorList& GetInstance(void);
  static const std::string& ClassName( void );
  ~DetectorList(){};

  bool Initialize(const char* FileName=NULL);
  bool Initialize(const std::string& filename);
  
  unsigned int GetCID(const TString &name) const;
  TString GetName    (const unsigned int &cid) const { return GetString(NAME,cid); }
  TString GetMaterial(const unsigned int &cid) const { return GetString(MATERIAL,cid); }
  TString GetColor   (const unsigned int &cid) const { return GetString(COLOR,cid); }
  bool GetFlag(const unsigned int &cid, int num) const;
  int  GetType(const unsigned int &cid) const { return GetData(TYPE,cid); }
  int  GetNum1(const unsigned int &cid) const { return GetData(NUM1,cid); }
  int  GetNum2(const unsigned int &cid) const { return GetData(NUM2,cid); }
  int  GetNum3(const unsigned int &cid) const { return GetData(NUM3,cid); }
  bool IsChamber (const unsigned int &cid) const { return IsList(CType_Chamber,cid); }
  bool IsCounter (const unsigned int &cid) const { return IsList(CType_Counter,cid); }
  bool IsMTDC    (const unsigned int &cid) const { return IsList(CType_MTDC,cid); }
  bool IsDetector(const unsigned int &cid) const { return !IsList(CType_NODET,cid); }
  
  int GetNlayers (const unsigned int &cid) const { return IsChamber(cid) ? GetNum1(cid) : 0 ; } 
  int GetNwires  (const unsigned int &cid) const { return IsChamber(cid) ? GetNum2(cid) : 0 ; } 
  int GetNsegs   (const unsigned int &cid) const { return !IsChamber(cid) ? GetNum1(cid) : 0 ; } 
  int GetNsensors(const unsigned int &cid) const { return IsCounter(cid) ? GetNum2(cid) : 0 ; } 

  int CounterID(int i) const { return CIDList.at(CType_Counter)[i]; }
  int ChamberID(int i) const { return CIDList.at(CType_Chamber)[i]; }

  void SetMaterial(const unsigned int &cid, std::string &mat)   { StringContainer.at(MATERIAL)[cid]=mat; }
  void SetMaterial(const unsigned int &cid, const TString &mat) { StringContainer.at(MATERIAL)[cid]=mat; }
  
  unsigned int nLists()    const { return CIDContainer.size(); }
  unsigned int nChambers() const { return CIDList.at(CType_Chamber).size(); }
  unsigned int nCounters() const { return CIDList.at(CType_Counter).size(); }
  unsigned int nMTDCs()    const { return CIDList.at(CType_MTDC).size(); }
  unsigned int nNODETs()   const { return CIDList.at(CType_NODET).size(); }

  std::map<TString,unsigned int> const GetList() const { return CIDContainer; }
  
  //  void Dump();

 private:
  DetectorList(){}
  DetectorList(const DetectorList& rhs);
  DetectorList& operator=(const DetectorList& rhs);
  
  void Clear();
  enum gCounterType{ CType_Chamber=0,
		     CType_Counter=1,
		     CType_MTDC=2,
		     CType_NODET=3
  };
  static const int ntype=4;
  enum gStringType{ NAME=0,
		    MATERIAL=1,
		    COLOR=2
  };

  enum gNumberType{ TYPE=0,
		    NUM1=1,
		    NUM2=2,
		    NUM3=3
  };
  
  std::map<TString,unsigned int> CIDContainer;
  std::map<int, std::map<unsigned int,TString> > StringContainer;
  std::map<int, std::map<unsigned int,int> > DataContainer;
  std::map<int, std::map<unsigned int, int> > FlagContainer;
  std::map<int, std::vector<unsigned int> > CIDList;
  
  bool IsList(const int &type,const unsigned int &cid) const;
  int GetData(const int &type,const unsigned int &cid) const;
  TString GetString(const int &type,const unsigned int &cid) const;
};
inline DetectorList&
DetectorList::GetInstance( void )
{
  static DetectorList g_instance;
  return g_instance;
}
inline const std::string&
DetectorList::ClassName( void )
{
  static std::string g_name("DetectorList");
  return g_name;
}

#endif
