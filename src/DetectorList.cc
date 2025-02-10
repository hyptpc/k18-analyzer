#include "DetectorList.hh"
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <std_ostream.hh>

bool DetectorList::Initialize(const std::string &filename){
  return Initialize(filename.data());
}
bool DetectorList::Initialize(const char* FileName){
  static const TString func_name = "DetectorList::Initialize";
  std::cout << "[" << func_name << "] Initialization start ..."<<std::endl;

  Clear();
  std::ifstream ifs( FileName );
  if( !ifs.is_open() ){
    hddaq::cerr << "#E " << func_name << " file open fail : "
		<< FileName << std::endl;
    return false;
  }
  std::string line;
  std::string tmpname,matname,color;
  int cid,type,num1,num2,num3,flag1,flag2;  
  while( ifs.good() && std::getline( ifs, line ) ){
    if( line[0]=='#' || line.empty() ) continue;
    if( TString(line).Contains("cid",TString::ECaseCompare::kIgnoreCase) ) continue;
    std::istringstream input_line( line );
    std::cout<<line<<std::endl;
    if( input_line 
	>> cid >> tmpname >> type
	>> num1 >> num2 >> num3
	>> flag1 >> flag2
	>> matname >> color )
      {
	CIDContainer[tmpname]=cid;
	StringContainer[NAME][cid]=tmpname;
	StringContainer[MATERIAL][cid]=matname;
	StringContainer[COLOR][cid]=color;
	FlagContainer[0][cid]=flag1;
	FlagContainer[1][cid]=flag2;
	DataContainer[TYPE][cid]=type;
	DataContainer[NUM1][cid]=num1;
	DataContainer[NUM2][cid]=num2;
	DataContainer[NUM3][cid]=num3;
	if(type>=0&&type<ntype)
	  CIDList[type].push_back(cid);
      }
  }
  return true;
}

void DetectorList::Clear()
{
  CIDContainer.clear();
  for(int i=0;i<3;i++) StringContainer[i].clear();
  for(int i=0;i<3;i++) DataContainer[i].clear();
  for(int i=0;i<3;i++) FlagContainer[i].clear();
  for(int i=0;i<ntype;i++)   CIDList[i].clear();
}

bool DetectorList::IsList(const int &type,const unsigned int &cid) const 
{ 
  if(type>=0&&type<ntype){
    if( find( CIDList.at(type).begin(), 
	      CIDList.at(type).end(), 
	      cid ) == CIDList.at(type).end() ) 
      return false;
    else return true;
  }
  else 
    return false;
}

int DetectorList::GetData(const int &type,const unsigned int &cid) const
{
  if(type<0||type>=ntype) return -1;
  std::map<unsigned int,int>::const_iterator it = DataContainer.at(type).find(cid);
  if(it != DataContainer.at(type).end()) return it->second;
  else return -1;
}

bool DetectorList::GetFlag(const unsigned int &cid,int type) const
{
  if(type<0&&type>1) return false;
  std::map<unsigned int,int>::const_iterator it = FlagContainer.at(type).find(cid);
  if(it != FlagContainer.at(type).end()) return it->second;
  else return false;
}

TString DetectorList::GetString(const int &type,const unsigned int &cid) const 
{
  std::map<unsigned int,TString>::const_iterator it = StringContainer.at(type).find(cid);
  if(it != StringContainer.at(type).end()) return it->second;
  else return "NULL";
}
unsigned DetectorList::GetCID(const TString  &name) const
{
  std::map<TString,unsigned int>::const_iterator it = CIDContainer.find(name);
  if(it != CIDContainer.end()) return it->second;
  else return -1;
}
