// DCTimeCorrMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>
#include <cmath>

#include "DCTimeCorrMan.hh"
#include "DetectorID.hh"

namespace
{
const unsigned int KEYMASK  = 0x000F;
// |0111|1111|0001|1111|0000|1111|1111|0011|
const unsigned int WMASK    = 0x00FF;      /* Wire Mask 8 Bits (0-255) */
const unsigned int LMASK    = 0x001F;      /* Layer Mask 5 Bits (0-31) */
const unsigned int CMASK    = 0x007F;      /* CID Mask 7 Bits (0-127) */
const int          WSHIFT   =  4;
const int          LSHIFT   = 16;
const int          CSHIFT   = 24;
const unsigned int KEYFLAG  = 0x0003;
inline int KEY( int cid, int layer, int wire){
  return ( ( (cid&CMASK)<<CSHIFT) |
           ( (layer&LMASK)<<LSHIFT) |
           ( (wire&WMASK)<<WSHIFT) |
           KEYFLAG );
}
const int MAXCHAR = 144;
const int MaxParam = 4;

const int nchamber=6;
// TString chmname[nchamber]={"BLC1a","BLC1b","BLC2a","BLC2b","BPC1","BPC2"};
// int chmid[nchamber]={DetIdBLC1a,DetIdBLC1b,DetIdBLC2a,DetIdBLC2b,DetIdBPC1,DetIdBPC2};
TString DefaultFileName="default";
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
DCTimeCorrMan::DCTimeCorrMan()
  :m_isready(false)
{
  FileNameCDC  = "default";
  FileNameBLDC = "default";
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
DCTimeCorrMan::DCTimeCorrMan( const DCTimeCorrMan &right )
{
  FileNameCDC = right.FileNameCDC;
  FileNameBLDC = right.FileNameBLDC;
  for( DCTimeCorrHistContainer::const_iterator i=right.dctimecorrContainer.begin();
       i!=right.dctimecorrContainer.end(); i++ ){
    dctimecorrContainer[i->first] = i->second;
  }
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
DCTimeCorrMan::~DCTimeCorrMan()
{
  dctimecorrContainer.clear();
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void DCTimeCorrMan::SetFileNames( const TString & filename1, const TString & filename2 )
{
  FileNameBLDC = filename1;
  FileNameCDC = filename2;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool DCTimeCorrMan::Initialize( const char *filename1, const char *filename2 )
{
  FileNameBLDC=filename1;
  FileNameCDC=filename2;
  return Initialize();
}
bool DCTimeCorrMan::Initialize( const std::string& filename1, const std::string& filename2 )
{
  FileNameBLDC=filename1;
  FileNameCDC=filename2;
  return Initialize();
}
bool DCTimeCorrMan::Initialize( const char *filename1 )
{
  FileNameBLDC=filename1;
  return Initialize();
}
bool DCTimeCorrMan::Initialize( const std::string& filename1 )
{
  FileNameBLDC=filename1;
  return Initialize();
}
bool DCTimeCorrMan::Initialize()
{
  static const TString funcname = "DCTimeCorrMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";

  int key;
  dctimecorrContainer.clear();
  TFile *tmpfile=gFile;
  TFile *f=new TFile(FileNameBLDC.Data());
  if(!f->IsOpen()){
    std::cout<<"!!!!!\t"<<FileNameBLDC<<" does not exist !!!!!"<<std::endl;
    FileNameBLDC=DefaultFileName;
    return true;
  }
  TGraph* gr;
  for(int ic=0;ic<nchamber;ic++){
    for(int lay=0;lay<8;lay++){
      // gr=(TGraph*)f->Get(Form("DCTimeCorr_%s_%d",NameDC[ic].Data(),lay+1));
      // if(!gr) continue;
      // //      std::cout<<Form("DCTimeCorr_%s_%d",NameDC[ic].Data(),lay+1)<<std::endl;
      // key = KEY(DetIdDC[ic],lay,0);
      // dctimecorrContainer[key] = *gr;
      // gROOT->Append(gr);
    }
  }
  f->Close();
  delete f;
  tmpfile->cd();
  std::cout << " finish." << std::endl;
  m_isready=true;
  return true;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double DCTimeCorrMan::CalcCValue( const int &cid, const int &layer, const int &wire,
				  const double &timemean, const double &timesub ) const
{
  static const TString funcname = "DCTimeCorrMan::CalcCValue";
  double ctime = timemean;
  int key;
  key = KEY(cid,layer,wire);
  DCTimeCorrHistContainer::const_iterator is;
  if( (is=(dctimecorrContainer.find(key))) != dctimecorrContainer.end() ){
    ctime = timemean - ((is->second).Eval(timesub));
    // if(abs(timesub)>100){
    //    std::cout<<cid<<"\t"<<(is->second).GetName()<<std::endl;
    //    std::cout<<"sub,mean,corr,ctime\t"<<timesub<<"\t"<<timemean<<"\t"<<(is->second).Eval(timesub)<<"\t"<<ctime<<std::endl;
    // }
    // }else if(cid==DetIdFDC){
    //   ctime=timemean;
  }else{
    std::cout << " cannot find parameters "
	      << " cid:" << cid << " layer:" << layer << " wire:"<< wire <<" timemean:" << timemean << " timesub:" << timesub
	      << std::endl;
  }
  return ctime;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
double DCTimeCorrMan::CalcDATValue( const int &cid, const int &seg, const int &ud,
				    const double &ctime, const double &de ) const
{
  static const TString funcname = "DCTimeCorrMan::CalcDATValue";
  double time = ctime;
  return time;
}
