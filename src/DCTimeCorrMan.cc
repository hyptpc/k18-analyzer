// DCTimeCorrMan.cpp

// Debug flag: set to 1 to enable debug output, 0 to disable
#define DC_TIME_CORR_DEBUG 0

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <new>
#include <string>
#include <utility>
#include <vector>

#include <TFile.h>
#include <TROOT.h>

#include "DCTimeCorrMan.hh"
#include "DetectorID.hh"

namespace
{
const UInt_t KEYMASK  = 0x000F;
// |0111|1111|0001|1111|0000|1111|1111|0011|
const UInt_t WMASK    = 0x00FF;      /* Wire Mask 8 Bits (0-255) */
const UInt_t LMASK    = 0x001F;      /* Layer Mask 5 Bits (0-31) */
const UInt_t CMASK    = 0x007F;      /* CID Mask 7 Bits (0-127) */
const Int_t  WSHIFT   =  4;
const Int_t  LSHIFT   = 16;
const Int_t  CSHIFT   = 24;
const UInt_t KEYFLAG  = 0x0003;
inline Int_t KEY( Int_t cid, Int_t layer, Int_t wire){
  return ( ( (cid&CMASK)<<CSHIFT) |
           ( (layer&LMASK)<<LSHIFT) |
           ( (wire&WMASK)<<WSHIFT) |
           KEYFLAG );
}
const Int_t MAXCHAR = 144;
const Int_t MaxParam = 4;

const std::vector<std::pair<TString, Int_t>> DCChambers = {
  {"BLC1a", DetIdBLC1a},
  {"BLC1b", DetIdBLC1b},
  {"BLC2a", DetIdBLC2a},
  {"BLC2b", DetIdBLC2b}
};
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
Bool_t DCTimeCorrMan::Initialize( const char *filename1, const char *filename2 )
{
  FileNameBLDC=filename1;
  FileNameCDC=filename2;
  return Initialize();
}
Bool_t DCTimeCorrMan::Initialize( const TString& filename1, const TString& filename2 )
{
  FileNameBLDC=filename1;
  FileNameCDC=filename2;
  return Initialize();
}
Bool_t DCTimeCorrMan::Initialize( const char *filename1 )
{
  FileNameBLDC=filename1;
  return Initialize();
}
Bool_t DCTimeCorrMan::Initialize( const TString& filename1 )
{
  FileNameBLDC=filename1;
  return Initialize();
}
Bool_t DCTimeCorrMan::Initialize()
{
  static const TString funcname = "DCTimeCorrMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ...";


  dctimecorrContainer.clear();
  TFile *tmpfile=gFile;
  TFile *f=new TFile(FileNameBLDC.Data());
  if(!f->IsOpen()){
    std::cout<<"!!!!!\t"<<FileNameBLDC<<" does not exist !!!!!"<<std::endl;
    FileNameBLDC=DefaultFileName;
    return true;
  }
  for(const auto& chamber : DCChambers){
    (void)chamber;
    for(Int_t lay=0;lay<8;lay++){
#if DC_TIME_CORR_DEBUG
      TGraph* gr = (TGraph*)f->Get(Form("DCTimeCorr_%s_%d",chamber.first.Data(),lay+1));
      if(!gr) continue;
      std::cout << "["<<funcname<<"] " 
                << Form("DCTimeCorr_%s_%d",chamber.first.Data(),lay+1) << std::endl;
      Int_t key = KEY(chamber.second,lay,0);
      dctimecorrContainer[key] = *gr;
      gROOT->Append(gr);
#endif
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
Double_t DCTimeCorrMan::CalcCValue( const Int_t &cid, const Int_t &layer, const Int_t &wire,
				  const Double_t &timemean, const Double_t &timesub ) const
{
  static const TString funcname = "DCTimeCorrMan::CalcCValue";
  Double_t ctime = timemean;
  Int_t key;
  key = KEY(cid,layer,wire);
  DCTimeCorrHistContainer::const_iterator is;
  if( (is=(dctimecorrContainer.find(key))) != dctimecorrContainer.end() ){
    ctime = timemean - ((is->second).Eval(timesub));
#if DC_TIME_CORR_DEBUG
    if(abs(timesub)>100){
       std::cout<<cid<<"\t"<<(is->second).GetName()<<std::endl;
       std::cout<<"sub,mean,corr,ctime\t"<<timesub<<"\t"<<timemean<<"\t"<<(is->second).Eval(timesub)<<"\t"<<ctime<<std::endl;
    }
#endif
  }else{
    std::cout << " cannot find parameters "
              << " cid:" << cid << " layer:" << layer << " wire:"<< wire <<" timemean:" << timemean << " timesub:" << timesub
              << std::endl;
  }
  return ctime;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
Double_t DCTimeCorrMan::CalcDATValue( const Int_t &cid, const Int_t &seg, const Int_t &ud,
				    const Double_t &ctime, const Double_t &de ) const
{
  static const TString funcname = "DCTimeCorrMan::CalcDATValue";
  Double_t time = ctime;
  return time;
}
