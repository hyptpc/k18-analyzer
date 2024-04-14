// HistTools.h

#ifndef HistTools_h
#define HistTools_h 1

#include <string>
#include <vector>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TMacro.h"
#include "TSystem.h"
#include "TString.h"

namespace hist
{
  inline void H1(TString name,  double val,
		 int nbin, double lbin, double ubin);
  inline void H1(TString name,  double val,
		 int nbin, double lbin, double ubin, std::vector<TString> add);
  inline void H1(TString name,  double val,
		 int nbin, double lbin, double ubin, double weight);
  inline void H1(TString name,  double val, double bins[3]);
  inline void H1(TString name,  double val, double bins[3], std::vector<TString> add);
  inline void H2(TString name,  double val1, double val2,
		 int nbinx, double lbinx, double ubinx,
		 int nbiny, double lbiny, double ubiny);
  inline void H2(TString name,  double val1, double val2,
		 int nbinx, double lbinx, double ubinx,
		 int nbiny, double lbiny, double ubiny, std::vector<TString> add);
  inline void H2(TString name,  double val1, double val2,double bins1[3],double bins2[3]);
  inline void H2(TString name,  double val1, double val2,double bins1[3],double bins2[3], std::vector<TString> add);
  inline void H2(TString name,  double val1, double val2,double bins[6]);
  inline void H2(TString name,  double val1, double val2,double bins[6], std::vector<TString> add);

  // inline void WriteFile(TString afile,TString direntry);
  // inline void WriteFile(TString afile);
  // inline void SaveFile(TString filename);
  // inline TFile* OpenFile(TString filename);
  // inline void Print(double* param,int n, bool RETURN=true);
};
// inline void hist::Print(double* param, int n, bool RETURN){
//   for(int i=0;i<n;i++)
//     std::cout<<param[i]<<"  ";
//   if(RETURN)   std::cout<<std::endl;
// }

// inline TFile* hist::OpenFile(TString filename){
//   TFile *f = new TFile(filename);
//   if(!f->IsOpen()){
//     std::cout<<" failed to open " <<filename<< "  !!!"<<std::endl;
//     return 0;
//     //    exit(false);
//   }
//   if(f->IsZombie()){
//     std::cout<<filename<< " is Zombie !!!"<<std::endl;
//     return 0;
//     //    exit(false);
//   }
//   return f;
// }

// inline void hist::WriteFile(TString afile,TString direntry){
//   if(afile==DefaultFileName) return;
//   //  std::cout<<afile<<std::endl;
//   TMacro *m = new TMacro(afile);
//   m->Write(direntry);
//   delete m;
//   return;
// }

// inline void hist::WriteFile(TString afile){
//   if(afile==DefaultFileName) return;
//   //  std::cout<<afile<<std::endl;
//   TMacro *m = new TMacro(afile);
//   m->Write(afile);
//   delete m;
//   return;
// }

// inline void hist::SaveFile(TString filename){
//   std::cout<<filename<<std::endl;
//   const char *dirname="file";
//   char *slash = (char*)strrchr(dirname,'/');
//   char *locdir;
//   if (slash) locdir = slash+1;
//   else       locdir = (char*)dirname;
//   TDirectory *savdir = gDirectory;
//   TDirectory *adir = savdir->GetDirectory(locdir);
//   if(!adir) adir= savdir->mkdir(locdir);
//   if(adir)  adir->cd();
//   else{
//     std::cout<<"!!!"<<std::endl;
//     return;
//   }
//   //  gSystem->cd(locdir);
//   //  void *dirp = gSystem->OpenDirectory(dirname);
//   //  if (!dirp) return;
//   //  TString afile = Form("%s/%s",dirname,filename.Data());
//   hist::WriteFile( filename);
//   //  gSystem->FreeDirectory(dirp);
//   savdir->cd();
//   return;
// }

inline void hist::H1(TString name, double val,
		     int nbinx, double lbinx, double ubinx)
{
  if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(name.Data())) ){
    //    std::cout<<"exist "<<name<<std::endl;
    h1->Fill(val);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1 = new TH1F(name,name, nbinx,lbinx,ubinx);
    h1->Fill(val);
  }
}
inline void hist::H1(TString name, double val,
		      int nbinx, double lbinx, double ubinx, std::vector<TString> add)
{
  for(int i=0;i<(int)add.size();i++){
    TString tmpname=name+add[i];
    if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(tmpname.Data())) ){
      //    std::cout<<"exist "<<name<<std::endl;
      h1->Fill(val);
    }else{
      //    std::cout<<"new "<<name<<std::endl;
      h1 = new TH1F(tmpname,tmpname, nbinx,lbinx,ubinx);
      h1->Fill(val);
    }
  }
}
inline void hist::H1(TString name, double val,
		      int nbinx, double lbinx, double ubinx,
		      double weight)
{
  if( TH1F* h1 = dynamic_cast<TH1F *>(gFile->Get(name.Data())) ){
    //    std::cout<<"exist "<<name<<std::endl;
    h1->Fill(val,weight);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1 = new TH1F(name,name, nbinx,lbinx,ubinx);
    h1->Fill(val,weight);
  }
}
inline void hist::H1(TString name, double val,double bins[3])
{
  H1(name,val,int(bins[0]),bins[1],bins[2]);
}
inline void hist::H1(TString name, double val,double bins[3], std::vector<TString> add)
{
  H1(name,val,int(bins[0]),bins[1],bins[2], add);
}
inline void hist::H2(TString name, double val1, double val2,
		      int nbinx, double lbinx, double ubinx,
		      int nbiny, double lbiny, double ubiny)
{
  if( TH2F* h1 = dynamic_cast<TH2F *>(gFile->Get(name.Data())) ){
    h1->Fill(val1,val2);
  }else{
    //    std::cout<<"new "<<name<<std::endl;
    h1= new TH2F(name,name, nbinx,lbinx,ubinx, nbiny,lbiny,ubiny); 
    h1->Fill(val1,val2);
  }  
}
inline void hist::H2(TString name, double val1, double val2,
		      int nbinx, double lbinx, double ubinx,
		      int nbiny, double lbiny, double ubiny, std::vector<TString> add)
{
  for(int i=0;i<(int)add.size();i++){
    TString tmpname=name+add[i];
    if( TH2F* h1 = dynamic_cast<TH2F *>(gFile->Get(tmpname.Data())) ){
      h1->Fill(val1,val2);
    }else{
      //    std::cout<<"new "<<name<<std::endl;
      h1= new TH2F(tmpname,tmpname, nbinx,lbinx,ubinx, nbiny,lbiny,ubiny); 
      h1->Fill(val1,val2);
    }  
  }
}
inline void hist::H2(TString name, double val1, double val2, double bins[6])
{
  H2(name,val1,val2,int(bins[0]),bins[1],bins[2],int(bins[3]),bins[4],bins[5]);
}
inline void hist::H2(TString name, double val1, double val2, double bins[6], std::vector<TString> add)
{
  H2(name,val1,val2,int(bins[0]),bins[1],bins[2],int(bins[3]),bins[4],bins[5], add);
}
inline void hist::H2(TString name, double val1, double val2, double bins1[3],double bins2[3])
{
  H2(name,val1,val2,int(bins1[0]),bins1[1],bins1[2],int(bins2[0]),bins2[1],bins2[2]);
}
inline void hist::H2(TString name, double val1, double val2, double bins1[3],double bins2[3], std::vector<TString> add)
{
  H2(name,val1,val2,int(bins1[0]),bins1[1],bins1[2],int(bins2[0]),bins2[1],bins2[2], add);
}
#endif
