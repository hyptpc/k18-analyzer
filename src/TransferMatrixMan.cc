// TransferMatrixMan.cpp

#include <string>
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <new>
#include "TransferMatrixMan.hh"

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TransferMatrixMan::TransferMatrixMan()
  :m_isready(false)
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TransferMatrixMan::TransferMatrixMan( const TString & filename )
{
  FileName = filename;
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TransferMatrixMan::TransferMatrixMan( const TransferMatrixMan &right )
{
  FileName=right.FileName;
  CentralMomentum=right.CentralMomentum;
}
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
TransferMatrixMan::~TransferMatrixMan()
{
  Clear();
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void TransferMatrixMan::Clear()
{
  CentralMomentum=0.;
  for(int i=0;i<36;i++){
    D5Matrix1st[i]=0;
    BLC1VIele[i]=0;
  }
  BLC2BLC1Matrix1st=0;
  for(int i=0;i<6;i++){
    for(int j=0;j<36;j++)
      D5Matrix2nd[i][j]=0.;
    BLC2BLC1Matrix2nd[i]=0;
  }
  BLC2VIMatrix=0;
  BLC1VIMatrix=0;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool TransferMatrixMan::Initialize( const TString & filename )
{
  SetFileName(filename);
  return Initialize();
}
const int MAXCHAR = 144;
// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
bool TransferMatrixMan::Initialize()
{
  static const TString funcname = "TransferMatrixMan::Initialize";
  std::cout << "[" << funcname << "] Initialization start ..."/* << std::endl*/;

  int nd;
  double temp;
  int tmp;
  double mat[36];
  char str[MAXCHAR];
  //  int ar[8];
  FILE *fp;

  if( (fp=fopen(FileName.Data(), "r"))==0 ){
    std::cerr << " File open fail. [" << FileName << "]" << std::endl;
    exit(-1);
    return false;
  }

  Clear();
  while( fgets(str,MAXCHAR,fp)!=0 ){
    if( str[0]=='#' ) continue;
    if( (nd=sscanf(str,"CentralMomentum: %lf", &temp)) == 1 ) {
      CentralMomentum=temp;
    }
    //    else if( str=="D5Matrix: start\n") {
    //      while( fgets(str,MAXCHAR,fp)!=0 || str=="D5Matrix: end\n" || i > 5){
    else if( (nd=sscanf(str,"D5Matrix: %d", &tmp)) == 1 ) {
      int iraw=0;
      while( fgets(str,MAXCHAR,fp)!=0 && iraw<6)
	if( (nd=sscanf(str,"%lf %lf %lf %lf %lf %lf", &mat[6*iraw+0],&mat[6*iraw+1],&mat[6*iraw+2],&mat[6*iraw+3],&mat[6*iraw+4],&mat[6*iraw+5])) == 6 ){	  
	  iraw++;
	}
      if(tmp==1){
	for(int i=0;i<36;i++) D5Matrix1st[i]=mat[i];
	D5Matrix.Use(6,6,D5Matrix1st);
      }
      else if(tmp>20&&tmp<26)
	for(int i=0;i<36;i++) D5Matrix2nd[tmp-21][i]=mat[i];
    }    
    else if( (nd=sscanf(str,"BHDBLC2Matrix: %d", &tmp)) == 1 ) {
      int iraw=0;
      while( fgets(str,MAXCHAR,fp)!=0 && iraw<6)
	if( (nd=sscanf(str,"%lf %lf %lf %lf %lf %lf", &mat[6*iraw+0],&mat[6*iraw+1],&mat[6*iraw+2],&mat[6*iraw+3],&mat[6*iraw+4],&mat[6*iraw+5])) == 6 ){	  
	  iraw++;
	}
      if(tmp==1){
	BLC2VIMatrix.Use(6,6,mat);
	BLC2VIMatrix.Invert();
      }    
    }
    else if( (nd=sscanf(str,"BHDBLC1Matrix: %d", &tmp)) == 1 ) {
      int iraw=0;
      while( fgets(str,MAXCHAR,fp)!=0 && iraw<6)
	if( (nd=sscanf(str,"%lf %lf %lf %lf %lf %lf", &mat[6*iraw+0],&mat[6*iraw+1],&mat[6*iraw+2],&mat[6*iraw+3],&mat[6*iraw+4],&mat[6*iraw+5])) == 6 ){	  
	  iraw++;
	}
      if(tmp==1){
	for(int i=0;i<36;i++) BLC1VIele[i]=mat[i];
      }    
    }
    else{
      std::cerr << "[" << funcname << "]: Invalid data format" << std::endl;
      std::cerr << std::string(str) << std::endl;
    }
  }
  BLC1VIMatrix.Use(6,6,BLC1VIele);
  BLC1VIMatrix.Invert();
  //  std::cout<<"VI to BLC1 matrix"<<std::endl;
  //  BLC1VIMatrix.Print();
  //  BLC1VIMatrix.Print();
  fclose(fp);
  //  std::cout <</* "[" << funcname << "] Initialization*/" finish." << std::endl;
  m_isready=true;
  return true;
}

// + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + -- + //
void TransferMatrixMan::CalcBLC1toBLC2(double *parblc1, double *parblc2,const int &order)
{
  //par x, x', y, y', mom
  double parin[6]={parblc1[0]/10.,parblc1[1],parblc1[2]/10.,parblc1[3],0.,parblc1[4]};
  //		   (parblc1[4]-CentralMomentum)/CentralMomentum*100.};
  TMatrixD in;
  in.Use(6,1,parin);
  TMatrixD out(6,1);
  out.Mult(D5Matrix,in);
  parblc2[0]=out[0][0]*10;
  parblc2[1]=out[1][0];
  parblc2[2]=out[2][0]*10;
  parblc2[3]=out[3][0];
  parblc2[4]=out[4][0]*10;
}
void TransferMatrixMan::CalcBLC2toBLC1(double *parblc2, double *parblc1,const int &order)
{}

void TransferMatrixMan::CalcBLC1toVI(double *parblc1, double *parvi,const int &order)
{
  //par x, x', y, y', mom
  //  double parin[6]={parblc1[0],parblc1[1],parblc1[2],parblc1[3],0.,(parblc1[4]-CentralMomentum)/CentralMomentum*100.};
  double parin[6]={parblc1[0],parblc1[1],parblc1[2],parblc1[3],0.,parblc1[4]};
  TMatrixD in;
  in.Use(6,1,parin);
  //  in.Print();
  TMatrixD out(6,1);
  //  BLC1VIMatrix.Print();
  out.Mult(BLC1VIMatrix,in);
  //  out.Print();
  //  parvi=out.GetMatrixArray();
  parvi[0]=out[0][0];
  parvi[1]=out[1][0];
  parvi[2]=out[2][0];
  parvi[3]=out[3][0];
  parvi[4]=out[4][0];
}

void TransferMatrixMan::CalcBLC2toVI(double *parblc2, double *parvi,const int &order)
{
  //par x, x', y, y', mom
  double parin[6]={parblc2[0],parblc2[1],parblc2[2],parblc2[3],0.,parblc2[4]};
  TMatrixD in;
  in.Use(6,1,parin);
  TMatrixD out(6,1);
  out.Mult(BLC2VIMatrix,in);
  parvi[0]=out[0][0];
  parvi[1]=out[1][0];
  parvi[2]=out[2][0];
  parvi[3]=out[3][0];
  parvi[4]=out[4][0];
}
