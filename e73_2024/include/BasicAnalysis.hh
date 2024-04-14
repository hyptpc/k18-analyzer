#ifndef BASICANALYSIS_HH
#define BASICANALYSIS_HH 1

#include <vector>

#include "RawData.hh"
#include "UserParamMan.hh"
#include "HodoAnalyzer.hh"
#include "DCAnalyzer.hh"
#include "CDCAnalyzer.hh"
#include "MTDCAnalyzer.hh"
#include "Particle.hh"

#include "setup.hh"
#include "Event.hh"

#include "TString.h"

using namespace e73_2024;
class BasicAnalysis
{
public:
  BasicAnalysis();
  ~BasicAnalysis() {};

private:
  int EventNumber;
  double Time0;
  double CTime0;
  double Time1;
  double CTime1;
  double TOF_BHTT0;
  double cTOF_BHTT0;
  double PbF2_AdcSum;
  double PbF2_EnergySum;
  double PbF2_EnergySumwT;
  double PbF2_EnergySum2;
  double PbF2_EnergySum2wT;
  double PbF2_EnergySum3;
  double PbF2_EnergySum4;
  double PbF2_EnergyMax;
  int    PbF2_Seg;
  double PbF2_Time;
  int sign_beam;
  bool TOFK;
  bool TOFP;
  bool TOFPi;
  bool ACHIT;
  bool isMC;
  bool isYamaga;
  bool TriggerFlag[32];
  int TriggerMode;
  TVector3 position[kNumHodo];
  std::vector<HodoEvent> hodo_container[kNumHodo];
  std::vector<int> good_trackid[kNumChm];
  std::vector<int> neutralCDH;

public:
  double time0()  { return Time0; }
  double ctime0() { return CTime0; }
  double time1()  { return Time1; }
  double ctime1() { return CTime1; }
  double tof_bhtt0() { return TOF_BHTT0; }
  double ctof_bhtt0() { return cTOF_BHTT0; }
  double pbf2_adcsum() { return PbF2_AdcSum; }
  double pbf2_enesum() { return PbF2_EnergySum; }
  double pbf2_enesumwt() { return PbF2_EnergySumwT; }
  double pbf2_enesum2() { return PbF2_EnergySum2; }
  double pbf2_enesum2wt() { return PbF2_EnergySum2wT; }
  double pbf2_enesum3() { return PbF2_EnergySum3; }
  double pbf2_enesum4() { return PbF2_EnergySum4; }
  double pbf2_emax() { return PbF2_EnergyMax; }
  int    pbf2_seg()    { return PbF2_Seg; }
  double pbf2_time()   { return PbF2_Time; }
  bool   tofk()   { return TOFK; }
  bool   tofp()   { return TOFP; }
  bool   tofpi()  { return TOFPi; }
  bool   achit()  { return ACHIT; }
  TVector3 pos(int khodo) { return position[khodo]; }
  bool   trig(int i)   { return TriggerFlag[i]; }
  int    nhit(int khodo) { return hodo_container[khodo].size(); }
  int    hodoseg(  int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].seg : -1; }
  double hodotime( int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].time: -999; }
  double hodoctime(int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].ctime: -999;}
  double hododeu(  int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].deu : -999; } 
  double hododed(  int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].ded : -999; }
  double hodode(   int khodo,int i=0) { return nhit(khodo)>i ? sqrt(hododeu(khodo,i)*hododed(khodo,i)) : -999; }

  double hodotu(  int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].tu : -999; } 
  double hodotd(  int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].td : -999; }
  double hodoctu(  int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].ctu : -999; } 
  double hodoctd(  int khodo,int i=0) { return nhit(khodo)>i ? hodo_container[khodo][i].ctd : -999; }

  double ngood(  int kchm          ){ return good_trackid[kchm].size(); }
  double goodid( int kchm, int i=0 ){ return good_trackid[kchm][i]; }
  double ncdh_neutral(){ return neutralCDH.size(); }
 
  void SetEvNum(int evnum){ EventNumber=evnum; }

  bool AnaInit(MTDCAnalyzer* mtdc, HodoAnalyzer *hodo,bool mc=false,bool slewing=false,bool yamaga=true);
  bool AnaFlag(MTDCAnalyzer* mtdc);
  bool CalcTime0(HodoAnalyzer* hodo);
  bool CheckBHTT0TOF(HodoAnalyzer* hodo, bool slewing=false);
  bool FillHodoRaw(RawData* raw,int kHodo);
  bool FillHodoDecoded(HodoAnalyzer* hodo,int kHodo,std::vector<TString> add);
  bool FillHodo1Decoded(HodoAnalyzer* hodo,int kHodo,std::vector<TString> add);
  bool FillBLDC(DCAnalyzer* raw){ return true; }
  bool FillCDC(RawData* raw){ return true; }
  bool AnaAC(RawData* raw);
  bool AnaDC(DCAnalyzer* dc, int kDC, double tmin=-30, double tmax=100, double chi2=100);
  bool AnaHodo(HodoAnalyzer* hodo, int kHodo);
  bool AnaBHT(HodoAnalyzer* hodo, int kHodo=kBHT);
  bool AnaHodo1(HodoAnalyzer* hodo, int kHodo);
  bool AnaPbF2(HodoAnalyzer* hodo,std::vector<TString> add);
  void SetPos(int khodo, TVector3 tmppos) { position[khodo]=tmppos; } 
  bool AnaCDHneutral(Particle* particle);

  SlewEvent get_slewevent();  
  DeuteronEvent get_devent(HodoAnalyzer* hodo);  
  CDHEvent get_cdhevent();  
};
#endif
