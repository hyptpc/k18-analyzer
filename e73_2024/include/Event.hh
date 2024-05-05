#ifndef EVENT_HH
#define EVENT_HH 1

#include <iostream>
#include <vector>

#include <TNamed.h>
#include <TObject.h>
#include <TString.h>

//_____________________________________________________________________________
struct EventTriggerFlag : TNamed
{
  EventTriggerFlag(const TString& name);
  Double_t trigger_flag;
};


class HitWire : public TObject
{
public:
  HitWire(){}
  ~HitWire(){}
  HitWire(const int &l,const int &w,const int &n):
    Layer(l),Wire(w),Nth(n)
  {
  }

private:
  Short_t Layer;
  Short_t Wire;
  Short_t Nth;

public:
  Int_t layer() const { return Layer; }
  Int_t wire()  const { return Wire; }
  Int_t nth()   const { return Nth; }
  void Print()   { std::cout<<"layer,wire,nth: "<<layer()<<", "<<wire()<<", "<<nth()<<std::endl; }
  ClassDef(HitWire,1);
};

class CDCTrackHits : public TObject
{
public:
  CDCTrackHits(){}
  ~CDCTrackHits(){ Clear(); }

private:
  std::vector<HitWire> hit;
  Float_t Param[5];

public:
  void Clear() { hit.clear(); }
  void addhit(const int &layer,const int &wire,const int &nth)
  { hit.push_back(HitWire(layer,wire,nth)); }
  void set_param(double param[5]) { for(int i=0;i<5;i++) Param[i]=param[i]; }
  void get_param(double param[5]) { for(int i=0;i<5;i++) param[i]=Param[i]; }
  int nhits() { return hit.size(); }
  HitWire get_hit(const int &i) { return i<nhits() ? hit[i] : HitWire(-1,-1,-1); }
  void Print() {
    std::cout<<"---------Track with nhits = "<<nhits()<<std::endl;
    for(int i=0;i<nhits();i++) get_hit(i).Print();
  }
  ClassDef(CDCTrackHits,1);
};

class CDCTrackContainer : public TObject
{
private:
  std::vector<CDCTrackHits> trackContainer;
public:
  CDCTrackContainer(){}
  ~CDCTrackContainer(){ Clear(); }

  void Clear(){ trackContainer.clear(); }
  void add_track(const CDCTrackHits& tra){ trackContainer.push_back(tra); }
  int ntracks() { return trackContainer.size(); }
  CDCTrackHits get_track(const int &i){ return i<ntracks() ? trackContainer[i]: CDCTrackHits(); }
  void Print() {
    std::cout<<"==== Event with ntracks = "<<ntracks()<<std::endl;
    for(int i=0;i<ntracks();i++) get_track(i).Print();
  }
  ClassDef(CDCTrackContainer,1);
};

struct HodoEvent
{
  int cid;
  int seg;
  double time;
  double ctime;
  double deu;
  double ded;
  double tu;
  double td;
  double ctu;
  double ctd;
};

namespace event_tree{
  inline TString get_slewleaf();
  inline TString get_dleaf();
  inline TString get_cdhleaf();
}

struct SlewEvent
{
  Float_t  mt_t0new,  mt_bht,  mt_t0,  mt_veto,  mt_def,  mt_btc;
  Float_t cmt_t0new, cmt_bht, cmt_t0, cmt_veto, cmt_def, cmt_btc;
  Float_t deu_t0new, deu_bht, deu_t0, deu_veto, deu_def, deu_btc;
  Float_t ded_t0new, ded_bht, ded_t0, ded_veto, ded_def, ded_btc;
  Short_t seg_t0new, seg_bht, seg_t0, seg_veto, seg_def, seg_btc;
};

inline TString event_tree::get_slewleaf()
{
  TString slewleaf="";
  slewleaf+="mt_t0new/F:mt_bht:mt_t0:mt_veto:mt_def:mt_btc";
  slewleaf+=":cmt_t0new/F:cmt_bht:cmt_t0:cmt_veto:cmt_def:cmt_btc";
  slewleaf+=":deu_t0new/F:deu_bht:deu_t0:deu_veto:deu_def:deu_btc";
  slewleaf+=":ded_t0new/F:ded_bht:ded_t0:ded_veto:ded_def:ded_btc";
  slewleaf+=":seg_t0new/S:seg_bht:seg_t0:seg_veto:seg_def:seg_btc";
  return slewleaf;
}


struct DeuteronEvent
{
  Float_t time_t0new,  time_bht,  time_t0,  time_veto,  time_def,  time_btc;
  Float_t de_t0new, de_bht, de_t0, de_veto, de_def, de_btc;
  Float_t de_rc[8];
  Short_t seg_t0new, seg_bht, seg_t0, seg_veto, seg_def, seg_btc;
};

inline TString event_tree::get_dleaf()
{
  TString dleaf="";
  dleaf+="time_t0new/F:time_bht:time_t0:time_veto:time_def:time_btc";
  dleaf+=":de_t0new/F:de_bht:de_t0:de_veto:de_def:de_btc";
  dleaf+=":de_rc[8]";
  dleaf+=":seg_t0new/S:seg_bht:seg_t0:seg_veto:seg_def:seg_btc";
  return dleaf;
}

struct CDHEvent
{
  Float_t  mt_t0,  mt_cdh,  mt_def;
  Float_t cmt_t0, cmt_cdh, cmt_def;
  Float_t tu_cdh, td_cdh;
  Float_t ctu_cdh, ctd_cdh;
  Float_t deu_t0, deu_cdh, deu_def;
  Float_t ded_t0, ded_cdh, ded_def;
  Float_t tof_bhtt0, ctof_bhtt0, fl_beam, fl_cdh;
  Float_t mom_cdc,fl_cdc,de_cdc,dt_cdc,dt_cdc2,cdt_cdc;
  Float_t dca_cdc,beta_cdc,mass2_cdc,chi2_cdc;
  Float_t vcdc_x, vcdc_y, vcdc_z;
  Float_t vcdh_phi_in, vcdh_z_in;
  Float_t vcdh_phi_out, vcdh_z_out;
  Float_t vbpc_x, vbpc_y, vbpc_z;
  Short_t seg_t0, seg_cdh, seg_def, pid_cdc, flag_fiducial;
};

inline TString event_tree::get_cdhleaf()
{
  TString leaf="";
  leaf+="mt_t0/F:mt_cdh:mt_def";
  leaf+=":cmt_t0:cmt_cdh:cmt_def";
  leaf+=":tu_cdh:td_cdh";
  leaf+=":ctu_cdh:ctd_cdh";
  leaf+=":deu_t0:deu_cdh:deu_def";
  leaf+=":ded_t0:ded_cdh:ded_def";
  leaf+=":tof_bhtt0:ctof_bhtt0:fl_beam:fl_cdh";
  leaf+=":mom_cdc:fl_cdc:de_cdc:dt_cdc:dt_cdc2:cdt_cdc";
  leaf+=":dca_cdc:beta_cdc:mass2_cdc:chi2_cdc";
  leaf+=":vcdc_x:vcdc_y:vcdc_z";
  leaf+=":vcdh_phi_in:vcdh_z_in";
  leaf+=":vcdh_phi_out:vcdh_z_out";
  leaf+=":vbpc_x:vbpc_y:vbpc_z";
  leaf+=":seg_t0/S:seg_cdh:seg_def:pid_cdc:flag_fiducial";
  return leaf;
}

#endif
